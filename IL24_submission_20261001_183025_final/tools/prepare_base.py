#!/usr/bin/env python3
"""Prepare the pinned official ESM-2 base; verify against the independently frozen A reference."""
import argparse, gc, hashlib, importlib.metadata, json, shutil, struct, sys, tempfile
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / 'models/base_manifest.json'


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'), ensure_ascii=True).encode()


def file_sha(path):
    result = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            result.update(block)
    return result.hexdigest()


def manifest_and_reference():
    manifest = json.loads(MANIFEST.read_text())
    reference_path = ROOT / manifest['reference']['file']
    if file_sha(reference_path) != manifest['reference']['file_sha256']:
        raise ValueError('Frozen reference file checksum mismatch; never regenerate it from a candidate model')
    reference = json.loads(reference_path.read_text())
    for key in ('tensor_sha256', 'config_sha256', 'tokenizer_sha256'):
        if reference[key] != manifest['reference'][key]:
            raise ValueError('Manifest/reference disagreement: ' + key)
    return manifest, reference


def tensor_record(state):
    import torch
    if sys.byteorder != 'little':
        raise RuntimeError('Fingerprint v1 requires a little-endian runtime')
    digest = hashlib.sha256(b'IL24-ESM2-TENSORS-v1\0')
    entries = []
    for name, tensor in sorted(state.items()):
        metadata = {'name': name, 'dtype': str(tensor.dtype), 'shape': list(tensor.shape)}
        header = canonical(metadata)
        raw = tensor.detach().cpu().contiguous().view(torch.uint8).numpy().tobytes()
        digest.update(struct.pack('>Q', len(header)))
        digest.update(header)
        digest.update(struct.pack('>Q', len(raw)))
        digest.update(raw)
        entries.append({**metadata, 'tensor_sha256': hashlib.sha256(raw).hexdigest(), 'bytes': len(raw)})
    return digest.hexdigest(), entries


def inspect_loaded(model, tokenizer):
    # Include every state_dict key and all registered buffers (including nonpersistent position IDs).
    # Aliased keys are hashed separately by name; no floating-point casts or missing-key exclusions.
    state = dict(model.state_dict())
    state.update(dict(model.named_buffers(remove_duplicate=False)))
    tensor_sha, entries = tensor_record(state)
    aliases = {}
    for name, tensor in state.items():
        aliases.setdefault(tensor.untyped_storage().data_ptr(), []).append(name)
    config = model.config.to_dict()
    config.pop('_name_or_path', None)  # Location only; all other loaded configuration fields are retained.
    tokens = {'vocab': tokenizer.get_vocab(), 'special_tokens_map': tokenizer.special_tokens_map,
              'all_special_ids': tokenizer.all_special_ids, 'model_max_length': tokenizer.model_max_length,
              'padding_side': tokenizer.padding_side, 'truncation_side': tokenizer.truncation_side,
              'model_input_names': tokenizer.model_input_names}
    return {'tensor_sha256': tensor_sha, 'tensors': entries,
            'shared_storage_groups': sorted(sorted(v) for v in aliases.values() if len(v) > 1),
            'config': config, 'config_sha256': hashlib.sha256(canonical(config)).hexdigest(),
            'tokenizer': tokens, 'tokenizer_sha256': hashlib.sha256(canonical(tokens)).hexdigest()}


def assert_reference(actual, reference):
    if actual['tensors'] != reference['tensors']:
        expected = {item['name']: item for item in reference['tensors']}
        found = {item['name']: item for item in actual['tensors']}
        changed = [name for name in sorted(set(expected) | set(found)) if expected.get(name) != found.get(name)]
        raise ValueError('Base tensor keys/dtype/shape/values differ from frozen A: ' + ', '.join(changed[:12]))
    for key in ('tensor_sha256', 'config_sha256', 'tokenizer_sha256', 'shared_storage_groups'):
        if actual[key] != reference[key]:
            raise ValueError('Base differs from frozen A: ' + key)


def verify_converted(directory):
    import torch
    from transformers import EsmModel, EsmTokenizer, AutoTokenizer
    _, reference = manifest_and_reference()
    torch.set_num_threads(8)
    directory = Path(directory).resolve()
    model, info = EsmModel.from_pretrained(directory, add_pooling_layer=False,
                                         local_files_only=True, output_loading_info=True)
    if any(info.get(key) for key in ('missing_keys', 'unexpected_keys', 'mismatched_keys', 'error_msgs')):
        raise ValueError('Incomplete or incompatible base checkpoint: ' + str(info))
    model.eval()
    tokenizer = AutoTokenizer.from_pretrained(directory, local_files_only=True)
    if type(tokenizer) is not EsmTokenizer:
        raise ValueError('Actual AutoTokenizer loader must select the historical EsmTokenizer class')
    actual = inspect_loaded(model, tokenizer)
    assert_reference(actual, reference)
    del model, tokenizer
    gc.collect()
    return {'model_dir': str(directory), 'base_reference': {key: actual[key] for key in
            ('tensor_sha256', 'config_sha256', 'tokenizer_sha256')}}


def normalize_identity(identity):
    """Bind legacy project checkpoints to A without claiming B has A's serialized file SHA."""
    _, reference = manifest_and_reference()
    expected = {key: reference[key] for key in ('tensor_sha256', 'config_sha256', 'tokenizer_sha256')}
    if not identity:
        raise ValueError('Missing source model identity')
    if 'base_reference' in identity:
        if identity['base_reference'] != expected:
            raise ValueError('Unknown base reference in checkpoint identity')
        return expected
    legacy = {key: value for key, value in identity.items() if key != 'model_dir'}
    if legacy != reference['legacy_file_identity']:
        raise ValueError('Legacy checkpoint is not bound to the independently verified historical base')
    return expected


def source_files(directory, manifest, download=False):
    directory = Path(directory).resolve()
    if download:
        directory.mkdir(parents=True, exist_ok=True)
    for asset in manifest['source_files']:
        path = directory / asset['filename']
        if not path.exists() and download:
            temporary = path.with_suffix(path.suffix + '.partial')
            print('Downloading ' + asset['filename'], flush=True)
            with urllib.request.urlopen(asset['url'], timeout=120) as response, temporary.open('wb') as output:
                shutil.copyfileobj(response, output, length=8 * 1024 * 1024)
            if temporary.stat().st_size != asset['bytes'] or file_sha(temporary) != asset['sha256']:
                temporary.unlink()
                raise ValueError('Official download checksum mismatch: ' + asset['filename'])
            temporary.replace(path)
        if not path.is_file() or path.stat().st_size != asset['bytes'] or file_sha(path) != asset['sha256']:
            raise ValueError('Official source checksum mismatch or missing file: ' + str(path))
    return directory


def convert_native(source, target, manifest):
    import torch, esm
    from transformers import EsmConfig, EsmModel, EsmTokenizer
    for name, expected_version in manifest['conversion']['dependencies'].items():
        if importlib.metadata.version(name) != expected_version:
            raise RuntimeError(f'Conversion requires {name}=={expected_version}')
    torch.set_num_threads(8)
    torch.manual_seed(42)
    # Unpickle only the two official assets after their pinned SHA-256 checks have passed.
    native, alphabet = esm.pretrained.load_model_and_alphabet_core(
        manifest['model_name'],
        torch.load(source / manifest['source_files'][0]['filename'], map_location='cpu', weights_only=False),
        torch.load(source / manifest['source_files'][1]['filename'], map_location='cpu', weights_only=False))
    config_dir = ROOT / 'models/base_config'
    for asset in manifest['configuration_files']:
        if file_sha(config_dir / asset['filename']) != asset['sha256']:
            raise ValueError('Frozen configuration file changed: ' + asset['filename'])
    model = EsmModel(EsmConfig.from_pretrained(config_dir, local_files_only=True), add_pooling_layer=False)
    mapping = {'embed_tokens.weight': 'embeddings.word_embeddings.weight',
               'emb_layer_norm_after.weight': 'encoder.emb_layer_norm_after.weight',
               'emb_layer_norm_after.bias': 'encoder.emb_layer_norm_after.bias',
               'contact_head.regression.weight': 'contact_head.regression.weight',
               'contact_head.regression.bias': 'contact_head.regression.bias'}
    parts = {'self_attn.q_proj': 'attention.self.query', 'self_attn.k_proj': 'attention.self.key',
             'self_attn.v_proj': 'attention.self.value', 'self_attn.out_proj': 'attention.output.dense',
             'self_attn_layer_norm': 'attention.LayerNorm', 'final_layer_norm': 'LayerNorm',
             'fc1': 'intermediate.dense', 'fc2': 'output.dense',
             'self_attn.rot_emb.inv_freq': 'attention.self.rotary_embeddings.inv_freq'}
    expected = model.state_dict()
    weights, excluded = {}, []
    for key, value in native.state_dict().items():
        destination = mapping.get(key)
        if key.startswith('layers.'):
            _, index, suffix = key.split('.', 2)
            for old, new in parts.items():
                if suffix == old or suffix.startswith(old + '.'):
                    destination = 'encoder.layer.' + index + '.' + suffix.replace(old, new, 1)
                    break
        if destination in expected:
            if value.dtype != expected[destination].dtype or value.shape != expected[destination].shape:
                raise ValueError('Unexpected native dtype/shape: ' + key)
            weights[destination] = value
        elif key.startswith('lm_head.'):
            excluded.append(key)  # Original EsmModel conversion has no masked-language-model head.
        else:
            raise ValueError('Unmapped native key: ' + key)
    weights['embeddings.position_embeddings.weight'] = torch.zeros_like(expected['embeddings.position_embeddings.weight'])
    if set(weights) != set(expected):
        raise ValueError('Missing/extra converted state keys: ' + str(set(weights) ^ set(expected)))
    model.load_state_dict(weights, strict=True)
    model.eval()
    target.mkdir()
    model.save_pretrained(target, safe_serialization=True)
    for asset in manifest['configuration_files']:
        shutil.copy2(config_dir / asset['filename'], target / asset['filename'])
    tokenizer = EsmTokenizer.from_pretrained(target, local_files_only=True)
    if tokenizer.get_vocab() != {token: index for index, token in enumerate(alphabet.all_toks)}:
        raise ValueError('Official alphabet and frozen tokenizer differ')
    return {'excluded_official_lm_head_keys': excluded, 'conversion': 'Original key mapping; float32 retained; unused absolute-position table zero-filled; no training/quantization'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument('--official-dir', type=Path, help='Offline original Meta .pt plus contact-regression .pt; source hashes checked')
    mode.add_argument('--converted-dir', type=Path, help='Offline already-converted EsmModel directory; tensors/config/tokenizer checked')
    parser.add_argument('--cache-dir', type=Path, default=Path.home() / '.cache/il24-esm2-official')
    parser.add_argument('--output', type=Path, help='New converted model directory; default external_models/esm2_650m_hf')
    parser.add_argument('--verify-only', action='store_true', help='With --converted-dir: verify without writing anything')
    args = parser.parse_args()
    manifest, _ = manifest_and_reference()
    if args.converted_dir:
        if args.output:
            parser.error('--converted-dir reuses that directory; do not supply --output')
        report = verify_converted(args.converted_dir)
    else:
        if args.verify_only:
            parser.error('--verify-only requires --converted-dir')
        output = (args.output or ROOT / 'external_models/esm2_650m_hf').resolve()
        source = source_files(args.official_dir or args.cache_dir, manifest, download=args.official_dir is None)
        if output.exists():
            report = verify_converted(output)
        else:
            output.parent.mkdir(parents=True, exist_ok=True)
            with tempfile.TemporaryDirectory(prefix='.esm2-prepare-', dir=output.parent) as work:
                target = Path(work) / 'model'
                details = convert_native(source, target, manifest)
                report = verify_converted(target)
                report.update(details)
                target.rename(output)
            report['model_dir'] = str(output)
            (output / 'preparation_record.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({'status': 'verified', **report}, indent=2), flush=True)


if __name__ == '__main__':
    main()
