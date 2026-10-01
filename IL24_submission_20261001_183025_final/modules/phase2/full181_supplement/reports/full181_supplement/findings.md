# Findings

- The compact release documents Boltz software 2.2.1 with boltz1 model, Python 3.10.12, torch 2.6.0+cu124, CUDA 12.4, RTX A6000 and PRODIGY 2.4.0.
- Legacy formal_msa launcher forces CUDA_VISIBLE_DEVICES=0; supplement requires isolated GPU assignment.
- Legacy launchers skip on confidence JSON existence only; supplement must also validate input/parameters, PDB sequences, confidence and PAE.
- Legacy input generator writes both old constructs and glyco inputs on import; do not import or execute it for the supplement.
- Existing MSA generator creates paired and unpaired searches for antigen/VH/VL together; preserve this batch boundary. Its manifest computes hashes, so a supplement entry point must honor current no-hash instruction.
- Historical stored SHA256 values may be cited only as historical metadata, never as freshly verified. Use exact sequence equality, exact command and path/size/mtime metadata for current provenance; no hash integrity claim.

## Final findings supersede earlier resource blocks
Root authorization made original resources and60 legacy originals available. Full empty/formal_msa are complete15/15 each. N-terminal contact participation7/15 and1/15; neither condition has N-terminal consensus. Common27–160 Full/Native consensus Jaccard0.708333/0.857143; Full/Immunogen0.03125/0.571429. A→H/A→L chain-pair confidence means are only0.1750/0.1712(empty) and0.1979/0.1912(MSA), so global confidence cannot establish antigen-Fv binding reliability. AllPRODIGY values are computational proxies.
