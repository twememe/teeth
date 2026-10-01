---
dataset_info:
  features:
  - name: vh
    dtype: string
  - name: label
    dtype:
      class_label:
        names:
          '0': non-binder
          '1': binder
  splits:
  - name: train
    num_bytes: 177739
    num_examples: 1310
  - name: eval
    num_bytes: 22104
    num_examples: 163
  - name: test
    num_bytes: 22462
    num_examples: 163
  download_size: 62689
  dataset_size: 222305
---
# Human IL-6 binding dataset

Nanobodies binding IL-6 were obtained from the [Github repo](https://github.com/cognano/AVIDa-hIL6) for [Tsuruta et al. (2023)](https://arxiv.org/abs/2306.03329). Labels for antibody sequences were provided from the Github repo as-is. 

Briefly, we first removed any nanobody sequence having lower than 75%; human germline sequence identity was determined using ANARCI. ![image](il6-binding-seqid.png)

Among the remaining 232084 sequences, we only use antibodies that have confirmed binding to one IL-6 variant or has no binding to any IL-6 variant, leading to 211920 sequences.
This was de-duplicated, leading to a total of 14,467 sequences. Non-binder sequences were randomly under-sampled using [imbalanced-learn](https://github.com/scikit-learn-contrib/imbalanced-learn).

During stratification, we also ensure that any two sequences in the training/validation/test sets have a minimum Levenshtein distance of 1 across the CDR3 region.

In total there are
* 1310 sequences in training
* 163 sequences in validation
* 163 sequences in test.

| vh_full | label |
| ------- | ----- |
| EVQ...  | 1     |
| QVQ...  | 0     |
| EVQ...  | 1     |
| EVK...  | 0     |

## License

The data is distributed under a CC-BY-NC 4.0 license.

## References

* [AVIDA-hIL6 publication](https://arxiv.org/abs/2306.03329)
* [ANARCI](https://academic.oup.com/bioinformatics/article/32/2/298/1743894)
