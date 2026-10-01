# Measured validation comparison

All rows use the same validation partition unless explicitly marked test. Labels represent numerical high scores, not a unified biological positive direction.

| Model | N | MCC | F1 | Spearman (probability, label) |
|---|---:|---:|---:|---:|
| random | 20000 | -0.0033540315163079733 | 0.4806498344370861 | -0.006402301789296079 |
| frozen | 20000 | 0.3858637125275409 | 0.7100100603621731 | 0.4850030852910464 |
| lora | 20000 | 0.49908265893815706 | 0.7543091899615969 | 0.5971478136869867 |

Theoretical balanced random reference: F1 = 0.5, MCC = 0 (population expectations, not measurements).
See each metrics JSON for source and mutated-partner subgroups.
