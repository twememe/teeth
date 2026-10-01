# Measured validation comparison

All rows use the same validation partition unless explicitly marked test. Labels represent numerical high scores, not a unified biological positive direction.

| Model | N | MCC | F1 | Spearman (probability, label) |
|---|---:|---:|---:|---:|
| lora | 8 | 0.4472135954999579 | 0.6666666666666666 | 0.39440531887330776 |

Theoretical balanced random reference: F1 = 0.5, MCC = 0 (population expectations, not measurements).
See each metrics JSON for source and mutated-partner subgroups.
