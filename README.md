# CCR data parser for Biothings Hub

This repo stems from the **biothings-data-parser-sample** code, and is used to parse **CCR** data for Biothings Hub.

# Data format

```json
{
    "_id": "chr10:g.1046704C>T",
    "primate_ai": {
        "chrom": "1",
        "start": 865645,
        "end": 865660,
        "ccr_pct": 91.316325156,
        "gene": "AL645608.1",
        "ranges": ["865645-865660","865665-865682"],
        "varflag": [false, false],
        "syn_density": 0.062,
        "cpg": 0.000,
        "cov_score": 30.580,
        "resid": 3.430,
        "resid_pctile": 11.335235448,
        "unique_key": 47070
    }
}
```

# Acknowledgement

All data acquired from Quinlan Lab at [https://github.com/quinlan-lab/ccr](https://github.com/quinlan-lab/ccr).
