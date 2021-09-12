# CCR data parser for Biothings Hub

This repo stems from the **biothings-data-parser-sample** code, and is used to parse **CCR** data for Biothings Hub.

# Data format

```json
{
    "_id": "chr1:g.865645_865660",
    "ccr": {
        "chrom": "1",
        "start": 865645,
        "end": 865660,
        "scores": [
          {
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
        ]
    }
}
```

Detailed field annotation can be found at [https://github.com/quinlan-lab/ccr#bed-file-columns](https://github.com/quinlan-lab/ccr#bed-file-columns).

# Results

If uploaded correctly, the BioThings studio should report 8,289,870 documents for the current data source (v2.20180420).

# Acknowledgement

All data are downloaded from Quinlan Lab at [https://github.com/quinlan-lab/ccr](https://github.com/quinlan-lab/ccr).
