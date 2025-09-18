# Data folder

This folder is reserved for raw and processed input data.

⚠️ Raw data are **NOT included** in the repository due to size and privacy constraints.

## Expected structure
Each sample should be placed in its own subdirectory (e.g., FFA1/, CTRL1/), containing the standard 10x Genomics output:
- barcodes.tsv.gz
- features.tsv.gz
- matrix.mtx.gz

## Example
data/
├── FFA1/filtered_feature_bc_matrix/
│   ├── barcodes.tsv.gz
│   ├── features.tsv.gz
│   └── matrix.mtx.gz
...

## Access
The dataset can be obtained from GEO accession [to be specified] or from your own sequencing runs.
