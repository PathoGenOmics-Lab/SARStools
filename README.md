# SARStools
A collection of useful scripts for SARS-CoV-2 data 
# SARS-CoV-2 Mutation Analyzer

This script processes mutations in SARS-CoV-2 sequences, determines the effects of these mutations on amino acid sequences, and outputs the results in a tabular format.

## Requirements

- Python 3.x
- pandas
- requests
- Biopython

You can install the required Python packages using pip:

```bash
pip install pandas requests biopython
```

## Input Example

``` tsv
Sample  Pos     Alt
seq1    241     G
seq1    3037    A
seq1    4321    A
seq1    9424    A
seq1    10198   A
seq1    10447   T
seq2    12880   A
seq2    15714   A
seq2    20055   A
seq2    25000   A
seq2    25584   A
seq2    26858   A
seq2    27259   A
```

## Running the Script

You can run the script from the command line as follows:
``` bash
python3 get_mutations.py -i input_file.tsv -o output_file.tsv
```