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

Go to a tsv file for the example: 
```bash
wget https://raw.githubusercontent.com/PathoGenOmics-Lab/SARStools/main/in.file
```

## Running the Script

You can run the script from the command line as follows:
```bash
python3 get_mutations.py -i input_file.tsv -o output_file.tsv
```
