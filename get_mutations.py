#!/usr/bin/env python
import argparse
import os
import pandas as pd
import requests
from Bio.Seq import Seq

def getsars_ref():
    url = "https://www.ebi.ac.uk/ena/browser/api/fasta/MN908947.3"
    response = requests.get(url)
    
    if response.status_code == 200:
        sars_ref_content = response.text
        return sars_ref_content
    else:
        print(f"Error in the EBI server: {response.status_code} {response.reason}")
        return None

def get_annotation():
    ranges = {
    'NSP1': (266, 805), 'NSP2': (806, 2719),'NSP3': (2720, 8554),
    'NSP4': (8555, 10054),'NSP5': (10055, 10972),'NSP6': (10973, 11842),
    'NSP7': (11843, 12091),'NSP8': (12092, 12685),'NSP9': (12686, 13024),
    'NSP10': (13025, 13441),'NSP11': (13442, 13468),'NSP12': (13468, 16236),
    'NSP13': (16237, 18039),'NSP14': (18040, 19620),'NSP15': (19621, 20658),
    'NSP16': (20659, 21552),'ORF1ab': (266, 21555),'S': (21563, 25384),
    'ORF3a': (25393, 26220),'E': (26245, 26472),'M': (26523, 27191),
    'ORF6': (27202, 27387),'ORF7a': (27394, 27759),'ORF7b': (27756, 27887),
    'ORF8': (27894, 28259),'N': (28274, 29533),'ORF10': (29558, 29674),
    'ORF1a': (266, 13468),'ORF1b': (13468, 21555)
    }
    return ranges

def read_tsv(file_path):
    df = pd.read_csv(file_path, sep='\t')
    df.columns = df.columns.str.upper()
    return df

def aa_annotation(ranges,ref_seq):
    for key, value in ranges.items():
        start = value[0]-1
        end = value[1]
        aa_seq = ref_seq[start:end]
        print(f"{key}: {aa_seq}")

def main():
    parser = argparse.ArgumentParser(description='Get mutations from SARS-CoV-2 sequences')
    parser.add_argument('-i', '--input', help='', required=True)
    #parser.add_argument('-o', '--output', help='', required=True)
    args = parser.parse_args()
    ref_seq = "".join(getsars_ref().split('\n')[1:])
    check_mut = read_tsv(args.input)

if __name__ == '__main__':
    main()