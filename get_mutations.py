#!/usr/bin/env python
import argparse
import os
import pandas as pd
import requests

def getsars_ref():
    url = "https://www.ebi.ac.uk/ena/browser/api/fasta/MN908947.3"
    response = requests.get(url)
    
    if response.status_code == 200:
        sars_ref_content = response.text
        return sars_ref_content
    else:
        print(f"Error in the EBI server: {response.status_code} {response.reason}")
        return None

def main():
    parser = argparse.ArgumentParser(description='Get mutations from SARS-CoV-2 sequences')
    #parser.add_argument('-i', '--input', help='', required=True)
    #parser.add_argument('-o', '--output', help='', required=True)
    args = parser.parse_args()

    ref_seq = "".join(getsars_ref().split('\n')[1:])
if __name__ == '__main__':
    main()