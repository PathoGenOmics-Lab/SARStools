#!/usr/bin/env python
import argparse
import pandas as pd
import requests
from Bio.Seq import Seq

def get_sars_ref():
    url = "https://www.ebi.ac.uk/ena/browser/api/fasta/MN908947.3"
    response = requests.get(url)
    
    if response.status_code == 200:
        sars_ref_content = response.text
        return "".join(sars_ref_content.split('\n')[1:])
    else:
        print(f"Error on the EBI server: {response.status_code} {response.reason}")
        return None

def get_annotation():
    ranges = {
        'NSP1': (266, 805), 'NSP2': (806, 2719), 'NSP3': (2720, 8554),
        'NSP4': (8555, 10054), 'NSP5': (10055, 10972), 'NSP6': (10973, 11842),
        'NSP7': (11843, 12091), 'NSP8': (12092, 12685), 'NSP9': (12686, 13024),
        'NSP10': (13025, 13441), 'NSP11': (13442, 13468), 'NSP12': (13468, 16236),
        'NSP13': (16237, 18039), 'NSP14': (18040, 19620), 'NSP15': (19621, 20658),
        'NSP16': (20659, 21552), 'ORF1ab': (266, 21555), 'S': (21563, 25384),
        'ORF3a': (25393, 26220), 'E': (26245, 26472), 'M': (26523, 27191),
        'ORF6': (27202, 27387), 'ORF7a': (27394, 27759), 'ORF7b': (27756, 27887),
        'ORF8': (27894, 28259), 'N': (28274, 29533), 'ORF10': (29558, 29674),
        'ORF1a': (266, 13468), 'ORF1b': (13468, 21555)
    }
    return ranges

def read_tsv(file_path):
    df = pd.read_csv(file_path, sep='\t')
    df.columns = df.columns.str.upper()  
    return df

def determine_amino_acids(mutations, ref_seq, annotation_ranges):
    effects = []
    for _, row in mutations.iterrows():
        pos_nucleotide = row['POS']
        alt_nucleotide = row['ALT']
        for gene, (start, end) in annotation_ranges.items():
            if pos_nucleotide >= start and pos_nucleotide <= end:
                gene_seq = Seq(ref_seq[start - 1:end])
                aa_index = (pos_nucleotide - start) // 3
                original_codon_start = aa_index * 3
                original_codon = gene_seq[original_codon_start:original_codon_start+3]
                
                aa_ref = original_codon.translate()
                mutated_codon = list(original_codon)
                mutated_codon[(pos_nucleotide - start) % 3] = alt_nucleotide
                mutated_codon = Seq("".join(mutated_codon))
                
                aa_mut = mutated_codon.translate()
                aa_pos = aa_index + 1

                effects.append((pos_nucleotide, gene, aa_ref, aa_mut, aa_pos))
    return effects

def main():
    parser = argparse.ArgumentParser(description='Get mutations from SARS-CoV-2 sequences')
    parser.add_argument('-i', '--input', help='Path to input TSV file', required=True)
    args = parser.parse_args()

    ref_seq = get_sars_ref()
    if ref_seq is None:
        return

    mutation_data = read_tsv(args.input)
    annotation_ranges = get_annotation()

    samples = mutation_data['SAMPLE'].unique()

    for sample in samples:
        sample_mutations = mutation_data[mutation_data['SAMPLE'] == sample]

        print(f"Sample {sample} mutations:")
        effects = determine_amino_acids(sample_mutations, ref_seq, annotation_ranges)
        for pos, gene, aa_ref, aa_mut, aa_pos in effects:
            print(f"Position {pos} in gene {gene}: Reference AA: {aa_ref}, Mutated AA: {aa_mut}, Amino Acid Position: {aa_pos}")

if __name__ == '__main__':
    main()
