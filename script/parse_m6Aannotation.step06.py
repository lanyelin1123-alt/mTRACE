import pandas as pd
import argparse

def prioritize_feature(features):
    priority = ["five_prime_utr", "CDS", "three_prime_utr", "exon", "intron", "TSS_upstream500","start_codon","stop_codon", "Intergenic"]
    features = set(features.split(","))
    
    mRNA_exons = {"five_prime_utr", "CDS", "three_prime_utr"}
    
    if mRNA_exons.issubset(features):
        return "mRNA_exons_mix"
    elif mRNA_exons & features and "exon" in features:
        return "exons_mix"
    else:
        for f in priority:
            if f in features:
                return f
    return "unknown"

def process_bed(input_file, output_file):
    df = pd.read_csv(input_file, sep="\t", header=None)
    df[4] = df[4].apply(prioritize_feature)
    df.to_csv(output_file, sep="\t", index=False, header=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter BED gene features based on priority rules.")
    parser.add_argument("-i", "--input", required=True, help="Input BED file")
    parser.add_argument("-o", "--output", required=True, help="Output BED file")
    args = parser.parse_args()
    
    process_bed(args.input, args.output)

