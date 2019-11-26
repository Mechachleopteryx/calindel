import re
import sys
import argparse
from tqdm import tqdm
import argparse
import json

def reverse_complement(dna):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    return ''.join([complement[base] for base in dna[::-1]])

def read_fasta(fn):
    reads = {}

    with open(fn) as f:
        readName = ""
        for line in f:
            if line[0] == ">":
                readName = line.replace(">", "").replace("@", "").rstrip().split(" ")[0]
                reads[readName] = ""
                continue

            reads[readName] += line.upper().rstrip()

    return reads

def get_num_lines(file_path):
    with open(file_path) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

def main():
    parser = argparse.ArgumentParser(description='Calculate stats from a PAF alignment (must have CIGAR)')

    parser.add_argument('paf', type=str,
                        help='PAF alignment file (required)')
    parser.add_argument('readfile', type=str,
                        help='FASTA file containing the mapped reads (required)')
    parser.add_argument('ref', type=str,
                        help='Reference genome (required)')
    parser.add_argument('--out', type=str,
                        help='Save output to a file')

    args = parser.parse_args()

    bases = ["A", "C", "T", "G"]

    tqdm.write("[INFO] Reading reads file...", file=sys.stderr)
    reads = read_fasta(args.readfile)

    tqdm.write("[INFO] Reading reference file...", file=sys.stderr)
    refs = read_fasta(args.ref)

    insertions = {
        "A": 0,
        "C": 0,
        "T": 0,
        "G": 0
    }
    
    deletions = {
        "A": 0,
        "C": 0,
        "T": 0,
        "G": 0
    }

    substitutions = {}
    
    for a in bases:
        for b in bases:
            if a != b:
                substitutions[a + "-" + b] = 0

    matches = 0
    total_bases = 0

    with open(args.paf) as f:
        pbar = tqdm(f, total=get_num_lines(args.paf))
        pbar.set_description_str("[INFO] Parsing PAF...")

        for line in pbar:
            info = line.split()
            
            readID = info[0]
            strand = info[4]
            refID = info[5]
            queryPos = int(info[2])
            queryEnd = int(info[3])
            targetPos = int(info[7])
            targetEnd = int(info[8])
            cigar = info[-1].split("cg:Z:")[1].split()[0] #cigar = "3M7N4M"

            if refID == "*":
                continue

            read = reads[readID].rstrip()
            readLen = len(read)

            ref = refs[refID].rstrip()

            #### temp fix for duplicated read ids
            if readLen != int(info[1]):
                continue

            read = read[queryPos:queryEnd]
            ref = ref[targetPos:targetEnd]

            if strand == "-":
                read = reverse_complement(read)

            cigarInfo = [(int(v[0]), v[1]) for v in re.findall(r"([0-9]+)([A-Z=]+)", cigar, re.I)]

            readPos = 0
            refPos = 0

            for pair in cigarInfo:
                if pair[1] == "D":
                    for i in range(pair[0]):
                        deletions[ref[refPos]] += 1
                        refPos += 1
                        total_bases += 1
                elif pair[1] == "I":
                    for i in range(pair[0]):
                        insertions[read[readPos]] += 1
                        readPos += 1
                        total_bases += 1
                elif pair[1] == "M":
                    for i in range(pair[0]):
                        if ref[refPos] == read[readPos]:
                            matches += 1
                        else:
                            substitutions[ref[refPos] + "-" + read[readPos]] += 1

                        readPos += 1
                        refPos += 1
                        total_bases += 1

                elif pair[1] == "N":
                    refPos += pair[0]


            # print((info, readID, refID, matches, subs))

    res = {
        "insertions": insertions,
        "deletions": deletions,
        "matches": matches,
        "substitutions": substitutions,
        "total_bases": total_bases
    }
    
    tqdm.write("[INFO] Done!", file=sys.stderr)
    if args.out:
        with open(args.out, 'w') as outfile:  
            tqdm.write("[INFO] Writing output file...", file=sys.stderr)
            json.dump(res, outfile)
    else:
        json.dump(res, sys.stdout)

if __name__ == "__main__":
    main()