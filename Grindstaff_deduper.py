#!/usr/bin/env python

import argparse
import sys
import re

def get_args():
    parser = argparse.ArgumentParser(description='Takes sorted SAM containing uniquely mapped single-end data and returns SAM with PCR duplicates removed.')
    parser.add_argument('-f', '--file', help='Absolute path to sorted SAM', required=True)
    parser.add_argument('-p', '--paired', help='Denotes paired end data; default is single end. This script currently only deals with single-end data so will return error message if -p is specified.', action='store_true', required=False)
    parser.add_argument('-u', '--umi', help='Text file containing list of UMI sequences used to construct library, one sequence per line. Do not set if randomers were used. This script does not currently support randomer usage and will return error message if -u is not specified.', required=False)
    return parser.parse_args()

args = get_args()

#script is not set up to run with paired data or randomers
if args.paired:
    print('Error: Grindstaff_deduper.py cannot handle paired-end data at this time.', file = sys.stderr)
    sys.exit(1)
if args.umi == None:
    print('Error: Grindstaff_deduper.py cannot handle libraries constructed with randomers at this time. If your library was constructed with UMIs, please add a file containing UMI sequences with -u --umi.', file = sys.stderr)
    sys.exit(1)

def strand_checker(flag: int) -> str:
    '''Takes bitwise SAM flag and returns '+' if the read maps to the + strand or - if the read maps to the - strand.'''
    if ((flag & 16) == 16):
        return '-'
    else:
        return '+'

def adj_start(pos: int, cigar: str, strand: str) -> int:
    '''Takes reported start position, CIGAR string, and strand ('+' or '-') and returns actual 5' start position of read.
    Adjusts for soft clipping at 5' of reads on + strand, and accounts for query-consuming sequences for reads on - strand.'''
    if strand == '+': #for + strand, check for soft clipping on leftmost side of CIGAR string and subtract from start pos
        if re.search('^[0-9]+S', cigar):
            return pos - int(re.search('(^[0-9]+)S', cigar).group(0)[:-1])
        else:
            return pos
    else: #for - strand, add up deletions (D), skipped regions from reference (N), matches and mismatches (M), and soft clipping on rightmost side of CIGAR string (S$).
        adj_pos = pos #adjusted position starts at reported position
        mdn = re.findall('([0-9]+)[MDN]', cigar)
        for i in mdn:
            adj_pos += int(i) #add all instances of D, N, and M
        right_s = re.search('[0-9]+S$', cigar) #if there is soft clipping on the rightmost side of CIGAR, add it to total
        if right_s:
            return adj_pos + int(right_s.group(0)[:-1]) - 1
        else:
            return adj_pos - 1


#create UMI set from UMI file
umi_list = []
with open(args.umi, 'r') as fh:
    for line in fh:
        line = line.strip()
        umi_list.append(line)
    umi_set = set(umi_list)

output = args.file + '_deduped'

infile = open(args.file, 'rt')
outfile = open(output, 'wt')

main_dict = {}

#main loop
for line in infile:
    line=line.strip()
    if line[0] == '@': #if line is header, write it to output file
        outfile.write(line + '\n')
    else:
        line = line.split('\t')
        umi = re.search('[ATCGN]+$', line[0]).group(0)
        if umi not in umi_set: #check that umi is in umi_set, and don't write out read if it's not
            continue
        else: #extract information about chrom, strand, and start position
            chrom = line[2]
            strand = strand_checker(int(line[1]))
            real_start = adj_start(int(line[3]), line[5], strand)
            if umi not in main_dict:
                outfile.write('\t'.join(line) + '\n')
                main_dict[umi] = {}
                main_dict[umi][chrom] = {}
                main_dict[umi][chrom][real_start] = [strand]
            elif chrom not in main_dict[umi]:
                outfile.write('\t'.join(line) + '\n')
                main_dict[umi][chrom] = {}
                main_dict[umi][chrom][real_start] = [strand]
            elif real_start not in main_dict[umi][chrom]:
                outfile.write('\t'.join(line) + '\n')
                main_dict[umi][chrom][real_start] = [strand]
            elif strand not in main_dict[umi][chrom][real_start]:
                outfile.write('\t'.join(line) + '\n')
                main_dict[umi][chrom][real_start] += [strand]
            else:
                continue
            
infile.close()
outfile.close()

exit
