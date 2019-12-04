#! /usr/bin/env python3
"""
Author: Jared Galloway

This file is an executable python script
which takes in a SAM formatted alignment file 
and writes out another SAM file with all PCR duplicates
removed. 

A PCR duplicate is defined as an aligned read which has:

* the same starting position (leftmost for positive 5` -> 3`)
and rightmost for the opposing orientation mapping

* The same unique molecular identifier (UMI)
"""

# IMPORTS

import argparse
from collections import defaultdict
import sys
import numpy as np
import re

# ARGPARSE

# TODO Clean this up
parser = argparse.ArgumentParser(description='This file is python code for a sam \
    file \
    de-duplication of PCR duplicates. It will run through a SORTED sam file and \
    output another file which only conatins uniquly mapped reads.')
parser.add_argument('-file', type=str,
    help='required arg, file path to the sam file you would like to de-duplicate')
parser.add_argument('-out', type=str, help='the file path to the output file \
    which will be the de-duplicated file will be written to.')
parser.add_argument('-paired', type=str, help='optional arg, designates \
    file is paired end (not single-end)')
parser.add_argument('-umi', default = None, type=str, help='a file containing \
    UMIs used in \
    the reads that are acceptable.')
args = parser.parse_args()

# SCRIPT

if __name__ == "__main__":

    # the data structure which will hold all of
    # the unique reads which I find as a buffer.
    # This buffer will be 'flushed' after each chromosome.
    # {"alignemt starting position" : ["full record", "umi"]}
    unique_chrom_set = set()
    
    # a list of unique UMI's.
    umi_list = []
    
    if args.umi != None:
        umi_fp = open(args.umi,"r")
        for line in umi_fp:
            umi_list.append(line.strip().split()[0])

    # actual same file pointer    
    sam_file_pointer = open(args.file, "r")
    
    # de-dup'd sam file
    de_dup_sam_file_pointer = open(args.out, "w")
    
    # now immediatly write out the header to our output
    header_line = sam_file_pointer.readline()
    while header_line[0] == '@':
        de_dup_sam_file_pointer.write(header_line)
        header_line = sam_file_pointer.readline()

    # now, get to the actual de-duping.
    # the first thing we would like to do is keep track
    # of the chromosome we are one. 
    alignment_record = header_line.strip().split()
    current_chromosome = alignment_record[2]
    for raw_record in sam_file_pointer:
        alignment_record = raw_record.strip().split()

        # need a valid umi
        umi = alignment_record[0].split(":")[-1]
        if umi in umi_list:

            # Parse the CIGAR
            record = raw_record.strip().split()
            flag = int(record[1])

            # if it's not mapped
            if ((flag & 4) == 4):
                continue

            start_position = int(record[3])
            position = int(record[3])
            cigar_string = record[5]

            # sugar, the re is useful for CIGAR parsing into 
            # a list of tuples so 
            # "15S20I" -> [(15,'S'),(20,'I')]
            is_positive = True if ((flag & 16) == 16) else False
            matches = re.findall(r'(\d+)([A-Z]{1})', cigar_string)

            # talked with Thomas Biondi on all these rules for CIGAR
            # if positive, then we just to need to adjust for leading
            # leading soft clipping
            if is_positive:
                if matches[0][1] == 'S':
                    position -= int(matches[0][0])
            
            # for negative strand:
            else:

                # ignore *leading* soft clipping 
                if matches[0][1] == 'S':
                    matches = matches[1:]

                # other wise we would like to simply
                # sum everything that is not an insersion
                # or X/= whatever those are.
                for match in matches:
                    if match[1] not in ['I','X','=']:
                        position += int(match[0])

            # after weve corrected for position and got valid umi
            # create a unique key that can be looked up in O(log n)
            # if this unique key exists inside our 
            unique_key = f"{position}_{umi}_{is_positive}"
            if unique_key not in unique_chrom_set:
                de_dup_sam_file_pointer.write(raw_record)
                unique_chrom_set.add(unique_key) 
        
        # if the alignment record chromosome is different,
        # it's time to flush the Buff!
        if alignment_record[2] != current_chromosome:
            unique_chrom_set = set()
            current_chromosome = alignment_record[2]
    
    # remember to close files for good practice :)
    sam_file_pointer.close()
    de_dup_sam_file_pointer.close()  






















