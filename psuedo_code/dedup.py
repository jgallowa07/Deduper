#! /usr/bin/env python3
"""
Author: Jared Galloway

Had helpful conversations with Thomas Biondi, Seve Villarruel, Nick Wagner, Stacey Coonrod, Pete Batzel

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
from helpers import *

# ARGPARSE

# TODO Clean this up
parser = argparse.ArgumentParser(description='This file is python code for a sam file \
    de-duplication of PCR duplicates. It will run through a SORTED sam file and \
    output another file which only conatins uniquly mapped reads.')
parser.add_argument('-file', type=str,
    help='required arg, file path to the sam file you would like to de-duplicate')
parser.add_argument('-out', type=str, help='the file path to the output file \
    which will be the de-duplicated file will be written to.')
parser.add_argument('-paired', type=str, help='optional arg, designates \
    file is paired end (not single-end)')
parser.add_argument('-umi', default = None, type=str, help='a file containing UMIs used in \
    the reads that are acceptable.')
args = parser.parse_args()

# SCRIPT
"""
if __name__ == "__main__":

    # the data structure which will hold all of
    # the unique reads which I find as a buffer.
    # This buffer will be 'flushed' after each chromosome.
    # {"alignemt starting position" : ["full record", "umi"]}
    unique_chrom_align = {}
    
    # a list of unique UMI's.
    umi_list = []
    
    if args.umi != None:
        # populate umi list!

    # actual same file pointer    
    sam_file_pointer = open(args.file, "r")
    
    # de-dup'd sam file
    de_dup_sam_file_pointer = open(args.out, "w")
    
    # now immediatly write out the header to our output
    header_line = sam_file_pointer.readline()
    while header_line[0] == '@':
        # write your stuff out immediately

    # now, get to the actual de-duping.
    # the first thing we would like to do is keep track
    # of the chromosome we are one. 
    
    current_chromosome = get starting chrom
    for raw_record in sam_file_pointer:
        # split the alignment recrd to a more useful list
        alignment_record = raw_record.strip().split()

        if umi in umi_list:
            # the helper function will process and add to dictionary 
            # if necessary
            #process_record(raw_record, umi, unique_chrom_align)
        
        # if the alignment record chromosome is different,
        # it's time to flush the Buff!
        # theres a better way to do this using a set rather than 
        # flughing the buffer. So I'm, going to try both ways and compare.
    
            # this function will write out the contents of our 
            # buffer to the output file

    # flush the buff one last time for the last chromosome.
    
    # remember to close files for good practice :)
"""






















