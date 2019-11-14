"""
Author: Jared Galloway

This file contains helper functions
for the script dedup.py
"""

import sys
import io
import re 

def process_record(raw_record, umi, unique_align_buffer):
    """
    this function will 'process' one record,
    by parsing the sam record, checking if a PCR duplicate 
    exists in the unique_alin_buffer, and adding it if not

    Returns None
    """
    record = raw_record.strip().split()
    flag = int(record[1])
    start_position = int(record[3])
    position = int(record[3])
    cigar_string = record[5]

    is_positive = True if ((flag & 16) == 16) else False
    matches = re.findall(r'(\d+)([A-Z]{1})', cigar_string)

    if is_positive:
        if matches[0][1] == 'S':
            position -= int(matches[0][0])
    
    else:
        for match in matches:
            if match[1] != 'I':
                position += int(match[0])

    unique_key = f"{position}_{umi}_{is_positive}"
    if unique_key not in unique_align_buffer:
        unique_align_buffer[unique_key] = [raw_record]
     
    return None


def flush_buffer(output_fp, unique_align_buffer):
    """
    this function takes in a unique align buffer,
    and writes all record to the fp given by output_fp
    """
    
    for key in unique_align_buffer:
        output_fp.write(unique_align_buffer[key][0])
    del unique_align_buffer

    return None
