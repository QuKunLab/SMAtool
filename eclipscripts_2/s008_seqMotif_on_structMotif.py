#!/usr/bin/python

import os
import sys

source_dir_path = sys.argv[1] + '/'    # for example, DMS_seq_on_struct
output_dir_path = sys.argv[2] + '/'    # for example, DMS_seqMotif_on_struct
try:
    os.popen('mkdir '+output_dir_path)
except:
    print output_dir_path+' folder exists'
#
source_files = os.listdir(source_dir_path)
seq_files = [x for x in source_files if x[-5:]=='fasta']
#
for seq in seq_files:
    input_path = source_dir_path + '/' + seq
    output_path = output_dir_path + '/' + seq[:-6]
    os.popen('meme ' + input_path + ' -rna -oc ' + output_path +
        ' -time 160000 -maxsize 1000000 -mod zoops -nmotifs 5 -minw 5 -maxw 6 -nostatus &')
