#!/usr/bin/python
#
import numpy
import os
import pandas
from optparse import OptionParser
#
opts = OptionParser()
usage = "run s4\nusage:%prog -i input folder -f fasta file -o output folder"
opts=OptionParser(usage=usage)
opts.add_option("-i", help="binding_sites and motif_start file location")
opts.add_option("-f", help="structure 's fasta and annotation file")
opts.add_option("-o", help="output files path")
options,arguments=opts.parse_args()
#
input_folder = options.i
BindingSites_file = input_folder + '/binding_sites_on_trans_25.csv'
BeginBase_file = input_folder + '/motif_start_base_25.csv'
info_file = options.f
output_dir = options.o  
width = 25
#
def ExtractFasta():
	motifs_df = pandas.DataFrame.from_csv(BindingSites_file, sep='\t')
	begin_df = pandas.DataFrame.from_csv(BeginBase_file, sep='\t')
	for protein in motifs_df.index.values:
    		binding_trans = motifs_df.ix[protein].values[0]
    		binding_trans = binding_trans.split(',')
    		begin_base = begin_df.ix[protein].values[0]
    		begin_base = begin_base.split(',')
   		info_csv = pandas.DataFrame.from_csv(info_file, sep='\t')
    		seqs = info_csv.ix[binding_trans].values
   		with open(output_dir + "/" + protein + '.fasta', 'w') as output:
        		for itran,tran in enumerate(binding_trans):
            			start, end = int(begin_base[itran]), int(begin_base[itran])+width
            			seq = seqs[itran][0]
            			print >> output, '> ', tran
 			        print >> output, seq[start:end]
def meme_SeqMotif():
	source_files = os.listdir(output_dir)
	seq_files = [x for x in source_files if x[-5:]=='fasta']
	for seq in seq_files:
    	    input_path = output_dir + '/' + seq
    	    output_path = output_dir + '/' + seq[:-6]
            os.popen('meme ' + input_path + ' -rna -oc ' + output_path +
                     ' -time 160000 -maxsize 1000000 -mod zoops -nmotifs 5 -minw 5 -maxw 6 -nostatus ')
if __name__=="__main__":
	ExtractFasta()
	meme_SeqMotif()

