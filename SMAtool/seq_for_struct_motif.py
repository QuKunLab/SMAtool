#!/usr/bin/python
#
import numpy
import os
import pandas
import sys
#
struct_df_dir = sys.argv[1] + '/'   
binding_sites_on_motif = struct_df_dir + '/binding_sites_on_trans.csv'
begin_base = struct_df_dir + '/motif_start_base.csv'
fasta_dir = sys.argv[2] + '/'
output_dir = sys.argv[3] + '/' 
try:
    os.popen('mkdir '+output_dir)
except:
    print output_dir+' folder exists'
width = int(sys.argv[5])
nmax = int(sys.argv[4])
#
motifs_df = pandas.DataFrame.from_csv(binding_sites_on_motif, sep='\t')
begin_df = pandas.DataFrame.from_csv(begin_base, sep='\t')
for protein in motifs_df.index.values:
    binding_trans = motifs_df.ix[protein].values[0]
    binding_trans = binding_trans.split(',')
    begin_base = begin_df.ix[protein].values[0]
    begin_base = begin_base.split(',')
    info_file = fasta_dir + protein.split('_')[0] + '_info.csv'
    info_csv = pandas.DataFrame.from_csv(info_file, sep='\t')
    seqs_df = info_csv.loc[binding_trans]
#
    if len(binding_trans)>=nmax:
        order = seqs_df['pvalue'].values.argsort()[::-1]
        begin_base = numpy.array(begin_base)[order[:nmax]]
        binding_trans = numpy.array(binding_trans)[order[:nmax]]
        seqs = seqs_df.iloc[order[:nmax]].values
    else:
        seqs = seqs_df.values
    with open(output_dir+protein+'.fasta', 'w') as output:
        for itran,tran in enumerate(binding_trans):
            start, end = int(begin_base[itran]), int(begin_base[itran])+width
            seq = seqs[itran][0]
            print >> output, '> ', tran+'-'+str(start)
            print >> output, seq[start:end]

