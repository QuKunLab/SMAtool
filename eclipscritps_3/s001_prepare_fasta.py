#!/usr/bin/python
#
import os,commands
import numpy
import pandas
import json
import copy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn
import sys
#
nmax_site=int(sys.argv[5])
domainMotif_width=sys.argv[6]
files = sys.argv[1]
out = sys.argv[2]
protein=sys.argv[3]
alphabet=sys.argv[4]
mRNA_ref=sys.argv[7]
try:
    os.popen('mkdir '+out)
except:
    print out+' folder exists'
#
#
def read_file(df_file, rna_list, nmax_site, protein, outpath):
    df = pandas.DataFrame.from_csv(df_file, sep='\t', header=None)
#    select_rna = list(set(df.index.values).intersection(rna_list))
#    print len(select_rna), len(list(set(df.index.values))), len(df.index.values)
#    df = df.loc[select_rna]
    print protein, len(df.index.values)
    df.index = [x+'_'+str(df.ix[i,1]) for i,x in enumerate(df.index.values)]
    name, seq, loop, pvalue = df.index.values, df.ix[:,2].values, df.ix[:,3].values, df.ix[:,5].values
    if len(pvalue) < 50: return
#    if len(pvalue) > nmax_site :
    pvalue_sort_index = numpy.argsort(pvalue)[::-1]
    name, seq = name[pvalue_sort_index[:nmax_site]], seq[pvalue_sort_index[:nmax_site]]
    loop, pvalue = loop[pvalue_sort_index[:nmax_site]], pvalue[pvalue_sort_index[:nmax_site]]
    seq_fasta = outpath + '/' + protein + '_seq.fasta'
    struct_fasta = outpath + '/' + protein + '_struct.fasta'
    with open(seq_fasta, 'w') as seq_f, open(struct_fasta, 'w') as struct_f:
        for i in range(0, len(name)):
            print >> seq_f, '> ' + name[i]
            print >> seq_f, seq[i]
            print >> struct_f, '> ' + name[i]
            print >> struct_f, loop[i]
    matrix = numpy.vstack((seq, loop, pvalue))
    matrix_df = pandas.DataFrame(matrix.T, index=name, columns=['sequence', 'structure', 'pvalue'])
    matrix_df.to_csv(outpath + '/' + protein + '_info.csv', sep='\t')
    return
def meme_DomainMotif(protein,out,out1):
        input_file = out1 + '/' + protein + '_struct.fasta'
        output = out + "/" + protein + '_structMotif'
	DomainMotif_width = str(domainMotif_width)
        os.popen('meme ' + input_file + ' -alph ' + alphabet + ' -oc ' + output +
        ' -time 80000 -maxsize 1000000 -mod zoops -nmotifs 20 -w ' + DomainMotif_width + ' -nostatus ')
#
mRNAs = [x[:-1] for x in open(mRNA_ref).readlines()]
#
input_Site = commands.getstatusoutput("wc -l " + files)[1].split(" ")[0]

if int(nmax_site) >= int(input_Site):
	out1=out + '/' + protein + "_" + "allsites"
else:
	out1=out + '/' + protein + "_" + str(nmax_site) + 'sites'
try:
   os.popen("mkdir " + out1)
except :
   print out1 + " folder exists"
read_file(files, mRNAs, nmax_site, protein, out1)
meme_DomainMotif(protein,out,out1)
#
#
