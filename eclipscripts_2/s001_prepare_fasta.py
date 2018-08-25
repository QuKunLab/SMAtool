#!/usr/bin/python
#
import os
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
nmax_site = 9000000
files = sys.argv[1]
out = sys.argv[2]
protein=sys.argv[3]
alphabet=sys.argv[4]#"/home/caipengfei/my_file/eCLIP_analysis/ref/dom_structure.alphabet"
try:
    os.popen('mkdir '+out)
except:
    print out+' folder exists'
#DMS_bam = '/home/caipengfei/my_file/eCLIP_analysis2/bamfiles/'
#PARS_bam = '/home/jinyonghao/project/eclipData/bamfiles/'
#DMS_out, PARS_out = './DMS_K562_45bp/', './PARS_K562_45bp/'
#
#
#def get_structAnot(bam):
#    dirs = [x for x in os.listdir(bam) if x[:4]=='K562']
#    dirs.sort()
#    proteins = [x[4:-6] for x in dirs]
#    files = []
#    for folder in dirs:
#        files.extend([bam+folder+'/'+x for x in os.listdir(bam+folder) if x[-10:]=='structAnot'])
#    return files, proteins
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
def meme_DomainMotif(protein,outpath):
        input_file = outpath + '/' + protein + '_struct.fasta'
        output = outpath + "/" + protein + '_structMotif'
	DomainMotif_width = str(25)
        os.popen('meme ' + input_file + ' -alph ' + alphabet + ' -oc ' + output +
        ' -time 80000 -maxsize 1000000 -mod zoops -nmotifs 20 -w ' + DomainMotif_width + ' -nostatus ')
#
mRNAs = [x[:-1] for x in open('grch38.mRNA.txt').readlines()]
#
read_file(files, mRNAs, nmax_site, protein, out)
meme_DomainMotif(protein,out)
#
#
