#!/usr/bin/python
#
import numpy
import pandas
import sys
import os
import commands
import pandas
#
#
folder_path = sys.argv[1]  # for example, 'DMS_struct_motif'
meme_path = sys.argv[2] # DMS_structMotif_fimo
struct_path = sys.argv[3] 
protein_name=sys.argv[4]
#
try:
    os.popen('mkdir '+meme_path)
except:
    print meme_path+' folder exists'
#
max_width = 25
#
def distribute(meme_txt, order, max_width):
    motif_order = str(order)
    lines = open(meme_txt).readlines()
    n_site = 0
    matrix, name = [], []
    for iline, line in enumerate(lines):
        words = line.split()
        if (len(words)>=4):
            if words[3] == 'maxsites=':
                nmax_site = int(words[4])
            if (words[0] == 'MOTIF') & (words[1] == motif_order):
                e_value, n_site, width = float(words[-1]), int(words[8]), int(words[5])
        if ' '.join(words) == 'Motif ' + str(order) + ' block diagrams' :
            nStart_name, nEnd_name = iline + 4, iline + 4 + n_site
        if ' '.join(words) == 'Motif ' + str(order) + ' position-specific probability matrix' :
            nStart_matrix, nEnd_matrix = iline + 3, iline + 3 + width
        if ' '.join(words) == 'Motif ' + str(order) + ' sites sorted by position p-value':
            nStart_begin, nEnd_begin = iline + 4, iline + 4 + n_site
    name = [x.split()[0] for x in lines[nStart_name:nEnd_name]]
    matrix = [x.split() for x in lines[nStart_matrix:nEnd_matrix]]
    matrix = numpy.asarray(matrix)
    begin = [int(x.split()[1])-1 for x in lines[nStart_begin:nEnd_begin]]
    return e_value, matrix, width, name, begin, n_site
#
bases = ['base_'+str(x) for x in numpy.arange(0, max_width)]
s_matrix, m_matrix, h_matrix, i_matrix = [], [], [], []
proteins, tran_name, motif_begin = [], [], []
def motif_search():
    meme_files = os.listdir(folder_path)
    meme_txt = [x for x in meme_files if x[-3:]=='txt']
    memetxt_path = folder_path + '/' + meme_txt[0]
    for i in range(1, 21):
        (status, output) = commands.getstatusoutput('grep "Motif ' + str(i) + '" ' + memetxt_path)
        if len(output) > 0:
            evalue, matrix, width, name, begin, n_site = distribute(memetxt_path, i, max_width)
            if (n_site>=1) & (evalue<1e-2):
                proteins.append(protein_name + '_' + str(i))
                s_matrix.append(matrix[:, 4])
                m_matrix.append(matrix[:, 3])
                h_matrix.append(matrix[:, 1])
                i_matrix.append(matrix[:, 2])
                tran_name.append(','.join(name))
                begin = list(map(str, begin))
                motif_begin.append(','.join(begin))
            else:
                print protein_name, i, evalue, n_site
motif_search()
s_matrix, m_matrix = numpy.asarray(s_matrix), numpy.asarray(m_matrix)
h_matrix, i_matrix = numpy.asarray(h_matrix), numpy.asarray(i_matrix)
#
s_df = pandas.DataFrame(s_matrix, index=proteins, columns=bases)
m_df = pandas.DataFrame(m_matrix, index=proteins, columns=bases)
h_df = pandas.DataFrame(h_matrix, index=proteins, columns=bases)
i_df = pandas.DataFrame(i_matrix, index=proteins, columns=bases)
#
#
for protein in s_df.index.values:
    s_prob = s_df.loc[protein].values
    m_prob = m_df.loc[protein].values
    h_prob = h_df.loc[protein].values
    i_prob = i_df.loc[protein].values
    matrix = numpy.vstack((s_prob, m_prob, h_prob, i_prob))
    with open(meme_path+'/'+ protein +'.meme', 'w') as output:
        print >> output, 'MEME version 4'
        print >> output
        print >> output, 'ALPHABET'
        print >> output, 'h'
        print >> output, 'i'
        print >> output, 'm'
        print >> output, 's'
        print >> output, 'END ALPHABET'
        print >> output
        print >> output, 'strands: + -'
        print >> output
        print >> output, 'Background letter frequencies'
        print >> output, 'h 0.200 i 0.150 m 0.150 s 0.500'
        print >> output
        print >> output, 'MOTIF ' + protein
        ns = '25'
        print >> output, 'letter-probability matrix: alength= 4 w= '+ns+' nsites= 50 E= 1.0e-2'
        for i in range(0, int(ns)): 
            print >> output, matrix[2,i], matrix[3,i], matrix[1,i], matrix[0,i]
def fimo_structMotif():
	files = [x for x in os.listdir(meme_path) if x[-5:]=='.meme']
	files.sort()
	for f in files:
	    f_path = meme_path + '/' + f
            print 'fimo --thresh 1e-10 -oc ' + f_path[:-5] + ' ' + f_path + ' ' + struct_path + '/' + f.split('_')[0] + '_struct.fasta'
            os.popen('fimo --thresh 1e-10 -oc ' + f_path[:-5] + ' ' + f_path + ' ' + struct_path + '/' + f.split('_')[0] + '_struct.fasta' + ' &')
fimo_structMotif()	
#
#
#
#
