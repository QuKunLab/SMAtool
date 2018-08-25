#!/usr/bin/python

import os
import commands
import numpy as np
import pandas
import sys


motif_path = sys.argv[1] + '/'   
output_path = sys.argv[2] + '/'
try:
    os.popen('mkdir '+output_path)
except:
    print output_path+' folder exists'
max_width = 6
#
domain_motif_dir = os.listdir(motif_path)
domain_motif_dir.sort()
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
    matrix = np.asarray(matrix)
    begin = [int(x.split()[1])-1 for x in lines[nStart_begin:nEnd_begin]]
    return e_value, matrix, width, name, begin, n_site
#
bases = ['base_'+str(x) for x in np.arange(0, max_width)]
A_matrix, C_matrix, G_matrix, U_matrix = [], [], [], []
proteins, tran_name, motif_begin, nsite = [], [], [], []
for folder in domain_motif_dir:
    folder_path = motif_path + folder
    meme_files = os.listdir(folder_path)
    meme_txt = [x for x in meme_files if x[-3:]=='txt']
    memetxt_path = folder_path + '/' + meme_txt[0]
    for i in range(1, 6):
        (status, output) = commands.getstatusoutput('grep "Motif ' + str(i) + '" ' + memetxt_path)
        if len(output) > 0:
            evalue, matrix, width, name, begin, n_site = distribute(memetxt_path, i, max_width)
            if matrix.shape[0] < max_width:
                matrix = np.row_stack((matrix, np.zeros(4)))
            if (n_site>=10) | (evalue<0.05):
                proteins.append(folder + '-' + str(i))
                A_matrix.append(matrix[:, 0])
                C_matrix.append(matrix[:, 1])
                G_matrix.append(matrix[:, 2])
                U_matrix.append(matrix[:, 3])
                tran_name.append(','.join(name))
                begin = list(map(str, begin))
                motif_begin.append(','.join(begin))
                nsite.append(n_site)
            else:
                print folder, i, evalue, n_site, width
A_matrix, C_matrix = np.asarray(A_matrix), np.asarray(C_matrix)
G_matrix, U_matrix = np.asarray(G_matrix), np.asarray(U_matrix)
print len(proteins), len(bases), A_matrix.shape
#
A_df = pandas.DataFrame(A_matrix, index=proteins, columns=bases)
C_df = pandas.DataFrame(C_matrix, index=proteins, columns=bases)
G_df = pandas.DataFrame(G_matrix, index=proteins, columns=bases)
U_df = pandas.DataFrame(U_matrix, index=proteins, columns=bases)
tran_df = pandas.DataFrame(tran_name, index=proteins, columns=['binding_sites_on_transcripts'])
begin_df = pandas.DataFrame(motif_begin, index=proteins, columns=['motif_start_base'])
nsite_df = pandas.DataFrame(nsite, index=proteins, columns=['n_site'])
A_df.to_csv(output_path + 'A.csv', sep='\t')
C_df.to_csv(output_path + 'C.csv', sep='\t')
G_df.to_csv(output_path + 'G.csv', sep='\t')
U_df.to_csv(output_path + 'U.csv', sep='\t')
tran_df.to_csv(output_path + 'binding_sites_on_trans.csv', sep='\t')
begin_df.to_csv(output_path + 'motif_start_base.csv', sep='\t')
nsite_df.to_csv(output_path + 'Nsite.csv', sep='\t')
#
#
