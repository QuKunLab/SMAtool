#!/usr/bin/python

import os
import commands
import numpy
import pandas
import sys
#
#
fimo_path = sys.argv[1] + '/' 
output_path = sys.argv[2] + '/' 
try:
    os.popen('mkdir '+output_path)
except:
    print output_path+' folder exists'
#
#
folders = [x for x in os.listdir(fimo_path) if x[-5:]!='.meme']
folders.sort()
s_matrix, m_matrix, h_matrix, i_matrix = [], [], [], []
tran_name, motif_begin, proteins = [], [], []
for folder in folders:
    print folder
    fimo = fimo_path + folder + '/' + 'fimo.txt'
    fimo_df = pandas.DataFrame.from_csv(fimo, sep='\t')
    if len(fimo_df.index.values)>=50:
        sites = ','.join(fimo_df['sequence name'].values)
        tran_name.append(sites)
        proteins.append(folder)
        motif_begin.append(','.join(map(str, fimo_df['start'].values)))
        sequences = fimo_df['matched sequence'].values
        s_mat, m_mat, h_mat, i_mat = numpy.zeros(25), numpy.zeros(25), numpy.zeros(25), numpy.zeros(25)
        for seq in sequences:
            s_mat = s_mat + numpy.array([1 if x=='s' else 0 for x in seq]) / float(len(sequences))
            m_mat = m_mat + numpy.array([1 if x=='m' else 0 for x in seq]) / float(len(sequences))
            h_mat = h_mat + numpy.array([1 if x=='h' else 0 for x in seq]) / float(len(sequences))
            i_mat = i_mat + numpy.array([1 if x=='i' else 0 for x in seq]) / float(len(sequences))
        s_matrix.append(s_mat)
        m_matrix.append(m_mat)
        h_matrix.append(h_mat)
        i_matrix.append(i_mat)
s_matrix, m_matrix = numpy.array(s_matrix), numpy.array(m_matrix)
h_matrix, i_matrix = numpy.array(h_matrix), numpy.array(i_matrix)
#
#
bases = numpy.arange(1, 26)
s_df = pandas.DataFrame(s_matrix, index=proteins, columns=bases)
m_df = pandas.DataFrame(m_matrix, index=proteins, columns=bases)
h_df = pandas.DataFrame(h_matrix, index=proteins, columns=bases)
i_df = pandas.DataFrame(i_matrix, index=proteins, columns=bases)
tran_df = pandas.DataFrame(tran_name, index=proteins, columns=['binding_sites_on_transcripts'])
begin_df = pandas.DataFrame(motif_begin, index=proteins, columns=['motif_start_base'])
s_df.to_csv(output_path + 'stem' + '.csv', sep='\t')
m_df.to_csv(output_path + 'multiloop' + '.csv', sep='\t')
h_df.to_csv(output_path + 'hairpin' + '.csv', sep='\t')
i_df.to_csv(output_path + 'interior' + '.csv', sep='\t')
tran_df.to_csv(output_path + 'binding_sites_on_trans' + '.csv', sep='\t')
begin_df.to_csv(output_path + 'motif_start_base' + '.csv', sep='\t')
#
#
