#!/usr/bin/python
#
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy
import pandas
import sys
import seaborn
from scipy.interpolate import spline
structMotif_DataFrame=sys.argv[1]
seqMotif_DataFrame=sys.argv[2]
plot_file=sys.argv[3]
struct_path=sys.argv[4]
seqMotif_onStruct_folder=sys.argv[5]
output_path=sys.argv[6]
#
input_df = pandas.DataFrame.from_csv(plot_file, sep='\t')
proteins = input_df.index.values
labels = input_df['label'].values
domain_trans_df = pandas.DataFrame.from_csv(structMotif_DataFrame + '/binding_sites_on_trans.csv', sep='\t')
domain_start_df = pandas.DataFrame.from_csv(structMotif_DataFrame + '/motif_start_base.csv', sep='\t')
seq_trans_df = pandas.DataFrame.from_csv(seqMotif_DataFrame + '/binding_sites_on_trans.csv', sep='\t')
seq_start_df = pandas.DataFrame.from_csv(seqMotif_DataFrame + '/motif_start_base.csv', sep='\t')
hairpin_df = pandas.DataFrame.from_csv(structMotif_DataFrame + '/hairpin.csv', sep='\t')
stem_df = pandas.DataFrame.from_csv(structMotif_DataFrame + '/stem.csv', sep='\t')
interior_df = pandas.DataFrame.from_csv(structMotif_DataFrame + '/interior.csv', sep='\t')
multiloop_df = pandas.DataFrame.from_csv(structMotif_DataFrame + '/multiloop.csv', sep='\t')
#
def draw(protein, nprotein, iprotein, label):
    domain_trans = domain_trans_df.loc[protein.split('-')[0]].values[0].split(',')
    domain_start = domain_start_df.loc[protein.split('-')[0]].values[0].split(',')
    domain_start = map(int, domain_start)
    d_start = {x:domain_start[ix] for ix,x in enumerate(domain_trans)}
    seq_trans = seq_trans_df.loc[protein].values[0].split(',')
    seq_start = seq_start_df.loc[protein].values[0].split(',')
    seq_start = map(int, seq_start)
    s_start = {x:seq_start[ix] for ix,x in enumerate(seq_trans)}
#
    info_df = pandas.DataFrame.from_csv(struct_path + "/" + protein.split('_')[0]+'_info.csv', sep='\t')
    hh, ss, ii, mm = numpy.zeros(25), numpy.zeros(25), numpy.zeros(25), numpy.zeros(25)
#    hv, sv, iv, mv = [], [], [], []
    for x in s_start.keys():
        seq_struct = info_df.ix[x, 'structure']
        letters = numpy.asarray([y for y in seq_struct[d_start[x]:d_start[x]+25]])
        h_struct, s_struct, i_struct, m_struct = numpy.zeros(25), numpy.zeros(25), numpy.zeros(25), numpy.zeros(25)
        h_struct[numpy.where(letters=='h')] = 1
        s_struct[numpy.where(letters=='s')] = 1
        i_struct[numpy.where(letters=='i')] = 1
        m_struct[numpy.where(letters=='m')] = 1
        hh[s_start[x]:s_start[x]+6] = hh[s_start[x]:s_start[x]+6] + h_struct[s_start[x]:s_start[x]+6]
        ss[s_start[x]:s_start[x]+6] = ss[s_start[x]:s_start[x]+6] + s_struct[s_start[x]:s_start[x]+6]
        ii[s_start[x]:s_start[x]+6] = ii[s_start[x]:s_start[x]+6] + i_struct[s_start[x]:s_start[x]+6]
        mm[s_start[x]:s_start[x]+6] = mm[s_start[x]:s_start[x]+6] + m_struct[s_start[x]:s_start[x]+6]
    n_avail = len(s_start)
    hh, ss, ii, mm = hh/float(n_avail), ss/float(n_avail), ii/float(n_avail), mm/float(n_avail)
#
    ax1 = plt.subplot2grid((60, 40*nprotein), (0, 40*iprotein+19), colspan=19, rowspan=49)
    xx = numpy.arange(0,25)
    xnew = numpy.linspace(xx.min(), xx.max(), 1000)
    h2, s2, i2, m2 = spline(xx, hh, xnew), spline(xx, ss, xnew), spline(xx, ii, xnew), spline(xx, mm, xnew)
    y1, y2, y3, y4 = m2+i2+s2+h2, m2+i2+s2, m2+i2, m2
    ax1.plot(y1, xnew, color='dimgray', linewidth=3)
#    ax1.plot(y2, xnew, color='lightblue')
#    ax1.plot(y3, xnew, color='palegoldenrod')
#    ax1.plot(y4, xnew, color='palegreen')
    ax1.fill_betweenx(xnew, 0, y1, facecolor='lightcoral')
    ax1.fill_betweenx(xnew, 0, y2, color='lightblue')
    ax1.fill_betweenx(xnew, 0, y3, color='palegoldenrod')
    ax1.fill_betweenx(xnew, 0, y4, color='palegreen')
    all_max = max(y1.max(), y2.max(), y3.max(), y4.max())
    ax1.set_xlim([0, all_max+0.02])
    ax1.set_ylim([24,0])
    plt.setp(ax1, yticks=[], xticks=[], xlabel='Seq-motif', title=label)
#
    hairpin = hairpin_df.ix[protein.split('-')[0]].values
    stem = stem_df.ix[protein.split('-')[0]].values
    interior = interior_df.ix[protein.split('-')[0]].values
    multiloop = multiloop_df.ix[protein.split('-')[0]].values
    bases = numpy.arange(1, 26)
    struct_matrix = numpy.column_stack((bases, multiloop+interior+stem+hairpin, multiloop+interior+stem, multiloop+interior, multiloop))
    struct_df = pandas.DataFrame(struct_matrix, index=xrange(25), columns=['bases','hairpin','stem','interior','multiloop'])
#
    ax2 = plt.subplot2grid((60, 40*nprotein), (0, 40*iprotein), colspan=19, rowspan=49)
    seaborn.barplot(x='hairpin', y='bases', data=struct_df, label='hairpin', color='lightcoral', orient='h', ax=ax2)
    seaborn.barplot(x='stem', y='bases', data=struct_df, label='stem', color='lightblue', orient='h', ax=ax2)
    seaborn.barplot(x='interior', y='bases', data=struct_df, label='interior', color='palegoldenrod', orient='h', ax=ax2)
    seaborn.barplot(x='multiloop', y='bases', data=struct_df, label='multiloop', color='palegreen', orient='h', ax=ax2)
    ax2.set_xlim([0,1])
    if iprotein==0:
        ax2.legend(bbox_to_anchor=(0, 1), ncol=1, fontsize=18)
        plt.setp(ax2, yticks=[], xticks=[], xlabel='Struct-motif', ylabel='BASE')#, title=label)
    else:
        plt.setp(ax2, yticks=[], xticks=[], xlabel='Struct-motif', ylabel='')#, title=label)
#
    image_file = seqMotif_onStruct_folder +'/' + protein.split('-')[0] + '/logo' + protein.split('-')[1] + '.png'
    image = plt.imread(image_file)
    ax3 = plt.subplot2grid((60, 4*nprotein), (52, 4*iprotein+1), colspan=3, rowspan=8)
    ax3.imshow(image, aspect='auto')
    ax3.axis('off')
#
    return n_avail
#
seaborn.set_context('paper', font_scale=1.8)
fig = plt.figure(figsize=(19, 11))
for i,protein in enumerate(proteins):
    print protein
    n_avail = draw(protein, len(proteins), i, labels[i])
    print 'number of available sites:', n_avail
fig.savefig( output_path + "/figure.pdf")
#plt.show()
#
