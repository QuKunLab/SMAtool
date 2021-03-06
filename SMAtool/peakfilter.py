#!/usr/bin/python
#
import os,sys
from optparse import OptionParser
#
opts = OptionParser()
usage = "run peakfilter\nusage:%prog -p p-value -r numbers of reads -f foldchange -i input bam.peak -o output filtered bam.peak"
opts = OptionParser(usage=usage)
opts.add_option("-p", help="threshold of p-value for experiments ")
opts.add_option("-r", help="threshold of numbers of experimental normalized reads ")
opts.add_option("-f", help="threshold of foldchange value which experiments was divided by control ")
opts.add_option("-i", help="input peak.txt generated by last step ")
opts.add_option("-o", help="output peak file filtered ")
options, arguments = opts.parse_args()
#
p_value=options.p
threshold_reads=options.r
foldchange=options.f
input_peak=options.i
output_peak=options.o
#
def isNum(value):
    try:
        float(value)
    except ValueError:
        return False
    else:
        return True
def valid_line(line):
	line=line.rstrip().split("\t")
	type_list=[isNum(i) for i in line[1:-1]]
	if False in type_list:
		return False
        if len(line) < 6:
                return False
        if float(p_value) < float(line[3]):
                return False
        if float(threshold_reads) > float(line[4]):
                return False
	try:
		if float(foldchange) > float(line[5]):
       	        	return False
	except:
		return True
def st_en(st,en):
	ss=[i2+";"+en[i1] for i1,i2 in enumerate(st)]
	return ' '.join(ss)
def packPval(pval):
	return ' '.join(map(str, pval))
def packTranscriptPeaks(name, st, en, pval, times):
	return '\t'.join([name, st_en(st,en), packPval(pval), packPval(times)])
name_list=[]
peak_data=open(input_peak)
out=open(output_peak, 'w')
for iline,line in enumerate(peak_data):
	if valid_line(line)==False:
		continue
	line=line.rstrip().split("\t")
	name=line[0]
	if name not in name_list:
		if len(name_list)==0:
			pass
		if len(name_list) > 0:
			print >> out, packTranscriptPeaks(name_list[-1], st, en, pval, times)	
		name_list.append(name)
		st,en,pval,times=[],[],[],[]
		en.append(line[2])
		st.append(line[1])
		pval.append(line[3])
		times.append(line[4])
	else:
		en.append(line[2])
                st.append(line[1])
                pval.append(line[3])
                times.append(line[4])
print >> out, packTranscriptPeaks(name, st, en, pval, times)		
