#!/usr/bin/env python
#
import os,sys
from optparse import OptionParser
#
opts = OptionParser()
usage = "run eCLIP data pipeline ,version 0.1.0\nusage:%prog -s scripts_dir -m bamfile_1 -n bamfile_2 -c mockfile -e number of extending --pv p_value --rt number of reads --fc fold change" \
	+ "-o output_dir -a annotation file -p bamfile_pars "
opts = OptionParser(usage=usage)
opts.add_option("-s", help="scripts_folder path")
opts.add_option("-m", help="one input .bam file")
opts.add_option("-n", help="another input .bam file")
opts.add_option("-c", help="control mock .bamfile")
opts.add_option("--pv", help="threshold of p_value")
opts.add_option("--rt", help="reads threshold of experiments")
opts.add_option("--fc", help="threshold of experiment's fold change which is computed with control")
opts.add_option("-e", help="number to extend of central position")
opts.add_option("-o", help="output file path")
opts.add_option("-a", help="genome annotation file")
opts.add_option("-p", help="pars bamfile data as ref")
options, arguments = opts.parse_args()
#
halfLength=2
minimum=500
annotation=options.a
scripts_folder = options.s
output_folder=options.o
input_bam1=options.m
input_bam2=options.n
input_mock=options.c
input_bam1_name=os.path.basename(input_bam1)
input_bam2_name=os.path.basename(input_bam2)
#input_bam1_rt=input_bam1+".rt"
input_bam1_rt=output_folder+"/"+input_bam1_name+".rt"
#input_bam2_rt=input_bam2+".rt"
input_bam2_rt=output_folder+"/"+input_bam2_name+".rt"
mocktab=output_folder +"/mock.rt"
#
#
def rtcount(input_bam,input_bam_rt,minimumNumber):
	 os.popen("samtools view "+ input_bam +" > "+ input_bam +".sam")
	 os.popen("python " + scripts_folder +"/rtCounts.py -l " + str(halfLength) + " -a " + annotation +" -m "+ str(minimumNumber) +" -i "+ input_bam +".sam"+" -o "+input_bam_rt )
	 os.popen("rm "+ input_bam +".sam")
#
#input_merged=input_bam1+".merged"
input_merged=output_folder+"/bam.merged"
def merge():
	os.popen("python " + scripts_folder +"/merge.py --t1 " + input_bam1_rt +" --t2 "+ input_bam2_rt + " -o " + input_merged )
#
ws = 5
thre = 5
sampleTime = 1000
#input_peak=input_bam1+".peak"
input_peak=output_folder+"/bam.peak"
def callPeak():
	os.popen("python " + scripts_folder + "/peak.py --ws " + str(ws) + " -n " + str(thre) + " -s " + str(sampleTime) + " -i " + input_merged + " -o " + input_peak )
#
peaktxt=output_folder + "/peak.txt"
def getRelEnrich():
	os.popen("python " + scripts_folder + "/mockenrichment.py -p " + input_peak + " -i " + input_merged + " -m " + mocktab + " -o " + peaktxt )
#
p_value=options.pv
reads_threshold=options.rt
foldchange=options.fc
filter_peaktxt=output_folder + "/peak_filtered.txt"
def peakfilter():
	os.popen("python " + scripts_folder + "/peakfilter.py -p " + str(p_value) + " -r " + str(reads_threshold) + " -f " + str(foldchange) + " -i " + peaktxt + " -o " + filter_peaktxt)
#
extend=options.e
structure = options.p
#structOutput1 = input_bam1+".struct"  
#structOutput2 = input_bam1+".structAnot"
structOutput1 = output_folder+"/bam.struct"
structOutput2 = output_folder+"/bam.structAnot"
def getStru():
	os.popen("python " + scripts_folder + "/getstructure.py -e " + str(extend) + " -i " + filter_peaktxt + " -s " + structure + " -o " + structOutput1 + " -p " + structOutput2)
#
if __name__ == "__main__": 
	rtcount(input_bam1,input_bam1_rt,str(minimum))
	rtcount(input_bam2,input_bam2_rt,str(minimum))
	rtcount(input_mock,mocktab,"0")
	print "RTCounts done\n"
	merge()
	print "merge done\n"
	callPeak()
	print "callpeak done\n"
	getRelEnrich()
	print "getRelEnrich done\n"
	peakfilter()
	print "peak filtered done\n"
	getStru()
	print "get structure done\n "
