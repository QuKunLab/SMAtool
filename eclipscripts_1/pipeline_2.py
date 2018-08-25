#!/usr/bin/python
#
import os,sys,time
from optparse import OptionParser
#
opts = OptionParser()
usage = "run eCLIP data pipeline ,version 0.1.0\nusage:%prog -s scripts_folder -n protein's name -a alphabet -e number of extending -w domainMotif_width -i structure annotation -o output_dir "
opts = OptionParser(usage=usage)
opts.add_option("-s", help="scripts_folder path")
opts.add_option("-n", help="protein's name or redefine itself")
opts.add_option("-a", help="alphabet of structure annotation(s,h,i,m)")
opts.add_option("-e", help="number to extend of central position")
opts.add_option("-w", help="domain motif width")
opts.add_option("-o", help="output file path")
opts.add_option("-i", help="input file by running pipeline_1")
options, arguments = opts.parse_args()
#
scripts_folder=options.s
protein_name=options.n
structAnot=options.i
alphabet=options.a
output_dir=options.o
extend=options.e
domainMotif_width=options.w
#
def str2fa_memeMotif():
        os.popen("python " + scripts_folder + "/s1-2_str2fa_memeMotif.py -n " + protein_name + " -i " + structAnot + " -a " + alphabet + " -e " + str(extend) + " -w " + domainMotif_width + " -o " + output_dir)
#
#memetxt = output_dir + "/" + protein_name + "/meme.txt"
AARSdata_folder = output_dir + "/" + protein_name + "_data"
def Motif2DataFrame():
	os.popen("python " + scripts_folder + "./s3_motif_to_dataframe.py -n " + protein_name + " -i " + output_dir + "/" + protein_name + " -o " + AARSdata_folder)
#
info_file = output_dir + "/" + protein_name + "_info.csv"
def SeqMotif():
	os.popen("python " + scripts_folder + "./s4_seq_for_domain_motif.py -i " + AARSdata_folder + " -f " + info_file + " -o " + AARSdata_folder)
#
if __name__=="__main__":
	str2fa_memeMotif()
	#while os.path.exists(memetxt)==False:
        #	time.sleep(30)
        #	if os.path.exists(memetxt):
        #        	break
	Motif2DataFrame()
	SeqMotif()
