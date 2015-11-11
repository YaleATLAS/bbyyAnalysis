#! /usr/bin/env python
import glob
import re
import subprocess

patterns = { "data15":"eos:h008/data_25ns/", "mc15":"eos:h008/mc_25ns/.*hh_yybb.*", "mc15_bkg":"eos:h008/mc_25ns/MGPy8_y[bjy][bj][bj].*" }

for list_pattern, dataset_pattern in sorted( patterns.items() ) :
  print "Working on",list_pattern
  list_name = glob.glob( "{0}*list".format(list_pattern) )[0]
  if "eos" in dataset_pattern :
    folder = "/".join( dataset_pattern.split(":")[1].split("/")[:-1] )+"/"
    file_pattern = dataset_pattern.split(":")[1].split("/")[-1]
    if len(file_pattern) == 0 : file_pattern = "."
    all_files = subprocess.Popen( [ "/afs/cern.ch/project/eos/installation/0.3.84-aquamarine/bin/eos.select", "ls", "root://eosatlas.cern.ch//eos/atlas/atlascerngroupdisk/phys-higgs/HSG1/MxAOD/{0}".format(folder) ], stdout=subprocess.PIPE).communicate()[0]
    file_names = [ file_name for file_name in all_files.rstrip().split("\n") if re.match(file_pattern,file_name) ]
    dataset_names = [ "root://eosatlas.cern.ch//eos/atlas/atlascerngroupdisk/phys-higgs/HSG1/MxAOD/{0}{1}".format(folder,file_name) for file_name in file_names if "root" in file_name ]
  else :
    dataset_names = glob.glob( "/afs/cern.ch/user/j/jrobinso/work/public/xAOD/mc15_13TeV/mc15*{0}*/*root*".format(dataset_pattern) )
  with open( list_name, "wb" ) as f_list : [ f_list.write( dataset_name+"\n" ) for dataset_name in sorted(dataset_names) ]
