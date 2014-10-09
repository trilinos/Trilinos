#!/usr/bin/env python

include_files_to_remove = [
  "AddSubdirectories",
  "TribitsAddAdvancedTest",
  "TribitsAddExecutable",
  "TribitsAddExecutableAndTest",
  "TribitsAddOptionAndDefine",
  "TribitsAddTest",
  "TribitsCopyFilesToBinaryDir",
  "TribitsLibraryMacros",
  "TribitsListHelpers",
  "TribitsPackageMacros",
  "TribitsSubPackageMacros",
  "TribitsTplDeclareLibraries"
  ]

import sys
import os

file_in = sys.argv[1]
file_out = sys .argv[2]

if file_in.find("/tribits/") != -1 or file_in.find("/TriBITS/") != -1:
  # Don't replace in a TriBITS system file
  sys.exit(0)

if os.path.islink(file_in):
  print "Ignoring soft link: ", file_in
  sys.exit(0)

file_in_str = open(file_in, 'r').read()
file_out_str = ""

file_in_array = file_in_str.split("\n")
if file_in_array[-1] == "":
  file_in_array = file_in_array[0:-1]

line_num = 0
for line in file_in_array:

  #print "line = '"+line+"'"
  trimed_line = line.strip()

  found_incl_match = False
  for incl_file in include_files_to_remove:
    if not found_incl_match:
      idx_uc_incl = trimed_line.find("INCLUDE("+incl_file)
      if idx_uc_incl == 0:
        found_incl_match = True
    if not found_incl_match:
      idx_lc_incl = trimed_line.find("include("+incl_file)
      if idx_lc_incl == 0:
        found_incl_match = True
    if found_incl_match:
      break

  if found_incl_match:
   print "Remove include line: "+file_in+":"+str(line_num)+": '"+line+"'"
  else:
    file_out_str += line + "\n"

  line_num += 1

open(file_out, 'w').write(file_out_str)
