#!/bin/python

import sys
import os
import subprocess
import commands


from optparse import OptionParser


#
# A) Get commandline options
#

clp = OptionParser()
  
clp.add_option(
  "--with-export-makefile", dest="exportMakefile", type="string",
  default=""
  )
  
clp.add_option(
  "--src-dir", dest="srcDir", type="string",
  default=""
  )
  
clp.add_option(
  "--build-dir", dest="buildDir", type="string",
  default=""
  )
  
(options, args) = clp.parse_args(sys.argv)

#
# B) Generate the makefile in the build tree!
#

generatedMakefile = \
"include "+options.exportMakefile+"\n" \
"\n" \
"all: libexternal_func.a\n" \
"\n" \
"external_func.o: "+options.srcDir+"/external_func.hpp "+options.srcDir+"/external_func.cpp\n" \
"\t$(TribitsExProj_CXX_COMPILER) $(TribitsExProj_CXX_COMPILER_FLAGS) -I. $(TribitsExProj_INCLUDE_DIRS) -o external_func.o -c "+options.srcDir+"/external_func.cpp\n" \
"\n" \
"libexternal_func.a: external_func.o\n" \
"\t$(TribitsExProj_AR) cr libexternal_func.a external_func.o\n" \
"\n"

open(options.buildDir+"/Makefile", 'w').write(generatedMakefile)


