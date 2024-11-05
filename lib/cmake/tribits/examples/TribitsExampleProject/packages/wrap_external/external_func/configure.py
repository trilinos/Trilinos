#!/usr/bin/env python3

import sys
import os
import subprocess

from optparse import OptionParser


#
# A) Get commandline options
#

clp = OptionParser()

clp.add_option(
  "--cxx", dest="cxx", type="string",
  default=""
  )

clp.add_option(
  "--cxx-flags", dest="cxxFlags", type="string",
  default=""
  )

clp.add_option(
  "--ar", dest="ar", type="string",
  default=""
  )

clp.add_option(
  "--include-dirs", dest="includeDirs", type="string",
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
"all: libexternal_func.a\n" \
"\n" \
"external_func.o: "+options.srcDir+"/external_func.hpp "+options.srcDir+"/external_func.cpp\n" \
"\t"+options.cxx+" "+options.cxxFlags+"" \
  " -I. "+options.includeDirs+" -o external_func.o -c "+options.srcDir+"/external_func.cpp\n" \
"\n" \
"libexternal_func.a: external_func.o\n" \
"\t"+options.ar+" cr libexternal_func.a external_func.o\n" \
"\n"

open(options.buildDir+"/Makefile", 'w').write(generatedMakefile)
