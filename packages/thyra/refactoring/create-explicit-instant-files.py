#!/usr/bin/env python



import sys
import os
import traceback

scriptsDir = os.path.abspath(os.path.dirname(sys.argv[0]))+"/../../../cmake/python"
sys.path.insert(0, scriptsDir)

from GeneralScriptSupport import *

#
# Read in the command-line arguments
#

usageHelp = r"""create-explicit-instant-files.py [OPTIONS]
"""

from optparse import OptionParser

clp = OptionParser(usage=usageHelp)

clp.add_option(
  "--file-base", dest="fileBase", type="string", default="",
  help="Base file name BaseFileName[_decl,_def].hpp." )

clp.add_option(
  "--class-name", dest="className", type="string", default="",
  help="Class name (defaults to <FILEBASE>)." )

clp.add_option(
  "--macro-name", dest="macroName", type="string", default="",
  help="Macro name." )

clp.add_option(
  "--namespace", dest="namespace", type="string",
  default="Thyra",
  help="Namespace." )

clp.add_option(
  "--explicit-instant-Define", dest="explicitInstantDefine", type="string",
  default="HAVE_THYRA_EXPLICIT_INSTANTIATION",
  help="Name of the explicit instantiation macro." )

clp.add_option(
  "--move-current-files", dest="moveCurrentFiles", action="store_true",
  help="Move the current existing files.",
  default=False )

(options, args) = clp.parse_args()

#
# Check/adjust the input arguments
#

if not options.fileBase:
  raise Exception("Error, must set --file-base")

if not options.className:
  options.className = options.fileBase

#
# Executble code
#

baseFileName = options.namespace+"_"+options.fileBase
declFile = baseFileName+"_decl.hpp"
defFile = baseFileName+"_def.hpp"
9
# A) Move the current files

if options.moveCurrentFiles:
  print "Move the current files ..."
  echoRunSysCmnd("eg mv "+baseFileName+"Decl.hpp "+declFile)
  echoRunSysCmnd("eg mv "+baseFileName+".hpp "+defFile)

# B) Create the XXX.hpp file

Namespace_FileBase_hpp = \
  "#include \""+declFile+"\"\n"+ \
  "#ifndef "+options.explicitInstantDefine+"\n"+ \
  "#  include \""+defFile+"\"\n"+ \
  "#endif\n"

writeStrToFile(baseFileName+".hpp", Namespace_FileBase_hpp)

echoRunSysCmnd("eg add "+baseFileName+".hpp")

# C) Create the XXX.cpp file

Namespace_FileBase_cpp = \
  "#include \""+declFile+"\"\n" + \
  "\n" + \
  "#ifdef "+options.explicitInstantDefine+"\n" + \
  "\n" + \
  "#include \""+defFile+"\"\n" + \
  "#include \"Teuchos_ExplicitInstantiationHelpers.hpp\"\n" + \
  "\n" + \
  "namespace "+options.namespace+" {\n" + \
  "\n"

if options.macroName:
  raise Exception("Do not support --macro-name yet!")
else:
  Namespace_FileBase_cpp += \
    "TEUCHOS_CLASS_TEMPLATE_INSTANT_SCALAR_TYPES("+options.className+")\n"

Namespace_FileBase_cpp += \
  "\n" + \
  "} // namespace "+options.namespace+"\n" + \
  "\n" + \
  "#endif // "+options.explicitInstantDefine+"\n"

writeStrToFile(baseFileName+".cpp", Namespace_FileBase_cpp)

echoRunSysCmnd("eg add "+baseFileName+".cpp")

# D) Print further instructions!

if options.moveCurrentFiles:
  print "TODO: Fix the include guards in "+defFile+"!"
  print "TODO: Fix the include path for the decl file in "+defFile+" and other files!"


