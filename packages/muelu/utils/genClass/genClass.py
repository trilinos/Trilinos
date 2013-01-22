#!/usr/bin/python

#import os.path
import os
import sys
import math

def getstatusoutput(cmd):
    """Return (status, output) of executing cmd in a shell."""
    pipe = os.popen(cmd + ' 2>&1', 'r')
    text = pipe.read()
    sts = pipe.close()
    if sts is None: sts = 0
    if text[-1:] == '\n': text = text[:-1]
    return sts, text


def deleteDir(path):
    """deletes the path entirely"""
    cmd = "rm -rf "+path
    result = getstatusoutput(cmd)
    if(result[0]!=0):
        raise RuntimeError(result[1])
      
def createDir(path):
    """deletes the path entirely"""
    cmd = "mkdir "+path
    result = getstatusoutput(cmd)
    if(result[0]!=0):
        raise RuntimeError(result[1])      

def runCommand(cmd):
    """deletes the path entirely"""
    result = getstatusoutput(cmd)
    #if(result[0]!=0):
    #    raise RuntimeError(result[1]) 
    return result[1]
        
###########
# MAIN routine
def main(argv=None):
    
    
  className = "MyTestFactory"
  
  baseClass = "SingleLevelFactoryBase"  # or "TwoLevelFactoryBase"
  
  # the template type can be
  # LO-GO-NO-LMO
  # SC-LO-GO-NO-LMO
  # SC-LO-GO
  # SC-LO
  templateType = "SC-LO-GO-NO-LMO" 

  # process input data
  classNameUpper = className.upper()
  
  if baseClass == "SingleLevelFactoryBase":
    DeclareInputName = "DeclareInput(Level &currentLevel)"
    BuildName        = "Build(Level &currentLevel)"
  elif baseClass == "TwoLevelFactoryBase":
    DeclareInputName = "DeclareInput(Level &fineLevel, Level &coarseLevel)"
    BuildName        = "Build(Level &fineLevel, Level &coarseLevel)"
  
  if templateType == "SC-LO-GO-NO-LMO":
    templateDefinition = "template <class Scalar = double, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps>"
    templateDefShort   = "template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>"
    templateParameters = "<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>"
  elif templateType == "SC-LO-GO":
    templateDefinition = "template <class Scalar = double, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal>"
    templateDefShort   = "template <class Scalar, class LocalOrdinal, class GlobalOrdinal>"
    templateParameters = "<Scalar, LocalOrdinal, GlobalOrdinal>"
  elif templateType == "SC-LO":
    templateDefinition = "template <class Scalar = double, class LocalOrdinal = int>"
    templateDefShort   = "template <class Scalar, class LocalOrdinal>"
    templateParameters = "<Scalar, LocalOrdinal>"
  elif templateType == "LO-GO-NO-LMO":
    templateDefinition = "template <class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps>"
    templateDefShort   = "template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>"
    templateParameters = "<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>"    
  else:
    print "unknown templateType. Must be either LO-GO-NO-LMO, SC-LO-GO-NO-LMO, SC-LO-GO or SC-LO"
    exit(-1)
    
  # define filenames
  decl_filename = "MueLu_" + className + "_decl.hpp"
  def_filename  = "MueLu_" + className + "_def.hpp"
  fwd_filename  = "MueLu_" + className + "_fwd.hpp"
  cpp_filename  = "MueLu_" + className + ".cpp"
  
  # generate declaration file
  if os.path.isfile(decl_filename):
    print "Error: " + decl_filename + " already exists..."
    sys.exit(-1)
    
  o = open(decl_filename,"a")
  for line in open("MueLu_DemoFactory_decl.hpp_tmpl"):
    line = line.replace("$CLASSNAMEUPPER", classNameUpper)
    line = line.replace("$CLASSNAME", className)
    line = line.replace("$BASECLASS", baseClass)
    line = line.replace("$TEMPLATEDEFINITION", templateDefinition)
    line = line.replace("$DECLAREINPUT", DeclareInputName)
    line = line.replace("$BUILD", BuildName)
    o.write(line)
  o.close()

  # generate definition file
  if os.path.isfile(def_filename):
    print "Error: " + def_filename + " already exists..."
    sys.exit(-1)
    
  o = open(def_filename,"a")
  for line in open("MueLu_DemoFactory_def.hpp_tmpl"):
    line = line.replace("$CLASSNAMEUPPER", classNameUpper)
    line = line.replace("$CLASSNAME", className)
    line = line.replace("$BASECLASS", baseClass)
    line = line.replace("$SHORTTEMPLATEDEFINITION", templateDefShort)
    line = line.replace("$TEMPLATEDEFINITION", templateDefinition)
    line = line.replace("$TEMPLATEPARAMETERS", templateParameters)
    line = line.replace("$DECLAREINPUT", DeclareInputName)
    line = line.replace("$BUILD", BuildName)
    o.write(line)
  o.close()	
if __name__ == "__main__":
  sys.exit(main())

