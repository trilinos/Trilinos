#!/usr/bin/env python

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
 

  
  jobfilename = "runJobs.sh"

  # open jobfile
  if os.path.isfile(jobfilename):
    os.remove(jobfilename)
  jobfile = open(jobfilename,"a")

  jobfile.write("#!/bin/bash\n")

  genFile("crada50x50x30-","muelu","driver_fullAMG.xml", 4, 1000,  1, 1.0 , "true", 1, 1.25, "false", "SGS", 3, 0.8, "ILU", 2, 8, jobfile)                                                                 
  genFile("crada50x50x30-","muelu","driver_fullAMG.xml", 4, 1000,  1, 1.1 , "true", 1, 1.25, "false", "SGS", 3, 0.8, "ILU", 2, 8, jobfile)                                                                 
  genFile("crada50x50x30-","muelu","driver_fullAMG.xml", 4, 1000,  1, 1.15, "true", 1, 1.25, "false", "SGS", 3, 0.8, "ILU", 2, 8, jobfile)                                                                 
  genFile("crada50x50x30-","muelu","driver_fullAMG.xml", 4, 1000,  1, 1.2 , "true", 1, 1.25, "false", "SGS", 3, 0.8, "ILU", 2, 8, jobfile)                                                                 
  genFile("crada50x50x30-","muelu","driver_fullAMG.xml", 4, 1000,  1, 1.25, "true", 1, 1.25, "false", "SGS", 3, 0.8, "ILU", 2, 8, jobfile)
  genFile("crada50x50x30-","muelu","driver_fullAMG.xml", 4, 1000,  1, 1.3 , "true", 1, 1.25, "false", "SGS", 3, 0.8, "ILU", 2, 8, jobfile)                                                                 


  genFile("crada50x50x30-","muelu","driver_fullAMG.xml", 4, 1000,  1, 1.25, "true", 2, 1.25, "false", "SGS", 3, 0.8, "ILU", 2, 8, jobfile)
  genFile("crada50x50x30-","muelu","driver_fullAMG.xml", 4, 1000,  2, 1.25, "true", 2, 1.25, "false", "SGS", 3, 0.8, "ILU", 2, 8, jobfile)


  genFile("crada50x50x30-","muelu","driver_fullAMG.xml", 4, 1000,  1, 1.2 , "true", 1, 1.25, "false", "SGS", 1, 0.8, "ILU", 2, 8, jobfile)                                                                 
  genFile("crada50x50x30-","muelu","driver_fullAMG.xml", 4, 1000,  1, 1.25, "true", 1, 1.25, "false", "SGS", 2, 0.8, "ILU", 2, 8, jobfile)
  genFile("crada50x50x30-","muelu","driver_fullAMG.xml", 4, 1000,  1, 1.3 , "true", 1, 1.25, "false", "SGS", 3, 0.8, "ILU", 2, 8, jobfile)                                                                 
  
  genFile("crada50x50x30-","muelu","driver_fullAMG.xml", 4, 1000,  1, 1.0 , "true", 1, 1.25, "false", "SGS", 3, 0.8, "ILU", 1, 8, jobfile)                                                                 
  genFile("crada50x50x30-","muelu","driver_fullAMG.xml", 4, 1000,  1, 1.1 , "true", 1, 1.25, "false", "SGS", 3, 0.8, "ILU", 1, 8, jobfile)                                                                 
  genFile("crada50x50x30-","muelu","driver_fullAMG.xml", 4, 1000,  1, 1.15, "true", 1, 1.25, "false", "SGS", 3, 0.8, "ILU", 1, 8, jobfile)                                                                 
  genFile("crada50x50x30-","muelu","driver_fullAMG.xml", 4, 1000,  1, 1.2 , "true", 1, 1.25, "false", "SGS", 3, 0.8, "ILU", 1, 8, jobfile)                                                                 
  genFile("crada50x50x30-","muelu","driver_fullAMG.xml", 4, 1000,  1, 1.25, "true", 1, 1.25, "false", "SGS", 3, 0.8, "ILU", 1, 8, jobfile)
  genFile("crada50x50x30-","muelu","driver_fullAMG.xml", 4, 1000,  1, 1.3 , "true", 1, 1.25, "false", "SGS", 3, 0.8, "ILU", 1, 8, jobfile)                                                                 
  
  
  # close the jobfile
  jobfile.close()

     
def genFile(seriesname,simname,tmplfile,maxLevel,maxCoarseSize, simplesweeps, simplerelax, usesimplec, coarsesimplesweeps, coarsesimplerelax, coarseusesimplec, relsmootype, relsmooiter, relsmoodamp, schursmootype, schursmoofill, procs, jobfile):
     
  # generate filenames
  filename_out = seriesname + simname + "_" + str(simplesweeps) + "Simple" + str(simplerelax) + "_" + str(coarsesimplesweeps) + str(usesimplec) + "Simple" + str(coarsesimplerelax) + str(coarseusesimplec) + "_" + str(relsmootype) + str(relsmooiter) + str(relsmoodamp) + "_" + str(schursmootype) + str(schursmoofill) + "_maxL" + str(maxLevel) + "_maxCoarseSize" + str(maxCoarseSize) + "_" + str(procs) + "proc"  

  
  if (relsmootype == "SGS"):
    relsmootype = "Symmetric Gauss-Seidel"
  elif (relsmootype == "GS"):
    relsmootype = "Gauss-Seidel"
  elif (relsmootype == "JAC"):
    relsmootype = "Jacobi"
    
    
  # generate filenames
  filename_log  = filename_out + ".log"     
  filename_xml  = filename_out + ".xml"

 
  # generate XML file for pre_exodus
  if os.path.isfile(filename_xml):
    os.remove(filename_xml)
  o = open(filename_xml,"a")
  for line in open(tmplfile):

    line = line.replace("$MAXLEVEL", str(maxLevel))
    line = line.replace("$MAXCOARSESIZE", str(maxCoarseSize))
    
    line = line.replace("$RELSMOOTYPE", str(relsmootype))
    line = line.replace("$RELSMOOITER", str(relsmooiter))
    line = line.replace("$RELSMOODAMP", str(relsmoodamp))
    
    line = line.replace("$SCHURSMOOTYPE", str(schursmootype))
    line = line.replace("$SCHURSMOOFILL", str(schursmoofill))

    
    line = line.replace("$COARSESIMPLESWEEPS", str(coarsesimplesweeps))
    line = line.replace("$COARSESIMPLERELAX", str(coarsesimplerelax))
    line = line.replace("$COARSEUSESIMPLEC", str(coarseusesimplec))
        
    line = line.replace("$SIMPLESWEEPS", str(simplesweeps))
    line = line.replace("$SIMPLERELAX", str(simplesweeps))
    line = line.replace("$USESIMPLEC", str(usesimplec))
    o.write(line)
  o.close()

  jobfile.write("mpirun -np " + str(procs) + " ./MueLu_CradaDriver.exe --xml=" + filename_xml + " --split=1 --thyra=0  --matrixfile=grid_50_50_30/A0.m --rhsfile=grid_50_50_30/b0.m --coordinatesfile=grid_50_50_30/coordinates.m --specialfile=grid_50_50_30/special.m --nspfile= --linAlgebra=Epetra --tol=1e-10 --maxits=150 --output=" + filename_out + " | tee " + filename_out + ".txt\n\n\n")
   
     
if __name__ == "__main__":
  sys.exit(main())
