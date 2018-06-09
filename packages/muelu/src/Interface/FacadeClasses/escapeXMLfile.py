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

  filename_input = "ic-test.xml"
  filename_output = "ic-test.out"
  


  # generate XML file for pre_exodus
  if os.path.isfile(filename_output):
    os.remove(filename_output)
  o = open(filename_output,"a")
  for line in open(filename_input):
    line = line.replace("\"", "\\\"")
    line = line.replace("<", "\"&lt;")
    line = line.replace(">", "&gt;\"")
    o.write(line)
  o.close()
     
if __name__ == "__main__":
  sys.exit(main())

