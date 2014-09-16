#! /usr/bin/env python
import os
import sys
import math
import subprocess

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

def clearWindow():
    os.system('cls' if os.name == 'nt' else 'clear')

def waitForKey():
    os.system('read -s -n 1 -p "Press any key to continue..."')
    print
 
def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

# some colors
class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    OKDARKGREEN = '\033[32m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'

def disable(self):
    self.HEADER = ''
    self.OKBLUE = ''
    self.OKGREEN = ''
    self.OKDARKGREEN = ''
    self.WARNING = ''
    self.FAIL = ''
    self.ENDC = ''

def print_part_output(filenum,filenameformat,line):
    filename = filenameformat.replace('#',str(filenum) )
    fout=open(filename,'w')
    fout.write(line)
    fout.close

def splitfileintofragments(inputfile,splittingtxt):
    
  filenameformat = inputfile + '_#.fragment'  
  file = open( inputfile )
  lines=file.read().split(splittingtxt)

  for i in range(0,len( lines ) ):
    print_part_output( i+1 , filenameformat, lines[i] )
    
  # move fragment files to this folder
  runCommand("mv ../src/*.fragment .")
    
################################# MAIN routine
if __name__ == '__main__':


  bAllDataPrepared = True  # set to false if something is missing to generate pdf file

  # split cpp file
  splitfileintofragments('../src/laplace2d.cpp','// TUTORIALSPLIT ===========================================================')
  
  # run simple easy input xml deck test
  cmd = "../src/MueLu_tutorial_laplace2d.exe --xml=../src/xml/s1_easy.xml > s1_easy.txt"
  runCommand(cmd)
  
  if os.path.isfile("s1_easy.txt") and os.access("s1_easy.txt", os.R_OK):
    print bcolors.WARNING+"Success"+bcolors.ENDC
    splitfileintofragments('s1_easy.txt','========================================================')
    splitfileintofragments('s1_easy.txt_3.fragment','--------------------------------------------------------------------------------')
    #runCommand("rm s1_easy.txt")
  else:
    print bcolors.FAIL+"Failure"+bcolors.ENDC
    bAllDataPrepared = False  # some data are missing
    runCommand("rm s1_easy.txt")
    
  # generate pdf for tutorial
  if bAllDataPrepared == True:
    runCommand("pdflatex main.tex")
    if os.path.isfile("main.pdf") and os.access("main.pdf", os.R_OK):
      print bcolors.WARNING+"Success"+bcolors.ENDC 
    else:
      print bcolors.FAIL+"Failure"+bcolors.ENDC
      
  else:
    print bcolors.FAIL+"Cannot generate pdf file due to missing data."+bcolors.ENDC
