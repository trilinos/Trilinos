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

  bLeadingLines = True
  
  for i in range(0,len( lines ) ):
    line = lines[i]
    if not line.strip():
      continue
    else:     
      print_part_output( i+1 , filenameformat, line )
    
  # move fragment files to this folder
  runCommand("mv ../src/*.fragment .")
    
    
def runXMLfile(xmlfile, executable):
  # run simple easy input xml deck test (10 multigrid levels)
  cmd = "mpirun -np 2 ../src/"+executable+" --xml=../src/xml/"+xmlfile+ ".xml > "+xmlfile+".txt"
  runCommand(cmd)
  if os.path.isfile(xmlfile+".txt") and os.access(xmlfile+".txt", os.R_OK):
    splitfileintofragments(xmlfile+'.txt','========================================================')
    splitfileintofragments(xmlfile+'.txt_3.fragment','--------------------------------------------------------------------------------')
    print xmlfile + bcolors.OKGREEN+" OK"+bcolors.ENDC
    #runCommand("rm s1_easy.txt")
  else:
    print xmlfile + bcolors.FAIL+" Failure"+bcolors.ENDC
    bAllDataPrepared = False  # some data are missing
    runCommand("rm " + xmlfile + ".txt")  
  
################################# MAIN routine
if __name__ == '__main__':

  bAllDataPrepared = True  # set to false if something is missing to generate pdf file

  print bcolors.OKDARKGREEN + "Prepare files... " + bcolors.ENDC
  
  if not os.path.isfile("../src/MueLu_Challenge_XML.exe") or not os.access("../src/MueLu_Challenge_XML.exe", os.R_OK):
    print bcolors.FAIL+"Failure 1: Could not copy executables from trilinos repository"+bcolors.ENDC
    print "You have to run the script in the binary folder. Make sure that the executables in the doc/Tutorial/src folder are built."
    bAllDataPrepared = False  # some data are missing    
    
  if not os.path.isfile("../src/MueLu_tutorial_laplace2d.exe") or not os.access("../src/MueLu_tutorial_laplace2d.exe", os.R_OK):
    print bcolors.FAIL+"Failure 2: Could not copy executables from trilinos repository"+bcolors.ENDC
    bAllDataPrepared = False  # some data are missing    
    
  if not os.path.isfile("../src/MueLu_tutorial_recirc2d.exe") or not os.access("../src/MueLu_tutorial_recirc2d.exe", os.R_OK):
    print bcolors.FAIL+"Failure 3: Could not copy executables from trilinos repository"+bcolors.ENDC
    bAllDataPrepared = False  # some data are missing    

  # here starts the preparation script

  print bcolors.OKDARKGREEN + "Extract version number from git repository... " + bcolors.ENDC
  # create version file
  cmd = "rm version.txt"
  runCommand(cmd)
  cmd = "git log --pretty=format:'%h' -n 1 > version.txt"
  runCommand(cmd)

  print bcolors.OKDARKGREEN + "Split source files in src folder for inclusion in pdf... " + bcolors.ENDC
  # split cpp file for first example
  splitfileintofragments('../src/laplace2d.cpp','// TUTORIALSPLIT ===========================================================')
  
  # split cpp file for ML interface example
  splitfileintofragments('../src/MLParameterList.cpp','// TUTORIALSPLIT ===========================================================')
  
  splitfileintofragments('../src/ScalingTestParamList.cpp','//============================================ SPLIT')      
  splitfileintofragments('../src/ScalingTest.cpp','// USER GUIDE ') 
  
  print bcolors.OKDARKGREEN + "Run test examples to include results in pdf... " + bcolors.ENDC
  # run simple easy input xml deck test
  cmd = "mpirun -np 2 ../src/MueLu_tutorial_laplace2d.exe --nx=50 --ny=50 --xml=../src/xml/s1_easy.xml > s1_easy.txt"
  runCommand(cmd)
  
  if os.path.isfile("s1_easy.txt") and os.access("s1_easy.txt", os.R_OK):
    splitfileintofragments('s1_easy.txt','========================================================')
    splitfileintofragments('s1_easy.txt_3.fragment','--------------------------------------------------------------------------------')
    print "s1_easy.txt" + bcolors.OKGREEN+" OK"+bcolors.ENDC
    #runCommand("rm s1_easy.txt")
  else:
    print "s1_easy.txt" + bcolors.FAIL+" Failure"+bcolors.ENDC
    bAllDataPrepared = False  # some data are missing
    runCommand("rm s1_easy.txt")
    
  runXMLfile("s1_easy_10levels", "MueLu_tutorial_laplace2d.exe")
  runXMLfile("s1_easy_3levels_unsmoothed", "MueLu_tutorial_laplace2d.exe")
  runXMLfile("s1_easy_3levels_smoothed", "MueLu_tutorial_laplace2d.exe")
  
  runXMLfile("s1_easy_jacobi", "MueLu_tutorial_laplace2d.exe")
  runXMLfile("s1_easy_jacobi2", "MueLu_tutorial_laplace2d.exe")
  runXMLfile("s1_easy_exercise", "MueLu_tutorial_laplace2d.exe")
  
  runXMLfile("s2_adv_b", "MueLu_tutorial_laplace2d.exe")
  runXMLfile("s2_adv_c", "MueLu_tutorial_laplace2d.exe")  

  runXMLfile("s3a", "MueLu_tutorial_recirc2d.exe") 
  runXMLfile("s3b", "MueLu_tutorial_recirc2d.exe") 
  runXMLfile("s3b1", "MueLu_tutorial_recirc2d.exe") 
  runXMLfile("s3b2", "MueLu_tutorial_recirc2d.exe") 
  runXMLfile("s3b3", "MueLu_tutorial_recirc2d.exe") 
  
  runXMLfile("s5a", "MueLu_tutorial_laplace2d.exe")
   
  # run rebalancing example
  cmd = "mpirun -np 4 ../src/MueLu_tutorial_laplace2d.exe --nx=300 --ny=300 --xml=../src/xml/s5a.xml > s5a.txt"
  runCommand(cmd)
  if os.path.isfile("s5a.txt") and os.access("s5a.txt", os.R_OK):
    splitfileintofragments('s5a.txt','========================================================')
    splitfileintofragments('s5a.txt_3.fragment','--------------------------------------------------------------------------------')
    print "s5a.txt " +  bcolors.OKGREEN+"OK"+bcolors.ENDC    
  else:
    print "s5a.txt " + bcolors.FAIL+"Failure"+bcolors.ENDC
    bAllDataPrepared = False  # some data are missing
    runCommand("rm s5a.txt")
    
  print bcolors.OKDARKGREEN + "Run LaTeX to generate PDF... " + bcolors.ENDC
  print bcolors.WARNING + "If the script stops here you can skip the step by pressing CTRL+C and run \"pdflatex main.tex\" by hand to fix the errors " + bcolors.ENDC
  # generate pdf for tutorial
  if bAllDataPrepared == True:
    runCommand("pdflatex main.tex")
    if os.path.isfile("main.pdf") and os.access("main.pdf", os.R_OK):
      print bcolors.WARNING+"Success"+bcolors.ENDC 
    else:
      print bcolors.FAIL+"Failure"+bcolors.ENDC
      
  else:
    print bcolors.FAIL+"Cannot generate pdf file due to missing data."+bcolors.ENDC
    
  # clean up
  print bcolors.OKDARKGREEN + "Clean up files... " + bcolors.ENDC
  #runCommand("rm *.fragment *.out *.txt")  
  #runCommand("mv main.pdf muelu_tutorial.pdf")

  print bcolors.OKDARKGREEN + "Finished. " + bcolors.ENDC
