##
## Read in the LastTest.log file, which is written when the Zoltan tests
## are run under CTest.  Extract the directory name and command for each
## test.
##
## Using this information, re-run each test, timing it.
##
## This is so we can run the same tests in Release 3.1 (pre CTest) and
## the trunk (post CTest) and compare the timings.
##
## We don't check answers.  We are just timing.  We could parse the
## LB_Eval output and say that tests pass if LB_Eval looks good.
##

import time
import subprocess
import os
import sys

testSuiteStart = time.clock()

f = open("LastTest.log")
lines = f.readlines()
f.close();

testlist = {}
ntests = 0

dirname="."

for line in lines:
  line.strip()
  if "Directory" in line:
    tokens = line.split(":")
    dirname = tokens[1].strip()
    tokens = dirname.split("/")
    dirname = tokens[-1]
    
  elif "mpiexec" in line and "-np" in line:
    tokens = line.split(":")
    if len(tokens) > 1:
      tokens = tokens[1].split("|")
      command = tokens[0].strip()
      if dirname in testlist:
        testlist[dirname].append(command+"\n")
      else:
        testlist[dirname] = [command+"\n"]

      ntests = ntests + 1

lines = []

print ntests," TESTS"

f = open("testOutput.txt","w")

def runtest(cmd):
  timeStart = time.clock()

  try:
    lcmd = cmd.strip().split()
    p = subprocess.Popen(lcmd, stdout = f, stderr = f, bufsize = 1)
  except:
    f.write("SOME FAILURE\n")
    pass
  else:
    p.wait()
    elapsedTime = time.clock() - timeStart
    f.write("Elapsed time: "+str(elapsedTime)+" "+dirname+" "+lcmd[1]+" "+lcmd[2]+" "+lcmd[4]+"\n")

for dirname in testlist:
  #
  # This assumes we are running in the zoltan/test directory
  #
  os.chdir(dirname)

  cmds = testlist[dirname]

  print len(cmds)," TESTS IN ",dirname

  for cmd in cmds:
    print "  RUNNING ",cmd
    runtest(cmd)

  os.chdir("..")

testSuiteEnd = time.clock() - testSuiteStart
f.write("\nTEST SUITE TOTAL TIME "+str(testSuiteEnd))
f.close()
