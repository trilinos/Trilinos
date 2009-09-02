import subprocess
import os
import sys

testSuiteStart = os.times()

f = open("LastTest.log")
lines = f.readlines()
f.close();

testlist = {}
ntests = 0

dirname="."
shortDirName=""

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
ntests = 0

f = open("testOutput.txt","w")

def runtest(cmd):
  timeStart = os.times()

  try:
    lcmd = cmd.strip().split()
    p = subprocess.Popen(lcmd, stdout = f, stderr = f, bufsize = 1)
  except:
    f.write("SOME FAILURE\n")
    pass
  else:
    p.wait()
    timeEnd = os.times()
    elapsedTime = timeEnd[4] - timeStart[4]
    resultText = "Elapsed time: "+str(elapsedTime)+" "+dirname+" "+lcmd[1]+" "+lcmd[2]+" "+lcmd[4]
    f.write(resultText+"\n")
    print ntests," ",resultText

for dirname in testlist:
  #
  # This assumes we are running in the zoltan/test directory
  #
  os.chdir(dirname)

  cmds = testlist[dirname]

  print len(cmds)," TESTS IN ",dirname

  for cmd in cmds:
    ntests = ntests + 1
    runtest(cmd)

  os.chdir("..")

testSuiteEnd = os.times()
elapsedTime = testSuiteEnd[4] - testSuiteStart[4]
f.write("\nTEST SUITE TOTAL TIME "+str(elapsedTime))
f.close()
