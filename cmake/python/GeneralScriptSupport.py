"""
Python module containing general support functions for creating scripts
"""

#
# Check the python version
#

import sys
if sys.version < '2.3':
  print "Error, Python version is "+sys.version+" < 2.3!"
  sys.exit(1)

#
# Import commands
#

import sys
import os
import re
import time


#
# Determine what system we are on:
#

rePlatformName = re.compile(r"^[a-zA-Z]+")
platformStr = rePlatformName.findall(sys.platform)[0]
#print "\nplatformStr =", platformStr


#
# Functions
#


def getScriptsDir():
  """Get the base directory where a file exists"""
  return os.path.abspath(os.path.dirname(sys.argv[0]))


def createDir( dirName ):
  """Create a directory if it does not exist"""
  if os.path.exists(dirName):
    if not os.path.isdir(dirName):
      print "\nError the path", dirName, "already exists but it is not a directory!"
      raise RuntimeError("Directory path error!")
    print "\nThe directory", dirName, "already exists!"
  else:
    print "\nCreating directory "+dirName+" ..."
    os.mkdir(dirName)


def createDirsFromPath(path):
  #print "\npath =", path
  pathList = path.split("/")
  #print "\npathList =", pathList
  if not pathList[0]:
    currDir = "/"
  for dir in pathList:
    currDir = os.path.join(currDir, dir)
    if currDir and not os.path.exists(currDir):
      #print "\ncurrDir =", currDir
      createDir(currDir)


def echoChDir(dirName):
  print "\nChanging current directory to \'"+dirName+"\'"
  os.chdir(dirName)
  print "\nCurrent directory is \'"+os.getcwd()+"\'\n"
  

def runSysCmnd(cmnd, throwExcept=True):
  """Run system command and optionally throw on failure"""
  sys.stdout.flush()
  sys.stderr.flush()
  try:
    errCode = os.system(cmnd)
  except OSError, e:
    errCode = 1 # Just some error code != 0 please!
  if errCode != 0 and throwExcept:
    raise RuntimeError('Error, the command \'%s\' failed with error code %d' \
                       % (cmnd,errCode) )
  return errCode


def echoRunSysCmnd(cmnd, throwExcept=True, msg=None, timeCmnd=False,
  verbose=True
  ):
  """Echo command to be run and run command with runSysCmnd()"""
  if verbose:
    print "\nRunning: "+cmnd+"\n"
  if msg and verbose:
    print "  "+msg+"\n"
  t1 = time.time()
  try:
    rtn = runSysCmnd(cmnd, throwExcept)
  finally:
    if timeCmnd:
      t2 = time.time()
      if verbose:
        print "\n  Runtime for command = %f minutes" % ((t2-t1)/60.0)
  return rtn


def getCmndOutput(cmnd, stripTrailingSpaces=False, throwOnError=True):
  """Run a shell command and return its output"""
  child = os.popen(cmnd)
  data = child.read()
  err = child.close()
  if err and throwOnError:
    raise RuntimeError, '%s failed w/ exit code %d' % (cmnd, err)
  if stripTrailingSpaces:
    return data.rstrip()
  return data
