# @HEADER
# ************************************************************************
#
#            TriBITS: Tribal Build, Integrate, and Test System
#                    Copyright 2013 Sandia Corporation
#
# Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
# the U.S. Government retains certain rights in this software.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
# 1. Redistributions of source code must retain the above copyright
# notice, this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright
# notice, this list of conditions and the following disclaimer in the
# documentation and/or other materials provided with the distribution.
#
# 3. Neither the name of the Corporation nor the names of the
# contributors may be used to endorse or promote products derived from
# this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
# EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
# PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# ************************************************************************
# @HEADER

"""
Python module containing general support functions for creating scripts
"""

#
# Check the python version
#

import sys
if sys.version < '2.6':
   print("Error, Python version is " + sys.version + " < 2.6!")
   sys.exit(1)

#
# Import commands
#
import os
import re
import math
import subprocess
import time
import datetime
import optparse
import traceback

#
# Byte array / string / unicode support across Python 2 & 3
#
# Note that the str class in Python 2 is an ASCII string (byte) array and in
# Python 3 it is a Unicode object. For Python 3 code that is backward compatible
# with Python 2, we sometimes need version-specific conversion functions to give
# us the data type we desire. These functions are:
#
#     b(x)    return a byte array of str x, much like b'<string const>' in
#             Python 3
#     s(x)    return a version-specific str object equivalent to x
#     u(x)    return a unicode object equivalent to x, much like
#             u'<string const>' in Python 2
#
if sys.version_info < (3,):
   # Python 2
   def b(x): return x
   def s(x): return x
   def u(x): return unicode(x)
else:
   # Python 3
   import codecs
   def b(x): return codecs.latin_1_encode(x)[0]
   def s(x):
      try:
         return x.decode("utf-8")
      except AttributeError:
         return x
   def u(x): return x

verboseDebug = False


#
# Determine what system we are on:
#

rePlatformName = re.compile(r"^[a-zA-Z]+")
platformStr = rePlatformName.findall(sys.platform)[0]
#print("\nplatformStr = " + platformStr)


######################################
# Script location functions
######################################


def getScriptBaseDir():
  return os.path.dirname(os.path.realpath(os.path.abspath(sys.argv[0])))


def getScriptName():
  return os.path.basename(os.path.dirname(sys.argv[0]))

def getCompleteFileDirname(filename):
  return os.path.dirname(os.path.realpath(os.path.abspath(filename)))

######################################
# List helper functions
######################################


def findInSequence(seq, item):
  for i in range(0, len(seq)):
    if seq[i] == item:
      return i
  return -1


def removeItemsFromList(list, items):
  numItemsRemoved = 0
  for item in items:
    if item in list:
      idx = list.index(item)
      del list[idx]
      numItemsRemoved = numItemsRemoved + 1
  return numItemsRemoved


######################################
# String helper functions
######################################


def stripWhitespaceFromStringList(stringList):
  return [x.strip() for x in stringList]


def isSubstrInMultiLineString(inputStr, findStr):
  return inputStr.find(findStr) >= 0


def getStrUnderlineStr(width):
  underlineStr = ""
  for i in range(width):
    underlineStr += "-"
  return underlineStr


def arrayToFormattedString(array_in, offsetStr = ""):
  sout = ""
  sout += offsetStr + "[\n"
  for i in range(0, len(array_in)):
    if i != len(array_in)-1:
      commaChar = ","
    else:
      commaChar = ""
    sout += (offsetStr + "  \'" + str(array_in[i]) + "\'"+commaChar+"\n")
  sout += offsetStr + "]\n"
  return sout


def extractLinesAfterRegex(string_in, regex_in):
  #print("regex_in = " + regex_in)
  reMatch = re.compile(regex_in)
  linesExtracted = ""
  foundRegex = False
  for line in string_in.strip().splitlines():
    #print("line = '" + line + "'")
    if not foundRegex:
      matchObj = reMatch.match(line)
      #print("matchObj = " + matchObj)
      if matchObj:
        foundRegex = True
    if foundRegex:
      linesExtracted += line + "\n"
  return linesExtracted


def extractLinesMatchingRegex(string_in, regex_in):
  #print("regex_in = " + regex_in)
  string_in = s(string_in)
  reMatch = re.compile(regex_in)
  linesExtracted = ""
  for line in string_in.strip().splitlines():
    #print("line = '" + line + "'")
    matchObj = reMatch.match(line)
    #print("matchObj = " + matchObj)
    if matchObj:
      linesExtracted += line + "\n"
  return linesExtracted
# NOTE: Above is *NOT* using tested!


def extractLinesMatchingSubstr(string_in, substr_in):
  #print("substr_in = '" + substr_in + "'")
  string_in = s(string_in)
  linesExtracted = ""
  for line in string_in.strip().splitlines():
    #print("line = '" + line + "'")
    if substr_in in line:
      #print("matched '" + substr_in + "'")
      linesExtracted += line + "\n"
  return linesExtracted
# NOTE: Above is *NOT* unit tested!


# Convert a dictionary to a string, using a sorted set of keys.
#
# This is needed to provide a portable string representation across various
# versions of Python and platforms (see TriBITS GitHub Issue #119).
def sorted_dict_str(d):
  items = []
  keys = list(d.keys())
  keys.sort()
  for key in keys:
    if isinstance(d[key], dict):
      value = sorted_dict_str(d[key])
    else:
      value = repr(d[key])
    items.append(repr(key) + ": " + value)
  return "{" + ", ".join(items) + "}"



##############################################
# System command unit testing utiltities
##############################################


class InterceptedCmndStruct:

  def __init__(self, cmndRegex, cmndReturn, cmndOutput):
    self.cmndRegex = cmndRegex
    self.cmndReturn = cmndReturn
    self.cmndOutput = cmndOutput

  def __str__(self):
    return "{cmndRegex='"+self.cmndRegex+"'," \
     " cmndReturn="+str(self.cmndReturn)+"," \
     " cmndOutput='"+str(self.cmndOutput)+"'}"

#
# Class that is used to record a set of commands that will be used to
# intercept commands
#

class SysCmndInterceptor:

  def __init__(self):
    self.__fallThroughCmndRegexList = []
    self.__interceptedCmndStructList = []
    self.__allowExtraCmnds = True

  def setFallThroughCmndRegex(self, cmndRegex):
    self.__fallThroughCmndRegexList.append(cmndRegex)

  def setInterceptedCmnd(self, cmndRegex, cmndReturn, cmndOutput=None):
    self.__interceptedCmndStructList.append(
       InterceptedCmndStruct(cmndRegex, cmndReturn, cmndOutput) )

  def setAllowExtraCmnds(self, allowExtraCmnds):
    self.__allowExtraCmnds = allowExtraCmnds

  def hasInterceptedCmnds(self):
     return len(self.__interceptedCmndStructList) > 0

  def getFallThroughCmndRegexList(self):
    return self.__fallThroughCmndRegexList[:]

  def getInterceptedCmndStructList(self):
    return self.__interceptedCmndStructList[:]

  def doProcessInterceptedCmnd(self, cmnd):
    #print("doProcessInterceptedCmnd(): cmnd='" + cmnd + "'")
    if self.isFallThroughCmnd(cmnd):
      return False
    if len(self.__interceptedCmndStructList) > 0:
      return True
    if not self.__allowExtraCmnds:
      return True
    return False

  def isFallThroughCmnd(self, cmnd):
     for interceptCmndStruct in self.__interceptedCmndStructList:
       if re.match(interceptCmndStruct.cmndRegex, cmnd):
         return False
     for cmndRegex in self.__fallThroughCmndRegexList:
       if re.match(cmndRegex, cmnd):
         return True
     return False

  def nextInterceptedCmndStruct(self, cmnd):
    assert(not self.isFallThroughCmnd(cmnd))
    if len(self.__interceptedCmndStructList) == 0:
      raise Exception("Error, cmnd='"+cmnd+"' is past the last expected command!")
    ics = self.__interceptedCmndStructList[0]
    if not re.match(ics.cmndRegex, cmnd):
      raise Exception("Error, cmnd='" + cmnd + "' did not match the" \
                      " expected regex='" + ics.cmndRegex + "'!")
    self.__interceptedCmndStructList.pop(0)
    return (ics.cmndReturn, ics.cmndOutput)

  def clear(self):
    self.__fallThroughCmndRegexList = []
    self.__interceptedCmndStructList = []
    self.__allowExtraCmnds = True

  def readCommandsFromStr(self, cmndsStr):
    lines = cmndsStr.splitlines()
    for line in lines:
      #print("line: '" + line + "'")
      if line == "":
        continue
      splitArray = line.split(':')
      (tag, entry) = (splitArray[0], ':'.join(splitArray[1:]))
      #(tag, entry) = line.split(':')
      #print("(tag, entry) = " + str((tag, entry)))
      if tag == "FT":
        self.__fallThroughCmndRegexList.append(entry.strip())
      elif tag == "IT":
        entryArray = entry.split(';')
        if len(entryArray) < 3:
          raise Exception("Error, invalid line {"+line+"}")
        cmndRegex = entryArray[0]
        cmndReturn = entryArray[1]
        cmndOutput = ""
        for cmndOutputEntry in entryArray[2:]:
          #print("cmndOutputEntry = {" + cmndOutputEntry + "}")
          cmndOutput += cmndOutputEntry.strip()[1:-1]+"\n"
        #print("(cmndRegex, cmndReturn, cmndOutput) = " +
        #      str((cmndRegex, cmndReturn, cmndOutput)))
        self.__interceptedCmndStructList.append(
          InterceptedCmndStruct(cmndRegex.strip(), int(cmndReturn), cmndOutput)
          )
      else:
        raise Exception("Error, invalid tag = '"+tag+"'!")

  def assertAllCommandsRun(self):
    if len(self.__interceptedCmndStructList) > 0:
      raise Exception("Error, all of the commands have not been run starting with" \
                      " the command " + str(self.__interceptedCmndStructList[0])
                      + "!")


g_sysCmndInterceptor = SysCmndInterceptor()


# Read the command interepts from a file?
cmndInterceptsFile = os.environ.get(
  "GENERAL_SCRIPT_SUPPORT_CMND_INTERCEPTS_FILE","")
if cmndInterceptsFile:
  cmndInterceptsFileStr = open(cmndInterceptsFile, 'r').read()
  print("\nReading system command intercepts from file '" +
        cmndInterceptsFile                                +
        "' with contents:\n"                              +
        "-----------------------------------\n"           +
        cmndInterceptsFileStr                             +
        "-----------------------------------\n")
  g_sysCmndInterceptor.readCommandsFromStr(cmndInterceptsFileStr)
  g_sysCmndInterceptor.setAllowExtraCmnds(False)


# Dump all commands being performed?
g_dumpAllSysCmnds = "GENERAL_SCRIPT_SUPPORT_DUMD_COMMANDS" in os.environ


def runSysCmndInterface(cmnd, outFile=None, rtnOutput=False, extraEnv=None, \
  workingDir="", getStdErr=False \
  ):
  if g_dumpAllSysCmnds:
    print("\nDUMP SYS CMND: " + cmnd + "\n")
  if outFile!=None and rtnOutput==True:
    raise Exception("Error, both outFile and rtnOutput can not be true!") 
  if g_sysCmndInterceptor.doProcessInterceptedCmnd(cmnd):
    (cmndReturn, cmndOutput) = g_sysCmndInterceptor.nextInterceptedCmndStruct(cmnd)
    if rtnOutput:
      if cmndOutput==None:
        raise Exception("Error, the command '"+cmnd+"' gave None output when" \
                        " non-null output was expected!")
      return (cmndOutput, cmndReturn)
    if outFile:
      writeStrToFile(outFile, cmndOutput)  
    return cmndReturn
  # Else, fall through
  if extraEnv:
    fullEnv = os.environ.copy()
    fullEnv.update(extraEnv)
  else:
    fullEnv = None
  pwd = None
  if workingDir:
    pwd = os.getcwd()
    os.chdir(workingDir)
  rtnObject = None
  try:
    if rtnOutput:
      if getStdErr:
        child = subprocess.Popen(cmnd, shell=True, stdout=subprocess.PIPE,
          stderr = subprocess.STDOUT, env=fullEnv)
      else:
        child = subprocess.Popen(cmnd, shell=True, stdout=subprocess.PIPE,
          env=fullEnv)
      data = child.stdout.read()
      #print("data = '" + str(data) + "'")
      child.wait()
      rtnCode = child.returncode
      #print("rtnCode = '" + str(rtnCode) + "'")
      rtnObject = (data, rtnCode)
    else:
      outFileHandle = None
      if outFile:
        outFileHandle = open(outFile, 'w')
      rtnCode = subprocess.call(cmnd, shell=True, stderr=subprocess.STDOUT,
        stdout=outFileHandle, env=fullEnv)
      rtnObject = rtnCode
  finally:
    if pwd: os.chdir(pwd)
  return rtnObject


######################################
# System interaction utilties
######################################


def runSysCmnd(cmnd, throwExcept=True, outFile=None, workingDir="",
  extraEnv=None \
  ):
  """Run system command and optionally throw on failure"""
  sys.stdout.flush()
  sys.stderr.flush()
  try:
    outFileHandle = None
    rtnCode = runSysCmndInterface(cmnd, outFile=outFile, extraEnv=extraEnv,
      workingDir=workingDir)
  except OSError as e:
    rtnCode = 1 # Just some error code != 0 please!
  if rtnCode != 0 and throwExcept:
    raise RuntimeError("Error, the command '%s' failed with error code %d"
                       % (cmnd, rtnCode))
  return rtnCode


def echoRunSysCmnd(cmnd, throwExcept=True, outFile=None, msg=None,
  timeCmnd=False, verbose=True, workingDir="", returnTimeCmnd=False,
  extraEnv=None
  ):
  """Echo command to be run and run command with runSysCmnd()"""
  if verbose:
    print("\nRunning: " + cmnd + "\n")
    if workingDir:
      print("  Running in working directory: " + workingDir + " ...\n")
    if extraEnv:
      print("  Appending environment:" + sorted_dict_str(extraEnv) + "\n")
    if outFile:
      print("  Writing console output to file " + outFile + " ...")
  if msg and verbose:
    print("  " + msg + "\n")
  t1 = time.time()
  totalTimeMin = -1.0
  try:
    rtn = runSysCmnd(cmnd, throwExcept, outFile, workingDir, extraEnv)
  finally:
    if timeCmnd:
      t2 = time.time()
      totalTimeMin = (t2-t1)/60.0
      if verbose:
        print("\n  Runtime for command = %f minutes" % totalTimeMin)
  if returnTimeCmnd:
    return (rtn, totalTimeMin)
  return rtn


def getCmndOutput(cmnd, stripTrailingSpaces=False, throwOnError=True, workingDir="", \
  getStdErr=False, rtnCode=False \
  ):
  """Run a shell command and return its output"""
  (data, errCode) = runSysCmndInterface(cmnd, rtnOutput=True, workingDir=workingDir,
    getStdErr=getStdErr)
  if errCode != 0:
    if throwOnError:
      raise RuntimeError('%s failed w/ exit code %d:\n\n%s' % (cmnd, errCode, data))
  dataToReturn = data
  if stripTrailingSpaces:
    dataToReturn = data.rstrip()
  if rtnCode:
    return (dataToReturn, errCode)
  return dataToReturn


def pidStillRunning(pid):
  #print("\npid = '" + pid + "'")
  cmnd = "kill -s 0 "+pid
  cmndReturn = runSysCmnd(cmnd, False)
  #print("\ncmndReturn = " + cmndReturn)
  return cmndReturn == 0


######################################
# File/path utilities 
######################################
     

def getFilePathArray(filePathStr):
  return filePathStr.split('/')


def joinDirs(dirArray):
  """
  Join directories.

  2009/06/09: rabartl: We should be able to just use os.path.join(...) but I
  found when used in at least on context that it resulted in not joining the
  elements but instead just returning the array.
  """
  dirPath = ""
  for dir in dirArray:
    if not dirPath:
      dirPath = dir
    else:
      dirPath = dirPath + "/" + dir
  return dirPath


def downDirsArray(numDownDirs):
  downDirsPathArray = []
  for i in range(0, numDownDirs):
    downDirsPathArray.append("..")
  #print("\ndownDirsPathArray = " + downDirsPathArray)
  return downDirsPathArray


def normalizePath(path):
  return os.path.normpath(path)


def echoChDir(dirName, verbose=True):
  if verbose:
    print("\nChanging current directory to \'" + dirName + "\'")
  if not os.path.isdir(dirName):
    raise OSError("Error, the directory \'"+dirName+"\' does not exist in the" \
      + " base directory \'"+os.getcwd()+"\"!" )
  os.chdir(dirName)
  if verbose:
    print("\nCurrent directory is \'" + os.getcwd() + "\'\n")


def createDir(dirName, cdIntoDir=False, verbose=False):
  """Create a directory if it does not exist"""
  if os.path.exists(dirName):
    if not os.path.isdir(dirName):
      errMsg = "\nError the path '" + dirName + \
               "'already exists but it is not a directory!"
      if verbose: print(errMsg)
      raise RuntimeError(errMsg)
    if verbose: print("\nThe directory " + dirName + "already exists!")
  else:
    if verbose: print("\nCreating directory " + dirName + " ...")
    os.mkdir(dirName)
  if cdIntoDir:
    echoChDir(dirName, verbose=verbose)


def createDirsFromPath(path):
  #print("\npath = " + path)
  pathList = path.split("/")
  #print("\npathList = " + pathList)
  if not pathList[0]:
    currDir = "/"
  for dir in pathList:
    currDir = os.path.join(currDir, dir)
    if currDir and not os.path.exists(currDir):
      #print("\ncurrDir = " + currDir)
      createDir(currDir)


def expandDirsDict(trilinosDirsDict_inout):

  for dir in list(trilinosDirsDict_inout):
    subdirsList = dir.split("/")
    #print("\nsubdirsList = " + subdirsList)
    for i in range(len(subdirsList)):
      trilinosDirsDict_inout.update({joinDirs(subdirsList[:i+1]) : 0})


def removeIfExists(fileName):
  if os.path.exists(fileName):
    echoRunSysCmnd("rm "+fileName)


def removeDirIfExists(dirName, verbose=False):
  if os.path.exists(dirName):
    if verbose:
      print("Removing existing directory '" + dirName + "' ...")
    echoRunSysCmnd("rm -rf "+dirName)


def writeStrToFile(fileName, fileBodyStr):
  open(fileName, 'w').write(fileBodyStr)


def readStrFromFile(fileName):
  return open(fileName, 'r').read()


def getFileNamesWithFileTag( baseDir, fileTag ):
  """Get a list of file names with a given date stamp"""
  return getCmndOutput(
    'cd %s && ls *%s*' % (baseDir, fileTag),
    throwOnError=False
    ).split()


def getFileNameFromGlob( baseDir, fileNameGlob ):
  return getCmndOutput("cd "+baseDir+" && ls "+fileNameGlob, True, False)


def isEmptyDir( absDir ):
  return (len(os.listdir(absDir)) == 0)


def getDirSizeInGb(dir):
  sizeIn1024Kb = int(getCmndOutput("du -s "+dir).split('\t')[0])
  #print("\nsizeIn1024Kb = " + str(sizeIn1024Kb))
  return float(sizeIn1024Kb)/1e+6 # Size in Gb!


def isPathChar(char):
  return (char.isalnum() or char == '/') and (not char == ' ')


# 2008/07/08: rabartl: This silly function is needed because on the sun
# machine (i.e. sass8000), the 'which' command returns some non-printable
# chars from the beginning of the output string that don't form a valid path.
# This was *very* hard to debug but this function is able to select the valid
# path string.  This has been tested at least on linux and the sun and should
# work anywhere.
def cleanBadPath(inputPath):
  cleanPath = ""
  for i in range(len(inputPath)-1, -1, -1):
    char = inputPath[i]
    if not isPathChar(char):
      break
    cleanPath = char + cleanPath
  return cleanPath


def getRelativePathFrom1to2(absPath1, absPath2):
  #print("\nabsPath1 =" + absPath1)
  #print("\nabsPath2 =" + absPath2)
  absPath1_array = absPath1.split('/')
  absPath2_array = absPath2.split('/')
  #print("\nabsPath1_array =" + absPath1_array)
  #print("\nabsPath2_array =" + absPath2_array)
  absPath1_array_len = len(absPath1_array)
  absPath2_array_len = len(absPath2_array)
  maxBaseDirDepth = min(absPath1_array_len, absPath2_array_len) 
  baseDirDepth = 0
  for dirIdx in range(0, maxBaseDirDepth):
    dir1 = absPath1_array[dirIdx]
    dir2 = absPath2_array[dirIdx]
    if dir1 != dir2:
      break
    baseDirDepth = baseDirDepth + 1
  #print("\nbaseDirDepth = %d" % baseDirDepth)
  numDownDirs = absPath1_array_len - baseDirDepth
  #print("\nnumDownDirs = %d" % numDownDirs)
  if numDownDirs > 0:
    downDirPath = joinDirs(downDirsArray(numDownDirs))
  else:
    downDirPath = "."
  #print("\ndownDirPath = '" + downDirPath + "'")
  if baseDirDepth == absPath2_array_len:
    upDirPath = "."
  else:
    upDirPath = joinDirs(absPath2_array[baseDirDepth:])
  #print("\nupDirPath = "  + upDirPath      )
  #print("\nupDirPath = '" + upDirPath + "'")
  relPath = os.path.join(downDirPath, upDirPath)
  if relPath == "./.":
    return "."
  return relPath


def getExecBaseDir(execName):
  whichOutput = getCmndOutput("type -p "+execName, True, False)
  # Handle the outpue 'execName is execFullPath' output
  execFullPath = whichOutput.split(' ')[-1]
  #print("\nexecFullPath = " + execFullPath)
  execNameMatchRe = r"^(.+)/"+execName
  execNameGroups = re.findall(execNameMatchRe, execFullPath)
  #print("\nexecNameGroups = " + execNameGroups)
  if not execNameGroups:
    return None
  execBaseDir = execNameGroups[0]
  #print("\nexecBaseDir = \"" + execBaseDir + "\"")
  #execBaseDir = cleanBadPath(execBaseDir)
  #print("\nexecBaseDir = \"" + execBaseDir + "\"")
  return execBaseDir


def extractAppendUniqueDirsDictFromFileNames(filenamesArray, dirsDict):
  for filename in filenamesArray:
    dirsDict.update( { normalizePath(os.path.dirname(filename)) : 0 } )


def copyFileAndReplaceTokens( scriptsDir, inputFile, tokenReplacementList,
  outputFile ):

  """Copies an input stub file and then does a set of token replacements"""

  echoRunSysCmnd("cp "+inputFile+" "+outputFile, verbose=verboseDebug)

  for tokenReplacementPair in tokenReplacementList:
    oldToken = tokenReplacementPair[0]
    newToken = tokenReplacementPair[1]
    echoRunSysCmnd( scriptsDir+"/token-replace.pl "+oldToken+" "+newToken\
      +" "+outputFile+" "+outputFile, verbose=verboseDebug );
    # ToDo: Replace above with native re commands


class TeeOutput(object):
  """
  Object that directs all calls to its write method to stdout as well
  as a file. This is to be used as a simple replacement for the Unix
  tee command.
  """
  def __init__(self, outputfile):
    """ Constructor takes a file-like object to write output to."""
    self._realstdout = sys.stdout
    self._outputfile = outputfile

  def _safe_outputfile_method(self, methodname, *args):
    """
    Calls the method specified by methodname with the given args on
    the internal file object if it is non-null.
    """
    if self._outputfile is not None:
      if hasattr(self._outputfile, methodname):
        method = getattr(self._outputfile, methodname)
        if method and callable(method):
          method(*args)

  def write(self, data):
    """
    Write the given data to stdout and to the log file.
    """
    self._realstdout.write(data)
    self._safe_outputfile_method('write', data)

  def flush(self):
    """
    Flush the internal file buffers.
    """
    self._realstdout.flush()
    self._safe_outputfile_method('flush')


######################################
# Shell argument helpers
######################################


reCmndLineArg = re.compile(r"(--.+=)(.+)")


def requoteCmndLineArgs(inArgs):
  argsStr = ""
  for arg in inArgs:
    splitArg = arg.split("=")
    newArg = None
    if len(splitArg) == 1:
      newArg = arg
    else:
      newArg = splitArg[0]+"=\""+'='.join(splitArg[1:])+"\""
    #print("\nnewArg = " + newArg)
    argsStr = argsStr+" "+newArg
  return argsStr


def commandLineOptionsToList(stringOptions):
  """
  Convert a string of space separated command line options to a python
  list of the individual optionstrings.
  TODO: Handle shell quoting.
  """
  return stringOptions.split()


class ConfigurableOptionParser(optparse.OptionParser):
  """
  OptionParser that accepts a python dictionary as a configuration
  file to provide default value overrides for the options.
  """
  def __init__(self, configuration, **kwargs):
    """
    Constructor accepts a configuration dictionary with default values
    for arguments and all of the OptionParser arguments as well.
    """
    optparse.OptionParser.__init__(self, **kwargs)
    self._configuration = configuration

  def add_option(self, *args, **kwargs):
    """
    Checks for a default override in the configuration dictionary and
    modifies the default and help arguments before dispatching them to
    the base class implementation.
    """
    if 'default' in kwargs:
      for arg in args:
        kwargs['default'] = self._configuration.get(arg, kwargs['default'])
    optparse.OptionParser.add_option(self, *args, **kwargs)


######################################
# Debugging support
######################################


def printStackTrace():
  sys.stdout.flush()
  traceback.print_exc()
    

class ErrorCaptureOptionParser(optparse.OptionParser):
  __sawError = None
  def __init__(self, usage="%prog [options]", version=None):
    optparse.OptionParser.__init__(self, usage, version)
    __sawError = False
  def error(self, msg):
    raise Exception("Received error message: " + msg)


######################################
# HTML directory browsing
######################################


def createIndexHtmlBrowserList(baseDir, fileDirList = None):
  htmlList = ""
  # Get the fileDirList if empty
  if not fileDirList:
    fileDirList = os.listdir(baseDir)
    fileDirList.sort()
  # Create the HTML header
  htmlList += "" \
    + "<ul>\n" \
    + "<li><a href=\"..\">..</a></li>\n"
  # Fill in links to directories first
  for fd in fileDirList:
    absFd = baseDir+"/"+fd
    #print("isfile(" + fd + ") = " + str(os.path.isfile(absFd)))
    #print("isdir("  + fd + ") = " + str(os.path.isdir(absFd) ))
    if not os.path.isfile(absFd):
      htmlList = htmlList \
                 +"<li>dir: <a href=\""+fd+"\">"+fd+"</a></li>\n"
  # Fill in links to regular files second
  for fd in fileDirList:
    absFd = baseDir+"/"+fd
    if os.path.isfile(absFd):
      if fd != 'index.html':
        htmlList = htmlList \
                 +"<li>file: <a href=\""+fd+"\">"+fd+"</a></li>\n"
  # Write the footer
  htmlList = htmlList \
    + "</ul>\n"
  return htmlList



def createIndexHtmlBrowserFile(baseDir, fileDirList):
  """Creates an HTML browser file as a returned string."""
  htmlFile = "" \
    + "<html>\n" \
    + "<head>\n" \
    + "<title>"+baseDir+"</title>\n" \
    + "</head>\n" \
    + "<body>\n" \
    + "<b>"+baseDir+"</b>\n" \
    + createIndexHtmlBrowserList(baseDir, fileDirList) \
    + "</body>\n" \
    + "</html>\n"
  return htmlFile


def createHtmlBrowserFiles(absBaseDir, depth, verbose=False):

  """Create a hierarchy of index.html files that will build a directory/file
  browser for a web server that will not allow directory/file browsing."""

  #print("\nEntering createHtmlBrowserFiles(" + absBaseDir + ",%d" % (depth) + ")")

  # Get the list of all of the files/directories in absBaseDir
  fileDirList = os.listdir(absBaseDir)
  fileDirList.sort()
  #print("\nfileDirList = " + str(fileDirList)
  #sys.stdout.flush()

  # Get the index.html file HTML
  indexHtml = createIndexHtmlBrowserFile(absBaseDir, fileDirList)
  #print("\nindexHtml:\n" + indexHtml)

  # Write the index.html file
  indexFileName = absBaseDir+"/index.html"
  if verbose:
    print("\nWriting " + indexFileName)
  open(indexFileName,'w').write(indexHtml)

  # Loop through all of the directories and recursively call this function
  if depth > 0:
    for fd in fileDirList:
      absFd = absBaseDir+"/"+fd
      if os.path.isdir(absFd):
        subDir = absFd
        #print("\nCalling createHtmlBrowserFiles(" + subDir + ",%d" % (depth-1)
        #      + ")")
        createHtmlBrowserFiles(absBaseDir+"/"+fd,depth-1)

  #print("\nLeaving createHtmlBrowserFiles(" + absBaseDir + ",%d" % (depth) +
  #      ")") 
