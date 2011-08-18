# @HEADER
# ************************************************************************
#
#            Trilinos: An Object-Oriented Solver Framework
#                 Copyright (2001) Sandia Corporation
#
#
# Copyright (2001) Sandia Corporation. Under the terms of Contract
# DE-AC04-94AL85000, there is a non-exclusive license for use of this
# work by or on behalf of the U.S. Government.  Export of this program
# may require a license from the United States Government.
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
# NOTICE:  The United States Government is granted for itself and others
# acting on its behalf a paid-up, nonexclusive, irrevocable worldwide
# license in this data to reproduce, prepare derivative works, and
# perform publicly and display publicly.  Beginning five (5) years from
# July 25, 2001, the United States Government is granted for itself and
# others acting on its behalf a paid-up, nonexclusive, irrevocable
# worldwide license in this data to reproduce, prepare derivative works,
# distribute copies to the public, perform publicly and display
# publicly, and to permit others to do so.
#
# NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED STATES DEPARTMENT
# OF ENERGY, NOR SANDIA CORPORATION, NOR ANY OF THEIR EMPLOYEES, MAKES
# ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR
# RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY
# INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS
# THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.
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
if sys.version < '2.4':
   print "Error, Python version is "+sys.version+" < 2.4!"
   sys.exit(1)

#
# Import commands
#

import sys
import os
import re
import math
import subprocess
import time
import datetime
import optparse
import traceback


verboseDebug = False


#
# Determine what system we are on:
#

rePlatformName = re.compile(r"^[a-zA-Z]+")
platformStr = rePlatformName.findall(sys.platform)[0]
#print "\nplatformStr =", platformStr


######################################
# Script location functions
######################################


def getScriptBaseDir():
  return os.path.abspath(os.path.dirname(sys.argv[0]))


def getScriptName():
  return os.path.basename(os.path.dirname(sys.argv[0]))


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
  #print "regex_in =", regex_in
  reMatch = re.compile(regex_in)
  linesExtracted = ""
  foundRegex = False
  for line in string_in.strip().split("\n"):
    #print "line = '" + line + "'"
    if not foundRegex:
      matchObj = reMatch.match(line)
      #print "matchObj =", matchObj
      if matchObj:
        foundRegex = True
    if foundRegex:
      linesExtracted += line + "\n"
  return linesExtracted
  

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
      raise Exception("Error, cmnd='"+cmnd+"' did not match the" \
        " expected regex='"+ics.cmndRegex+"'!")
    self.__interceptedCmndStructList.pop(0)
    return (ics.cmndReturn, ics.cmndOutput)

  def clear(self):
    self.__fallThroughCmndRegexList = []
    self.__interceptedCmndStructList = []
    self.__allowExtraCmnds = True

  def readCommandsFromStr(self, cmndsStr):
    lines = cmndsStr.split('\n')
    for line in lines:
      if line == "":
        continue
      splitArray = line.split(':')
      (tag, entry) = (splitArray[0], ':'.join(splitArray[1:]))
      #(tag, entry) = line.split(':')
      #print "(tag, entry) =", (tag, entry)
      if tag == "FT":
        self.__fallThroughCmndRegexList.append(entry.strip())
      elif tag == "IT":
        (cmndRegex, cmndReturn, cmndOutput) = entry.split(';')
        self.__interceptedCmndStructList.append(
          InterceptedCmndStruct(cmndRegex.strip(), int(cmndReturn),
            cmndOutput.strip()[1:-1] )
          )
      else:
        raise Exception("Error, invalid tag = '"+tag+"'!")

  def assertAllCommandsRun(self):
    if len(self.__interceptedCmndStructList) > 0:
      raise Exception("Error, all of the commands have not been run starting with" \
        " the command "+str(self.__interceptedCmndStructList[0])+"!")


g_sysCmndInterceptor = SysCmndInterceptor()


# Read the command interepts from a file?
cmndInterceptsFile = os.environ.get(
  "GENERAL_SCRIPT_SUPPORT_CMND_INTERCEPTS_FILE","")
if cmndInterceptsFile:
  cmndInterceptsFileStr = open(cmndInterceptsFile, 'r').read()
  print "\nReading system command intercepts from file '"+cmndInterceptsFile+"' with contents:\n" \
    "-----------------------------------\n" \
    +cmndInterceptsFileStr+ \
    "-----------------------------------\n"
  g_sysCmndInterceptor.readCommandsFromStr(cmndInterceptsFileStr)
  g_sysCmndInterceptor.setAllowExtraCmnds(False)


# Dump all commands being performed?
g_dumpAllSysCmnds = os.environ.has_key("GENERAL_SCRIPT_SUPPORT_DUMD_COMMANDS")


def runSysCmndInterface(cmnd, outFile=None, rtnOutput=False):
  if g_dumpAllSysCmnds:
    print "\nDUMP SYS CMND: " + cmnd + "\n"
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
  if rtnOutput:
    child = os.popen(cmnd)
    data = child.read()
    rtnCode = child.close()
    return (data, rtnCode)
  else:
    outFileHandle = None
    if outFile:
      outFileHandle = open(outFile, 'w')
    rtnCode = subprocess.call(cmnd, shell=True, stderr=subprocess.STDOUT,
      stdout=outFileHandle)
    return rtnCode


######################################
# System interaction utilties
######################################


def runSysCmnd(cmnd, throwExcept=True, outFile=None, workingDir=""):
  """Run system command and optionally throw on failure"""
  sys.stdout.flush()
  sys.stderr.flush()
  try:
    if workingDir:
      pwd = os.getcwd()
      os.chdir(workingDir)
    #rtnCode = subprocess.call(cmnd, shell=True)
    outFileHandle = None
    rtnCode = runSysCmndInterface(cmnd, outFile)
#    if outFile:
#       outFileHandle = open(outFile, 'w')
#    rtnCode = subprocess.call(cmnd, shell=True, stderr=subprocess.STDOUT, stdout=outFileHandle)
  except OSError, e:
    rtnCode = 1 # Just some error code != 0 please!
  if workingDir:
    os.chdir(pwd)
  if rtnCode != 0 and throwExcept:
    raise RuntimeError('Error, the command \'%s\' failed with error code %d' \
                       % (cmnd,rtnCode) )
  return rtnCode


def echoRunSysCmnd(cmnd, throwExcept=True, outFile=None, msg=None,
  timeCmnd=False, verbose=True, workingDir="", returnTimeCmnd=False
  ):
  """Echo command to be run and run command with runSysCmnd()"""
  if verbose:
    print "\nRunning: "+cmnd+"\n"
    if workingDir:
      print "  Running in working directory: "+workingDir+" ...\n"
    if outFile:
      print "  Writing console output to file "+outFile+" ..."
  if msg and verbose:
    print "  "+msg+"\n"
  t1 = time.time()
  totalTimeMin = -1.0
  try:
    rtn = runSysCmnd(cmnd, throwExcept, outFile, workingDir)
  finally:
    if timeCmnd:
      t2 = time.time()
      totalTimeMin = (t2-t1)/60.0
      if verbose:
        print "\n  Runtime for command = %f minutes" % totalTimeMin
  if returnTimeCmnd:
    return (rtn, totalTimeMin)
  return rtn


def getCmndOutput(cmnd, stripTrailingSpaces=False, throwOnError=True, workingDir=""):
  """Run a shell command and return its output"""
  pwd = None
  if workingDir:
    pwd = os.getcwd()
    os.chdir(workingDir)
  try:
    (data, err) = runSysCmndInterface(cmnd, rtnOutput=True)
    if err:
      if throwOnError:
        raise RuntimeError, '%s failed w/ exit code %d' % (cmnd, err)
    if stripTrailingSpaces:
      return data.rstrip()
    return data
  finally:
    if pwd: os.chdir(pwd)


def pidStillRunning(pid):
  #print "\npid = '"+pid+"'"
  cmnd = "kill -s 0 "+pid
  cmndReturn = runSysCmnd(cmnd, False)
  #print "\ncmndReturn =", cmndReturn
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
  #print "\ndownDirsPathArray =", downDirsPathArray
  return downDirsPathArray


def normalizePath(path):
  return os.path.normpath(path)


def echoChDir(dirName, verbose=True):
  if verbose:
    print "\nChanging current directory to \'"+dirName+"\'"
  if not os.path.isdir(dirName):
    raise OSError("Error, the directory \'"+dirName+"\' does not exist in the" \
      + " base directory \'"+os.getcwd()+"\"!" )
  os.chdir(dirName)
  if verbose:
    print "\nCurrent directory is \'"+os.getcwd()+"\'\n"


def createDir(dirName, cdIntoDir=False, verbose=False):
  """Create a directory if it does not exist"""
  if os.path.exists(dirName):
    if not os.path.isdir(dirName):
      errMsg = "\nError the path '"+dirName+"'already exists but it is not a directory!"
      if verbose: print errMsg
      raise RuntimeError(errMsg)
    if verbose: print "\nThe directory", dirName, "already exists!"
  else:
    if verbose: print "\nCreating directory "+dirName+" ..."
    os.mkdir(dirName)
  if cdIntoDir:
    echoChDir(dirName, verbose=verbose)


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


def expandDirsDict(trilinosDirsDict_inout):

  for dir in trilinosDirsDict_inout.keys():
    subdirsList = dir.split("/")
    #print "\nsubdirsList =", subdirsList
    for i in range(0, len(subdirsList)):
      trilinosDirsDict_inout.update(
        { joinDirs(subdirsList[:i+1]) : 0 }
        )


def removeIfExists(fileName):
  if os.path.exists(fileName):
    echoRunSysCmnd("rm "+fileName)


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
  #print "\nsizeIn1024Kb =", sizeIn1024Kb
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
  #print "\nabsPath1 =", absPath1
  #print "\nabsPath2 =", absPath2
  absPath1_array = absPath1.split('/')
  absPath2_array = absPath2.split('/')
  #print "\nabsPath1_array =", absPath1_array
  #print "\nabsPath2_array =", absPath2_array
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
  #print "\nbaseDirDepth =", baseDirDepth
  numDownDirs = absPath1_array_len - baseDirDepth
  #print "\nnumDownDirs =", numDownDirs
  if numDownDirs > 0:
    downDirPath = joinDirs(downDirsArray(numDownDirs))
  else:
    downDirPath = "."
  #print "\ndownDirPath = '" + downDirPath + "'"
  if baseDirDepth == absPath2_array_len:
    upDirPath = "."
  else:
    upDirPath = joinDirs(absPath2_array[baseDirDepth:])
  #print "\nupDirPath =", upDirPath
  #print "\nupDirPath = '" + upDirPath + "'"
  relPath = os.path.join(downDirPath, upDirPath)
  if relPath == "./.":
    return "."
  return relPath


def getExecBaseDir(execName):
  whichOutput = getCmndOutput("type -p "+execName, True, False)
  # Handle the outpue 'execName is execFullPath' output
  execFullPath = whichOutput.split(' ')[-1]
  #print "\nexecFullPath =", execFullPath
  execNameMatchRe = r"^(.+)/"+execName
  execNameGroups = re.findall(execNameMatchRe, execFullPath)
  #print "\nexecNameGroups =", execNameGroups
  if not execNameGroups:
    return None
  execBaseDir = execNameGroups[0]
  #print "\nexecBaseDir = \""+execBaseDir+"\""
  #execBaseDir = cleanBadPath(execBaseDir)
  #print "\nexecBaseDir = \""+execBaseDir+"\""
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


######################################
# Shell argument helpers
######################################


reCmndLineArg = re.compile(r"(--.+=)(.+)")


def addOptionParserChoiceOption(
  optionName,
  optionDest,
  choiceOptions,
  defaultChoiceIndex,
  helpStr,
  optionParser
  ):
  """ Add a general choice option to a optparse.OptionParser object"""
  defaultOptionValue = choiceOptions[defaultChoiceIndex]
  optionParser.add_option(
    optionName,
    dest=optionDest,
    type="choice",
    choices=choiceOptions,
    default=defaultOptionValue,
    help='%s Choices = (\'%s\').  [default = \'%s\']'
    % (helpStr, '\', \''.join(choiceOptions), defaultOptionValue)
    )


def requoteCmndLineArgs(inArgs):
  argsStr = ""
  for arg in inArgs:
    splitArg = arg.split("=")
    newArg = None
    if len(splitArg) == 1:
      newArg = arg
    else:
      newArg = splitArg[0]+"=\""+'='.join(splitArg[1:])+"\""
    #print "\nnewArg =", newArg
    argsStr = argsStr+" "+newArg
  return argsStr


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
    raise "Received error message: "+msg


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
    #print "isfile("+fd+") =", os.path.isfile(absFd)
    #print "isdir("+fd+") =", os.path.isdir(absFd)
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

  #print "\nEntering createHtmlBrowserFiles("+absBaseDir+",%d"%(depth)+")"

  # Get the list of all of the files/directories in absBaseDir
  fileDirList = os.listdir(absBaseDir)
  fileDirList.sort()
  #print "\nfileDirList =", fileDirList
  #sys.stdout.flush()

  # Get the index.html file HTML
  indexHtml = createIndexHtmlBrowserFile(absBaseDir, fileDirList)
  #print "\nindexHtml:\n", indexHtml

  # Write the index.html file
  indexFileName = absBaseDir+"/index.html"
  if verbose:
    print "\nWriting "+indexFileName
  open(indexFileName,'w').write(indexHtml)

  # Loop through all of the directories and recursively call this function
  if depth > 0:
    for fd in fileDirList:
      absFd = absBaseDir+"/"+fd
      if os.path.isdir(absFd):
        subDir = absFd
        #print "\nCalling createHtmlBrowserFiles("+subDir+",%d"%(depth-1)+")"
        createHtmlBrowserFiles(absBaseDir+"/"+fd,depth-1)

  #print "\nLeaving createHtmlBrowserFiles("+absBaseDir+",%d"%(depth)+")"
