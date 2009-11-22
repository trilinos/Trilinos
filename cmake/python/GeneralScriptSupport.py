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


#
# Functions
#


def getScriptBaseDir():
  return os.path.abspath(os.path.dirname(sys.argv[0]))


def getScriptName():
  return os.path.basename(os.path.dirname(sys.argv[0]))


def findInSequence(seq, item):
 for i in range(0, len(seq)):
   if seq[i] == item:
     return i
 return -1
   

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


#########################################
# System command unit testing utiltities
#########################################


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
     for cmndRegex in self.__fallThroughCmndRegexList:
       if re.match(cmndRegex, cmnd):
         return True
     return False

  def nextInterceptedCmndStruct(self, cmnd):
    assert(not self.isFallThroughCmnd(cmnd))
    if len(self.__interceptedCmndStructList) == 0:
      raise Exception("Error, cmnd='"+cmnd+"' is past the last expected command!")
    ics = self.__interceptedCmndStructList.pop(0)
    if re.match(ics.cmndRegex, cmnd):
      return (ics.cmndReturn, ics.cmndOutput)
    raise Exception("Error, cmnd='"+cmnd+"' did not match the" \
      " expected regex='"+ics.cmndRegex+"'!")

  def clear(self):
    self.__fallThroughCmndRegexList = []
    self.__interceptedCmndStructList = []
    self.__allowExtraCmnds = True

  def readCommandsFromStr(self, cmndsStr):
    lines = cmndsStr.split('\n')
    for line in lines:
      if line == "":
        continue
      (tag, entry) = line.split(':')
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
      writeStrToFile(cmndOutput, outFile)  
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


###############
# File helpers
###############


def removeIfExists(fileName):
  if os.path.exists(fileName):
    echoRunSysCmnd("rm "+fileName)


def writeStrToFile(fileBodyStr, fileName):
  open(fileName, 'w').write(fileBodyStr)


def readStrFromFile(fileName):
  return open(fileName, 'r').read()


#
# System interaction utilties
#


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
      else:
        return ""
    if stripTrailingSpaces:
      return data.rstrip()
    return data
  finally:
    if pwd: os.chdir(pwd)


#####################
# Other Stuff
#####################


def printStackTrace():
  sys.stdout.flush()
  traceback.print_exc()


def getFileNamesWithFileTag( baseDir, fileTag ):
  """Get a list of file names with a given date stamp"""
  return getCmndOutput(
    'cd %s && ls *%s*' % (baseDir, fileTag),
    throwOnError=False
    ).split()


def pidStillRunning(pid):
  #print "\npid = '"+pid+"'"
  cmnd = "kill -s 0 "+pid
  cmndReturn = runSysCmnd(cmnd, False)
  #print "\ncmndReturn =", cmndReturn
  return cmndReturn == 0


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
    

class ErrorCaptureOptionParser(optparse.OptionParser):
  __sawError = None
  def __init__(self, usage="%prog [options]", version=None):
    optparse.OptionParser.__init__(self, usage, version)
    __sawError = False
  def error(self, msg):
    raise "Received error message: "+msg


def stripWhitespaceFromStringList(stringList):
  return [x.strip() for x in stringList]


def removeItemsFromList(list, items):
  itemsRemoved = 0
  for item in items:
    if item in list:
      idx = list.index(item)
      del list[idx]
      itemsRemoved = itemsRemoved + 1
  return itemsRemoved


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
    #print "\nnewArg =", newArg
    argsStr = argsStr+" "+newArg
  return argsStr


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


def getNewDateStr():
  return getCmndOutput(r"date '+%F-%H-%M'",True)


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


reDateDirName = re.compile(r"[0-9]{4,4}-[0-9]{2,2}-[0-9]{2,2}-[0-9]{2,2}-[0-9]{2,2}")


def isNameDateDir(name):
  if (reDateDirName.match(name)):
    return True
  return False


def getDateDirs(absBaseDir):

  fileDirList = os.listdir(absBaseDir)
  fileDirList.sort()
  #print "\nfileDirList =", fileDirList
  #sys.stdout.flush()

  dateDirs = []

  for fd in fileDirList:
    absFd = absBaseDir+"/"+fd
    if os.path.isdir(absFd) and isNameDateDir(fd):
      dateDirs.append(fd)

  dateDirs.sort(reverse=True)
  return dateDirs


def getBuildDirs(absBaseDir):

  fileDirList = os.listdir(absBaseDir)
  #print "\nfileDirList =", fileDirList
  #sys.stdout.flush()

  buildDirs = []

  for fd in fileDirList:
    absFd = absBaseDir+"/"+fd
    if os.path.isdir(absFd):
      buildDirs.append(fd)

  return buildDirs


def getDateTimeFromDateStr(dateStr):
  dta = dateStr.split("-")
  return datetime.datetime(
    year=int(dta[0]), month=int(dta[1]), day=int(dta[2]),
    hour=int(dta[3]), minute=int(dta[4]) )


relativeTimeDeltaFormatHelp = \
  "The format for relative time delta is wW+dD+hH+mM for w weeks,"\
  +" d days, h hours, and m minutes (all integer values)."


def convertTimeDeltaStrToTimeDelta(timeDeltaStr):
  """Converts time delta strings of the form wW+dD+hH+mM to timedelta object"""
  #print "\ntimeDeltaStr =", timeDeltaStr
  if timeDeltaStr:
    tdsa = timeDeltaStr.split("+")
  else:
    tdsa = [" "]
  #print "\ntdsa =", tdsa
  weeks = 0
  days = 0
  hours = 0
  minutes = 0
  for dt in tdsa:
    valStr = dt[:-1]
    if dt[-1] == 'W':
      weeks = int(valStr)
    elif dt[-1] == 'D':
      days = int(valStr)
    elif dt[-1] == 'H':
      hours = int(valStr)
    elif dt[-1] == 'M':
      minutes = int(valStr)
    else:
      raise "Error, the time-delta string \""+timeDeltaStr+"\" is not a valid"\
        +" format.  "+relativeTimeDeltaFormatHelp
  return datetime.timedelta(
    weeks=weeks, days=days, hours=hours, minutes=minutes )


def thinOutFilesAndDirsInBaseDir(maxFileKb, absBaseDir, localBaseDir,
  protectedFilesAndDirs, noOp, verbose
  ):

  #localVerboseDebug = True
  localVerboseDebug = False

  if localVerboseDebug:
    print "\nEntering thinOutFilesAndDirsInBaseDir(",maxFileKb, ",", absBaseDir,\
     ",", localBaseDir, "... )"

  absLocalBaseDir = os.path.join(absBaseDir, localBaseDir)

  fileDirList = os.listdir(absLocalBaseDir)
  fileDirList.sort()
  #print "\nfileDirList =", fileDirList
  #sys.stdout.flush()

  for fd in fileDirList:

    # Get the relative file/directory path to the base directory
    localFd = os.path.join(localBaseDir, fd)
    if localVerboseDebug: print "\nlocalFd =", localFd

    # Get the absolute file or directory path
    absFd = os.path.join(absBaseDir, localFd)
    if localVerboseDebug: print "\nabsFd =", absFd

    # Determine if this is a protected file/directory
    isProtectedFd = (localFd in protectedFilesAndDirs)

    # Thin out or remove file/directory if it is not protected
    if not isProtectedFd:
      
      # If this is a directory, thin it out and remove it if needed/requested
      if os.path.isdir(absFd):
        thinOutFilesAndDirsInBaseDir(maxFileKb, absBaseDir,
          os.path.join(localBaseDir, fd), protectedFilesAndDirs, noOp, verbose )
        if isEmptyDir(absFd) and maxFileKb < 0 and not noOp:
          runSysCmnd("rm -rf "+absFd)

      # If this is a regular file, then copy over it or delete if needed/requested
      elif os.path.isfile(absFd) and not noOp:
        if maxFileKb > 0 and int(float(os.path.getsize(absFd))/1000.0) > maxFileKb:
          open(absFd, 'w').write(
           'Removed file because it was bigger than %d KB\n' % maxFileKb )
        elif maxFileKb <= 0:
          if localVerboseDebug:
            echoRunSysCmnd("rm -f "+absFd, False)
          else:
            runSysCmnd("rm -f "+absFd, False)

    else:

      if verbose:
        print "\nLeaving untouched \""+localFd+"\" since it is protected!"

  if localVerboseDebug:
    print "\nLeaving thinOutFilesAndDirsInBaseDir(",maxFileKb, ",", absBaseDir,\
     ",", localBaseDir, ", ... )"


def removeLargerFilesInBaseDir(maxFileKb, absBaseDir, protectedFilesAndDirs,
  noOp, verbose
  ):
  thinOutFilesAndDirsInBaseDir(maxFileKb, absBaseDir, "", protectedFilesAndDirs,
    noOp, verbose)


def deleteAllButProtectedInBaseDir(absBaseDir, protectedFilesAndDirs,
  noOp, verbose
  ):
  thinOutFilesAndDirsInBaseDir(-1, absBaseDir, "", protectedFilesAndDirs,
    noOp, verbose)
