"""Set the python search path to include the library directory created by
distutils."""

# System includes
import os
import sys

# Get the build information
uname        = os.uname()
sysName      = uname[0].lower()
sysNameLen   = len(sysName)
vInfo        = sys.version_info
pyVersion    = str(vInfo[0]) + "." + str(vInfo[1])
pyVersionLen = len(pyVersion)

# Search the build directory for this machine
build   = os.path.join("..", "build")
dirList = os.listdir(build)
libDir  = None
for dir in dirList:
    if ((dir[:4]             == "lib."   ) and
        (dir[4:sysNameLen+4] == sysName  ) and
        (dir[-pyVersionLen:] == pyVersion)):
        libDir = os.path.join(build, dir)
        break

# Insert the library directory name at the beginning of
# the python search path
if libDir:
    sys.path.insert(0,libDir)
