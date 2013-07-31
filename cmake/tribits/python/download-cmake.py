#!/usr/bin/env python

# @HEADER
# ************************************************************************
#
#            TriBITS: Tribial Build, Integrate, and Test System
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

#
# General scripting support
#
# NOTE: Included first to check the version of python!
#

from GeneralScriptSupport import *

from CMakeBinaries import *

import glob
import re
import shutil
import string
import subprocess
import sys
import tarfile
import urllib
import zipfile
import socket


#
# Discover domain (detect if in .sandia.gov) if possible.
# Use it to set default_http_proxy value appropriately.
#

default_http_proxy = ""
hostname = ""

hostname = socket.getfqdn()

if hostname != '':
  hostname = hostname.strip()
  hostname = hostname.lower()

  if hostname[-11:] == ".sandia.gov":
    default_http_proxy = "http://wwwproxy.sandia.gov:80/"


#
# Detailed documentation
#

usageHelp = r"""download-cmake.py [--install-dir=$HOME/cmake ...]

Script that downloads and installs one of the verified versions of CMake for
building Trilinos. The verified versions of CMake and their available
installers are indicated by the files CMakeVersions.py and CMakeBinaries.py.

There are multiple verified versions: you can choose which one you want to use
by passing the --installer-type argument to this script.

Valid values for --installer-type are:

  'min' - minimum required version of CMake that can build Trilinos
  'release' - default / recommended version, may be same as min
  'rc' - latest release candidate, may be same as release
  'dev' - latest development build

If your computer is in the sandia.gov domain, this script should automatically
detect that and use the correct value for --http-proxy to enable downloading
files from the web with http. If you need to pass a proxy in manually, pass
the fully qualified form, as in: --http-proxy="http://wwwproxy.sandia.gov:80/"
If you do not specify an --http-proxy and one is not automatically set due to
sandia.gov domain detection, then the http_proxy environment variable is
honored.

By default, if you just type:

   $ SOME_DIR/download-cmake.py

then the directory download_area will get created in the local working
directory. It will contain the extracted install tree of a pre-built
binary CMake suitable for this platform. That extracted install tree may
be used directly from the download_area directory, or you may move or
copy it to any location of your choosing.

This script will also pull down the list of available installers from the CMake
web site and detect the latest available installer for the v2.8 and development
builds. After running the detection phase of the script, a new CMakeVersions.py
is written in the download_area directory. Periodically as needed, this
generated CMakeVersions.py should be committed as the official CMakeVersions.py
file in the same directory with this script.

You can control various parts of the process with the options (see below). Call
this script with -h or --help for a list of possible options.

Example using the --install-dir option:

  $ SOME_DIR/download-cmake.py --install-dir=$HOME/cmake

This usage would install CMake and the other executables in $HOME/cmake/bin.
NOTE: You will have to update your PATH variable to include whatever directory
you choose to install CMake in.

NOTE: If you need to use sudo to install in some place that requires root
privileges, or if you need to install to a directory that already exists, do:

  $ ./download-cmake.py --skip-install [other options]
  $ cd download_area
  $ sudo cp -r cmake-2.8.0-Linux-i386/ /usr/local

Alternatively, if you have sudo access you can do the install in one shot as:

  $ sudo ./download-cmake.py --install-dir=/usr/local

After you have done a successful install, you may remove the downloaded files:

  $ rm -r download_area

Enjoy CMake!

"""


#
# Read in the commandline arguments
#

from optparse import OptionParser

clp = OptionParser(usage=usageHelp)

if sys.platform == 'darwin':
  defInstallDir = "/Applications"
  defSymlinksDir = "/usr/local/bin"
else:
  defInstallDir = "/usr/local"
  defSymlinksDir = ""

clp.add_option(
  "--all-platforms", dest="allPlatforms", action="store_true", default=False,
  help="Download and/or extract tarballs for all platforms (default = just this platform)" )

clp.add_option(
  "--http-proxy", dest="httpProxy", type="string", default=default_http_proxy,
  help="Proxy in the form 'http://server:port/' - use if you are behind a firewall with respect to downloading from http://www.cmake.org (default = \""+default_http_proxy+"\")." )

clp.add_option(
  "--install-dir", dest="installDir", type="string", default=defInstallDir,
  help="The install directory for CMake (default = "+defInstallDir+"." )

clp.add_option(
  "--installer-type", dest="installerType", type="string", default="release",
  help="Which CMake installer: min, release, rc, or dev (default = release)." )

clp.add_option(
  "--skip-detect", dest="skipDetect", action="store_true", default=False,
  help="Skip detecting the latest available builds" )

clp.add_option(
  "--skip-download", dest="skipDownload", action="store_true", default=False,
  help="Skip the download step" )

clp.add_option(
  "--skip-extract", dest="skipExtract", action="store_true", default=False,
  help="Skip the extract step" )

clp.add_option(
  "--skip-install", dest="skipInstall", action="store_true", default=False,
  help="Skip the install step" )

clp.add_option(
  "--symlinks", dest="symlinks", action="store_true", default=False,
  help="Create symlinks to the installed CMake executables in --symlinks-dir" )

clp.add_option(
  "--symlinks-dir", dest="symlinksDir", type="string", default=defSymlinksDir,
  help="The directory in which to create symlinks to CMake executables (default = "+defSymlinksDir+")." )

(options, args) = clp.parse_args()


# "Latest available CMake build" is defined as the most recent build on the
# server in vdir with pre-built binary tarballs available for Mac, Linux and
# Windows.
#
def DetectLatestCMakeBuilds(basedir, baseurl, vdir):
  url = ''.join([baseurl, "/", vdir, "/?C=M;O=D"])
  filename = ''.join([basedir, "/CMake_", vdir, ".html"])

  try:
    createDir(basedir)
  except:
    if not os.path.exists(basedir):
      raise

  print "Querying " + url + "..."

  proxies = None # if None, use proxy from env var http_proxy
  if not options.httpProxy == "":
    proxies = {'http': options.httpProxy}

  opener = urllib.FancyURLopener(proxies=proxies)
  opener.retrieve(url, filename)

  print "Detecting ..."

  lines = []
  regex = re.compile(
    "href.*cmake-[0-9.]+(-rc[0-9]+)*(-g[0-9a-fA-F]+)*-(Darwin-universal.tar.gz|Linux-i386.tar.gz|win32-x86.zip)"
  )

  f = open(filename)
  alllines = f.readlines()
  for line in alllines:
    if regex.search(line):
      lines.append(line)

  count = 0
  found = 0
  version_iterator = ""

  versionRegEx = re.compile(r'.*-([0-9.]+-rc[0-9]+|[0-9.]+-g[0-9a-fA-F]+|[0-9.]+)-.*')
  dateRegEx = re.compile(r'^[0-9].[0-9].([0-9.]+-rc[0-9]+|[0-9.]+-g[0-9a-fA-F]+|[0-9.]+)$')
  hrefRegEx = re.compile(r'^.*href="([^"]+)".*$')

  for line in lines:
    version = versionRegEx.match(line).group(1)

    if version == "" or version == line:
      print "error: line does not match version extraction regex"
      print " line: [" + line + "]"
      sys.exit(1)

    date = dateRegEx.match(version).group(1)

    # l, m, w == found an installer for Linux, Mac, Windows respectively
    # When encountering a new version, reset back to zeroes...
    #
    if(version_iterator != version):
      version_iterator = version
      l = 0
      m = 0
      w = 0
      lhref = ""
      mhref = ""
      whref = ""

    href = hrefRegEx.match(line).group(1)

    if re.search('Linux', line) != None:
      lhref = href
      l = 1
    elif re.search('Darwin', line) != None:
      mhref = href
      m = 1
    elif re.search('win32', line) != None:
      whref = href
      w = 1
    else:
      print "error: unexpected non-matching line"
      sys.exit(1)

    count = count + 1

    if l == 1 and m == 1 and w == 1:
      found = 1
      print "Detected latest available CMake " + vdir + " build: " + version
      break

  if not found:
    print "error: could not find a " + vdir + " version with all 3 platforms available"
    return ()

  return (('linux2', lhref, version), ('darwin', mhref, version), ('win32', whref, version))


def Download(basedir, url):
  cmps = url.rsplit("/", 1)
  href = cmps[1]
  filename = ''.join([basedir, "/", href])

  print 'Downloading ' + href + '...'

  try:
    createDir(basedir)
  except:
    if not os.path.exists(basedir):
      raise

  proxies = None # if None, use proxy from env var http_proxy
  if not options.httpProxy == "":
    proxies = {'http': options.httpProxy}

  opener = urllib.FancyURLopener(proxies=proxies)
  opener.retrieve(url, filename)


def Extract(basedir, url):
  cmps = url.rsplit("/", 1)
  href = cmps[1]
  filename = ''.join([basedir, "/", href])

  print 'Extracting ' + href + '...'

  if href[-4:] == ".zip":
    if sys.version < '2.6':
      if sys.platform == 'win32':
        print "error: cannot extract zip files on win32 with older python < 2.6"
      else:
        print "warning: avoiding zipfile.extractall on older python < 2.6"
        print "         skipping this extraction..."
    else:
      z = zipfile.ZipFile(filename)
      z.extractall(basedir)
      z.close()
  else:
    if sys.version < '2.6':
      if sys.platform == 'win32':
        print "error: cannot extract tar files on win32 with older python < 2.6"
      else:
        print "warning: avoiding tarfile.extractall on older python < 2.6"
        print "         trying command line tar instead..."
        origDir = os.getcwd()
        echoChDir(basedir)
        echoRunSysCmnd("tar -xzf " + href)
        echoChDir(origDir)
    else:
      t = tarfile.open(filename)
      t.extractall(basedir)
      t.close()


def Install(basedir, url):
  cmps = url.rsplit("/", 1)
  href = cmps[1]

  if href[-4:] == ".zip":
    href = href[:-4]
  elif href[-7:] == ".tar.gz":
    href = href[:-7]

  dirname = ''.join([basedir, "/", href])

  print 'Installing ' + href + '...'
  print '  src dir: [' + dirname + ']'
  print '  dst dir: [' + options.installDir + ']'

  if sys.platform == 'win32':
    if os.path.exists(options.installDir):
      print "error: --install-dir '" + options.installDir + "' already exists - remove it or rename it and try again -- or manually copy the source directory '" + dirname + "' to the final installation location..."
      sys.exit(1)

    shutil.copytree(dirname, options.installDir)
  else:
    # avoid the "copytree doesn't work if dir already exists" problem by using
    # the sys command "cp"
    try:
      createDir(options.installDir)
    except:
      if not os.path.exists(options.installDir):
        raise

    echoRunSysCmnd("cp -r " + dirname + "/* " + options.installDir)

  # After installing, create symlinks if requested.
  #   (This chunk needs work if ever it must work on Windows. For now, it's a
  #    Linux/Mac only chunk of python... Uses 'ln', no '.exe' suffix...)
  #
  if options.symlinks and options.symlinksDir != '':
    cmake_file = ""
    cmake_install_topdir = ""
    pre = ""

    cmake_files = glob.glob(dirname + "/*")
    if len(cmake_files) >= 1:
      cmake_file = cmake_files[0]

    cmps = cmake_file.rsplit("/", 1)
    if len(cmps) >= 1:
      cmake_install_topdir = cmps[1]

    if os.path.exists(options.installDir + "/" + cmake_install_topdir + "/Contents/bin/cmake"):
      pre = cmake_install_topdir + "/Contents/bin"
    elif os.path.exists(options.installDir + "/bin/cmake"):
      pre = "bin"

    if pre == '':
      print "error: could not determine CMake install tree structure - cannot create symlinks into unexpected directory structure"
      sys.exit(1)

    if not os.path.exists(options.symlinksDir):
      echoRunSysCmnd("mkdir -p \"" + options.symlinksDir + "\"")

    for exe in ('ccmake', 'cmake', 'cmake-gui', 'cmakexbuild', 'cpack', 'ctest'):
      if os.path.exists(options.installDir + "/" + pre + "/" + exe):
        print "Creating " + exe + " symlink..."
        echoRunSysCmnd("ln -fs \"" + options.installDir + "/" + pre + "/" + exe + "\" \"" + options.symlinksDir + "/" + exe + "\"")


def DownloadForPlatform(p):
  if options.allPlatforms:
    return True
  if p == sys.platform:
    return True
  return False


def PrintDetectedDownloads(detected):
  print ""
  print "Detected CMake downloads available:"

  sorted_keys = detected.keys()
  sorted_keys.sort()

  detected_urls = list()

  for k in sorted_keys:
    for v in detected[k]:
      if DownloadForPlatform(v[0]):
        detected_urls.append(cmake_baseurl + "/" + k + "/" + v[1])

  for u in detected_urls:
    print "[" + u + "]"


def PrintVerifiedDownloads():
  print ""
  print "Verified CMake downloads:"

  verified_urls = list()

  for v in cmake_min_binaries:
    if DownloadForPlatform(v[0]):
      verified_urls.append(v[1])

  for v in cmake_release_binaries:
    if DownloadForPlatform(v[0]):
      verified_urls.append(v[1])

  for v in cmake_rc_binaries:
    if DownloadForPlatform(v[0]):
      verified_urls.append(v[1])

  for v in cmake_dev_binaries:
    if DownloadForPlatform(v[0]):
      verified_urls.append(v[1])

  for u in verified_urls:
    print "[" + u + "]"


# Read file "CMakeVersions.py" from the same directory that this script lives
# in. Then write an auto_updated edition of CMakeVersions.py in download_dir
# based on the detected latest available build for each auto_update vdir entry.
#
# Each line that matches the "auto_update" regex will be updated using the
# latest data in the "detected" dictionary. Other lines will be written out
# verbatim.
#
def ReadWriteCMakeVersionsFile(download_dir, detected):
  rfname = os.path.dirname(os.path.abspath(__file__)) + "/CMakeVersions.py"

  fr = open(rfname)
  lines = fr.readlines()

  wfname = download_dir + "/CMakeVersions.py"
  fw = open(wfname, "w")

  regex = re.compile(
    "cmake_version_(.*) = \"(.*)\" # auto_update ([^ \t\n]*)"
  )

  for line in lines:
    if regex.search(line):
      #fw.write("# original: ")
      #fw.write(line)

      type = regex.match(line).group(1)
      vdir = regex.match(line).group(3)

      binaries = detected[vdir]
      binary_0 = binaries[0]
      version = binary_0[2]

      auto_update_line = "cmake_version_" + type + " = \"" + version + "\" # auto_update " + vdir + "\n"

      fw.write(auto_update_line)
    else:
      fw.write(line)

  print ""
  print "Wrote new '" + wfname + "' -- copy to '" + rfname + "' (if different) to use newly detected installers."


#
# The main script
#

print ""
print "**************************************************************************"
print "Script: download-cmake.py \\"

if options.allPlatforms:
  print "  --all-platforms \\"
print "  --http-proxy="+options.httpProxy+" \\"
print "  --install-dir="+options.installDir+" \\"
print "  --installer-type="+options.installerType+" \\"
if options.skipDetect:
  print "  --skip-detect \\"
if options.skipDownload:
  print "  --skip-download \\"
if options.skipExtract:
  print "  --skip-extract \\"
if options.skipInstall:
  print "  --skip-install \\"
if options.symlinks:
  print "  --symlinks \\"
if options.symlinksDir != '':
  print "  --symlinks-dir="+options.symlinksDir+" \\"

if not options.httpProxy and not default_http_proxy:
  print "\nWARNING: Could not detect default http proxy for '"+hostname+"'!"

download_dir = "download_area"

binaries = None
if options.installerType == 'min':
  binaries = cmake_min_binaries
if options.installerType == 'release':
  binaries = cmake_release_binaries
if options.installerType == 'rc':
  binaries = cmake_rc_binaries
if options.installerType == 'dev':
  binaries = cmake_dev_binaries
if binaries == None:
  print "error: unknown --installer-type: [" + options.installerType + "]"
  sys.exit(1)

print ""
print ""
print "A) Detect the latest available builds of CMake ..."
print "    (requires network access to www.cmake.org)"
print ""

if options.skipDetect:
  print "Skipping on request ..."
else:
  detected = dict()

  for vdir in cmake_vdirs:
    detected[vdir] = DetectLatestCMakeBuilds(download_dir, cmake_baseurl, vdir)

  PrintDetectedDownloads(detected)

  PrintVerifiedDownloads()

  ReadWriteCMakeVersionsFile(download_dir, detected)


print ""
print ""
print "B) Download CMake for --installer-type '" + options.installerType + "' ..."
print "    (requires network access to www.cmake.org)"
print ""

if options.skipDownload:
  print "Skipping on request ..."
else:
  for binary in binaries:
    if DownloadForPlatform(binary[0]):
      Download(download_dir, binary[1])


print ""
print ""
print "C) Extract the CMake install tree ..."
print ""

if options.skipExtract:
  print "Skipping on request ..."
else:
  for binary in binaries:
    if DownloadForPlatform(binary[0]):
      Extract(download_dir, binary[1])


print ""
print ""
print "D) Install (copy the CMake install tree) ..."
print ""

if options.skipInstall:
  print "Skipping on request ..."
else:
  for binary in binaries:
    if binary[0] == sys.platform:
      Install(download_dir, binary[1])
