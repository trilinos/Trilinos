#!/usr/bin/env python

#
# General scripting support
#
# NOTE: Included first to check the version of python!
#

from GeneralScriptSupport import *


#
# Read in the commandline arguments
#

usageHelp = r"""install-git.py [--install-dir=SOMEDIR ...]

Tool that checks out, configures, builds, and installs git/gitk/eg in one shot!

By default, if you just type:

   $ SOME_DIR/install-git.py --do-all

then the directory Git.base will get created in the local working directory
and it will contain a tarball for git and eg and the build files. NOTE: This
requires that you run this as root (or with an account that has root
privileges).  For not running as root, you need to specify --install-dir.

You can control various parts of the process with the options (see below).

The one option that you may need to change if you do not have root privileges
is the --install-dir option which is set to /usr/local/bin by default.  For
example, you might just type:

  $ SOME_DIR/install-git.py --install-dir=$HOME --do-all

and then it would install git/eg and the other executables in $HOME/bin.
NOTE: You will have to update your PATH variable to include whatever directory
you choose to install git/eg in.

NOTE: If you need to use sudo to install in /usr/local/bin or some other place
that needs root privileges, do:

  $ install-git.py --checkout --untar --configure --build [other options]
  $ cd Git.base/build
  $ sudo make install

After you have done a successful install, you might want to do:

  $ rm -r Git.base

in order to remove the source and build files.

Good luck with git/eg!

"""

gitVersion = "1.6.5.2"
egVersion = "66"

from optparse import OptionParser

clp = OptionParser(usage=usageHelp)

clp.add_option(
  "--checkout", dest="checkout", action="store_true", default=False,
  help="Do the checkout of the tarball" )

clp.add_option(
  "--cvs-command", dest="cvsCommand", type="string",
  default="cvs -d :ext:software.sandia.gov:/space/CVS co -d Git.base Trilinos3PL/Git.base",
  help="Command used to check out git tarball." )

clp.add_option(
  "--untar", dest="untar", action="store_true", default=False,
  help="Do the untar of the git sources" )

clp.add_option(
  "--configure", dest="configure", action="store_true", default=False,
  help="Configure git to build" )

clp.add_option(
  "--install-dir", dest="installDir", type="string",
  default="/usr/local",
  help="The install directory for git (default = /usr/local)." )

clp.add_option(
  "--build", dest="build", action="store_true", default=False,
  help="Build git and related executables" )

clp.add_option(
  "--make-options", dest="makeOptions", type="string",
  default="",
  help="The options to pass to make for git." )

clp.add_option(
  "--install", dest="install", action="store_true", default=False,
  help="Install the git and related executables along with eg (just a simple copy)." )

clp.add_option(
  "--do-all", dest="doAll", action="store_true", default=False,
  help="Same as: --checkout --untar --configure --build --install" )

(options, args) = clp.parse_args()

if options.doAll:
  options.checkout = True
  options.untar = True
  options.configure = True
  options.build = True
  options.install = True

#
# Execute the commands
#

baseDir = os.getcwd()

getBaseDir = baseDir+"/Git.base"
gitSrcBuildDir = getBaseDir+"/git-"+gitVersion

print ""
print "A) Checkout the tarball for Git ..."
print ""

if options.checkout:
  echoRunSysCmnd(options.cvsCommand)
else:
  print "Skipping on request ..."

print ""
print "B) Untar the tarball ..."
print ""

if options.untar:
  echoChDir(getBaseDir)
  echoRunSysCmnd("tar -xzvf git-"+gitVersion+".tar.gz")
else:
  print "Skipping on request ..."


print ""
print "C) Configure Git (only insource builds are supported) ..."
print ""


if options.configure:
  echoChDir(gitSrcBuildDir)
  echoRunSysCmnd("./configure --prefix="+options.installDir)
else:
  print "Skipping on request ..."


print ""
print "D) Build Git ..."
print ""

if options.build:
  echoChDir(gitSrcBuildDir)
  echoRunSysCmnd("make "+options.makeOptions)
else:
  print "Skipping on request ..."


print ""
print "E) Install Git and eg ..."
print ""

if options.install:
  echoChDir(gitSrcBuildDir)
  echoRunSysCmnd("make install")
  egInstallFile = options.installDir+"/bin/eg"
  echoRunSysCmnd("cp ../eg."+egVersion+".pl "+egInstallFile)
  echoRunSysCmnd("chmod a+rx "+egInstallFile)
else:
  print "Skipping on request ..."

print "[End]"
