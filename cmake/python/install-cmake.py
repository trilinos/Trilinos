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

usageHelp = r"""install-cmake.py [--install-dir=SOMEDIR ...]

Tool that checks out, configures, builds, and installs cmake in one shot!

By default, if you just type:

   $ SOME_DIR/install-cmake.py

then the directory CMake.base will get created in the local working directory
and it will contain the checked out sources for CMake and the build
files. NOTE: This requires that you run this as root (or with an account that
has root privileges).

You can control various parts of the process with the options (see below).

The one option that you may need to change if you do not have root privileges
is the --install-dir option which is set to /usr/local/bin by default.  For
example, you might just type:

  $ SOME_DIR/install-cmake.py --install-dir=$HOME

and then it would install cmake and the other executables in $HOME/bin.  NOTE:
You will have to update your PATH variable to include whatever directory you
choose to install CMake in.

NOTE: If you need to use sudo to install in /usr/local/bin or some other place
that needs root privileges, do:

  $ install-cmake.py --skip-install [other options]
  $ cd CMake.base/build
  $ sudo make install

After you have done a successful install, you might want to do:

  $ rm -r CMake.base

in order to remove the source and build files.

Good luck with CMake!

"""

from optparse import OptionParser

clp = OptionParser(usage=usageHelp)

clp.add_option(
  "--skip-checkout", dest="skipCheckout", action="store_true", default=False,
  help="Skip the checkout" )

clp.add_option(
  "--cvs-command", dest="cvsCommand", type="string",
  default="cvs -d :ext:software.sandia.gov:/space/CVS co -d CMake.base Trilinos3PL/CMakeDev/CMake.tar.gz",
  help="Command used to check out CMake tarball." )

clp.add_option(
  "--skip-untar", dest="skipUntar", action="store_true", default=False,
  help="Skip the untar step" )

clp.add_option(
  "--skip-configure", dest="skipConfigure", action="store_true", default=False,
  help="Skip the configure step" )

clp.add_option(
  "--install-dir", dest="installDir", type="string",
  default="/usr/local",
  help="The install directory for CMake (default = /usr/local)." )

clp.add_option(
  "--skip-build", dest="skipBuild", action="store_true", default=False,
  help="Skip the build step" )

clp.add_option(
  "--make-options", dest="makeOptions", type="string",
  default="",
  help="The options to pass to make." )

clp.add_option(
  "--skip-install", dest="skipInstall", action="store_true", default=False,
  help="Skip the install step" )

(options, args) = clp.parse_args()

#
# Execute the commands
#

print ""
print "A) Checkout the tarball for CMake ..."
print ""

if options.skipCheckout:
  print "Skipping on request ..."
else:
  echoRunSysCmnd(options.cvsCommand)

echoChDir("CMake.base")


print ""
print "B) Untar the tarball ..."
print ""

if options.skipUntar:
  print "Skipping on request ..."
else:
  echoRunSysCmnd("tar -xzvf CMake.tar.gz")


print ""
print "C) Configure CMake ..."
print ""

createDir("build")
echoChDir("build")

if options.skipConfigure:
  print "Skipping on request ..."
else:
  echoRunSysCmnd("../CMake/configure --prefix="+options.installDir)


print ""
print "D) Build CMake ..."
print ""

if options.skipBuild:
  print "Skipping on request ..."
else:
  echoRunSysCmnd("make "+options.makeOptions)


print ""
print "E) Install CMake ..."
print ""

if options.skipInstall:
  print "Skipping on request ..."
else:
  echoRunSysCmnd("make install")
