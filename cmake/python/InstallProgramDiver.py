

from GeneralScriptSupport import *
from optparse import OptionParser


class InstallProgramDriver:


  def __init__(self, installObj):
    self.installObj = installObj


  def runDriver(self):

    productName = self.installObj.getProductName()
    baseDirName = self.installObj.getBaseDirName()

    scriptName = getScriptName();


    #
    # 1) Set up the help text
    #
      
    usageHelp = scriptName+\
""" [OPTIONS] [--install-dir=<SOMEDIR> ...]

Tool that checks out, untars, configures, builds, and installs """+productName+""" in one shot!

By default, if you just type:

   $ SOME_DIR/"""+scriptName+""" --do-all

then the directory """+baseDirName+""" will get created in the local working directory
and it will contain a tarball for """+productName+""" and the build files. NOTE: This
requires that you run this as root (or with an account that has root
privileges).  For not running as root, you need to specify --install-dir.

You can control various parts of the process with the options (see below).

The one option that you may need to change if you do not have root privileges
is the --install-dir option which is set to /usr/local/bin by default.  For
example, you might just type:

  $ SOME_DIR/"""+scriptName+""" --install-dir=$HOME --do-all

and then it would install "++" and the other executables in $HOME/bin.
NOTE: You will have to update your PATH variable to include whatever directory
you choose to install """+productName+""" in.

NOTE: If you need to use sudo to install in /usr/local/bin or some other place
that needs root privileges, do:

  $ SOME_DIR/"""+scriptName+""" --install-dir=$HOME --checkout --untar --configure --build
  $ sudo SOME_DIR/"""+scriptName+""" --install-dir=$HOME --install

This appears to work on some systems.

After you have done a successful install, you might want to do:

  $ rm -r """+baseDirName+"""

in order to remove the source and build files.

Good luck with """+productName+"""!

"""

    #
    # 2) Parse the command-line
    #

    clp = OptionParser(usage=usageHelp)
    
    clp.add_option(
      "--checkout", dest="checkout", action="store_true", default=False,
      help="Do the checkout of the tarball" )
    
    clp.add_option(
      "--checkout-command", dest="checkoutCommand", type="string",
      default=self.installObj.getDefaultCheckoutCommand(),
      help="Command used to check out "+productName+" tarball(s)." )
    
    clp.add_option(
      "--untar", dest="untar", action="store_true", default=False,
      help="Do the untar of the "+productName+" sources" )
    
    clp.add_option(
      "--configure", dest="configure", action="store_true", default=False,
      help="Configure "+productName+" to build" )
    
    clp.add_option(
      "--install-dir", dest="installDir", type="string",
      default="/usr/local",
      help="The install directory for "+productName+" (default = /usr/local)." )
    
    clp.add_option(
      "--build", dest="build", action="store_true", default=False,
      help="Build "+productName+" and related executables" )
    
    clp.add_option(
      "--make-options", dest="makeOptions", type="string",
      default="",
      help="The options to pass to make for "+productName+"." )
    
    clp.add_option(
      "--install", dest="install", action="store_true", default=False,
      help="Install "+productName )
    
    clp.add_option(
      "--show-final", dest="showFinal", action="store_true", default=False,
      help="Show final instructions for using "+productName )
    
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
      options.showFinal = True
      

    #
    # 3) Echo the command-line options
    #


    
    #
    # 4) Execute the commands
    #
    
    baseDir = os.getcwd()
    
    productBaseDir = baseDir+"/"+baseDirName

    self.installObj.setup(options)

    print ""
    print "A) Checkout the tarball(s) for "+productName+" ..."
    print ""
    
    if options.checkout:
      self.installObj.doCheckout(options, baseDir, productBaseDir)
    else:
      print "Skipping on request ..."
    
    print ""
    print "B) Untar the tarball(s) and set up ready to configure ..."
    print ""
    
    if options.untar:
      self.installObj.doUntar(options, baseDir, productBaseDir)
    else:
      print "Skipping on request ..."
    
    
    print ""
    print "C) Configure "+ProductName+" ..."
    print ""
    
    
    if options.configure:
      self.installObj.doConfigure(options, baseDir, productBaseDir)
    else:
      print "Skipping on request ..."
    
    
    print ""
    print "D) Build "+ProductName+" ..."
    print ""
    
    if options.build:
      self.installObj.doBuild(options, baseDir, productBaseDir)
    else:
      print "Skipping on request ..."
    
    
    print ""
    print "E) Install "+ProductName+" ..."
    print ""
    
    if options.install:
      self.installObj.doInstall(options, baseDir, productBaseDir)
    else:
      print "Skipping on request ..."
    
    print "[End]"
