#########################################################
# Unit testing code for TrilinosPackageFilePathUtils.py #
######################################################### 


from InstallProgramDriver import *
import unittest


class DummyInstall:

  def __init__(self):
    self.dummy = None

  def getProductName(self):
    return "Dummy-0.1"

  def getScriptName(self):
    return "dummy-install.py"

  def getExtraHelpStr(self):
    return "\nThis is the extra help string!"

  def getBaseDirName(self):
    return "Dummy.base"

  def injectExtraCmndLineOptions(self, clp):
    clp.add_option(
      "--checkout-cmnd", dest="checkoutCmnd", type="string",
      default="git clone software.sandia.gov:/space/git/TrilinosToolset/Dummy "\
      +self.getBaseDirName()+"/Dummy",
      help="Command used to check out "+self.getProductName()+" tarball(s)." )

  def echoExtraCmndLineOptions(self, inOptions):
    cmndLine = ""
    cmndLine += "  --checkout-cmnd='" + inOptions.checkoutCmnd + "' \\\n"
    return cmndLine

  def setup(self, inOptions):
    self.options = inOptions


class test_InstallProgramDriver(unittest.TestCase):


  def test_default(self):
    ipd = InstallProgramDriver(DummyInstall())
    sys.argv = ["scriptName"]
    self.assertEqual(ipd.runDriver(), None)


#
# Main
#

if __name__ == '__main__':
  unittest.main()
