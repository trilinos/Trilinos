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

  def getDefaultCheckoutCommand(self):
    return "git clone software.sandia.gov:/space/git/TrilinosToolset/Dummy "\
      +self.getBaseDirName()+"/Dummy"

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
