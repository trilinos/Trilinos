#########################################################
# Unit testing code for TrilinosPackageFilePathUtils.py #
######################################################### 


from InstallProgramDriver import *
import unittest


class DummyInstall:

  def __init__(self):
    self.dummy = None

  def getProductName(self):
    return "Dummy"

  def getBaseDirName(self):
    return "Dummy.base"

  def getDefaultCheckoutCommand(self):
    return "git clone software.sandia.gov:/space/git/TrilinosToolset/Dummy "\
      +self.getBaseDirName()+"/Dummy"

  def setup(self, inOptions):
    self.options = inOptions


class test_InstallProgramDriver(unittest.TestCase):


  def test_help(self):
    ipd = InstallProgramDriver(DummyInstall())
    self.assertEqual(ipd.runDriver(), None)


#
# Main
#

if __name__ == '__main__':
  unittest.main()
