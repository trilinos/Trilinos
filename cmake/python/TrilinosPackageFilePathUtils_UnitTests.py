#########################################################
# Unit testing code for TrilinosPackageFilePathUtils.py #
######################################################### 


from TrilinosPackageFilePathUtils import *
import unittest


#
# Test isGlobalBuildFileRequiringGlobalRebuild
#

class test_isGlobalBuildFileRequiringGlobalRebuild(unittest.TestCase):


  def test_Trilinos_version_h(self):
    self.assertEqual( isGlobalBuildFileRequiringGlobalRebuild( 'Trilinos_version.h' ), True )


  def test_CMakeLists_txt(self):
    self.assertEqual( isGlobalBuildFileRequiringGlobalRebuild( 'CMakeLists.txt' ), True )


  def test_TrilinosPackages_cmake(self):
    self.assertEqual( isGlobalBuildFileRequiringGlobalRebuild( 'cmake/TrilinosPackages.cmake' ), False )


  def test_TrilinosTPLs_cmake(self):
    self.assertEqual( isGlobalBuildFileRequiringGlobalRebuild( 'cmake/TrilinosTPLs.cmake' ), False )


  def test_TrilinosCMakeQuickstart_txt(self):
    self.assertEqual( isGlobalBuildFileRequiringGlobalRebuild( 'cmake/TrilinosCMakeQuickstart.txt' ), False )


  def test_experimental_build_test_cmake(self):
    self.assertEqual( isGlobalBuildFileRequiringGlobalRebuild( 'cmake/ctest/experimental_build_test.cmake' ),
      False )


  def test_DependencyUnitTests(self):
    self.assertEqual( isGlobalBuildFileRequiringGlobalRebuild( 'cmake/DependencyUnitTests/blah' ),
      False )


  def test_FindTPLBLAS(self):
    self.assertEqual( isGlobalBuildFileRequiringGlobalRebuild( 'cmake/TPLs/FindTPLBLAS.cmake' ),
      True )


  def test_FindTPLLAPACK(self):
    self.assertEqual( isGlobalBuildFileRequiringGlobalRebuild( 'cmake/TPLs/FindTPLLAPACK.cmake' ),
      True )


  def test_FindTPLMPI(self):
    self.assertEqual( isGlobalBuildFileRequiringGlobalRebuild( 'cmake/TPLs/FindTPLMPI.cmake' ),
      True )


  def test_FindTPLDummy(self):
    self.assertEqual( isGlobalBuildFileRequiringGlobalRebuild( 'cmake/TPLs/FindTPLDummy.cmake' ),
      False )


  def test_SetNotFound(self):
    self.assertEqual( isGlobalBuildFileRequiringGlobalRebuild( 'cmake/utils/SetNotFound.cmake' ),
      True )


testingTrilinosDepsXmlInFile = getScriptBaseDir()+"/UnitTests/TrilinosPackageDependencies.gold.xml"
trilinosDependencies = getTrilinosDependenciesFromXmlFile(testingTrilinosDepsXmlInFile)
  
#print "\ntrilinosDependencies:\n", trilinosDependencies


updateOutputStr = """
? packages/triutils/doc/html
M CMakeLists.txt
M cmake/python/checkin-test.py
M cmake/python/dump-cdash-deps-xml-file.py
A packages/nox/src/dummy.C
P packages/stratimikos/dummy.blah
M packages/thyra/src/Thyra_ConfigDefs.hpp
M packages/thyra/CMakeLists.txt
M packages/ifpack2/CMakeLists.txt
M demos/FEApp/src/CMakeLists.txt
"""

updateOutputList = updateOutputStr.split("\n")


class testTrilinosPackageFilePathUtils(unittest.TestCase):


  def test_getPackageNameFromPath_01(self):
    self.assertEqual(
      getPackageNameFromPath( trilinosDependencies, 'packages/teuchos/CMakeLists.txt' ),
      'Teuchos' )


  def test_getPackageNameFromPath_02(self):
    self.assertEqual(
      getPackageNameFromPath( trilinosDependencies, 'packages/thyra/src/blob.cpp' ),
      'ThyraCoreLibs' )


  def test_getPackageNameFromPath_03(self):
    self.assertEqual(
      getPackageNameFromPath( trilinosDependencies, 'cmake/CMakeLists.txt' ),
      'TrilinosFramework' )


  def test_getPackageNameFromPath_04(self):
    self.assertEqual(
      getPackageNameFromPath( trilinosDependencies, 'cmake/CMakeLists.txt' ),
      'TrilinosFramework' )


  def test_getPackageNameFromPath_noMatch(self):
    self.assertEqual(
      getPackageNameFromPath( trilinosDependencies, 'packages/blob/blob' ), '' )


  def test_extractFilesListMatchingPattern_01(self):

    modifedFilesList = extractFilesListMatchingPattern( updateOutputList,
      re.compile(r"^[MA] (.+)$") )

    modifedFilesList_expected = \
      [
        "CMakeLists.txt",
        "cmake/python/checkin-test.py",
        "cmake/python/dump-cdash-deps-xml-file.py",
        "packages/nox/src/dummy.C",
        "packages/thyra/src/Thyra_ConfigDefs.hpp",
        "packages/thyra/CMakeLists.txt",
        "packages/ifpack2/CMakeLists.txt",
        "demos/FEApp/src/CMakeLists.txt",
      ]

    self.assertEqual( modifedFilesList, modifedFilesList_expected )


  def test_getPackagesListFromFilePathsList_01(self):

    filesList = extractFilesListMatchingPattern( updateOutputList,
      re.compile(r"^[AMP] (.+)$") )
    
    packagesList = getPackagesListFromFilePathsList( trilinosDependencies, filesList )

    packagesList_expected = \
      [u"TrilinosFramework", u"Stratimikos", u"ThyraCoreLibs", u"Thyra"]

    self.assertEqual( packagesList, packagesList_expected )


  def test_getPackageCheckinEmailAddressesListFromFilePathsList_01(self):

    filesList = extractFilesListMatchingPattern( updateOutputList,
      re.compile(r"^[AMP] (.+)$") )
    
    packagesList = getPackageCheckinEmailAddressesListFromFilePathsList(
      trilinosDependencies, filesList )

    packagesList_expected = [
      u"trilinos-checkins@software.sandia.gov",
      u"trilinosframework-checkins@software.sandia.gov",
      u"stratimikos-checkins@software.sandia.gov",
      u"thyra-checkins@software.sandia.gov",
      ]

    self.assertEqual( packagesList, packagesList_expected )


  def test_get_trilinos_packages_from_files_list_01(self):

    writeStrToFile( "modifiedFiles.txt",
      "CMakeLists.txt\n" \
      "cmake/python/checkin-test.py\n" \
      "cmake/python/dump-cdash-deps-xml-file.py\n" \
      "packages/thyra/src/Thyra_ConfigDefs.hpp\n" \
      "packages/thyra/CMakeLists.txt\n" \
      )

    self.assertEqual(
      getCmndOutput(getScriptBaseDir()+"/get-trilinos-packages-from-files-list.py" \
        " --files-list-file=modifiedFiles.txt", True),
      "ALL_PACKAGES;TrilinosFramework;Thyra"
      )

    self.assertEqual(
      getCmndOutput(getScriptBaseDir()+"/get-trilinos-packages-from-files-list.py" \
        " --files-list-file=modifiedFiles.txt --deps-xml-file="+testingTrilinosDepsXmlInFile,
        True),
      "ALL_PACKAGES;TrilinosFramework;ThyraCoreLibs;Thyra"
      )


class testFilterPackagesList(unittest.TestCase):


  def test_get_PS(self):
    self.assertEqual(
      getCmndOutput(getScriptBaseDir()+"/filter-packages-list.py" \
        " --deps-xml-file="+testingTrilinosDepsXmlInFile+"" \
        " --input-packages-list=Teuchos,Thyra,Phalanx,Stokhos --keep-types=PS",
        True),
      "Teuchos,Thyra"
      )


  def test_get_PS_SS(self):
    self.assertEqual(
      getCmndOutput(getScriptBaseDir()+"/filter-packages-list.py" \
        " --deps-xml-file="+testingTrilinosDepsXmlInFile+"" \
        " --input-packages-list=Teuchos,Thyra,Phalanx,Stokhos --keep-types=PS,SS",
        True),
      "Teuchos,Thyra,Phalanx"
      )


  def test_get_PS_SS_EX(self):
    self.assertEqual(
      getCmndOutput(getScriptBaseDir()+"/filter-packages-list.py" \
        " --deps-xml-file="+testingTrilinosDepsXmlInFile+"" \
        " --input-packages-list=Teuchos,Thyra,Phalanx,Stokhos --keep-types=PS,SS,EX",
        True),
      "Teuchos,Thyra,Phalanx,Stokhos"
      )


  def test_get_SS(self):
    self.assertEqual(
      getCmndOutput(getScriptBaseDir()+"/filter-packages-list.py" \
        " --deps-xml-file="+testingTrilinosDepsXmlInFile+"" \
        " --input-packages-list=Teuchos,Thyra,Phalanx,Stokhos --keep-types=SS",
        True),
      "Phalanx"
      )


  def test_get_PS_EX(self):
    self.assertEqual(
      getCmndOutput(getScriptBaseDir()+"/filter-packages-list.py" \
        " --deps-xml-file="+testingTrilinosDepsXmlInFile+"" \
        " --input-packages-list=Teuchos,Thyra,Phalanx,Stokhos --keep-types=PS,EX",
        True),
      "Teuchos,Thyra,Stokhos"
      )


def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(testTrilinosPackageFilePathUtils))
    return suite


if __name__ == '__main__':
  unittest.main()
