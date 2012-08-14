# @HEADER
# ************************************************************************
#
#            Trilinos: An Object-Oriented Solver Framework
#                 Copyright (2001) Sandia Corporation
#
#
# Copyright (2001) Sandia Corporation. Under the terms of Contract
# DE-AC04-94AL85000, there is a non-exclusive license for use of this
# work by or on behalf of the U.S. Government.  Export of this program
# may require a license from the United States Government.
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
# NOTICE:  The United States Government is granted for itself and others
# acting on its behalf a paid-up, nonexclusive, irrevocable worldwide
# license in this data to reproduce, prepare derivative works, and
# perform publicly and display publicly.  Beginning five (5) years from
# July 25, 2001, the United States Government is granted for itself and
# others acting on its behalf a paid-up, nonexclusive, irrevocable
# worldwide license in this data to reproduce, prepare derivative works,
# distribute copies to the public, perform publicly and display
# publicly, and to permit others to do so.
#
# NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED STATES DEPARTMENT
# OF ENERGY, NOR SANDIA CORPORATION, NOR ANY OF THEIR EMPLOYEES, MAKES
# ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR
# RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY
# INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS
# THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.
#
# ************************************************************************
# @HEADER

#########################################################
# Unit testing code for TribitsPackageFilePathUtils.py #
######################################################### 


from TribitsPackageFilePathUtils import *
import unittest


#
# Test isGlobalBuildFileRequiringGlobalRebuild
#

class test_isGlobalBuildFileRequiringGlobalRebuild(unittest.TestCase):


  def test_CMakeLists_txt(self):
    self.assertEqual( isGlobalBuildFileRequiringGlobalRebuild( 'CMakeLists.txt' ), True )


  def test_PackagesList_cmake(self):
    self.assertEqual( isGlobalBuildFileRequiringGlobalRebuild( 'PackagesList.cmake' ), False )


  def test_TPLsList_cmake(self):
    self.assertEqual( isGlobalBuildFileRequiringGlobalRebuild( 'TPLsList.cmake' ), False )


  def test_Version_cmake(self):
    self.assertEqual( isGlobalBuildFileRequiringGlobalRebuild( 'Version.cmake' ), True )


  def test_Anything_cmake(self):
    self.assertEqual( isGlobalBuildFileRequiringGlobalRebuild( 'Anything.cmake' ), True )


  def test_TrilinosCMakeQuickstart_txt(self):
    self.assertEqual( isGlobalBuildFileRequiringGlobalRebuild( 'cmake/TrilinosCMakeQuickstart.txt' ),
      False )


  def test_TPLsList_cmake(self):
    self.assertEqual( isGlobalBuildFileRequiringGlobalRebuild( 'cmake/ExtraRepositoriesList.cmake' ),
      False )


  def test_experimental_build_test_cmake(self):
    self.assertEqual( isGlobalBuildFileRequiringGlobalRebuild( 'cmake/ctest/experimental_build_test.cmake' ),
      False )


  def test_cmake_ctest_drivers_something(self):
    self.assertEqual( isGlobalBuildFileRequiringGlobalRebuild( 'something/cmake/ctest/drivers/machine/somefile.cmake' ),
      False )


  def test_something_cmake_ctest_drivers_something(self):
    self.assertEqual( isGlobalBuildFileRequiringGlobalRebuild( 'cmake/ctest/drivers/machine/somefile.cmake' ),
      False )


  def test_cmake_UnitTests(self):
    self.assertEqual( isGlobalBuildFileRequiringGlobalRebuild( 'cmake/anything/UnitTests/CMakeLists.txt' ),
      False )


  def test_FindTPLBLAS_cmake(self):
    self.assertEqual( isGlobalBuildFileRequiringGlobalRebuild( 'cmake/TPLs/FindTPLBLAS.cmake' ),
      False )


  def test_FindTPLLAPACK_cmake(self):
    self.assertEqual( isGlobalBuildFileRequiringGlobalRebuild( 'cmake/TPLs/FindTPLLAPACK.cmake' ),
      False )


  def test_FindTPLMPI_cmake(self):
    self.assertEqual( isGlobalBuildFileRequiringGlobalRebuild( 'cmake/TPLs/FindTPLMPI.cmake' ),
      False )


  def test_FindTPLDummy_cmake(self):
    self.assertEqual( isGlobalBuildFileRequiringGlobalRebuild( 'cmake/TPLs/FindTPLDummy.cmake' ),
      False )


  def test_SetNotFound_cmake(self):
    self.assertEqual( isGlobalBuildFileRequiringGlobalRebuild( 'cmake/utils/SetNotFound.cmake' ),
      True )


testingTrilinosDepsXmlInFile = getScriptBaseDir()+"/UnitTests/TrilinosPackageDependencies.gold.xml"
trilinosDependencies = getProjectDependenciesFromXmlFile(testingTrilinosDepsXmlInFile)
  
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


class testProjectPackageFilePathUtils(unittest.TestCase):


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


  def test_get_trilinos_packages_from_files_list_01(self):

    writeStrToFile( "modifiedFiles.txt",
      "CMakeLists.txt\n" \
      "cmake/python/checkin-test.py\n" \
      "cmake/python/dump-cdash-deps-xml-file.py\n" \
      "packages/thyra/src/Thyra_ConfigDefs.hpp\n" \
      "packages/thyra/CMakeLists.txt\n" \
      )

    self.assertEqual(
      getCmndOutput(getScriptBaseDir()+"/get-tribits-packages-from-files-list.py" \
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
