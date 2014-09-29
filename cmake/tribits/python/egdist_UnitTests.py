# @HEADER
# ************************************************************************
#
#            TriBITS: Tribal Build, Integrate, and Test System
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

################################
# Unit testing code for egdist #
################################

from GeneralScriptSupport import *
import unittest
import sys
import imp


utilsDir = getScriptBaseDir()+"/utils"
tribitsDir = os.path.abspath(getScriptBaseDir()+"/..")
commonToolsGitDir = tribitsDir+"/common_tools/git"
#print "commonToolsGitDir = ", commonToolsGitDir

sys.path.append(commonToolsGitDir)
from egdist import *


#
# Test data
#


repoVersionFile_withSummary_1 = """*** Base Git Repo: MockTrilinos
sha1_1 [Mon Sep 23 11:34:59 2013 -0400] <author_1@ornl.gov>
First summary message
*** Git Repo: extraTrilinosRepo
sha1_2 [Fri Aug 30 09:55:07 2013 -0400] <author_2@ornl.gov>
Second summary message
*** Git Repo: extraRepoOnePackage
sha1_3 [Thu Dec 1 23:34:06 2011 -0500] <author_3@ornl.gov>
Third summary message
"""

repoVersionFile_withoutSummary_1 = """*** Base Git Repo: MockTrilinos
sha1_1 [Mon Sep 23 11:34:59 2013 -0400] <author_1@ornl.gov>
*** Git Repo: extraRepoTwoPackages
sha1_2 [Fri Aug 30 09:55:07 2013 -0400] <author_2@ornl.gov>
*** Git Repo: extraRepoOnePackageThreeSubpackages
sha1_3 [Thu Dec 1 23:34:06 2011 -0500] <author_3@ornl.gov>
"""


#
# Unit tests for functions in egdist
#


class test_egdist_getRepoVersionDictFromRepoVersionFileString(unittest.TestCase):


  def setUp(self):
    None


  def test_repoVersionFile_withSummary_1(self):
    repoVersionDict = \
      getRepoVersionDictFromRepoVersionFileString(repoVersionFile_withSummary_1)
    expectedDict = {
      'MockTrilinos': 'sha1_1',
      'extraTrilinosRepo': 'sha1_2',
      'extraRepoOnePackage': 'sha1_3'
      }
    #print "repoVersionDict =\n", repoVersionDict
    self.assertEqual(repoVersionDict, expectedDict)


  def test_repoVersionFile_withoutSummary_1(self):
    repoVersionDict = \
      getRepoVersionDictFromRepoVersionFileString(repoVersionFile_withoutSummary_1)
    expectedDict = {
      'MockTrilinos': 'sha1_1',
      'extraRepoTwoPackages': 'sha1_2',
      'extraRepoOnePackageThreeSubpackages': 'sha1_3'
      }
    #print "repoVersionDict =\n", repoVersionDict
    self.assertEqual(repoVersionDict, expectedDict)


# ToDo: Add unit tests for requoteCmndLineArgsIntoArray!


#
# Test entire script egdist#


egdistPath = commonToolsGitDir+"/egdist"
egdistPathNoColor = egdistPath+" --dist-no-color"
egdistPathMock = egdistPathNoColor+" --with-eg-git=mockeg --dist-no-opt"

mockProjectDir = tribitsDir+"/package_arch/UnitTests/MockTrilinos"
unitTestDataDir = tribitsDir+"/python/UnitTests"

testBaseDir = os.getcwd()



def getCmndOutputInMockProjectDir(cmnd):
  os.chdir(mockProjectDir)
  cmndOut = getCmndOutput(cmnd)
  os.chdir(testBaseDir)
  return cmndOut


class test_egdist(unittest.TestCase):


  def setUp(self):
    None


  def test_default(self):
    cmndOut = getCmndOutput(egdistPathNoColor)
    cmndOut_expected = "Must specify eg/git command. See 'eg/git --help' for options."
    self.assertEqual(cmndOut, cmndOut_expected)


  def test_help(self):
    cmndOut = getCmndOutput(egdistPath+" --help")
    cmndOutList = cmndOut.split("\n")
    cmndOutFirstLine = cmndOutList[0] 
    cmndOutFirstLineAfterComma = cmndOutFirstLine.split(":")[1].strip() 
    cmndOutFirstLineAfterComma_expected = "egdist [egdist options] [OPTIONS]"
    self.assertEqual(cmndOutFirstLineAfterComma, cmndOutFirstLineAfterComma_expected)


  def test_noEgGit(self):
    cmndOut = getCmndOutput(egdistPathNoColor+" --with-eg-git= log")
    cmndOut_expected = "Can't find eg or git, please set --with-eg-dist"
    self.assertEqual(cmndOut, cmndOut_expected)


  def test_log_args(self):
    cmndOut = getCmndOutputInMockProjectDir(egdistPathMock+" log HEAD -1")
    cmndOut_expected = \
      "\n*** Base Git Repo: MockTrilinos\n" \
      "['mockeg', 'log', 'HEAD', '-1']\n"
    self.assertEqual(cmndOut, cmndOut_expected)


  def test_log_args_extra_repo_1(self):
    cmndOut = getCmndOutputInMockProjectDir(
      egdistPathMock+" --dist-extra-repos=extraTrilinosRepo log HEAD -1")
    cmndOut_expected = \
      "\n*** Base Git Repo: MockTrilinos\n" \
      "['mockeg', 'log', 'HEAD', '-1']\n\n" \
      "*** Git Repo: extraTrilinosRepo\n" \
      "['mockeg', 'log', 'HEAD', '-1']\n"
    self.assertEqual(cmndOut, cmndOut_expected)


  def test_log_args_extra_repo_2_not_first(self):
    cmndOut = getCmndOutputInMockProjectDir(
      egdistPathMock+\
        " --dist-extra-repos=extraTrilinosRepo,extraRepoOnePackage "+\
        " --dist-not-extra-repos=extraTrilinosRepo "+\
        " log HEAD -1"
      )
    cmndOut_expected = \
      "\n*** Base Git Repo: MockTrilinos\n" \
      "['mockeg', 'log', 'HEAD', '-1']\n\n" \
      "*** Git Repo: extraRepoOnePackage\n" \
      "['mockeg', 'log', 'HEAD', '-1']\n"
    self.assertEqual(cmndOut, cmndOut_expected)


  def test_log_args_extra_repo_2_not_second(self):
    cmndOut = getCmndOutputInMockProjectDir(
      egdistPathMock+\
        " --dist-extra-repos=extraTrilinosRepo,extraRepoOnePackage "+\
        " --dist-not-extra-repos=extraTrilinosRepo "+\
        " log HEAD -1"
      )
    cmndOut_expected = \
      "\n*** Base Git Repo: MockTrilinos\n" \
      "['mockeg', 'log', 'HEAD', '-1']\n\n" \
      "*** Git Repo: extraRepoOnePackage\n" \
      "['mockeg', 'log', 'HEAD', '-1']\n"
    self.assertEqual(cmndOut, cmndOut_expected)


  def test_log_args_extra_repo_1_not_base(self):
    cmndOut = getCmndOutputInMockProjectDir(
      egdistPathMock+\
        " --dist-extra-repos=extraTrilinosRepo "+\
        " --dist-not-base-repo "+\
        " log HEAD -1"
      )
    cmndOut_expected = \
      "\n*** Git Repo: extraTrilinosRepo\n" \
      "['mockeg', 'log', 'HEAD', '-1']\n"
    self.assertEqual(cmndOut, cmndOut_expected)


  def test_log_version_file(self):
    cmndOut = getCmndOutputInMockProjectDir(
      egdistPathMock+ \
      " --dist-version-file="+unitTestDataDir+"/versionFile_withSummary_1.txt"+\
      " log _VERSION_ --some -other args")
    cmndOut_expected = \
      "\n*** Base Git Repo: MockTrilinos\n" \
      "['mockeg', 'log', 'sha1_1', '--some', '-other', 'args']\n"
    self.assertEqual(cmndOut, cmndOut_expected)


  def test_log_version_file_extra_repo_1(self):
    cmndOut = getCmndOutputInMockProjectDir(
      egdistPathMock+ \
      " --dist-version-file="+unitTestDataDir+"/versionFile_withSummary_1.txt"+ \
      " --dist-extra-repos=extraTrilinosRepo"+ \
      " log _VERSION_")
    cmndOut_expected = \
      "\n*** Base Git Repo: MockTrilinos\n" \
      "['mockeg', 'log', 'sha1_1']\n" \
      "\n*** Git Repo: extraTrilinosRepo\n['mockeg', 'log', 'sha1_2']\n"
    self.assertEqual(cmndOut, cmndOut_expected)


  def test_log_version_file_extra_repo_2(self):
    cmndOut = getCmndOutputInMockProjectDir(
      egdistPathMock+ \
      " --dist-version-file="+unitTestDataDir+"/versionFile_withSummary_1.txt"+ \
      " --dist-extra-repos=extraRepoOnePackage,extraTrilinosRepo"+ \
      " log _VERSION_")
    cmndOut_expected = \
      "\n*** Base Git Repo: MockTrilinos\n" \
      "['mockeg', 'log', 'sha1_1']\n" \
      "\n*** Git Repo: extraRepoOnePackage\n['mockeg', 'log', 'sha1_3']\n" \
      "\n*** Git Repo: extraTrilinosRepo\n['mockeg', 'log', 'sha1_2']\n"
    self.assertEqual(cmndOut, cmndOut_expected)


  def test_log_HEAD_version_file_extra_repo_1(self):
    cmndOut = getCmndOutputInMockProjectDir(
      egdistPathMock+ \
      " --dist-version-file="+unitTestDataDir+"/versionFile_withSummary_1.txt"+ \
      " --dist-extra-repos=extraTrilinosRepo"+ \
      " log HEAD ^_VERSION_")
    cmndOut_expected = \
      "\n*** Base Git Repo: MockTrilinos\n" \
      "['mockeg', 'log', 'HEAD', '^sha1_1']\n" \
      "\n*** Git Repo: extraTrilinosRepo\n['mockeg', 'log', 'HEAD', '^sha1_2']\n"
    self.assertEqual(cmndOut, cmndOut_expected)


  def test_version_file_invalid_extra_repo(self):
    cmndOut = getCmndOutputInMockProjectDir(
      egdistPathMock+ \
      " --dist-version-file="+unitTestDataDir+"/versionFile_withSummary_1.txt"+ \
      " --dist-extra-repos=extraRepoTwoPackages"+ \
      " log _VERSION_")
    cmndOut_expected = \
      "\n*** Base Git Repo: MockTrilinos\n['mockeg', 'log', 'sha1_1']\n" \
      "\n*** Git Repo: extraRepoTwoPackages\nExtra repo 'extraRepoTwoPackages' is not in the list of extra repos ['extraTrilinosRepo', 'extraRepoOnePackage'] read in from version file."
    self.assertEqual(cmndOut, cmndOut_expected)


  def test_log_not_version_file_2(self):
    cmndOut = getCmndOutputInMockProjectDir(
      egdistPathMock+ \
      " --dist-version-file="+unitTestDataDir+"/versionFile_withSummary_1.txt"+ \
      " --dist-version-file2="+unitTestDataDir+"/versionFile_withSummary_1_2.txt"+ \
      " log _VERSION_ ^_VERSION2_")
    cmndOut_expected = \
      "\n*** Base Git Repo: MockTrilinos\n" \
      "['mockeg', 'log', 'sha1_1', '^sha1_1_2']\n"
    self.assertEqual(cmndOut, cmndOut_expected)


  def test_log_not_version_file_2_extra_repo_1(self):
    cmndOut = getCmndOutputInMockProjectDir(
      egdistPathMock+ \
      " --dist-version-file="+unitTestDataDir+"/versionFile_withSummary_1.txt"+ \
      " --dist-version-file2="+unitTestDataDir+"/versionFile_withSummary_1_2.txt"+ \
      " --dist-extra-repos=extraTrilinosRepo"+ \
      " log _VERSION_ ^_VERSION2_")
    cmndOut_expected = \
      "\n*** Base Git Repo: MockTrilinos\n" \
      "['mockeg', 'log', 'sha1_1', '^sha1_1_2']\n" \
      "\n*** Git Repo: extraTrilinosRepo\n['mockeg', 'log', 'sha1_2', '^sha1_2_2']\n"
    self.assertEqual(cmndOut, cmndOut_expected)


  def test_log_since_until_version_file_2_extra_repo_1(self):
    cmndOut = getCmndOutputInMockProjectDir(
      egdistPathMock+ \
      " --dist-version-file="+unitTestDataDir+"/versionFile_withSummary_1.txt"+ \
      " --dist-version-file2="+unitTestDataDir+"/versionFile_withSummary_1_2.txt"+ \
      " --dist-extra-repos=extraTrilinosRepo"+ \
      " log _VERSION2_.._VERSION_")
    cmndOut_expected = \
      "\n*** Base Git Repo: MockTrilinos\n" \
      "['mockeg', 'log', 'sha1_1_2..sha1_1']\n" \
      "\n*** Git Repo: extraTrilinosRepo\n['mockeg', 'log', 'sha1_2_2..sha1_2']\n"
    self.assertEqual(cmndOut, cmndOut_expected)
  # The above test ensures that it repalces the SHA1s for in the same cmndline args


if __name__ == '__main__':
  unittest.main()
