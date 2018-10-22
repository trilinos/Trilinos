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

#########################################################
# Unit testing code for TribitsPackageFilePathUtils.py #
######################################################### 

import os
import sys

from ProjectCiFileChangeLogic import *

import unittest

class test_TribitsExampleProject_ProjectCiFileChangeLogic(unittest.TestCase):

  def check(self, filePath, expectedResult):
    dpcl = ProjectCiFileChangeLogic()
    self.assertEqual(
      dpcl.isGlobalBuildFileRequiringGlobalRebuild(filePath),
      expectedResult )

  def test_CMakeLists_txt(self):
    self.check('CMakeLists.txt', True)

  def test_PackagesList_cmake(self):
    self.check('PackagesList.cmake', True)
    
  def test_TPLsList_cmake(self):
    self.check('TPLsList.cmake', True)

  def test_Version_cmake(self):
    self.check('Version.cmake', True)

  def test_Anything_cmake(self):
    self.check('Anything.cmake', True)

  def test_cmake_TrilinosCMakeQuickstart_txt(self):
    self.check('cmake/TrilinosCMakeQuickstart.txt', False)

  def test_cmake_TPLsList(self):
    self.check('cmake/ExtraRepositoriesList.cmake', False)

  def test_cmake_std(self):
    self.check('cmake/std/anything', True)

  def test_cmake_std_atdm_anything(self):
    self.check('cmake/std/atdm/anything', False)

  def test_cmake_std_atdm_cmake(self):
    self.check('cmake/std/atdm/ATDMDevEnvSettings.cmake', False)

  def test_cmake_std_atdm_system_tweaks_cmake(self):
    self.check('cmake/std/atdm/shiller/tweaks/INTEL-RELEASE-OPENMP-HSW.cmake', False)

  def test_cmake_std_sems_anything(self):
    self.check('cmake/std/sems/anything', False)

  def test_cmake_ctest(self):
    self.check('cmake/ctest/anything', True)

  def test_cmake_ctest_drivers(self):
    self.check('cmake/ctest/drivers/anything', True)

  def test_cmake_ctest_drivers_machinedir(self):
    self.check('cmake/ctest/drivers/some_machine/ctest_driver.cmake', False)

  def test_cmake_somewhere_else(self):
    self.check('cmake/somewhere_else/anthing.cmake', True)

  def test_cmake_TPLs(self):
    self.check('cmake/TPLs/anything', True)

  def test_doc(self):
    self.check('doc/anything', False)


if __name__ == '__main__':
  unittest.main()
