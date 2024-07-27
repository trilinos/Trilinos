# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
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

  def test_packages_something(self):
    self.check('packages/something', False)

  def test_packages_framework(self):
    self.check('packages/framework/something', True)

  def test_dotgithub_workflows(self):
    self.check('.github/workflows/something', True)


if __name__ == '__main__':
  unittest.main()
