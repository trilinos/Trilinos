#!/usr/bin/env python
# -*- coding: utf-8 -*-

################################################################################
# Unit testing code for create_trilinos_github_test_failure_issue.py           #
################################################################################


import sys
import shutil
import unittest

import create_trilinos_github_test_failure_issue as CTGTFI


class test_stripGentConfigBuildName(unittest.TestCase):

  def test_1(self):
    genConfigBuildName = CTGTFI.stripGentConfigBuildName(
      "PR-10472-test-ats2_cuda-10.1.243-gnu-8.3.1-blah-blah-package-enables-911")
    genConfigBuildName_expected = \
      "ats2_cuda-10.1.243-gnu-8.3.1-blah-blah-package-enables"
    self.assertEqual(
      genConfigBuildName, genConfigBuildName_expected)


class test_getUniqueGenConfigBuildNamesList(unittest.TestCase):

  def test_1(self):
    genConfigBuildNamesList = CTGTFI.getUniqueGenConfigBuildNamesList(
      [
        "PR-10472-test-ats2_cuda-10.1.243-gnu-8.3.1-blah-blah-package-enables-911",
        "PR-10472-test-ats2_cuda-10.1.243-gnu-8.3.1-blah-blah-package-enables-912",
        "PR-10472-test-sems-rhel7_blah-blah-124",
        "PR-10473-test-ats2_cuda-10.1.243-gnu-8.3.1-blah-blah-package-enables-913",
        "PR-10472-test-sems-rhel7_blah-blah-125",
        "PR-10473-test-sems-rhel7_blah-blah-126",
        ]
       )
    genConfigBuildNamesList_expected = [
      "ats2_cuda-10.1.243-gnu-8.3.1-blah-blah-package-enables",
      "sems-rhel7_blah-blah",
       ]
    self.assertEqual(
      genConfigBuildNamesList, genConfigBuildNamesList_expected)


if __name__ == '__main__':
  unittest.main()
