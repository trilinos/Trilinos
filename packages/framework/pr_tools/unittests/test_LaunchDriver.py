#!/usr/bin/env python3
# -*- mode: python; py-indent-offset: 4; py-continuation-offset: 4 -*-
'''
Tests for LaunchDriver
'''
from __future__ import print_function
import sys
sys.dont_write_bytecode = True

import os
sys.path.insert(1, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import unittest

import LaunchDriver as ld



class Test_LaunchDriver(unittest.TestCase):
    def setUp(self):
      pr_tools_path = os.path.dirname(os.path.realpath(__file__))
      self.build_name = "rhel7_stack0"
      args = ["--supported-systems="+pr_tools_path+"/supporting_files/supported-systems.ini", "--build-name="+self.build_name]
      self.args0 = args + ["--driver="+pr_tools_path+"/DriverTest0.sh"]
      self.args1 = args + ["--driver="+pr_tools_path+"/DriverTest1.sh"]


    ## Test LaunchDriver methods
    def testUnitGetLaunchEnv(self):
      env = ld.get_launch_env("dne")
      self.assertEqual(env, "")


    def testUnitGetLaunchCmd(self):
      cmd = ld.get_launch_cmd("dne")
      self.assertEqual(cmd, " ")


    def testUnitGetDriverArgs(self):
      args = ld.get_driver_args("dne1")
      self.assertEqual(args, " --on_dne1")


    ## Test LaunchDriver main
    def testIntegration0(self):
      with self.assertRaises(SystemExit) as se:
        ld.main(self.args0)
      self.assertEqual(se.exception.code, 0)


    def testIntegration1(self):
      with self.assertRaises(SystemExit) as se:
        ld.main(self.args1)
      self.assertEqual(se.exception.code, 1)


if __name__ == '__main__':
    unittest.main()  # pragma nocover
