#!/usr/bin/env python
# -*- mode: python; py-indent-offset: 4; py-continuation-offset: 4 -*-
"""
Custom PR Executor for Standard testing
"""

import os
import subprocess

from . import TrilinosPRConfigurationBase
from gen_config import GenConfig
from pathlib import Path


class TrilinosPRConfigurationStandard(TrilinosPRConfigurationBase):
    """
    Implements Standard mode Trilinos Pull Request Driver
    """
    def __init__(self, args):
        super(TrilinosPRConfigurationStandard, self).__init__(args)


    def execute_test(self):
        """
        Execute the test
        """
        self.message("+" + "="*78 + "+")
        self.message("|   E X E C U T E   S T A N D A R D   P U L L R E Q E S T   T E S T")
        self.message("+" + "="*78 + "+")

        #
        # Typically, we execute the test from $WORKSPACE/TFW_testing_single_configure_prototype
        # We'll skip it if we're doing a dry-run.
        #
        self.message("")
        #self.message("--- Change directory to {}".format(self.working_directory_ctest))
        #if not self.args.dry_run:
        #    os.chdir(self.working_directory_ctest)
        self.chdir_logged(self.arg_build_dir, create_if_missing=True)

        # Use GenConfig to write the configure script for cmake
        genconfig_arglist = ["-y",
                             "--force",
                             "--cmake-fragment",
                             os.path.join(self.arg_workspace_dir, self.config_script),
                             self.arg_pr_genconfig_job_name
                             ]

        genconfig_inifile = Path(self.arg_pr_gen_config_file)

        self.message( "--- GenConfig:")
        gc = GenConfig(genconfig_arglist, gen_config_ini_file=genconfig_inifile)

        if not self.args.dry_run:
            gc.write_cmake_fragment()

        # Execute the call to ctest.
        # - NOTE: simple_testing.cmake can be found in the TFW_single_configure_support_scripts
        #         repository.
        cmd = ['ctest',
                "-S", f"{self.arg_ctest_driver}",
               f"-Dsource_dir:PATH={self.arg_source_dir}",
               f"-Dbuild_dir:PATH={self.arg_build_dir}",
               f"-Dbuild_name:STRING={self.pullrequest_build_name}",
                "-Dskip_by_parts_submit:BOOL=OFF",
                "-Dskip_update_step:BOOL=ON",
                "-Ddashboard_model:STRING='Experimental'",
               f"-Ddashboard_track:STRING='{self.arg_pullrequest_cdash_track}'",
               f"-DPARALLEL_LEVEL:STRING={self.concurrency_build}",
               f"-DTEST_PARALLEL_LEVEL:STRING={self.concurrency_test}",
                "-Dconfigure_script:FILEPATH=" + os.path.join(self.arg_workspace_dir, self.config_script),
               f"-Dpackage_enables:FILEPATH={self.arg_filename_packageenables}",
               f"-Dsubprojects_file:FILEPATH={self.arg_filename_subprojects}"
             ]

        self.message( "--- ctest version:")
        if not self.args.dry_run:
            try:
                # Check the version of CTest
                subprocess.check_call(['ctest', '--version'], env=os.environ)
            except subprocess.CalledProcessError as err:
                self.message("-- ERROR: `ctest --version` failed.")
                self.message(err)
                return 1

        self.message( "--- ctest command:")
        self.message(f"--- pwd:\n{os.getcwd()}")
        # self.message( "--- cmd:\n{}".format(" \\\n   ".join(cmd)))
        self.message( "--- cmd:")
        for i in range(len(cmd)):
            self.message(f"  ARG {i} {cmd[i]}")
        self.message("")

        if not self.args.dry_run:
            try:
                subprocess.check_call(cmd, env=os.environ)
                # Note: check_call will throw an exception if there's a problem.
                self.message("--- OK")
            except subprocess.CalledProcessError as err:
                self.message("--- ctest command failed!")
                self.message("--- error:")
                self.message(err)
                return 1
        else:
            self.message("--- SKIPPED DUE TO DRYRUN")
        self.message("")

        return 0




