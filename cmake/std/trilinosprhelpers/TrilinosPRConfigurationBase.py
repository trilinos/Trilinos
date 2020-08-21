#!/usr/bin/env python
# -*- mode: python; py-indent-offset: 4; py-continuation-offset: 4 -*-
import configparser
import multiprocessing
import os
import re
import subprocess
import sys
from textwrap import dedent


from . import setenvironment
from . import sysinfo
from . import jenkinsenv



class TrilinosPRConfigurationBase(object):
    """
    Trilinos Pull Request configuration driver
    """
    def __init__(self, args):
        self.args                  = args
        self._config_data          = None
        self._trilinos_pr_env      = None
        self._mem_per_core         = None
        self._max_cores_allowed    = None
        self._max_test_parallelism = None
        self._concurrency_build    = None
        self._concurrency_test     = None


    @property
    def arg_github_pr_number(self):
        """
        Argument Wrapper: This property wraps the value provided in self.args
        to provide a convenient way to override this value if needed for some
        specialty reason or for a customized test.

        Returns:
            self.args.github_pr_number
        """
        return self.args.github_pr_number


    @property
    def arg_job_number(self):
        """
        Argument Wrapper: This property wraps the value provided in self.args
        to provide a convenient way to override this value if needed for some
        specialty reason or for a customized test.

        Returns:
            self.args.job_number
        """
        return self.args.job_number


    @property
    def arg_req_mem_per_core(self):
        """
        Argument Wrapper: This property wraps the value provided in self.args
        to provide a convenient way to override this value if needed for some
        specialty reason or for a customized test.

        Returns:
            float self.args.req_mem_per_core
        """
        return float(self.args.req_mem_per_core)


    @property
    def arg_max_cores_allowed(self):
        """
        Argument Wrapper: This property wraps the value provided in self.args
        to provide a convenient way to override this value if needed for some
        specialty reason or for a customized test.

        Returns:
            self.args.max_cores_allowed
        """
        return int(self.args.max_cores_allowed)


    @property
    def arg_num_concurrent_tests(self):
        """
        Argument Wrapper: This property wraps the value provided in self.args
        to provide a convenient way to override this value if needed for some
        specialty reason or for a customized test.

        Returns:
            self.args.num_concurrent_tests
        """
        return int(self.args.num_concurrent_tests)


    @property
    def subprojects_file(self):
        """
        This property generates the subprojects_list file to pass into CTest.
        """
        return os.path.join(self.arg_workspace_dir, "package_subproject_list.cmake")


    @property
    def arg_package_enables_file(self):
        """
        This property generates the package_enables filename. It's in a property so it can
        easily be changed class-wide. The default behavior of this property is to return
        the `package_enables` entry from the program arguments.
        """
        return self.args.package_enables


    @property
    def arg_workspace_dir(self):
        """
        Returns the Jenkins workspace directory for the PR.
        The default behavior of this property is to return the value of
        the `workspaceDir` argument from the program arguments.
        """
        return self.args.workspaceDir


    @property
    def arg_pr_config_file(self):
        """
        Returns the configuration file that we'd want to load.
        The default behavior is to return the value in args.configfile.
        """
        return self.args.configfile


    @property
    def arg_pr_jenkins_job_name(self):
        """
        The Jenkins job name that is executing this Pull Request test.
        Default is to use the value in args.job_base_name.
        """
        return self.args.job_base_name


    @property
    def working_directory_ctest(self):
        """
        Generate the working directory for where we should launch the CTest command from.
        For PR testing this should be in $WORKSPACE/TFW_testing_single_configure_prototype
        """
        return os.path.join(self.arg_workspace_dir, 'TFW_testing_single_configure_prototype')


    @property
    def trilinos_pr_env(self):
        if self._trilinos_pr_env is None:
            self._trilinos_pr_env = jenkinsenv.TrilinosJenkinsEnv()
        return self._trilinos_pr_env


    @property
    def config_data(self):
        if self._config_data is None:
            self._config_data = setenvironment.SetEnvironment(self.args.configfile, self.arg_pr_jenkins_job_name)
        return self._config_data


    @property
    def config_script(self):
        return self.get_property_from_config("CONFIG_SCRIPT_MAP", self.arg_pr_jenkins_job_name)


    @property
    def max_cores_allowed(self):
        """
        Calculate the true max number of cores we can use. This is the minimum of
        the available cores detected on the system and the max cores allowed by
        argument.

        If the max_cores_allowd is not > 0 then we use the # of cores on the system.
        """
        if self._max_cores_allowed is None:
            num_cores = multiprocessing.cpu_count()
            if self.arg_max_cores_allowed > 0:
                self._max_cores_allowed = min(num_cores, int(self.arg_max_cores_allowed))
            else:
                self._max_cores_allowed = num_cores
        return self._max_cores_allowed


    @property
    def max_test_parallelism(self):
        """
        Maximum parallelism of your tests. This is obtained from the configuration.ini
        file from the [PR_JOB_PARAMETERS] section key max-test-parallelism value.
        If this value is missing from the configuration file, we default to 1.
        """
        if self._max_test_parallelism is None:
            try:
                self._max_test_parallelism = int(self.get_property_from_config("PR_JOB_PARAMETERS","max-test-parallelism"))
            except:
                self._max_test_parallelism = 1
        return self._max_test_parallelism


    @property
    def concurrency_build(self):
        """
        Concurrency to use for building Trilinos
        This is equvalent to running `make -j <concurrency_build>` from the command line.
        """
        if self._concurrency_build is None:
            si = sysinfo.SysInfo()

            self._concurrency_build = si.compute_num_usable_cores(req_mem_gb_per_core = self.arg_req_mem_per_core,
                                                                  max_cores_allowed   = self.max_cores_allowed)

        return self._concurrency_build


    @property
    def concurrency_test(self):
        """
        Concurrency to use for running Trilinos tests
        This is equivalent to the command `ctest -j <concurrency_test>` if running ctest at the
        command line.

        If not overridden, we'll compute the cores / max_test_parallelism
        """
        if self._concurrency_test is None:
            if self.arg_num_concurrent_tests > 0:
                self._concurrency_test = self.arg_num_concurrent_tests

            num_cores = self.max_cores_allowed
            self._concurrency_test = max(1, int(num_cores / self.max_test_parallelism))
        return self._concurrency_test


    @property
    def pullrequest_build_name(self):
        """
        Generate the build name string to report back to CDash.
        PR-<PR Number>-test-<Jenkins Job Name>-<Job Number">
        """
        output = "PR-{}-test-{}-{}".format(self.arg_github_pr_number, self.arg_pr_jenkins_job_name, self.arg_job_number)
        return output


    @property
    def pullrequest_cdash_track(self):
        """
        Attempt to load the envvar PULLREQUEST_CDASH_TRACK, but if it's missing,
        we just use "Pull Request"
        """
        try:
            output = self.trilinos_pr_env.pullrequest_cdash_track
        except:
            print("WARNING: envvar PULLREQUEST_CDASH_TRACK missing. Using default 'Pull Request'")
            output = "Pull Request"
        return output


    def get_property_from_config(self, section, option, default=None):
        """
        Helper to load a property from the ConfigParser configuration data.

        Args:
            section (str) : The section name we're loading (i.e., [SECTION_NAME])
            option  (str) : The name of a property inside a section to load (i.e., property: value)
            default       : Default value to return if the request fails. Default: None

        Returns:
            str: The value of a property inside some section as a string. For example:
                 [SECTION_NAME]
                 property: <value>
        """
        output = default
        try:
            output = self.config_data.config.get(section, option)
        except configparser.NoSectionError:
            print("WARNING: Configuration section '{}' does not exist.".format(section))
            print("       : Returning default value: '{}'".format(output))
        except configparser.NoOptionError:
            print("WARNING: Configuration section '{}' has no key '{}'".format(section,option))
            print("       : Returning default value: {}".format(output))
        return output


    def get_multi_property_from_config(self, section, option, default=None, delimeter=','):
        """
        This works the same as get_property_from_config but it allows you to separate a key
        into multiple entries by adding a space then an extra string to distinguish entries.
        Results are merged together into a single delimited string.

        For example, if the section looks like this:

            [section]
            option_one A: foo,bar
            option_one B: baz

        Then the output of a search for `get_multi_property_from_config('section', 'option_one')`
        will return:  "foo,bar,baz"

        We have this because of a limitation in the configparser .ini files in that keys
        inside a section must be unique.
        """
        output = default
        try:
            for option_full in self.config_data.config.options(section):
                if option == option_full.split(" ")[0]:
                    if output is None:
                        output = self.config_data.config.get(section, option_full)
                    else:
                        output = delimeter.join([output, self.config_data.config.get(section, option_full)])
        except configparser.NoSectionError:
            print("WARNING: Configuration section '{}' does not exist.".format(section))
            print("       : Returning default value: {}".format(output))
        return output


    def create_package_enables_file(self, dryrun=False):
        """
        Generate the packageEnables.cmake file.
        If the ENABLE_MAP section in the properties.ini file contains
        a match to this job then we use that information to generate the
        packageEnables.cmake and package_subproject_list.cmake files.
        Otherwise, we generate the list of packages to build by looking
        at the git changelogs via the TriBiTS `get-changed-trilinos-packages.sh`
        script.
        """
        job_name = self.arg_pr_jenkins_job_name

        #enable_map_entry = self.get_property_from_config("ENABLE_MAP", job_name)
        enable_map_entry = self.get_multi_property_from_config("ENABLE_MAP", job_name, delimeter=" ")

        try:
            # Generate files using ATDM/TriBiTS Scripts
            if enable_map_entry is None:
                cmd = [os.path.join( self.arg_workspace_dir,
                                    'Trilinos',
                                    'commonTools',
                                    'framework',
                                    'get-changed-trilinos-packages.sh'),
                       os.path.join('origin', self.args.targetBranch),
                       'HEAD',
                       'packageEnables.cmake',
                       'package_subproject_list.cmake']

                print("")
                print("packageEnables Command: \n$ {}\n".format(" \\\n    ".join(cmd)))

                if not dryrun:
                    subprocess.check_call(cmd)
                else:
                    print("-\n- SKIPPED DUE TO DRYRUN\n-")
            else:
                # Use the values in the PACKAGE_ENABLES section of the .ini file
                with open('packageEnables.cmake',  'w') as f_out:
                    f_out.write(dedent('''\
                        MACRO(PR_ENABLE_BOOL  VAR_NAME  VAR_VAL)
                          MESSAGE("-- Setting ${VAR_NAME} = ${VAR_VAL}")
                          SET(${VAR_NAME} ${VAR_VAL} CACHE BOOL "Set in $CMAKE_PACKAGE_ENABLES_OUT")
                        ENDMACRO()
                        '''))

                    for entry in enable_map_entry.split(" "):
                        f_out.write("PR_ENABLE_BOOL(Trilinos_ENABLE_{} ON)\n".format(entry))
                    #    '''
                    #    PR_ENABLE_BOOL(Trilinos_ENABLE_''' + enable_map_entry + ''' ON)
                    #    '''))
                with open ('package_subproject_list.cmake', 'w') as f_out:
                    f_out.write(dedent('''\
                        set(CTEST_LABELS_FOR_SUBPROJECTS ''' + enable_map_entry + ''')
                        '''))

            print("")
            print("Enabled Packages:")
            cmd = ['cmake', '-P', 'packageEnables.cmake']
            cmake_rstring=None
            if not dryrun:
                cmake_rstring = subprocess.check_output(cmd, stderr=subprocess.STDOUT)
            else:
                print("-\n- SKIPPED DUE TO DRYRUN\n-")
                cmake_rstring = str.encode("")
            cmake_rstring = cmake_rstring.decode('utf-8')
            print(cmake_rstring)

        except subprocess.CalledProcessError as cpe:
            print('There was an issue generating packageEnables.cmake. '
                  'The error code was: {}'.format(cpe.returncode))
            raise cpe

        return 0


    def validate_branch_constraints(self):
        """
        Verify that the source branch is allowed.

        For the `master` branch, we only allow the source branch to be
        a protected branch named with the scheme `master_merge_YYYYMMDD_HHMMSS`
        """
        print("")
        print("Validate target branch constraints:")
        print("- Target branch is '{}'".format(self.args.targetBranch))

        re_master_merge_source = "master_merge_[0-9]{8}_[0-9]{6}"
        if "master" == self.args.targetBranch:
            print("- Target branch is 'master'. Checking source branch constraints...")
            if not re.match(re_master_merge_source, self.args.sourceBranch):
                message  = "---------------------------------------------------------------------------\n"
                message += "ERROR: Source branch is NOT trilinos/Trilinos::master_merge_YYYYMMDD_HHMMSS\n"
                message += "       This violates Trilinos policy for pull requests into the master\n"
                message += "       branch.\n"
                message += "       Source branch provided is {}\n".format(self.args.sourceBranch)
                message += "       Perhaps you forgot to set `develop` as the target in your PR?\n"
                message += "---------------------------------------------------------------------------"
                #print(message)
                sys.exit(message)

        print("- target branch constraints OK")
        print("")
        return 0


    def prepare_test(self):
        """
        Test preparation
        Things that we generally do to prepare the environment for the test
        can go here.
        """
        # Validate the branch constraints (i.e., if targetBranch is master, then
        # sourceBranch must be master_merge_YYYYMMDD_HHMMSS)
        self.validate_branch_constraints()

        print("Configuration Parameters")
        print("------------------------")
        print(">>> arg_pr_jenkins_job_name  = {}".format(self.arg_pr_jenkins_job_name))
        print(">>> arg_job_number           = {}".format(self.arg_job_number))
        print(">>> arg_github_pr_number     = {}".format(self.arg_github_pr_number))
        print(">>> arg_max_cores_allowed    = {}".format(self.arg_max_cores_allowed))
        print(">>> arg_req_mem_per_core     = {}".format(self.arg_req_mem_per_core))
        print(">>> arg_num_concurrent_tests = {}".format(self.arg_num_concurrent_tests))
        print(">>> arg_pr_config_file       = {}".format(self.arg_pr_config_file))
        print(">>> arg_package_enables_file = {}".format(self.arg_package_enables_file))
        print(">>> arg_workspace_dir        = {}".format(self.arg_workspace_dir))
        print(">>> config_script            = {}".format(self.config_script))
        print(">>> pullrequest_build_name   = {}".format(self.pullrequest_build_name))
        print(">>> pullrequest_cdash_track  = {}".format(self.pullrequest_cdash_track))
        print(">>> subprojects_file         = {}".format(self.subprojects_file))
        print(">>> concurrency_build        = {}".format(self.concurrency_build))
        print(">>> concurrency_test         = {}".format(self.concurrency_test))
        print(">>> max_cores_allowed        = {}".format(self.max_cores_allowed))
        print(">>> max_test_parallelism     = {}".format(self.max_test_parallelism))
        print(">>> working_directory_ctest  = {}".format(self.working_directory_ctest))
        #print(">>>         = {}".format(self.))
        print("")


        print("="*80)
        print("=  E N V I R O N M E N T   S E T   U P   S T A R T")
        print("="*80)
        tr_config = setenvironment.SetEnvironment(self.arg_pr_config_file, self.arg_pr_jenkins_job_name)

        rval = 0
        if not self.args.dry_run:
            rval = tr_config.apply(throw_on_error=True)
            print("rval = {}".format(rval))
        else:
            tr_config.pretty_print()
            print("\nNOTICE: ENVVARS not set due to dry-run flag.")

        if rval:
            msg = "ERROR: There was a problem configuring the environment."
            print(msg)
            raise Exception(msg)

        envvars = [
            "SEMS_",
            "TRILINOS_",
            "PULLREQUEST",
            "JENKINS",
            "KOKKOS",
            "BUILD",
            "JOB",
            "FORCE_CLEAN",
            "proxy",
            "PROXY",
            "PATH",
            "OMP_",
            "WORKSPACE"
            ]
        print("")
        tr_config.pretty_print_envvars(envvar_filter=envvars)

        print("="*80)
        print("=  E N V I R O N M E N T   S E T   U P   C O M P L E T E")
        print("="*80)
        print("")

        print("="*80)
        print("=  Create packageEnables.cmake")
        print("="*80)
        self.create_package_enables_file(dryrun=self.args.dry_run)
        print("")


        return 0


    def execute_test(self):
        """
        Executes the test. This method should be considered 'virtual' and should
        be overridden.
        """
        raise NotImplementedError("This method must be overridden.")




