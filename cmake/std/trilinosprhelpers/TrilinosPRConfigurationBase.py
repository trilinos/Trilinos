#!/usr/bin/env python
# -*- mode: python; py-indent-offset: 4; py-continuation-offset: 4 -*-
"""
This file contains the base class for the Pull Request test driver.
"""
import configparser
import multiprocessing
import os
import re
import subprocess
import sys
from textwrap import dedent

sys.dont_write_bytecode = True


from . import setenvironment
from . import sysinfo
#from . import jenkinsenv



class TrilinosPRConfigurationBase(object):
    """
    Trilinos Pull Request configuration driver. This should be
    treated as an Abstract Base Class because the `execute()`
    method is implemented as a stub.

    Note: Attributes / Properties prefixed with `arg_` come from the
        command line arguments that is passed into this class c'tor.
        We provide these properties to give a read-only access to the
        properties.

    Attributes:
        arg_pullrequest_number: PR Number on Github
        arg_pullrequest_cdash_track: The Pull Request track to post results to.
                                     Defaults to "Pull Request" in arguments.
        arg_jenkins_job_number: Job Number of the Jenkins job.
        arg_req_mem_per_core: Memory required per core (GB)
        arg_max_cores_allowed: Max cores allowed for building (i.e., make -j <num cores>)
        arg_num_concurrent_tests: Testing concurrency (i.e., ctest -j <num cores>)
        arg_filename_packageenables: Path to the `packageEnables.cmake` file.
        arg_workspace_dir: Path to the workspace (this would be the ${WORKSPACE}
            variable set by Jenkins.)
        arg_pr_config_file: The config.ini file that specifies the configuration to load.
        arg_pr_jenkins_job_name: The Jenkins Job Name.
        filename_subprojects: The subprojects file.
        working_directory_ctest: Gen. working dir where TFW_testing_single_configure_prototype
            is executed from.
        config_data: The setenvironment.SetEnvironment class instance containing
                     the parsed config.ini file data.
        config_script: Returns the configuration script from the configuration file.
        max_cores_allowed: Absolute maximum number of cores allowed.
            For example, if the job weight is 32 (even if we have an 80 core system)
            this would be set to 32.
        max_test_parallelism: Maximum test parallelism in any test. For example, if
            the highest individual test parallelism is 4 cores or threads then this
            would be set to 4.
        concurrency_build: Concurrency to use for building Trilinos
        concurrency_test: Concurrency to use for running Trilinos tests.
        pullrequest_build_name: PR build name reported to CDash.
    """
    def __init__(self, args):
        self.args                  = args
        self._config_data          = None
        self._mem_per_core         = None
        self._max_cores_allowed    = None
        self._max_test_parallelism = None
        self._concurrency_build    = None
        self._concurrency_test     = None


    # --------------------
    # A R G U M E N T S
    # --------------------

    @property
    def arg_pullrequest_number(self):
        """
        Argument Wrapper: This property wraps the value provided in self.args
        to provide a convenient way to override this value if needed for some
        specialty reason or for a customized test.

        Returns:
            self.args.pullrequest_number
        """
        return self.args.pullrequest_number


    @property
    def arg_pullrequest_cdash_track(self):
        """
        Argument Wrapper: This property wraps the value provided in self.args
        to provide a convenient way to override this value if needed for some
        specialty reason or for a customized test.

        Returns:
            self.args.pullrequest_cdash_track
        """
        return self.args.pullrequest_cdash_track


    @property
    def arg_jenkins_job_number(self):
        """
        Argument Wrapper: This property wraps the value provided in self.args
        to provide a convenient way to override this value if needed for some
        specialty reason or for a customized test.

        Returns:
            self.args.jenkins_job_number
        """
        return self.args.jenkins_job_number


    @property
    def arg_req_mem_per_core(self):
        """
        Argument Wrapper: This property wraps the value provided in self.args
        to provide a convenient way to override this value if needed for some
        specialty reason or for a customized test.

        This encodes the # of GB/core required.

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
    def arg_filename_packageenables(self):
        """
        This property generates the packageEnables filename. It's in a property so it can
        easily be changed class-wide. The default behavior of this property is to return
        the `filename_packageenables` entry from the program arguments.
        """
        return self.args.filename_packageenables


    @property
    def arg_workspace_dir(self):
        """
        Returns the Jenkins workspace directory for the PR.
        The default behavior of this property is to return the value of
        the `workspace_dir` argument from the program arguments.
        """
        return self.args.workspace_dir


    @property
    def arg_pr_config_file(self):
        """
        Returns the configuration file that we'd want to load.
        The default behavior is to return the value in args.pullrequest_config_file.
        """
        return self.args.pullrequest_config_file


    @property
    def arg_pr_jenkins_job_name(self):
        """
        The Jenkins job name that is executing this Pull Request test.
        Default is to use the value in args.pullrequest_build_name.
        """
        return self.args.pullrequest_build_name


    @property
    def arg_filename_subprojects(self):
        """
        This property generates the subprojects_list file to pass into CTest.
        """
        return self.args.filename_subprojects


    # --------------------
    # P R O P E R T I E S
    # --------------------

    @property
    def working_directory_ctest(self):
        """
        Generate the working directory for where we should launch the CTest command from.
        For PR testing this should be in $WORKSPACE/TFW_testing_single_configure_prototype
        """
        return os.path.join(self.arg_workspace_dir, 'TFW_testing_single_configure_prototype')


    @property
    def config_data(self):
        """
        Configuration data.

        This is a setenvironment.SetEnvironment class instance containing the information
        loaded from the config.ini file.
        """
        if self._config_data is None:
            self._config_data = setenvironment.SetEnvironment(self.arg_pr_config_file,
                                                              self.arg_pr_jenkins_job_name)
        return self._config_data


    @property
    def config_script(self):
        """
        Returns the configuration script from the configuration file.

        This function searches the [CONFIG_SCRIPT_MAP] section in the config.ini file
        and looks for the key matching the Jenkins job name and returns the
        value from that KV pair.  This would be the .cmake file that maps to the
        job.

        Returns:
            String containing the job-specific configuration script to load.
        """
        return self.get_property_from_config("CONFIG_SCRIPT_MAP", self.arg_pr_jenkins_job_name)


    @property
    def max_cores_allowed(self):
        """
        Calculate the true max number of cores we can use.

        This is the minimum of the available cores detected on the system
        and the max cores allowed by argument.

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

        If num_concurrent_tests is > 0 then we just use this value.
        Otherwise, we compute the cores using num_available_cores / max_test_parallelism.
        """
        if self._concurrency_test is None:
            # If >0 then we use the value provided by the user
            if self.arg_num_concurrent_tests > 0:
                self._concurrency_test = self.arg_num_concurrent_tests

            # otherwise we calculate based on number of allowable cores and max test parallelism.
            else:
                num_cores = self.max_cores_allowed
                self._concurrency_test = max(1, int(num_cores / self.max_test_parallelism))
        return self._concurrency_test


    @property
    def pullrequest_build_name(self):
        """
        Generate the build name string to report back to CDash.

        PR-<PR Number>-test-<Jenkins Job Name>-<Job Number">
        """
        output = "PR-{}-test-{}-{}".format(self.arg_pullrequest_number, self.arg_pr_jenkins_job_name, self.arg_jenkins_job_number)
        return output


    # --------------------
    # M E T H O D S
    # --------------------

    def get_property_from_config(self, section, option, default=None):
        """
        Helper to load a property from the ConfigParser configuration data
        that was loaded up by the SetEnvironment class. We use this to load
        information from free-form sections that we include int he .ini file
        which encode information that isn't used by the parser but is used
        to encode Trilinos PR configuration data.

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
                        output = delimeter.join(
                            [output, self.config_data.config.get(section, option_full)]
                        )
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
                       os.path.join('origin', self.args.target_branch_name),
                       'HEAD',
                       'packageEnables.cmake',
                       'package_subproject_list.cmake']

                print("")
                print("packageEnables Command: \n$ {}\n".format(" \\\n    ".join(cmd)))

                if not dryrun:
                    subprocess.check_call(cmd)
                else:
                    print("")
                    print("--- SKIPPED DUE TO DRYRUN")
                    print("")
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
                print("")
                print("--- SKIPPED DUE TO DRYRUN")
                print("")
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
        print("--- Target branch is '{}'".format(self.args.target_branch_name))

        re_master_merge_source = "master_merge_[0-9]{8}_[0-9]{6}"
        if "master" == self.args.target_branch_name:
            print("--- Target branch is 'master'. Checking source branch constraints...")
            if not re.match(re_master_merge_source, self.args.source_branch_name):
                message  = "+" + "="*78 + "+\n"
                message += "ERROR: Source branch is NOT trilinos/Trilinos::master_merge_YYYYMMDD_HHMMSS\n"
                message += "       This violates Trilinos policy for pull requests into the master\n"
                message += "       branch.\n"
                message += "       Source branch provided is {}\n".format(self.args.source_branch_name)
                message += "       Perhaps you forgot to set `develop` as the target in your PR?\n"
                message += "+" + "="*78 + "+\n"
                #print(message)
                sys.exit(message)

        print("--- target branch constraints OK")
        print("")
        return 0


    def prepare_test(self):
        """
        Prepares a test environment for exeution.

        This includes tasks like determining the # of cores to use, setting
        environment variables, loading environment modules, etc.
        """
        # Validate the branch constraints (i.e., if target_branch_name is master, then
        # source_branch_name must be master_merge_YYYYMMDD_HHMMSS)
        self.validate_branch_constraints()

        print("+" + "="*78 + "+")
        print("Configuration Parameters")
        print("+" + "="*78 + "+")
        print("--- arg_filename_packageenables = {}".format(self.arg_filename_packageenables))
        print("--- arg_filename_subprojects    = {}".format(self.arg_filename_subprojects))
        print("--- arg_jenkins_job_number      = {}".format(self.arg_jenkins_job_number))
        print("--- arg_max_cores_allowed       = {}".format(self.arg_max_cores_allowed))
        print("--- arg_num_concurrent_tests    = {}".format(self.arg_num_concurrent_tests))
        print("--- arg_pr_config_file          = {}".format(self.arg_pr_config_file))
        print("--- arg_pr_jenkins_job_name     = {}".format(self.arg_pr_jenkins_job_name))
        print("--- arg_pullrequest_number      = {}".format(self.arg_pullrequest_number))
        print("--- arg_pullrequest_cdash_track = {}".format(self.arg_pullrequest_cdash_track))
        print("--- arg_req_mem_per_core        = {}".format(self.arg_req_mem_per_core))
        print("--- arg_workspace_dir           = {}".format(self.arg_workspace_dir))
        print("")
        print("--- concurrency_build           = {}".format(self.concurrency_build))
        print("--- concurrency_test            = {}".format(self.concurrency_test))
        print("--- config_script               = {}".format(self.config_script))
        print("--- max_cores_allowed           = {}".format(self.max_cores_allowed))
        print("--- max_test_parallelism        = {}".format(self.max_test_parallelism))
        print("--- pullrequest_build_name      = {}".format(self.pullrequest_build_name))
        print("--- working_directory_ctest     = {}".format(self.working_directory_ctest))
        #print("---         = {}".format(self.))
        print("")


        print("+" + "="*68 + "+")
        print("|   E N V I R O N M E N T   S E T   U P   S T A R T")
        print("+" + "="*68 + "+")
        tr_config = setenvironment.SetEnvironment(self.arg_pr_config_file, self.arg_pr_jenkins_job_name)

        rval = 0
        if not self.args.dry_run:
            rval = tr_config.apply(throw_on_error=True)
            print("apply() rval: {}".format(rval))
        else:
            tr_config.pretty_print()
            print("")
            print("--- NOTICE: ENVVARS not set due to dry-run flag.")
            print("")

        if rval:
            msg = "ERROR: There was a problem configuring the environment."
            print(msg)
            raise Exception(msg)

        # Environment variables that we wish to print out.
        envvars_to_print = [
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
        tr_config.pretty_print_envvars(envvar_filter=envvars_to_print)

        print("+" + "="*68 + "+")
        print("|   E N V I R O N M E N T   S E T   U P   C O M P L E T E")
        print("+" + "="*68 + "+")

        print("--- Create packageEnables.cmake")
        self.create_package_enables_file(dryrun=self.args.dry_run)
        print("")

        return 0


    def execute_test(self):
        """
        Executes the test. This method must be overridden by a subclass.
        """
        raise NotImplementedError("This method must be overridden.")




