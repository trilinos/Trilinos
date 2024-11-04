#!/usr/bin/env python3
# -*- mode: python; py-indent-offset: 4; py-continuation-offset: 4 -*-
"""
This file contains the base class for the Pull Request test driver.
"""
import configparserenhanced
import inspect
import multiprocessing
import os
from pathlib import Path
import re
import subprocess
import sys
from textwrap import dedent

sys.dont_write_bytecode = True

from .sysinfo import SysInfo
from LoadEnv.load_env import LoadEnv
import setenvironment
from .sysinfo import gpu_utils



class TrilinosPRConfigurationBase(object):
    """
    Trilinos Pull Request configuration driver. This should be
    treated as an Abstract Base Class because the `execute_test()`
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
        arg_ccache_enable: Enable ccache.
        arg_dashboard_build_name: A shortened genconfig build name
                                  for posting to a testing dashboard.
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
        self.load_env_ini_file     = None
        self._config_data          = None
        self._mem_per_core         = None
        self._max_cores_allowed    = None
        self._max_test_parallelism = None
        self._concurrency_build    = None
        self._concurrency_test     = None
        self._debug_level          = 1
        self._arg_extra_configure_args = None


    # --------------------
    # A R G U M E N T S
    # --------------------

    @property
    def arg_extra_configure_args(self):
        """
        Argument Wrapper: This property wraps the value provided in self.args
        to provide a convenient way to override this value if needed for some
        specialty reason or for a customized test.

        This parameter stores extra configure arguments that will be passed
        to the cmake call when configuring Trilinos.

        Returns:
            self.args.extra_configure_args
        """
        if not self._arg_extra_configure_args:
            if gpu_utils.has_nvidia_gpus():
                self.message("-- REMARK: I see that I am running on a machine that has NVidia GPUs; I will feed TriBITS some data enabling GPU resource management")
                slots_per_gpu = 2
                gpu_indices = gpu_utils.list_nvidia_gpus()
                self.message(f"-- REMARK: Using {slots_per_gpu} slots per GPU")
                self.message(f"-- REMARK: Using GPUs {gpu_indices}")
                self._arg_extra_configure_args = f"-DTrilinos_AUTOGENERATE_TEST_RESOURCE_FILE:BOOL=ON;-DTrilinos_CUDA_NUM_GPUS:STRING={len(gpu_indices)};-DTrilinos_CUDA_SLOTS_PER_GPU:STRING={slots_per_gpu}" + (";" + self.args.extra_configure_args if self.args.extra_configure_args else "")
            else:
                self._arg_extra_configure_args = self.args.extra_configure_args
        return self._arg_extra_configure_args

    @property
    def arg_ctest_driver(self):
        """
        Argument Wrapper: This property wraps the value provided in self.args
        to provide a convenient way to override this value if needed for some
        specialty reason or for a customized test.

        This parameter stores the location of the CTest driver script that gets
        loaded by the -S argument.

        Returns:
            self.args.ctest_driver
        """
        return self.args.ctest_driver

    @property
    def arg_ctest_drop_site(self):
        """
        Argument Wrapper: This property wraps the value provided in self.args
        to provide a convenient way to override this value if needed for some
        specialty reason or for a customized test.

        This parameter stores the location of the CTest driver script that gets
        loaded by the -S argument.

        Returns:
            self.args.ctest_drop_site
        """
        return self.args.ctest_drop_site

    @property
    def arg_use_explicit_cachefile(self):
        """
        Argument Wrapper: This property wraps the value provided in self.args
        to provide a convenient way to override this value if needed for some
        specialty reason or for a customized test.

        This parameter stores whether or not an explicit cachefile directive
        will be passed (as opposed to using -C).

        Returns:
            self.args.use_explicit_cachefile
        """
        return self.args.use_explicit_cachefile

    @property
    def arg_build_dir(self):
        """
        Argument Wrapper: This property wraps the value provided in self.args
        to provide a convenient way to override this value if needed for some
        specialty reason or for a customized test.

        Returns:
            self.args.build_dir
        """
        return self.args.build_dir


    @property
    def arg_source_dir(self):
        """
        Argument Wrapper: This property wraps the value provided in self.args
        to provide a convenient way to override this value if needed for some
        specialty reason or for a customized test.

        Returns:
            self.args.source_dir
        """
        return self.args.source_dir


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
    def arg_pr_env_config_file(self):
        """
        Returns the configuration file that we'd want to load.
        The default behavior is to return the value in args.pullrequest_env_config_file.
        """
        return self.args.pullrequest_env_config_file


    @property
    def arg_pr_gen_config_file(self):
        """
        Returns the configuration file that we'd want to load.
        The default behavior is to return the value in args.pullrequest_gen_config_file.
        """
        return self.args.pullrequest_gen_config_file


    @property
    def arg_pr_jenkins_job_name(self):
        """
        The Jenkins job name that is executing this Pull Request test.
        Default is to use the value in args.pullrequest_build_name.
        """
        return self.args.pullrequest_build_name


    @property
    def arg_pr_genconfig_job_name(self):
        """
        The Jenkins job name that is executing this Pull Request test.
        Default is to use the value in args.pullrequest_build_name.
        """
        return self.args.genconfig_build_name

    @property
    def arg_dashboard_build_name(self):
        """
        The simplified genconfig build name containing only the
        special attributes of the full build name.
        Default is to use the value in args.dashboard_build_name.
        """
        return self.args.dashboard_build_name

    @property
    def arg_filename_subprojects(self):
        """
        This property generates the subprojects_list file to pass into CTest.
        """
        return self.args.filename_subprojects

    @property
    def arg_ccache_enable(self):
        """Is ccache enabled?"""
        return self.args.ccache_enable

    # --------------------
    # P R O P E R T I E S
    # --------------------

    @property
    def working_directory_ctest(self):
        """
        Generate the working directory for where we should launch the CTest command from.
        For PR testing this should be in $WORKSPACE/TFW_testing_single_configure_prototype

        This is hard-coded to the location where the ctest driver is and will be set to the
        working dir when CTest is called.

        DEPRECATION:
            - This may be deprecated with the parameter --ctest-driver
        """
        #return os.path.join(self.arg_workspace_dir, 'TFW_testing_single_configure_prototype')
        # Set to the location where ctest-driver.cmake lives.
        return os.path.join(self.arg_workspace_dir, "pr-ctest-framework", "cmake")


    @property
    def config_data(self):
        """
        Configuration data.

        This is a setenvironment.SetEnvironment class instance containing the information
        loaded from the config.ini file.
        """
        if self._config_data is None:
            self._load_env_config_data = configparserenhanced.ConfigParserEnhanced(
                Path(self.arg_pr_env_config_file)
                ).configparserenhanceddata

            if not self._load_env_config_data.has_section("load-env"):
                msg = f"'{self.load_env_ini_file}' must contain a 'load-env' section."
                raise ValueError(self.get_formatted_msg(msg))

            pr_specs_key = "pullrequest-specs"
            if not self._load_env_config_data.has_option("load-env",
                                                        pr_specs_key):
                raise ValueError(
                    f"'{self.load_env_ini_file}' must contain the "
                    "following in the 'load-env' section: "
                    "{key} : /path/to/{key}.ini".format(key=pr_specs_key)
                )

            # try a relative path first
            self._config_data_path = Path(
                os.path.join(Path(self.arg_pr_env_config_file).parent,
                self._load_env_config_data['load-env']['pullrequest-specs'])
                ).resolve()
            # if the file does not exist at the relative path try an
            # absolute path
            if not self._config_data_path.is_file():
                self._config_data_path = Path(
                self._load_env_config_data['load-env']['pullrequest-specs']
                ).resolve()

            self._config_data = configparserenhanced.ConfigParserEnhanced(
                self._config_data_path).configparserenhanceddata

        return self._config_data


    @property
    def config_script(self):
        """
        Returns the configuration script name

        This arbitrary name will be used for all runs until an override is established

        Returns:
            String containing the job-specific configuration script to load.
        """
        return "generatedPRFragment.cmake"


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
                self._max_test_parallelism = \
                    int(self.get_property_from_config("PR_JOB_PARAMETERS",
                                                      "max-test-parallelism"))
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
            si = SysInfo()

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

        PR-<PR Number>-test-<Jenkins Job Name>-<Job Number>
        """
        if "Pull Request" in self.arg_pullrequest_cdash_track:
            output = f"PR-{self.arg_pullrequest_number}-test-{self.arg_pr_genconfig_job_name}"
            if not self.arg_jenkins_job_number or "UNKNOWN" not in str(self.arg_jenkins_job_number):
                output = f"{output}-{self.arg_jenkins_job_number}"
        elif self.arg_dashboard_build_name != "__UNKNOWN__":
            output = self.arg_dashboard_build_name
        else:
            output = self.arg_pr_genconfig_job_name            
        return output


    @property
    def dashboard_model(self):
        """
        Generate the dashboard model for CDash

        Nightly, Continuous, Experimental
        """
        if self.arg_pullrequest_cdash_track in ["Pull Request", "Experimental"]:
            return "Experimental"
        return "Nightly"


    # --------------------
    # M E T H O D S
    # --------------------

    def get_formatted_msg(self, msg):
        pass

    def message(self, text, debug_level_override=None):
        """
        A simple wrapper to print out a message to the console output that
        can also provide (if debugging) a prefix that indicates what file,
        line, and function the call is coming from.

        Args:
            text (str): The text that will be printed.
        """
        rval = None
        debug_level = self._debug_level if debug_level_override is None else debug_level_override

        if debug_level == 0:
            rval = print(text)
        else:
            sframe   = inspect.stack()[1]
            filename = os.path.basename(sframe.filename)
            lineno   = sframe.lineno
            function = sframe.function
            rval = print(f"{filename}:{lineno} {function}()> {text}")
        sys.stdout.flush()
        sys.stderr.flush()
        return rval


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
            output = self.config_data.get(section, option)
        except KeyError:
            self.message("WARNING: Configuration section '{}' does not exist.".format(section))
            self.message("       : Returning default value: '{}'".format(output))
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
            for option_full in self.config_data.options(section):
                if option == option_full.split(" ")[0]:
                    if output is None:
                        output = self.config_data.get(section, option_full)
                    else:
                        output = delimeter.join(
                            [output, self.config_data.get(section, option_full)]
                        )
        except KeyError:
            self.message("WARNING: Configuration section '{}' does not exist.".format(section))
            self.message("       : Returning default value: {}".format(output))
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

        enable_map_entry = self.get_multi_property_from_config("ENABLE_MAP", job_name, delimeter=" ")

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
                   'package_subproject_list.cmake',
                   '2>&1']

            self.message("")
            self.message(f"packageEnables Command:")
            self.message("{}".format(" \\\n    ".join(cmd)), debug_level_override=1)
            #print("packageEnables Command: \n$ {}\n".format(" \\\n    ".join(cmd)))

            if not dryrun:
                try:
                    sys.stdout.flush()
                    sys.stderr.flush()
                    subprocess.check_call(cmd)

                except subprocess.CalledProcessError as cpe:
                    self.message("--- There was an issue generating `packageEnables.cmake`.")
                    self.message("--- The error code was: {}".format(cpe.returncode))
                    self.message("--- Console Output:\n{}".format(cpe.output))
                    sys.stdout.flush()
                    sys.stderr.flush()
                    raise cpe
            else:
                self.message("")
                self.message("--- SKIPPED DUE TO DRYRUN")
                self.message("")
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

        self.message("")
        self.message("Enabled Packages:")
        cmd = ['cmake', '-P', 'packageEnables.cmake', '2>&1']
        cmake_rstring=None

        if not dryrun:
            try:
                cmake_rstring = subprocess.check_output(cmd, stderr=subprocess.STDOUT)
            except subprocess.CalledProcessError as cpe:
                self.message("--- There was an issue generating `packageEnables.cmake`.")
                self.message("--- The error code was: {}\n".format(cpe.returncode))
                self.message("--- Console Output:\n{}".format(cpe.output))
                raise cpe
        else:
            self.message("")
            self.message("--- SKIPPED DUE TO DRYRUN")
            self.message("")
            cmake_rstring = str.encode("")

        cmake_rstring = cmake_rstring.decode('utf-8')
        self.message(cmake_rstring)

        return 0


    def prepare_test(self):
        """
        Prepares a test environment for exeution.

        This includes tasks like determining the # of cores to use, setting
        environment variables, loading environment modules, etc.
        """

        self.message("+" + "-"*78 + "+")
        self.message("Configuration Parameters")
        self.message("+" + "-"*78 + "+")
        self.message("--- arg_filename_packageenables = {}".format(self.arg_filename_packageenables))
        self.message("--- arg_filename_subprojects    = {}".format(self.arg_filename_subprojects))
        self.message("--- arg_jenkins_job_number      = {}".format(self.arg_jenkins_job_number))
        self.message("--- arg_max_cores_allowed       = {}".format(self.arg_max_cores_allowed))
        self.message("--- arg_num_concurrent_tests    = {}".format(self.arg_num_concurrent_tests))
        self.message("--- arg_pr_env_config_file      = {}".format(self.arg_pr_env_config_file))
        self.message("--- arg_pr_gen_config_file      = {}".format(self.arg_pr_gen_config_file))
        self.message("--- arg_pr_jenkins_job_name     = {}".format(self.arg_pr_jenkins_job_name))
        self.message("--- arg_pr_genconfig_job_name   = {}".format(self.arg_pr_genconfig_job_name))
        self.message("--- arg_dashboard_build_name    = {}".format(self.arg_dashboard_build_name))
        self.message("--- arg_pullrequest_number      = {}".format(self.arg_pullrequest_number))
        self.message("--- arg_pullrequest_cdash_track = {}".format(self.arg_pullrequest_cdash_track))
        self.message("--- arg_req_mem_per_core        = {}".format(self.arg_req_mem_per_core))
        self.message("--- arg_workspace_dir           = {}".format(self.arg_workspace_dir))
        self.message("--- arg_source_dir              = {}".format(self.arg_source_dir))
        self.message("--- arg_build_dir               = {}".format(self.arg_build_dir))
        self.message("--- arg_ctest_driver            = {}".format(self.arg_ctest_driver))
        self.message("--- arg_ctest_drop_site         = {}".format(self.arg_ctest_drop_site))
        self.message("--- arg_ccache_enable           = {}".format(self.arg_ccache_enable))
        self.message("")
        self.message("--- concurrency_build           = {}".format(self.concurrency_build))
        self.message("--- concurrency_test            = {}".format(self.concurrency_test))
        self.message("--- config_script               = {}".format(self.config_script))
        self.message("--- max_cores_allowed           = {}".format(self.max_cores_allowed))
        self.message("--- max_test_parallelism        = {}".format(self.max_test_parallelism))
        self.message("--- pullrequest_build_name      = {}".format(self.pullrequest_build_name))
        self.message("--- working_directory_ctest     = {}".format(self.working_directory_ctest))
        #print("---         = {}".format(self.))
        self.message("")


        self.message("+" + "-"*68 + "+")
        self.message("|   E N V I R O N M E N T   S E T   U P   S T A R T")
        self.message("+" + "-"*68 + "+")
        tr_env = LoadEnv([self.arg_pr_genconfig_job_name, "--force"],
                         load_env_ini_file=Path(self.arg_pr_env_config_file))
        tr_env.load_set_environment()

        rval = 0
        if not self.args.dry_run:
            rval = tr_env.apply_env()
            self.message("--- Environment setup completed ({})".format(rval))
        else:
            if tr_env.set_environment is None:
                tr_env.load_set_environment()
            tr_env.set_environment.pretty_print_actions(tr_env.parsed_env_name)
            self.message("")
            self.message("--- NOTICE: ENVVARS not set due to dry-run flag.")
            self.message("")

        if rval:
            msg = "ERROR: There was a problem configuring the environment."
            self.message(msg)
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
            "WORKSPACE",
            "CC",
            "CXX",
            "F77",
            "F90",
            "FC",
            "MODULESHOME"
            ]
        self.message("")
        tr_env.set_environment.pretty_print_envvars(envvar_filter=envvars_to_print)

        self.message("+" + "-"*68 + "+")
        self.message("|   E N V I R O N M E N T   S E T   U P   C O M P L E T E")
        self.message("+" + "-"*68 + "+")

        self.message("+" + "-"*68 + "+")
        self.message("|   G e n e r a t e   `packageEnables.cmake`   S T A R T I N G")
        self.message("+" + "-"*68 + "+")

        self.create_package_enables_file(dryrun=self.args.dry_run)

        self.message("+" + "-"*68 + "+")
        self.message("|   G e n e r a t e   `packageEnables.cmake`   C O M P L E T E D")
        self.message("+" + "-"*68 + "+")
        self.message("")

        return 0


    def execute_test(self):
        """
        Executes the test. This method must be overridden by a subclass.
        """
        raise NotImplementedError("This method must be overridden.")


    def chdir_logged(self, dest_dir, create_if_missing=False):
        """
        An extra verbose wrapper for chdir. This is helpful because it can provide
        a decent audit log of what paths have been changed. Optionally, you can also
        tell it to create a directory if the destination is missing.

        Args:
            dest_dir (str,path): The destination directory.
            create_if_missing (bool): If ``True`` then we attempt to recursively make
                the directory plus intermediate paths and then change to the dir.
        """
        self.message("--- CHDIR ---")
        self.message(f"current path: {os.getcwd()}")
        self.message(f"new path    : {dest_dir}")
        if not self.args.dry_run:
            try:
                os.chdir(self.arg_build_dir)
            except FileNotFoundError:
                self.message("WARNING: Path is missing!")
                self.message("         Will attempt to create the missing directory.")
                os.makedirs(dest_dir, exist_ok=True)
                os.chdir(dest_dir)
        else:
            self.message("Skipped (Dry-Run)")
        self.message("--- OK")
