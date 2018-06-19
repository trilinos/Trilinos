# This file allows project-level configuration of the checkin-test system to
# set project options that are required for all developers. Machine or package
# specific options should not be placed in this file.

# This is a dictionary that specifies project-specific options for the
# checkin-test script that should be used by all developers. This
# includes default command line arguments that must be passed to
# checkin-test as well as settings for specific builds.
configuration = {

    # The default command line arguments that should be used by all
    # developers.
    'defaults': {
        '--send-email-to-on-push': 'trilinos-checkin-tests@software.sandia.gov',
        '--no-rebase' : '',
        },

    # CMake options for various build configurations. All entries in
    # this dictionary should be Python lists of -D arguments to cmake.
    'cmake': {
        
        # Default options that are common to all builds.
        'common': [
            # For graceful disables, we want to turn this on.
            '-DTrilinos_DISABLE_ENABLED_FORWARD_DEP_PACKAGES=ON',
            ],

        # Setup for the builds that should be run by default for a
        # standard checkin. This is a list of tuples so a preference
        # for build order can be expressed (e.g. if a project's
        # developers prefer one case to fail earlier than another).
        'default-builds': [

            ('MPI_RELEASE_DEBUG_SHARED_PT_OPENMP', [
                '-DTrilinos_CONFIGURE_OPTIONS_FILE:STRING=cmake/std/GCC-4.8.4-OpenMPI-1.10.1-MpiReleaseDebugSharedPtOpenMP.cmake',
                ]),

            ], # default-builds

        }, # cmake

    } # configuration
