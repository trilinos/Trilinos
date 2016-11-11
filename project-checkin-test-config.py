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
        },

    # CMake options for various build configurations. All entries in
    # this dictionary should be Python lists of -D arguments to cmake.
    'cmake': {
        
        # Default options that are common to all builds.
        'common': [
            # Shared libs safes a ton of disk space and catches more errors
            '-DBUILD_SHARED_LIBS=ON',
            # For graceful disables, we want to turn this on
            '-DTrilinos_DISABLE_ENABLED_FORWARD_DEP_PACKAGES=ON',
            # We want to see tracing of added tests to help in debugging
            # problems.
            '-DTrilinos_TRACE_ADD_TEST=ON',
            ],

        # Setup for the builds that should be run by default for a
        # standard checkin. This is a list of tuples so a preference
        # for build order can be expressed (e.g. if a project's
        # developers prefer one case to fail earlier than another).
        'default-builds': [

            # Options for the MPI_DEBUG build.
            ('MPI_RELEASE_DEBUG_SHARED', [
                '-DTPL_ENABLE_MPI=ON',
                '-DCMAKE_BUILD_TYPE=RELEASE',
                '-DTrilinos_ENABLE_DEBUG=ON',
                '-DBUILD_SHARED_LIBS=ON',
                '-DTrilinos_ENABLE_DEBUG_SYMBOLS=ON',
                '-DTrilinos_ENABLE_CI_TEST_MODE=ON',
                '-DTrilinos_ENABLE_EXPLICIT_INSTANTIATION=ON',
                '-DTrilinos_ENABLE_SECONDARY_TESTED_CODE=OFF',
                '-DTrilinos_ENABLE_TESTS=ON',
                '-DTeuchos_ENABLE_DEFAULT_STACKTRACE=OFF',
                ]),

#            # Options for the SERIAL_RELEASE build.
#            ('SERIAL_RELEASE_SHARED', [
#                '-DTPL_ENABLE_MPI:BOOL=OFF',
#                '-DCMAKE_BUILD_TYPE:STRING=RELEASE',
#                '-DTrilinos_ENABLE_DEBUG:BOOL=OFF',
#                '-DBUILD_SHARED_LIBS=ON',
#                '-DTrilinos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=OFF',
#                ]),

            ], # default-builds

        }, # cmake

    } # configuration
