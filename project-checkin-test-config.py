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
        
        # Options that are common to all builds.
        'common': [
            '-DTPL_ENABLE_Pthread:BOOL=OFF',
            '-DTPL_ENABLE_BinUtils:BOOL=OFF',
            ],

        # Setup for the builds that should be run by default for a
        # standard checkin. This is a list of tuples so a preference
        # for build order can be expressed (e.g. if a project's
        # developers prefer one case to fail earlier than another).
        'default-builds': [

            # Options for the MPI_DEBUG build.
            ('MPI_DEBUG', [
                '-DTPL_ENABLE_MPI:BOOL=ON',
                '-DCMAKE_BUILD_TYPE:STRING=RELEASE',
                '-DTrilinos_ENABLE_DEBUG:BOOL=ON',
                '-DTrilinos_ENABLE_CHECKED_STL:BOOL=ON',
                '-DTrilinos_ENABLE_DEBUG_SYMBOLS:BOOL=ON',
                '-DTrilinos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON',
                '-DTeuchos_ENABLE_DEFAULT_STACKTRACE:BOOL=OFF',
                ]),

            # Options for the SERIAL_RELEASE build.
            ('SERIAL_RELEASE', [
                '-DTPL_ENABLE_MPI:BOOL=OFF',
                '-DCMAKE_BUILD_TYPE:STRING=RELEASE',
                '-DTrilinos_ENABLE_DEBUG:BOOL=OFF',
                '-DTrilinos_ENABLE_CHECKED_STL:BOOL=OFF',
                '-DTrilinos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=OFF',
                ]),

            ], # default-builds

        }, # cmake

    } # configuration

