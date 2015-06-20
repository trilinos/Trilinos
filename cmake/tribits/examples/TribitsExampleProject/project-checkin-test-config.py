#
# Define project-specific options for the checkin-test script for
# TribitsExampleProject.
#

configuration = {

    # Default command line arguments
    'defaults': {
        '--send-email-to-on-push': 'trilinos-checkin-tests@software.sandia.gov',
        },

    # CMake options (-DVAR:TYPE=VAL) cache variables.
    'cmake': {
        
        # Options that are common to all builds.
        'common': [],

        # Defines --default-builds, in order.
        'default-builds': [
            # Options for the MPI_DEBUG build.
            ('MPI_DEBUG', [
                '-DTPL_ENABLE_MPI:BOOL=ON',
                '-DCMAKE_BUILD_TYPE:STRING=RELEASE',
                '-DTribitsExProj_ENABLE_DEBUG:BOOL=ON',
                '-DTribitsExProj_ENABLE_CHECKED_STL:BOOL=ON',
                '-DTribitsExProj_ENABLE_DEBUG_SYMBOLS:BOOL=ON',
                ]),
            # Options for the SERIAL_RELEASE build.
            ('SERIAL_RELEASE', [
                '-DTPL_ENABLE_MPI:BOOL=OFF',
                '-DCMAKE_BUILD_TYPE:STRING=RELEASE',
                '-DTribitsExProj_ENABLE_DEBUG:BOOL=OFF',
                '-DTribitsExProj_ENABLE_CHECKED_STL:BOOL=OFF',
                ]),
            ], # default-builds

        }, # cmake

    } # configuration
