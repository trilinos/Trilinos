configuration = {

    # The default command line arguments
    'defaults': {
        '--send-email-to-on-push': 'trilinos-checkin-tests@software.sandia.gov',
        '--enable-packages': 'TriBITS',
        '--no-rebase': '',
        '--no-append-test-results': ''
        },
    # NOTE: The last two options --no-rebase --no-append-test-results ensure
    # that the checkin-test.py script does *not* modify the SHA1s of the local
    # commits when pushing to the master repo.  This is important when git
    # commit SHA1s are snaphsotted into other TriBITS projects, for example,
    # in their <projectDir>/cmake/tribits/ directory.

    'cmake': {
        
        'common': [
            ],

        'default-builds': [

            # Options for the MPI_DEBUG build.
            ('MPI_DEBUG', [
                '-DTPL_ENABLE_MPI:BOOL=ON',
                '-DCMAKE_BUILD_TYPE:STRING=DEBUG',
                '-DTriBITSProj_ENABLE_DEBUG:BOOL=ON',
                ]),

            # Options for the SERIAL_RELEASE build.
            ('SERIAL_RELEASE', [
                '-DTPL_ENABLE_MPI:BOOL=OFF',
                '-DCMAKE_BUILD_TYPE:STRING=RELEASE',
                '-DTriBITSProj_ENABLE_DEBUG:BOOL=OFF',
                ]),

            ], # default-builds

        }, # cmake

    } # configuration
