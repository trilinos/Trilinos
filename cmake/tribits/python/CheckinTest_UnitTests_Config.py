configuration = {

    'defaults': {
        '--send-email-to-on-push': 'trilinos-checkin-tests@software.sandia.gov',
        },

    'cmake': {
        
        'common': [
            '-DTPL_ENABLE_Pthread:BOOL=OFF',
            '-DTPL_ENABLE_BinUtils:BOOL=OFF',
            ],
        'default-builds': [    
          # This is a list of tuples so we can preserve order.
          ('MPI_DEBUG', [
            '-DTPL_ENABLE_MPI:BOOL=ON',
            '-DCMAKE_BUILD_TYPE:STRING=RELEASE',
            '-DTrilinos_ENABLE_DEBUG:BOOL=ON',
            '-DTrilinos_ENABLE_CHECKED_STL:BOOL=ON',
            '-DTrilinos_ENABLE_DEBUG_SYMBOLS:BOOL=ON',
            '-DTrilinos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON',
            '-DTeuchos_ENABLE_DEFAULT_STACKTRACE:BOOL=OFF',
            ]),
          # Order matters here for the expected test output.
          ('SERIAL_RELEASE', [
            '-DTPL_ENABLE_MPI:BOOL=OFF',
            '-DCMAKE_BUILD_TYPE:STRING=RELEASE',
            '-DTrilinos_ENABLE_DEBUG:BOOL=OFF',
            '-DTrilinos_ENABLE_CHECKED_STL:BOOL=OFF',
            '-DTrilinos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=OFF',
            ]),
          ]
    },

}

