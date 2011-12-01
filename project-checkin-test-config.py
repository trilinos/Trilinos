# This file allows project-level configuration of the checkin-test system to
# set project options that are required for all developers. Machine or package
# specific options should not be placed in this file.

# This is a dictionary with key-value pairs that map to the script's
# command line arguments. These are used to add project-specific
# defaults for the arguments that all developers should use.
configuration = {
    '--extra-cmake-options': '-DTPL_ENABLE_Pthread:BOOL=OFF -DTPL_ENABLE_BinUtils:BOOL=OFF',
    '--send-email-to-on-push': 'trilinos-checkin-tests@software.sandia.gov',
}

