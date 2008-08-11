#!/bin/bash

#
# Main driver script for all of of the automated cron jobs run on
# gabriel.sandia.gov
#

# Here is the crontab script that I am currently using with this:
#
#    # (Use to post in the top of your crontab)
#    # ----------------- minute (0 - 59)
#    # |  -------------- hour (0 - 23)
#    # |  |  ----------- day of month (1 - 31)
#    # |  |  |  -------- month (1 - 12)
#    # |  |  |  |  ----- day of week (0 - 7) (Sunday=0 or 7)
#    # |  |  |  |  |
#    # *  *  *  *  *  command to be executed
#     00 00  *  *  0,2,4,6 /mnt/disk2/rabartl/Trilinos.nightly-tests/Trilinos/commonTools/test/harness/drivers/gabriel-crontab-driver.sh serial
#     00 00  *  *  1,3,5 /mnt/disk2/rabartl/Trilinos.nightly-tests/Trilinos/commonTools/test/harness/drivers/gabriel-crontab-driver.sh mpi
#

# In this crontab setup, we run the serial debug on sunday, tuesday, thursday,
# and saturday and we run mpi optimized on monday, wednesday and friday.

TRILINOS_BUILD=$1; shift

SKIP_CHARON=$1; shift

if [ "$SKIP_CHARON" != "skip_charon" ]; then

  # Charon + Trilinos Integration tests
  /mnt/disk2/rabartl/Charon.nightly-tests-4/Utils/prj-management/asc-level-2-milestone/charon-tridev-nightly-tests-on-gabriel.sh \
    &> /mnt/disk2/rabartl/Charon.nightly-tests-4/charon-tridev-nightly-tests-on-gabriel.out

fi

# Trilinos specific tests
/mnt/disk2/rabartl/Trilinos.nightly-tests/Trilinos/commonTools/test/harness/drivers/gabriel-nightly-tests.sh $TRILINOS_BUILD \
  &> /mnt/disk2/rabartl/Trilinos.nightly-tests/gabriel-nightly-tests.out 

# ToDo: Build doxygen documentation for Trilinos given what is in:
#
#   /mnt/disk2/rabartl/Trilinos.nightly-tests/Trilinos
#
# and host that from gabriel.sandia.gov!
