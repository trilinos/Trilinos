In this directory, are scripts for running checkin-test.py for the developers of Zoltan2.  
The hope is that by using these scripts, the overhead of using checkin-test will be reduced.  
These scripts represent varying degrees of thoroughness (and thus varying amounts of time) 
in the testing process.


First steps:
------------

1. Create directory CHECKIN in Trilinos source directory (trunk) that you are modifying.
2. Copy files MPI_ss.config and SERIAL_ss.config to CHECKIN directory.  These are necessary
   since Zoltan2 is secondary stable and checkin-test.py won't build secondary stable packages
   by default.  These config files have -DTrilinos_ENABLE_SECONDARY_STABLE_CODE=ON set to turn
   on secondary stable.  
3. Modify MPI_ss.config (in CHECKIN) so the MPI path points to your MPI version.


Description of scripts:
----------------------

checkinZoltan2-fast.sh :
   
   Fastest version.  Only enables and tests Zoltan2.  Disables packages that depend on Zoltan.
   Appropriate if you are pushing very minor changes that can't possibly break the rest of Trilinos
   (e.g., if you are working in an experimental section of code that is protected by #ifdefs).
   Otherwise, it is more appropriate to use a different script that tests more thoroughly.		      

 
checkinZoltan2-medium.sh  
   Medium version: not the fastest, not the slowest.  Reasonable compromise between thoroughness
   of testing and the time to completion.



checkinZoltan2-complete.sh
   Slowest version.  Enables and tests Zoltan2 and all packages that depend (optional or required)
   on Zoltan2.  Appropriate perhaps if you are pushing substantial changes that could affect
   many packages that depend on Zoltan 2.  Otherwise, perhaps overkill.

  
checkinZoltan2-push.sh
   Script to push changes after checkin-test has been run.  This copies a log of the testing that
   has occurred into the git message log.



Suggested TPLs installed on common development machines:
--------------------------------------------------------

