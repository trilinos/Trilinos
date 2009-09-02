#!/bin/csh
# ************************************************************************
# 
#                                 AMESOS
#                 Copyright (2004) Sandia Corporation
# 
# Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
# license for use of this work by or on behalf of the U.S. Government.
# 
# This library is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation; either version 2.1 of the
# License, or (at your option) any later version.
#  
# This library is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#  
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
# USA
# Questions? Contact Jonathan Hu or Ray Tuminaro ({jhu,rstumin}@sandia.gov.
# 
# ************************************************************************
#
## NOTE: Those wishing to cusomize this script to run test exe's
## that have already been autotool'ed should read lines beginning with '##'

# $1 - Used only for automated testing.  No action required by script owner.
#      This parameter names the subdirectory from which the script is being
#      run.  This assists a developer in figuring out which tests failed.
# $2 - Indicates if the test is an automated nightly test.  No action required
#      by script owner.

set error = None
set AnError = False
set printexitvalue
if( "$2" == "True" ) then # $2 is an optional parameter indicating if 
			  # this is an automated test or not
    # file2 is the log that is created and put into a results email if 
    # errors occur.
    set file2 = ../logMpiErrors.txt
    rm -f $file2
    # Echo some important information into the log file to help developers
    # figure out which tests failed.
    #'file' is a shorter log that is retained even if all tests pass.
    set file = ../log`eval uname`.txt
    rm -f $file
## IMPORTANT: Specify the script owner(s) on the following line
## For one owner type "owner@abc.com", for multiple owners
## "owner1@abc.com, owner2@def.com"
    echo "amesos-regression@software.sandia.gov" >>& $file
    echo "Script owner(s) is listed on the previous line." >>& $file
## List the Trilinos package being tested on the following line
    echo "Package being tested: ML  " >>& $file
    echo "Name of subdirectory: " $1 >>& $file
else
    set file = log_mpi_`eval date +%d%b%Y_%H%M%S`
    rm -f $file
endif
echo $file
echo "Date: " `eval date` >>& $file
echo `uname -a` >>& $file
## Different directory structures will require different setups.
## This file assumes a structure like that of epetra - exe's live in 
## a direct subdirectory of 'epetra/test' 

## Keep in mind that file and file2-4 live in 'package_name/test'
## Also, 'package_name/test' is the current directory
## It is recommended that all directory changing be done relative to
## the current directory because scripts live in the source directory,
## but are invoked from various build directories

set mpigo = `printenv TRILINOS_TEST_HARNESS_MPIGO_COMMAND`

if ("$mpigo" == "") then
    set mpigo = "mpirun -np "
endif

## List the subdirectories of 'test' containing test exe's in the foreach loop
## if directory structure is like that of epetra.
cd ../example
touch fakefile.exe
set tt = (*.exe)
# get rid of spaces
set exefiles = `echo ${tt} | sed "s/ //g"`
# avoids problem of no executables in subdirectory
if ( ${exefiles} != 'fakefile.exe' ) then
      /bin/rm -f fakefile.exe
	  foreach g(*.exe)
                  echo "[TEST]" >>& $file
		  echo "[TEST] ############" $g "##############" >>& $file
                  echo "[TEST]" >>& $file
		  if( "$2" == "True" ) then
                    # ================== #
                    # run with 1 process #		    
                    # ================== #
		    $mpigo 1 ./$g >>& /dev/null
		     if( $status != 0 ) then
		        # A test failed.
			    set AnError = True
			    echo "[TEST] ******** Test w/ 1 proc failed ********" >>& $file
                            echo "[TEST] ******** now re-running with output" >>& $file
                            $mpigo 2 ./$g >>& $file
			    echo "[TEST] Errors for script " $g " are listed above." >>& $file2
		      else
		        # Tests passed
			    echo "[TEST] ******** Test w/ 1 proc passed ********" >>& $file
		      endif
                    # ==================== #
                    # run with 4 processes #		    
                    # ==================== #
		    $mpigo 4 ./$g >>& /dev/null
                    # run with 1 process		    
		     if( $status != 0 ) then
		        # A test failed.
			    set AnError = True
			    echo "[TEST] ******** Test w/ 4 proc failed ********" >>& $file
                            echo "[TEST] ******** now re-running with output" >>& $file
                            $mpigo 2 ./$g >>& $file
			    echo "[TEST] Errors for script " $g " are listed above." >>& $file2
		      else
		        # Tests passed
			    echo "[TEST] ******** Test w/ 4 proc passed ********" >>& $file
		      endif
		  else
		    # This is not an automated test.
		    $mpigo 1 ./$g >>& $file
		    $mpigo 4 ./$g >>& $file
		  endif
	  end
else
  echo "[TEST] *** no examples found ***" >>& $file
  echo "[TEST] *** Perhaps a package the test depends on didn't build? ****" >>& $file
endif
/bin/rm -f fakefile.exe

# copy Summary file and Error file to standard out
if ( "$3" == "True" ) then
    echo "@#@#@#@#  Summary file @#@#@#@#@"
    cat $file
    if( "$AnError" == "True" ) then
        echo "@#@#@#@# Error file - output from same tests run with verbose enabled @#@#@#@#@"
	cat $file2
    endif
endif

## At this point, it is assumed that the current directory is
## 'package_name/test'
if ( "$2" == "True" ) then
    rm $file
    rm -f $file2
endif

if ( "$AnError" == "True" ) then
    exit 1
else
    echo "End Result: TEST PASSED"
    exit 0
endif

