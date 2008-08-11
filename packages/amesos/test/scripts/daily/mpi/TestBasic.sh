#!/bin/csh
# ************************************************************************
# 
#                 Amesos: Direct Sparse Solver Package
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
# Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
# 
# ************************************************************************
#
#  The only change that I made to change this from dscpack to klu 
#  was to replace the single occurence of dscpack with klu
#  I also added HAVE_AMESOS_KLU to Makefile.am in this directory.
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
    echo "Package being tested: Amesos  " >>& $file
    echo "Name of subdirectory: " $1 >>& $file
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

cd Test_Basic
set TestRan = False
foreach g  ( TestBasic.csh )
                set TestRan = True
		echo "############" $g " ##############" >>& ../$file
		if( "$2" == "True" ) then
                    echo "############ KLU ##############" >>& ../$file
		    source $g KLU
		    if( $status != 0 ) then
		    # A test failed.
			set AnError = True
			echo "  ******** Test failed ********" >>& ../$file
			echo "Errors for script " $g " are listed above." >>& ../$file2
			cat SST.summary >> ../$file
		    else
		    # Tests passed
			echo "******** Test passed ********" >>& ../$file
		    endif
		endif
end
if ( $TestRan != "True" ) then
  echo "No executables were found  " 
  set AnError = True
endif

cd ..

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

