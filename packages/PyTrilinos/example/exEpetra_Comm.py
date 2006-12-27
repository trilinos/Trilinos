#! /usr/bin/env python

# @HEADER
# ************************************************************************
# 
#                PyTrilinos: Python Interface to Trilinos
#                   Copyright (2005) Sandia Corporation
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
# Questions? Contact Bill Spotz (wfspotz@sandia.gov) 
# 
# ************************************************************************
# @HEADER

# ============================================================================
# Shows the usage of Epetra.PyComm().
#
# The output of this example is trivial for serial runs. If you have
# configured Trilinos with MPI support, you can try something like:
# $ mpirun -np 4 python ./exComm.py
#
# \author Marzio Sala, 9214; Bill Spotz, 1433
#
# \date Last updated on 15-Oct-05
# ============================================================================

# "from PyTrilinos import ..." syntax.  Here, the setpath module adds the build
# directory, including "PyTrilinos", to the front of the search path.  We thus
# use "import ..." for Trilinos modules.  This prevents us from accidentally
# picking up a system-installed version and ensures that we are testing the
# build module.
from numpy    import *
from optparse import *

parser = OptionParser()
parser.add_option("-t", "--testharness", action="store_true",
                  dest="testharness", default=False,
                  help="test local build modules; prevent loading system-installed modules")
parser.add_option("-v", "--verbosity", type="int", dest="verbosity", default=2,
                  help="set the verbosity level [default 2]")
options,args = parser.parse_args()
if options.testharness:
    import setpath
    import Epetra
else:
    try:
        import setpath
        import Epetra
    except ImportError:
        from PyTrilinos import Epetra
        print "Using system-installed Epetra"

def main():

    # Defines a communicator, which will be an Epetra.SerialComm or
    # an Epetra.MpiComm depending on how Trilinos was configured
    Comm = Epetra.PyComm()
    base = Comm.MyPID() + 1
    # Defines here some source variable.
    source = [1.0*base, 2.0*base, 3.0*base]
    print "PE = ", base, ", source = ", source

    # get the mininum element
    #target = Comm.GlobalOp(Epetra.MINALL, source)
    target = Comm.MinAll(source)
    if Comm.MyPID() == 0:
        print "MINALL = ", target
    Comm.Barrier()

    # get the maximum element
    target = Comm.MaxAll(source)
    if Comm.MyPID() == 0:
        print "MAXALL = ", target
    Comm.Barrier()

    # sum all the elements
    target = Comm.SumAll(source)
    if Comm.MyPID() == 0:
        print "SUMALL = ", target
    Comm.Barrier()

    # scansum
    target = Comm.ScanSum(source)
    print "PE = ", base, ", SCANSUM = ", target
    Comm.Barrier()

    # broadcast from processor 0
    if Comm.MyPID() == 0:
        source = [10, 20]
    else:
        source = [0, 0]

    #target = Comm.GlobalOp(Epetra.BROADCAST, source, 0)
    aSource = array(source,dtype='i')
    Comm.Broadcast(aSource,0)
    print "PE = ", base, ", BROADCAST = ", aSource

    # synchronize processors
    Comm.Barrier()

    if Comm.MyPID() == 0: print "End Result: TEST PASSED"

# This is a standard Python construct.  Put the code to be executed in a
# function [typically main()] and then use the following logic to call the
# function if the script has been called as an executable from the UNIX command
# line.  This also allows, for example, this file to be imported from a python
# debugger and main() called from there.
if __name__ == "__main__":
    main()
