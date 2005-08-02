#! /usr/bin/env python

# ============================================================================
# Shows the usage of Epetra.PyComm().
# Since the pure virtual Epetra_Comm requires integer and double pointers
# (that are usually not defined within Python), you can replace the
# Epetra_Comm operations with the following:
#
# >>> Comm = PyComm()
# >>> target = Comm(<type>, source, root)
#
# where <type> is one of the following:
# - Epetra.MINALL
# - Epetra.MAXALL
# - Epetra.SUMALL
# - Epetra.SCANSUM
# - Epetra.BROADCAST
# and source is a Python's list.
# NOTE: all the elements in a list MUST have the same type. Also,
#       only integers and doubles are supported, since only these two
#       types are supported by Epetra_Comm.
#
# The root processor (specified by root) is required only by Epetra.BROADCAST,
# and ignored in all other cases.
#
#
# The output of this example is trivial for serial runs. If you have
# configured Trilinos with MPI support, you can try something like:
# $ mpirun -np 4 python ./exComm.py
#
# \author Marzio Sala, 9214
#
# \date Last updated on 01-Aug-05
# ============================================================================

# "from PyTrilinos import ..." syntax.  Here, the setpath module adds the build
# directory, including "PyTrilinos", to the front of the search path.  We thus
# use "import ..." for Trilinos modules.  This prevents us from accidentally
# picking up a system-installed version and ensures that we are testing the
# build module.
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
  base = Comm.MyPID()
  # Defines here some source variable.
  source = [1.0 * base, 2.0 * base, 3.0 * base]
  print "PE = ", base, ", source = ", source

  # get the mininum element
  target = Comm.GlobalOp(Epetra.MINALL, source)
  if Comm.MyPID() == 0:
    print "MINALL = ", target
  Comm.Barrier()
  
  # get the maximum element
  target = Comm.GlobalOp(Epetra.MAXALL, source)
  if Comm.MyPID() == 0:
    print "MAXALL = ", target
  Comm.Barrier()
  
  # sum all the elements
  target = Comm.GlobalOp(Epetra.SUMALL, source)
  if Comm.MyPID() == 0:
    print "SUMALL = ", target
  Comm.Barrier()
  
  # scansum
  target = Comm.GlobalOp(Epetra.SCANSUM, source)
  print "PE = ", base, ", SCANSUM = ", target
  Comm.Barrier()

  # broadcast from processor 0
  if Comm.MyPID() == 0:
    source = [10, 20]
  else:
    source = [0, 0]
  
  target = Comm.GlobalOp(Epetra.BROADCAST, source, 0)
  print "PE = ", base, ", BROADCAST = ", target

# This is a standard Python construct.  Put the code to be executed in a
# function [typically main()] and then use the following logic to call the
# function if the script has been called as an executable from the UNIX command
# line.  This also allows, for example, this file to be imported from a python
# debugger and main() called from there.
if __name__ == "__main__":
    main()
