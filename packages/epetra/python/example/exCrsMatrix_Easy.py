#! /usr/bin/env python

# Shows how to build a serial or distributed Epetra_CrsMatrix in python, 
# the easy, MATLAB-like way.

try:
  import setpath
  import Epetra
except:
  from PyTrilinos import Epetra
  print "Using system-installed Epetra"

def main():

  # define the communicator (Serial or parallel, depending on your configure
  # line), then initialize a distributed matrix of size 4. The matrix is empty,
  # `0' means to allocate for 0 elements on each row (better estimates make the
  # code faster). `NumMyElements' is the number of rows that are locally hosted 
  # by the calling processor; `MyGlobalElements' is the global ID of locally 
  # hosted rows.
  Comm              = Epetra.PyComm()
  NumGlobalElements = 4
  Map               = Epetra.Map(NumGlobalElements, 0, Comm)
  Matrix            = Epetra.CrsMatrix(Epetra.Copy, Map, 0)
  NumMyElements     = Map.NumMyElements()
  MyGlobalElements  = Map.MyGlobalElements()
  
  # Populates the matrix, element-by-element. Note that you *must* use global
  # IDs for rows and columns. However, you can only set *local* rows. When the
  # matrix is completed, call FillComplete() to freeze its structure. Matrix
  # elements cannot be added after FillComplete(); however, existing elements
  # can still be modified.
  for i in MyGlobalElements:
    if i > 0:
      Matrix[i, i - 1] = -1
    if i < NumGlobalElements - 1:
      Matrix[i, i + 1] = -1
    Matrix[i, i] = 2.
  
  Matrix.FillComplete()
  
  # Prints matrix diagonal elements. Two important notes:
  # - you must first call FillComplete()
  # - if you ask for a non-valid element, the returned value is 0.0
  # - only elements in local rows can be queried.
  for i in MyGlobalElements:
    print "PE%d: Matrix(%d, %d) = %e" %(Comm.MyPID(), i, i, Matrix[i, i])
  
  Comm.Barrier()
  
  # Prints the nonzero elements of the matrix, row-by-row
  # 
  for i in MyGlobalElements:
    Indices, Values = Matrix[i]
    for j in xrange(len(Indices)):
      print "PE%d: Matrix(%d, %d) = %e" % (Comm.MyPID(), i, Indices[j], Values[j])

  if Comm.MyPID() == 0: print "End Result: TEST PASSED"

# This is a standard Python construct.  Put the code to be executed in a
# function [typically main()] and then use the following logic to call the
# function if the script has been called as an executable from the UNIX
# command
# line.  This also allows, for example, this file to be imported from a python
# debugger and main() called from there.
if __name__ == "__main__":
  main()
