//@HEADER
// ************************************************************************
// 
//               Epetra: Linear Algebra Services Package 
//                 Copyright (2001) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER

// Epetra_Object Test routine

#include "Epetra_Object.h"
#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#include <mpi.h>
#else
#include "Epetra_SerialComm.h"
#endif

#include "../epetra_test_err.h"
int main(int argc, char *argv[]) {

  int ierr = 0;
#ifdef EPETRA_MPI

  // Initialize MPI

  MPI_Init(&argc,&argv);
  int size, rank; // Number of MPI processes, My process ID

  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  Epetra_MpiComm comm( MPI_COMM_WORLD );

#else

  int size = 1; // Serial case (not using MPI)
  int rank = 0;

  Epetra_SerialComm comm;

#endif

  bool verbose = false;

  // Check if we should print results to standard out
  if (argc>1) if (argv[1][0]=='-' && argv[1][1]=='v') verbose = true;


  // I'm alive !!!
  if (verbose) cout << comm <<endl;

  Epetra_Object obj;

  // Test Epetra_Object label and the method to get the label attribute
  char* ObjLabel = obj.Label();
  const char* ObjLabel1 = "Epetra::Object";
  if (verbose) cout << endl << endl << "This should say " << ObjLabel1 << ": " << ObjLabel << endl << endl << endl;
  EPETRA_TEST_ERR(strcmp(ObjLabel1,ObjLabel),ierr);

  // Test Epetra_Object SetLabel attribute set method
  const char* NewObjLabel = "New name for Epetra_Object";
  obj.SetLabel(NewObjLabel);
  char* NewObjLabel1 = obj.Label(); 
  if (verbose) cout << endl << "This should say " << NewObjLabel << ": " << NewObjLabel1 << endl << endl << endl;
  EPETRA_TEST_ERR(strcmp(NewObjLabel1,NewObjLabel),ierr);

  // Test GetRacebackMode and SetTracebackMode methods
  EPETRA_TEST_ERR(!(obj.GetTracebackMode()==DefaultTracebackMode),ierr);
  if (verbose) cout << endl <<"Default Traceback Mode value = " << obj.GetTracebackMode() << endl;

  obj.SetTracebackMode(DefaultTracebackMode-1);
  if (verbose) cout << "Set Traceback Mode value to one less than default = " << obj.GetTracebackMode() << endl << endl;
  Epetra_Object obj0;
  EPETRA_TEST_ERR(!(obj0.GetTracebackMode()==DefaultTracebackMode-1),ierr);

  // Test constructors other than the default
  Epetra_Object obj1(1); // pass only TracebackMode
  int TbM = obj1.GetTracebackMode();
  if (verbose) cout << endl << endl << "This should say 1: " << TbM << endl << endl;
  EPETRA_TEST_ERR(!(1==TbM),ierr);

  Epetra_Object obj2(NewObjLabel); // pass only a label
  char* NewObjLabel2 = obj2.Label();
  if (verbose) cout << endl << endl << "This should say " << NewObjLabel << ": " << NewObjLabel2 << endl << endl << endl;
  EPETRA_TEST_ERR(strcmp(NewObjLabel2,NewObjLabel),ierr);

  Epetra_Object obj3(NewObjLabel,1); // pass a label and a TracebackMode
  char* NewObjLabel3 = obj3.Label();
  int TbM1 = obj3.GetTracebackMode();
  if (verbose) cout << endl << "This should say " << NewObjLabel << "," << "1: " << NewObjLabel3 << "," << TbM1 << endl << endl << endl;
  EPETRA_TEST_ERR(strcmp(NewObjLabel3,NewObjLabel),ierr);
  EPETRA_TEST_ERR(!(1==TbM1),ierr);
  
  Epetra_Object obj4(obj3); // copy constructor
  char* NewObjLabel4 = obj4.Label();
  int TbM2 = obj4.GetTracebackMode();
  if (verbose) cout << endl << "This should say " << NewObjLabel << "," << "1: " << NewObjLabel4 << "," << TbM2 << endl << endl << endl;
  EPETRA_TEST_ERR(strcmp(NewObjLabel4,NewObjLabel),ierr);
  EPETRA_TEST_ERR(!(1==TbM2),ierr);
  

#ifdef EPETRA_MPI
  MPI_Finalize();
#endif
  return ierr;
}

/*
  end of file main.cc
*/
