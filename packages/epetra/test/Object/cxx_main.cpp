//@HEADER
// ************************************************************************
// 
//               Epetra: Linear Algebra Services Package 
//                 Copyright 2011 Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
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
#include "Epetra_Version.h"
                                            
int main(int argc, char *argv[]) {

  int ierr = 0;
#ifdef EPETRA_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm comm( MPI_COMM_WORLD );

#else
  Epetra_SerialComm comm;
#endif

  bool verbose = false;

  // Check if we should print results to standard out
  if (argc>1) if (argv[1][0]=='-' && argv[1][1]=='v') verbose = true;

  if (verbose && comm.MyPID()==0)
    cout << Epetra_Version() << endl << endl;

  // I'm alive !!!
  if (verbose) cout << comm <<endl;

  Epetra_Object obj;

  // Test Epetra_Object label and the method to get the label attribute
  const char* ObjLabel = obj.Label();
  const char* ObjLabel1 = "Epetra::Object";
  if (verbose) cout << endl << endl << "This should say " << ObjLabel1 << ": " << ObjLabel << endl << endl << endl;
  EPETRA_TEST_ERR(strcmp(ObjLabel1,ObjLabel),ierr);

  // Test Epetra_Object SetLabel attribute set method
  const char* NewObjLabel = "New name for Epetra_Object";
  obj.SetLabel(NewObjLabel);
  const char* NewObjLabel1 = obj.Label(); 
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
  const char* NewObjLabel2 = obj2.Label();
  if (verbose) cout << endl << endl << "This should say " << NewObjLabel << ": " << NewObjLabel2 << endl << endl << endl;
  EPETRA_TEST_ERR(strcmp(NewObjLabel2,NewObjLabel),ierr);

  Epetra_Object obj3(NewObjLabel,1); // pass a label and a TracebackMode
  const char* NewObjLabel3 = obj3.Label();
  int TbM1 = obj3.GetTracebackMode();
  if (verbose) cout << endl << "This should say " << NewObjLabel << "," << "1: " << NewObjLabel3 << "," << TbM1 << endl << endl << endl;
  EPETRA_TEST_ERR(strcmp(NewObjLabel3,NewObjLabel),ierr);
  EPETRA_TEST_ERR(!(1==TbM1),ierr);
  
  Epetra_Object obj4(obj3); // copy constructor
  const char* NewObjLabel4 = obj4.Label();
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
