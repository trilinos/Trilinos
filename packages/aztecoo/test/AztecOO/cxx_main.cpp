//@HEADER
// ***********************************************************************
//
//        AztecOO: An Object-Oriented Aztec Linear Solver Package
//                 Copyright (2002) Sandia Corporation
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
// ***********************************************************************
//@HEADER

#include "AztecOO.h"
#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#include "mpi.h"
#endif
#include "Epetra_SerialComm.h"
#include "Epetra_Map.h"
#include "Epetra_BlockMap.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_LinearProblem.h"
#include "AztecOO_StatusTestMaxIters.h"
#include "AztecOO_StatusTestResNorm.h"
#include "AztecOO_StatusTestCombo.h"
#include "AztecOO_Operator.h"
#include "AztecOO_Version.h"

bool argument_is_present(const char* argument,
                         int argc,
                         char** argv);

double resid2norm(Epetra_CrsMatrix& A,
                  Epetra_Vector& x,
                  Epetra_Vector& b);

int test_azoo_as_precond_op(Epetra_CrsMatrix& A,
                            Epetra_Vector& x,
                            Epetra_Vector& b,
                            bool verbose);

int test_azoo_with_ilut(Epetra_CrsMatrix& A,
                        Epetra_Vector& x,
                        Epetra_Vector& b,
                        bool verbose);

int main(int argc, char *argv[])
{
#ifdef EPETRA_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm comm;
#endif

  int numprocs = comm.NumProc();
  int localproc = comm.MyPID();

  ///////////////////////////////////////////////////
  //First figure out whether verbose output is on.
  //(and if so, only turn it on for proc 0)

  bool verbose = false;
  if (localproc == 0) {
    verbose = argument_is_present("-v", argc, argv);
  }

  ///////////////////////////////////////////////////

  int local_n = 20;
  int global_n = numprocs*local_n;

  Epetra_Map emap(global_n, 0, comm);

  Epetra_CrsMatrix A(Copy, emap, 3);
  Epetra_Vector x(emap), b(emap);

  x.PutScalar(1.0);

  int myFirstGlobalRow = localproc*local_n;
  int globalCols[3];
  double values[3];

  for(int i=0; i<local_n; ++i) {
    int globalRow = myFirstGlobalRow +i;

    int numcols = 0;
    if (globalRow > 0) {
      globalCols[numcols] = globalRow-1;
      values[numcols++] = -1.0;
    }

    globalCols[numcols] = globalRow;
    values[numcols++] = 4.0;

    if (globalRow < global_n-1) {
      globalCols[numcols] = globalRow+1;
      values[numcols++] = -1.0;
    }

    A.InsertGlobalValues(globalRow, numcols, values, globalCols);
  }

  A.FillComplete();

  A.Multiply(false, x, b);
  x.PutScalar(0.0);

  double initial_norm = resid2norm(A, x, b);

  if (verbose) {
    cout << "Initial 2-norm of b-A*x: "<<initial_norm<<endl;
  }

  int err = test_azoo_as_precond_op(A, x, b, verbose);
  if (err != 0) {
    cout << "test_azoo_as_precond_op err, test FAILED."<<endl;
    return(err);
  }

  err = test_azoo_with_ilut(A, x, b, verbose);
  if (err != 0) {
    cout << "test_azoo_with_ilut err, test FAILED."<<endl;
    return(err);
  }

  cout << "********* Test passed **********" << endl;

#ifdef EPETRA_MPI
  MPI_Finalize() ;
#endif

  return(0);
}

bool argument_is_present(const char* argument,
                         int argc,
                         char** argv)
{
  if (argument == NULL || argc < 1) return(false);

  for(int i=0; i<argc; ++i) {
    if (strcmp(argument, argv[i]) == 0) {
      return(true);
    }
  }
  return(false);
}

double resid2norm(Epetra_CrsMatrix& A,
                  Epetra_Vector& x,
                  Epetra_Vector& b)
{
  Epetra_Vector r(x);
  A.Multiply(false, x, r);

  r.Update(1.0, b, -1.0);//r  =  b - r  =  b - A*x

  double nrm2;
  r.Norm2(&nrm2);

  return(nrm2);
}

int test_azoo_as_precond_op(Epetra_CrsMatrix& A,
                            Epetra_Vector& x,
                            Epetra_Vector& b,
                            bool verbose)
{
  AztecOO* azoo0 = new AztecOO(&A, &x, &b);

  azoo0->SetAztecOption(AZ_solver, AZ_cg);
  azoo0->SetAztecOption(AZ_output, AZ_none);
  azoo0->SetAztecOption(AZ_conv, AZ_none);

  AztecOO_Operator azooOp(azoo0, 10);

  Epetra_LinearProblem elp(&A, &x, &b);

  x.PutScalar(0.0);
  AztecOO* azoo1 = new AztecOO(elp);
  azoo1->SetPrecOperator(&azooOp);
  azoo1->SetAztecOption(AZ_solver, AZ_gmres);
  azoo1->SetAztecOption(AZ_output, AZ_none);
  azoo1->SetAztecOption(AZ_conv, AZ_none);

  if (verbose) {
    cout << "testing recursive solve (AztecOO as"
    << " preconditioner for another AztecOO)."<<endl;
  }

  int maxiters = 100;
  double tolerance = 1.e-12;
  azoo1->Iterate(maxiters, tolerance);

  double resid1 = resid2norm(A, x, b);
  if (verbose) {
    cout << "residual 2-norm after recursive solve: "
        << resid1 << endl;
  }

  if (resid1 > 1.e-6) {
    return(-1);
  }

  if (verbose) {
    cout << "now make sure the precond AztecOO instance"
      << " hasn't been corrupted."<<endl;
  }

  x.PutScalar(0.0);
  azoo0->Iterate(&A, &x, &b, maxiters, tolerance);

  double resid0 = resid2norm(A, x, b);
  if (verbose) {
    cout << "residual 2-norm: " << resid0 << endl;
  }
  if (resid0 > 1.e-6) {
    return(-1);
  }

  delete azoo0;
  delete azoo1;

  return(0);
}

int test_azoo_with_ilut(Epetra_CrsMatrix& A,
                        Epetra_Vector& x,
                        Epetra_Vector& b,
                        bool verbose)
{
  x.PutScalar(0.0);

  AztecOO* azoo0 = new AztecOO(&A, &x, &b);

  azoo0->SetAztecOption(AZ_solver, AZ_gmres);
  azoo0->SetAztecOption(AZ_output, AZ_none);
  azoo0->SetAztecOption(AZ_conv, AZ_none);
  azoo0->SetAztecOption(AZ_precond, AZ_dom_decomp);
  azoo0->SetAztecOption(AZ_subdomain_solve, AZ_ilut);
  azoo0->SetAztecOption(AZ_keep_info, 1);

  if (verbose) {
    cout << "testing AztecOO with GMRES and ILUT, AZ_keep_info==1" << endl;
  }

  int maxiters = 100;
  double tolerance = 1.e-12;
  int err = azoo0->Iterate(maxiters, tolerance);
  if (err != 0) {
    if (verbose) cout << "AztecOO::Iterate return err="<<err<<endl;
    return(err);
  }

  double resid = resid2norm(A, x, b);
  if (verbose) {
    cout << "residual 2-norm after GMRES/ILUT solve: "
        << resid << endl;
  }

  if (resid > 1.e-6) {
    return(-1);
  }

  if (verbose) {
    cout << "solving with GMRES/ILUT again, AZ_pre_calc==AZ_reuse"
        << endl << "(will error out if factors weren't kept from"
         << " previous solve)"<<endl;
  }

  azoo0->SetAztecOption(AZ_pre_calc, AZ_reuse);
  x.PutScalar(0.0);
  err = azoo0->Iterate(maxiters, tolerance);
  if (err != 0) {
    if (verbose) cout << "AztecOO::Iterate return err="<<err<<endl;
    return(err);
  }

  double resid2 = resid2norm(A, x, b);
  if (verbose) {
    cout << "after second GMRES/ILUT solve, residual 2-norm: "
      << resid2 <<endl;
  }
  if (resid2 > 1.e-6) {
    return(-1);
  }

  delete azoo0;

  return(0);
}

