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

// Epetra_BlockMap Test routine

#include "EpetraExt_Version.h"

#include "Epetra_CrsMatrix.h"
#include "Epetra_VbrMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_LocalMap.h"
#include "Epetra_IntVector.h"
#include "Epetra_Map.h"

#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#include <mpi.h>
#else
#include "Epetra_SerialComm.h"
#endif

#include "Epetra_Time.h"
#include "EpetraExt_Transpose_RowMatrix.h"
#include "Trilinos_Util.h"

int checkResults(Epetra_RowMatrix * A, Epetra_CrsMatrix * transA, 
                 Epetra_Vector * xexact, bool verbose);

int main(int argc, char *argv[])
{
  int i;

#ifdef EPETRA_MPI
  // Initialize MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm comm;
#endif

  // Uncomment to debug in parallel int tmp; if (comm.MyPID()==0) cin >> tmp; comm.Barrier();

  bool verbose = false;
  bool veryVerbose = false;

  // Check if we should print results to standard out
  if (argc>1) if (argv[1][0]=='-' && argv[1][1]=='v') verbose = true;

  if (!verbose) comm.SetTracebackMode(0); // This should shut down any error traceback reporting

  if (verbose) cout << comm << endl << flush;

  bool verbose1 = verbose;
  if (verbose) verbose = (comm.MyPID()==0);

  if (verbose)
    cout << EpetraExt::EpetraExt_Version() << endl << endl;

  int nx = 128;
  int ny = comm.NumProc()*nx; // Scale y grid with number of processors

  // Create funky stencil to make sure the matrix is non-symmetric (transpose non-trivial):

  // (i-1,j-1) (i-1,j  )
  // (i  ,j-1) (i  ,j  ) (i  ,j+1)
  // (i+1,j-1) (i+1,j  )

  int npoints = 7;

  int xoff[] = {-1,  0,  1, -1,  0,  1,  0};
  int yoff[] = {-1, -1, -1,  0,  0,  0,  1};

  Epetra_Map * map;
  Epetra_CrsMatrix * A;
  Epetra_Vector * x, * b, * xexact;
	
  Trilinos_Util_GenerateCrsProblem(nx, ny, npoints, xoff, yoff, comm, map, A, x, b, xexact);

  if (nx<8)
  {
    cout << *A << endl;
    cout << "X exact = " << endl << *xexact << endl;
    cout << "B       = " << endl << *b << endl;
  }

  // Construct transposer 
  Epetra_Time timer(comm);

  double start = timer.ElapsedTime();

  bool IgnoreNonLocalCols = false;
  bool MakeDataContiguous = true;
  EpetraExt::RowMatrix_Transpose transposer( MakeDataContiguous );

  if (verbose) cout << "\nTime to construct transposer  = " << timer.ElapsedTime() - start << endl;
  
  Epetra_CrsMatrix & transA = dynamic_cast<Epetra_CrsMatrix&>(transposer(*A));

  start = timer.ElapsedTime();
  if (verbose) cout << "\nTime to create transpose matrix  = " << timer.ElapsedTime() - start << endl;
 	
  // Now test output of transposer by performing matvecs
  int ierr = 0;
  ierr += checkResults(A, &transA, xexact, verbose);


  // Now change values in original matrix and test update facility of transposer
  // Add 2 to the diagonal of each row
  double Value = 2.0;
  for (i=0; i< A->NumMyRows(); i++)
  A->SumIntoMyValues(i, 1, &Value, &i);

  start = timer.ElapsedTime();
  transposer.fwd();

  if (verbose) cout << "\nTime to update transpose matrix  = " << timer.ElapsedTime() - start << endl;
 	
  ierr += checkResults(A, &transA, xexact, verbose);

  delete A;
  delete b;
  delete x;
  delete xexact;
  delete map;

  if (verbose) cout << endl << "Checking transposer for VbrMatrix objects" << endl<< endl;

  int nsizes = 4;
  int sizes[] = {4, 6, 5, 3};

  Epetra_VbrMatrix * Avbr;
  Epetra_BlockMap * bmap;

  Trilinos_Util_GenerateVbrProblem(nx, ny, npoints, xoff, yoff, nsizes, sizes,
                                   comm, bmap, Avbr, x, b, xexact);

  if (nx<8)
  {
    cout << *Avbr << endl;
    cout << "X exact = " << endl << *xexact << endl;
    cout << "B       = " << endl << *b << endl;
  }

  start = timer.ElapsedTime();
  EpetraExt::RowMatrix_Transpose transposer1( MakeDataContiguous );

  Epetra_CrsMatrix & transA1 = dynamic_cast<Epetra_CrsMatrix&>(transposer1(*Avbr));
  if (verbose) cout << "\nTime to create transpose matrix  = " << timer.ElapsedTime() - start << endl;
 	
  // Now test output of transposer by performing matvecs
;
  ierr += checkResults(Avbr, &transA1, xexact, verbose);

  // Now change values in original matrix and test update facility of transposer
  // Scale matrix on the left by rowsums

  Epetra_Vector invRowSums(Avbr->RowMap());

  Avbr->InvRowSums(invRowSums);
  Avbr->LeftScale(invRowSums);

  start = timer.ElapsedTime();
  transposer1.fwd();
  if (verbose) cout << "\nTime to update transpose matrix  = " << timer.ElapsedTime() - start << endl;
 	
  ierr += checkResults(Avbr, &transA1, xexact, verbose);

  delete Avbr;
  delete b;
  delete x;
  delete xexact;
  delete bmap;

#ifdef EPETRA_MPI
  MPI_Finalize();
#endif

  return ierr;
}

int checkResults(Epetra_RowMatrix * A, Epetra_CrsMatrix * transA, 
                 Epetra_Vector * xexact, bool verbose) {

  int n = A->NumGlobalRows();

  if (n<100) cout << "A transpose = " << endl << *transA << endl;

  Epetra_Vector x1(View,A->OperatorDomainMap(), &((*xexact)[0]));
  Epetra_Vector b1(A->OperatorRangeMap());

  A->SetUseTranspose(true);

  Epetra_Time timer(A->Comm());
  double start = timer.ElapsedTime();
  A->Apply(x1, b1);
  if (verbose) cout << "\nTime to compute b1: matvec with original matrix using transpose flag  = " << timer.ElapsedTime() - start << endl;

  if (n<100) cout << "b1 = " << endl << b1 << endl;
  Epetra_Vector x2(View,transA->OperatorRangeMap(), &((*xexact)[0]));
  Epetra_Vector b2(transA->OperatorDomainMap());
  start = timer.ElapsedTime();
  transA->Multiply(false, x2, b2);
  if (verbose) cout << "\nTime to compute b2: matvec with transpose matrix                      = " << timer.ElapsedTime() - start << endl;

  if (n<100) cout << "b1 = " << endl << b1 << endl;

  double residual;
  Epetra_Vector resid(A->OperatorDomainMap());

  resid.Update(1.0, b1, -1.0, b2, 0.0);
  assert(resid.Norm2(&residual)==0);
  if (verbose) cout << "Norm of b1 - b2 = " << residual << endl;

  int ierr = 0;

  if (residual > 1.0e-10) ierr++;

  if (ierr!=0 && verbose) cerr << "Status: Test failed" << endl;
  else if (verbose) cerr << "Status: Test passed" << endl;

  return(ierr);
}
