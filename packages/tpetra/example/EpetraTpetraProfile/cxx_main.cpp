//@HEADER
// ************************************************************************
// 
//               Tpetra: Templated Linear Algebra Services Package 
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

#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_Time.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#ifdef EPETRA_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Trilinos_Util.h"

int main(int argc, char *argv[]) {

#ifdef EPETRA_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm comm (MPI_COMM_WORLD);
#else
  Epetra_SerialComm comm;
#endif

  cout << comm << endl;

  int MyPID = comm.MyPID();

  bool verbose = false;
  bool verbose1 = true;
  if (MyPID==0) verbose = true;

  if(argc < 2 && verbose) {
    cerr << "Usage: " << argv[0] << " HPC_filename " << endl << "where:" << endl
	 << "HB_filename        - filename and path of a Harwell-Boeing data set" << endl
	 << "Example:" << endl << argv[0] << " mymatrix.hb" << endl << endl;
    return(1);
  }
  if (comm.NumProc()!=1) {
    cerr << "This driver does not support more than one processor at this time."<< endl
	 << "The insertion of indices must be preceded by conversion to global index space first, if using more than one processor."
	 << endl;
    return(1);
  }

  // Uncomment the next three lines to debug in mpi mode
  //int tmp;
  //if (MyPID==0) cin >> tmp;
  //comm.Barrier();

  Epetra_Map * map;
  Epetra_CrsMatrix * A;
  Epetra_Vector * x; 
  Epetra_Vector * b;
  Epetra_Vector * xexact;

  Trilinos_Util_ReadHb2Epetra(argv[1], comm, map, A, x, b, xexact);

  bool smallProblem = false;
  int dim = map->NumGlobalElements();
  int nnz = A->NumGlobalNonzeros();
  if (dim<100) smallProblem = true;

  if (smallProblem)
    cout << "Original Matrix = " << endl << *A   << endl;

  x->PutScalar(0.0);

  double norm1a = A->NormInf();
  if (verbose)
    cout << "Problem Dimension        = " << dim << endl
	 << "Number of matrix entries = " << nnz << endl;
  
  Epetra_Time timer(comm);
  comm.Barrier();

  int numEntries;
  double * values;
  int * indices;

  double tstart = timer.ElapsedTime();
  Epetra_CrsMatrix Ae(Copy, *map, 0);
  for (int i=0; i< dim; i++) {
    A->ExtractMyRowView(i, numEntries, values, indices);
    Ae.InsertGlobalValues(i, numEntries, values, indices);
  }
  double epetraInsertTime = timer.ElapsedTime() - tstart;
  Ae.FillComplete();
  double epetraFillCompleteTime = timer.ElapsedTime() - epetraInsertTime;

  // Next, compute how many times we should call the Multiply method, assuming a rate of 100 MFLOPS and a desired time of 1 second total.
  int niters = (int) (100000000.0/((double) 2*nnz));

  Epetra_Flops counter;
  Epetra_Vector bcomp(*map);
  Ae.SetFlopCounter(counter);
  tstart = timer.ElapsedTime();
  for (int i=0; i<niters; i++) Ae.Multiply(false, *xexact, bcomp);
  double epetraMatvecTime = timer.ElapsedTime() - tstart;
  double numFlops = Ae.Flops(); // Total number of FLOPS in Multiplies
  if (verbose)
    cout << "\n\n****************************************************" << endl
	 << "    Epetra Insertion  time (sec)           = " << epetraInsertTime<< endl
	 << "    Epetra FillComplete time (sec)         = " << epetraFillCompleteTime << endl
	 << "    Epetra time for " << niters << " Matvecs (sec) = " << epetraMatvecTime << endl
	 << "    Epetra Total time (sec)           = " 
	 << epetraInsertTime+epetraFillCompleteTime+epetraMatvecTime  << endl
	 << "    Epetra MFLOPS = " << numFlops/epetraMatvecTime/1000000 << endl<<endl;

  
  if (smallProblem)
  cout << " X          = " << endl << *xexact << endl
       << " B expected = " << endl << *b << endl
       << " B computed = " << endl << bcomp << endl;

  Epetra_Vector resid(bcomp);

  resid.Update(1.0, *b, -1.0, bcomp, 0.0); // resid = xcomp - xexact
  double residual;
  resid.Norm2(&residual);
  double normb, normbexact;
  bcomp.Norm2(&normb);
  b->Norm2(&normbexact);

  if (verbose) 
    cout << "2-norm of computed RHS                               = " << normb << endl
	 << "2-norm of exact RHS                                  = " << normbexact << endl
	 << "2-norm of difference between computed and exact RHS  = " << residual << endl;
    
  // These objects were explicitly "new'ed" in ReadHb2Epetra
  delete A;
  delete x;
  delete b;
  delete xexact;
  delete map;

#ifdef EPETRA_MPI
  MPI_Finalize() ;
#endif

return 0 ;
}
