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

// Epetra includes
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
// TriUtils includes
#include "Trilinos_Util.h"
// Tpetra includes
#include "Tpetra_SerialPlatform.hpp"
#include "Tpetra_ElementSpace.hpp"
#include "Tpetra_VectorSpace.hpp"
#include "Tpetra_CisMatrix.hpp"

void test(Epetra_Comm& comm, Epetra_Map*& map, Epetra_CrsMatrix*& A, Epetra_Vector*& xexact,
					Epetra_Vector*& b, int dim, int nnz, bool verbose, bool smallProblem);

void tpetraAlone(Epetra_Comm& comm, Epetra_Map*& map, Epetra_CrsMatrix*& A, Epetra_Vector*& xexact,
								 Epetra_Vector*& b, int dim, int nnz, bool verbose, bool smallProblem);


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
  //bool verbose1 = true; // *WARNING-UNUSED*
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

  Epetra_Map* map;
  Epetra_CrsMatrix* A;
  Epetra_Vector* x; 
  Epetra_Vector* b;
  Epetra_Vector* xexact;

  Trilinos_Util_ReadHb2Epetra(argv[1], comm, map, A, x, b, xexact);

  bool smallProblem = false;
  int dim = map->NumGlobalElements();
  int nnz = A->NumGlobalNonzeros();
  if(dim < 100) 
		smallProblem = true;

  if(verbose && smallProblem)
    cout << "Original Matrix = " << endl << *A   << endl;

  x->PutScalar(0.0);

  //double norm1a = A->NormInf(); // *WARNING-UNUSED*
  if(verbose)
    cout << "Problem Dimension        = " << dim << endl
				 << "Number of matrix entries = " << nnz << endl;

	// ------------------------------------------------------------------
	// start of performance testing
	// ------------------------------------------------------------------

	test(comm, map, A, xexact, b, dim, nnz, verbose, smallProblem);
	//tpetraAlone(comm, map, A, xexact, b, dim, nnz, verbose, smallProblem);

	// ------------------------------------------------------------------
	// end of performance testing
	// ------------------------------------------------------------------
    
  // These objects were explicitly "new'ed" in ReadHb2Epetra
  delete A;
  delete x;
  delete b;
  delete xexact;
  delete map;

#ifdef EPETRA_MPI
  MPI_Finalize() ;
#endif

	return(0);
}

//=========================================================================================
// main testing function: does performance testing on both Epetra and Tpetra
//=========================================================================================
void test(Epetra_Comm& comm, Epetra_Map*& map, Epetra_CrsMatrix*& A, Epetra_Vector*& xexact, 
					Epetra_Vector*& b, int dim, int nnz, bool verbose, bool smallProblem) {
	// ------------------------------------------------------------------
	// create Tpetra versions of map, xexact, and b
	// ------------------------------------------------------------------

	// create Tpetra VectorSpace<int, double> , named vectorspace
	// should be compatible with map.
	if(!map->LinearMap())
		cerr << "*** Epetra_Map is not contiguous, can't create VectorSpace (yet). ***" << endl;
	Tpetra::SerialPlatform<int, double> platformV;
	Tpetra::SerialPlatform<int, int> platformE;
	Tpetra::ElementSpace<int> elementspace(map->NumGlobalElements(), map->NumMyElements(), map->IndexBase(), platformE);
	Tpetra::VectorSpace<int, double> vectorspace(elementspace, platformV);

	// create Tpetra Vector<int, double>, named xexact_t
	// should be identical to xexact
	Tpetra::Vector<int, double> xexact_t(xexact->Values(), xexact->GlobalLength(), vectorspace);

	// create Tpetra Vector<int, double>, named b_t
	// should be identical to b
	Tpetra::Vector<int, double> b_t(b->Values(), b->GlobalLength(), vectorspace);

	// ------------------------------------------------------------------
	// other initialization stuff
	// ------------------------------------------------------------------

  Epetra_Time timer(comm);
  comm.Barrier();

  int numEntries;
  double* values;
  int* indices;

	// ------------------------------------------------------------------
	// measure time to do creation and insertions
	// ------------------------------------------------------------------

  double tstart = timer.ElapsedTime();
  Epetra_CrsMatrix Ae(Copy, *map, 0);
  for (int i = 0; i < dim; i++) {
    A->ExtractMyRowView(i, numEntries, values, indices);
    Ae.InsertGlobalValues(i, numEntries, values, indices);
  }
  double epetraInsertTime = timer.ElapsedTime() - tstart;

	tstart = timer.ElapsedTime();
  Tpetra::CisMatrix<int, double> At(vectorspace);
  for (int i = 0; i < dim; i++) {
    A->ExtractMyRowView(i, numEntries, values, indices);
    At.submitEntries(Tpetra::Insert, i, numEntries, values, indices);
  }
  double tpetraInsertTime = timer.ElapsedTime() - tstart;

	// ------------------------------------------------------------------
	// measure time to do fillComplete
	// ------------------------------------------------------------------

  Ae.FillComplete();
  double epetraFillCompleteTime = timer.ElapsedTime() - tpetraInsertTime;

	At.fillComplete();
	double tpetraFillCompleteTime = timer.ElapsedTime() - epetraFillCompleteTime;

	// ------------------------------------------------------------------
	// measure time to do multiply/apply
	// ------------------------------------------------------------------

  // Next, compute how many times we should call the Multiply method, assuming a rate of 100 MFLOPS and a desired time of 1 second total.
  int niters = (int) (100000000.0/((double) 2*nnz));
  if (smallProblem) niters = 1;

  Epetra_Flops counter;
  Epetra_Vector bcomp_e(*map);
  Ae.SetFlopCounter(counter);
  tstart = timer.ElapsedTime();
  for(int i = 0; i < niters; i++) 
		Ae.Multiply(false, *xexact, bcomp_e);
  double epetraMatvecTime = timer.ElapsedTime() - tstart;
  double epetraNumFlops = Ae.Flops(); // Total number of Epetra FLOPS in Multiplies

  Teuchos::Flops flops;
  Tpetra::Vector<int, double> bcomp_t(vectorspace);
  At.setFlopCounter(flops);
  tstart = timer.ElapsedTime();
  for(int i = 0; i < niters; i++) 
		At.apply(xexact_t, bcomp_t); // At * xexact_t = bcomp_t
  double tpetraMatvecTime = timer.ElapsedTime() - tstart;
  double tpetraNumFlops = At.getFlops(); // Total number of Tpetra FLOPS in Multiplies

	/*if(smallProblem) {
		cout << "=============================================================" << endl; 
		cout << "bcomp_e:" << endl << bcomp_e << endl; 
		cout << "bcomp_t:" << endl << bcomp_t << endl; 
		cout << "xexact:" << endl << *xexact << endl; 
		cout << "xexact_t:" << endl << xexact_t << endl; 
		cout << "=============================================================" << endl; 
		}*/

	// ------------------------------------------------------------------
	// output results
	// ------------------------------------------------------------------

  if (verbose)
    cout << "\n\n****************************************************" << endl
				 << "    Epetra Insertion  time (sec)           = " << epetraInsertTime<< endl
				 << "    Epetra FillComplete time (sec)         = " << epetraFillCompleteTime << endl
				 << "    Epetra time for " << niters << " Matvecs (sec) = " << epetraMatvecTime << endl
				 << "    Epetra Total time (sec)           = " 
				 << epetraInsertTime+epetraFillCompleteTime+epetraMatvecTime  << endl
				 << "    Epetra Total flops = " << epetraNumFlops << endl
				 << "    Epetra MFLOPS = " << epetraNumFlops/epetraMatvecTime/1000000 << endl<<endl;

  if (verbose)
    cout << "\n\n****************************************************" << endl
				 << "    Tpetra Insertion  time (sec)           = " << tpetraInsertTime<< endl
				 << "    Tpetra FillComplete time (sec)         = " << tpetraFillCompleteTime << endl
				 << "    Tpetra time for " << niters << " Matvecs (sec) = " << tpetraMatvecTime << endl
				 << "    Tpetra Total time (sec)           = " 
				 << tpetraInsertTime+tpetraFillCompleteTime+tpetraMatvecTime  << endl
				 << "    Tpetra Total flops = " << tpetraNumFlops << endl
				 << "    Tpetra MFLOPS = " << tpetraNumFlops/tpetraMatvecTime/1000000 << endl<<endl;

  if (smallProblem)
		cout << " X          = " << endl << *xexact << endl
				 << " B expected = " << endl << *b << endl
				 << " B computed (Epetra) = " << endl << bcomp_e << endl
				 << " B computed (Tpetra) = " << endl << bcomp_t << endl;
	
	// ------------------------------------------------------------------
	// calculate & output residuals
	// ------------------------------------------------------------------

  Epetra_Vector resid_e(bcomp_e);
	//Tpetra::Vector<int, double> resid_t(bcomp_t);
	Tpetra::Vector<int, double> resid_t(bcomp_t.scalarPointer(), bcomp_t.getNumMyEntries(), bcomp_t.vectorSpace());

	/*cout << "=============================================================" << endl;
	cout << "bcomp_t addrs: " << &bcomp_t << " " << bcomp_t.scalarPointer() << endl; 
	cout << "resid_t addrs: " << &resid_t << " " << resid_t.scalarPointer() << endl; 
	cout << "bcomp_t contents: "; bcomp_t.printValues(cout); 
	cout << "resid_t contents: "; resid_t.printValues(cout); 
	cout << "=============================================================" << endl;*/

  resid_e.Update(1.0, *b, -1.0, bcomp_e, 0.0); // resid = xcomp - xexact
	resid_t.update(1.0, b_t, -1.0, bcomp_t, 0.0);
	/*if(smallProblem) {
		cout << "=============================================================" << endl; 
		cout << "resid_e:" << endl << resid_e << endl; 
		cout << "resid_t:" << endl << resid_t << endl; 
		cout << "bcomp_t:" << endl << bcomp_t << endl; 
		cout << "=============================================================" << endl; 
		}*/
  double residual_e, residual_t;
  resid_e.Norm2(&residual_e);   // residual_e = 2norm or resid_e
	residual_t = resid_t.norm2(); // residual_t = 2norm of resid_t
  double normb_e, normb_t, normb_exact;
  bcomp_e.Norm2(&normb_e);   // normb_e = 2norm of bcomp_e
	normb_t = bcomp_t.norm2(); // normb_t = 2norm of bcomp_t
  b->Norm2(&normb_exact);    // normb_exact = 2norm of b

  if (verbose) 
    cout << "2-norm of computed RHS (Epetra)                              = " << normb_e << endl
				 << "2-norm of computed RHS (Tpetra)                              = " << normb_t << endl
				 << "2-norm of exact RHS                                          = " << normb_exact << endl
				 << "2-norm of difference between computed and exact RHS (Epetra) = " << residual_e << endl
				 << "2-norm of difference between computed and exact RHS (Tpetra) = " << residual_t << endl;


}


//=========================================================================================
// Tpetra code from test, isolated and with timings removed
//=========================================================================================
void tpetraAlone(Epetra_Comm& comm, Epetra_Map*& map, Epetra_CrsMatrix*& A, Epetra_Vector*& xexact, 
								 Epetra_Vector*& b, int dim, int nnz, bool verbose, bool smallProblem) {
	// ------------------------------------------------------------------
	// create Tpetra versions of map, xexact, and b
	// ------------------------------------------------------------------

	// create Tpetra VectorSpace<int, double> , named vectorspace
	// should be compatible with map.
	if(!map->LinearMap())
		cerr << "*** Epetra_Map is not contiguous, can't create VectorSpace (yet). ***" << endl;
	Tpetra::SerialPlatform<int, double> platformV;
	Tpetra::SerialPlatform<int, int> platformE;
	Tpetra::ElementSpace<int> elementspace(map->NumGlobalElements(), map->NumMyElements(), map->IndexBase(), platformE);
	Tpetra::VectorSpace<int, double> vectorspace(elementspace, platformV);

	// create Tpetra Vector<int, double>, named xexact_t
	// should be identical to xexact
	Tpetra::Vector<int, double> xexact_t(xexact->Values(), xexact->GlobalLength(), vectorspace);

	// create Tpetra Vector<int, double>, named b_t
	// should be identical to b
	Tpetra::Vector<int, double> b_t(b->Values(), b->GlobalLength(), vectorspace);

	// ------------------------------------------------------------------
	// do creation, insertions, and fillComplete
	// ------------------------------------------------------------------
  int numEntries;
  double* values;
  int* indices;

  Tpetra::CisMatrix<int, double> At(vectorspace);
  for (int i = 0; i < dim; i++) {
    A->ExtractMyRowView(i, numEntries, values, indices);
    At.submitEntries(Tpetra::Insert, i, numEntries, values, indices);
  }

	At.fillComplete();

	cout << "===== [ETP] At, after fillComplete ================" << endl;
	cout << At;
	cout << "===================================================" << endl;

	// ------------------------------------------------------------------
	// do apply
	// ------------------------------------------------------------------

  int niters = (int) (100000000.0/((double) 2*nnz));
  if(smallProblem) niters = 1;
  Tpetra::Vector<int, double> bcomp_t(vectorspace);

	/*cout << "===== [ETP] before apply ==========================" << endl;
	cout << "xexact_t: \t" << &xexact_t << " \t"; xexact_t.printValues(cout);
	cout << "b_t: \t" << &b_t << " \t"; b_t.printValues(cout);
	cout << "bcomp_t: \t" << &bcomp_t << " \t"; bcomp_t.printValues(cout);
	cout << "===================================================" << endl;*/

  for(int i = 0; i < niters; i++) 
		At.apply(xexact_t, bcomp_t); // At * xexact_t = bcomp_t

	// ------------------------------------------------------------------
	// output results
	// ------------------------------------------------------------------

	/*cout << "===== [ETP] after apply ===========================" << endl;
	cout << "xexact_t: " << &xexact_t << " "; xexact_t.printValues(cout);
	cout << "b_t: " << &b_t << " "; b_t.printValues(cout);
	cout << "bcomp_t: " << &bcomp_t << " "; bcomp_t.printValues(cout);
	cout << "===================================================" << endl;*/
}
