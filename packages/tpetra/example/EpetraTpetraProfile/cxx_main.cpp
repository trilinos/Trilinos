// @HEADER
// ***********************************************************************
// 
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
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
// @HEADER

// Epetra includes
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_Time.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Export.h"
#include "Epetra_Version.h"
#ifdef EPETRA_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
// TriUtils includes
#include "Trilinos_Util.h"
// Tpetra includes
#include "Tpetra_ElementSpace.hpp"
#include "Tpetra_VectorSpace.hpp"
#include "Tpetra_CisMatrix.hpp"
#include "Tpetra_Version.hpp"
#ifdef TPETRA_MPI
#include "Tpetra_MpiPlatform.hpp"
#else
#include "Tpetra_SerialPlatform.hpp"
#endif
// Teuchos includes
#include "Teuchos_Array.hpp"

void test(Epetra_Comm& comm, Epetra_Map& map, Epetra_CrsMatrix& A, Epetra_Vector& xexact,
		  Epetra_Vector& b, int dim, bool verbose, bool smallProblem, std::string const name);

void outputResults(bool const verbose, int niters, std::string const name, int neq, int nnz,
				   double epetraInsertTime, double epetraFillCompleteTime, 
				   double epetraMatvecTime, double epetraNumFlops, 
				   double tpetraInsertTime, double tpetraFillCompleteTime,
				   double tpetraMatvecTime, double tpetraNumFlops);

int main(int argc, char *argv[]) {
	
#ifdef EPETRA_MPI
	MPI_Init(&argc,&argv);
	Epetra_MpiComm comm (MPI_COMM_WORLD);
#else
	Epetra_SerialComm comm;
#endif
    
	int myPID = comm.MyPID();
	
	bool verbose = false;
	if(myPID == 0) verbose = true; // verbose is true only on the root node
	
	if(verbose) {
		cout << "\n===========================================================================================\n";
		cout << "== Epetra Tpetra Profile" << endl;
		cout << "== " << Epetra_Version() << endl;
		cout << "== " << Tpetra::version() << endl;
		cout << "===========================================================================================\n";
	}
	
	if(argc < 2 && verbose) {
		cerr << "Usage: " << argv[0] << " HPC_filename " << endl 
		<< "where:" << endl 
		<< "HB_filename        - filename and path of a Harwell-Boeing data set" << endl
		<< "Example:" << endl 
		<< argv[0] << " mymatrix.hb" << endl << endl;
		return(1);
	}

#ifdef EPETRA_MPI
	char processor_name[MPI_MAX_PROCESSOR_NAME];
	int namelen;
	MPI_Get_processor_name(processor_name,&namelen);
	for(int i = 0; i < comm.NumProc(); i++) {
		if(i == myPID) {
			cout << "Image " << i << " is " << processor_name << endl;
		}
		comm.Barrier();
		comm.Barrier();
		comm.Barrier();
	}
	if(verbose) {
		cout << "===========================================================================================\n\n";
	}
#endif

	// convert filename to std::string for output
	std::string name(argv[1]); // construct string from null-terminated C-string
	std::string::size_type idx = name.rfind('/');
	if(idx != std::string::npos) { // we want just the filename, without the path to it
		name = name.substr(idx+1);
	}
	
	// ------------------------------------------------------------------
	// Use TrilUtil's ReadHb2Epetra to read in data file
	// ------------------------------------------------------------------
	
	Epetra_Map* readMap;
	Epetra_CrsMatrix* readA; 
	Epetra_Vector* readx; 
	Epetra_Vector* readb;
	Epetra_Vector* readxexact;
	
	Trilinos_Util_ReadHb2Epetra(argv[1], comm, readMap, readA, readx, readb, readxexact);
	
	// ------------------------------------------------------------------
	// Use an Epetra Import/Export to distribute the data globally
	// ------------------------------------------------------------------
	
	// Create uniform distributed map
	Epetra_Map map(readMap->NumGlobalElements(), 0, comm);
	
	// Create Exporter to distribute read-in matrix and vectors
	Epetra_Export exporter(*readMap, map);
	Epetra_CrsMatrix A(Copy, map, 0);
	Epetra_Vector b(map);
	Epetra_Vector xexact(map);
	
	// redistribute the vectors
	b.Export(*readb, exporter, Add);
	xexact.Export(*readxexact, exporter, Add);
	
	// redistribute the matrix
	A.Export(*readA, exporter, Add);
	
	// ------------------------------------------------------------------
	// Other initial stuff
	// ------------------------------------------------------------------
	
	bool smallProblem = false;
	int dim = map.NumGlobalElements();
	if(dim < 100) 
		smallProblem = true;
	
	/*if(smallProblem) {
		if(verbose)
			cout << "Original Matrix = " << endl;
		cout << A << endl;
	}
	
	if(verbose)
		cout << "Problem Dimension        = " << dim << endl
		     << "Number of matrix entries = " << nnz << endl;*/
	
	// ------------------------------------------------------------------
	// start of performance testing
	// ------------------------------------------------------------------
	
	// convert dim and nnz from global values to local values?
	dim = map.NumMyElements();
	
	test(comm, map, A, xexact, b, dim, verbose, smallProblem, name);
	
	// ------------------------------------------------------------------
	// end of performance testing
	// ------------------------------------------------------------------
    
	// These objects were explicitly "new'ed" in ReadHb2Epetra
	delete readMap;
	delete readA;
	delete readx;
	delete readb;
	delete readxexact;
	
	
#ifdef EPETRA_MPI
	MPI_Finalize() ;
#endif
	
	return(0);
}

//=========================================================================================
// main testing function: does performance testing on both Epetra and Tpetra
//=========================================================================================
void test(Epetra_Comm& comm, Epetra_Map& map, Epetra_CrsMatrix& A, Epetra_Vector& xexact, 
		  Epetra_Vector& b, int dim, bool verbose, bool smallProblem, std::string const name) {
	// ------------------------------------------------------------------
	// create Tpetra versions of map, xexact, and b
	// ------------------------------------------------------------------
	
	// create Tpetra VectorSpace<int, double> , named vectorspace
	// should be compatible with map.
	if(!map.LinearMap())
		cerr << "*** Epetra_Map is not contiguous, can't create VectorSpace (yet). ***" << endl;
#ifdef TPETRA_MPI
	Tpetra::MpiPlatform<int, double> platformV(MPI_COMM_WORLD);
	Tpetra::MpiPlatform<int, int> platformE(MPI_COMM_WORLD);
#else
	Tpetra::SerialPlatform<int, double> platformV;
	Tpetra::SerialPlatform<int, int> platformE;
#endif
	
	Tpetra::ElementSpace<int> elementspace(map.NumGlobalElements(), map.NumMyElements(), map.IndexBase(), platformE);
	Tpetra::VectorSpace<int, double> vectorspace(elementspace, platformV);

	// create Tpetra Vector<int, double>, named xexact_t
	// should be identical to xexact
	Tpetra::Vector<int, double> xexact_t(xexact.Values(), xexact.MyLength(), vectorspace);
	
	// create Tpetra Vector<int, double>, named b_t
	// should be identical to b
	Tpetra::Vector<int, double> b_t(b.Values(), b.MyLength(), vectorspace);

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
	
	double startTime = timer.ElapsedTime();
	Epetra_CrsMatrix Ae(Copy, map, 0);
	for(int i = 0; i < dim; i++) {
		int GIDi = A.GRID(i);
		A.ExtractGlobalRowView(GIDi, numEntries, values, indices);
		Ae.InsertGlobalValues(GIDi, numEntries, values, indices);
	}
	double epetraInsertTime = timer.ElapsedTime() - startTime;
	
	startTime = timer.ElapsedTime();
	Tpetra::CisMatrix<int, double> At(vectorspace);
	for(int i = 0; i < dim; i++) {
		int GIDi = A.GRID(i);
		A.ExtractGlobalRowView(GIDi, numEntries, values, indices);
		At.submitEntries(Tpetra::Insert, GIDi, numEntries, values, indices);
	}
	double tpetraInsertTime = timer.ElapsedTime() - startTime;

	// ------------------------------------------------------------------
	// measure time to do fillComplete
	// ------------------------------------------------------------------
	
	Ae.FillComplete();
	Ae.OptimizeStorage();
	double epetraFillCompleteTime = timer.ElapsedTime() - tpetraInsertTime;
	
	At.fillComplete();
	double tpetraFillCompleteTime = timer.ElapsedTime() - epetraFillCompleteTime;
	
	// ------------------------------------------------------------------
	// measure time to do multiply/apply
	// ------------------------------------------------------------------
	
	// First we call A.FillComplete() so that we can call A.NumGlobalNonzeros().
	// This will give us the nnz value we use. (We don't extract nnz from a call
	// to Ae.NumGlobalNonzeros() because we want to keep the states of Ae and At
	// as identical as possible.)
	A.FillComplete();
	int neq = A.NumGlobalRows();
	int nnz = A.NumGlobalNonzeros();
	
	// Next, compute how many times we should call the Multiply method, 
	// assuming a rate of 100 MFLOPS and a desired time of 1 second total.
	int niters = static_cast<int>(100000000.0 / static_cast<double>(2 * nnz));	
	if(smallProblem) // "small" problems are usually used for diagnostics, so only do it once
		niters = 1;
	
	Epetra_Vector bcomp_e(map);
	Epetra_Flops flops_e;
	Ae.SetFlopCounter(flops_e);
	startTime = timer.ElapsedTime();
	for(int i = 0; i < niters; i++) 
		Ae.Multiply(false, xexact, bcomp_e);
	double epetraMatvecTime = timer.ElapsedTime() - startTime;
	double epetraNumFlops = Ae.Flops(); // Total number of Epetra FLOPS in Multiplies
	
	Tpetra::Vector<int, double> bcomp_t(vectorspace);
	Teuchos::Flops flops_t;
	At.setFlopCounter(flops_t);
	startTime = timer.ElapsedTime();
	for(int i = 0; i < niters; i++) 
		At.apply(xexact_t, bcomp_t); // At * xexact_t = bcomp_t
	double tpetraMatvecTime = timer.ElapsedTime() - startTime;
	double tpetraNumFlops = At.getFlops(); // Total number of Tpetra FLOPS in Multiplies
	
	// ------------------------------------------------------------------
	// output results
	// ------------------------------------------------------------------
	
	outputResults(verbose, niters, name, neq, nnz,
				  epetraInsertTime, epetraFillCompleteTime, epetraMatvecTime, epetraNumFlops,
				  tpetraInsertTime, tpetraFillCompleteTime, tpetraMatvecTime, tpetraNumFlops);
	
	if(smallProblem) { // ** TODO ** This needs to be massaged for parallel output
		if(verbose) cout << "\n X          = " << endl;
		cout << xexact << endl;
		if(verbose) cout << " B expected = " << endl;
		cout << b << endl;
		if(verbose) cout << " B computed (Epetra) = " << endl;
		cout << bcomp_e << endl;
		if(verbose) cout << " B computed (Tpetra) = " << endl;
		cout << bcomp_t << endl;
	}
	
	// ------------------------------------------------------------------
	// calculate & output residuals
	// ------------------------------------------------------------------
	
	Epetra_Vector resid_e(bcomp_e);
	// make level 2 deep copy, Tpetra::Vector cpy ctr would only make level 1 deep copy
	Tpetra::Vector<int, double> resid_t(bcomp_t.scalarPointer(), bcomp_t.getNumMyEntries(), bcomp_t.vectorSpace());
	
	resid_e.Update(1.0, b, -1.0, bcomp_e, 0.0); // resid = b - bcomp
	resid_t.update(1.0, b_t, -1.0, bcomp_t, 0.0);
	double residual_e, residual_t;
	resid_e.Norm2(&residual_e);   // residual_e = 2norm or resid_e
	residual_t = resid_t.norm2(); // residual_t = 2norm of resid_t
	double normb_e, normb_t, normb_exact;
	bcomp_e.Norm2(&normb_e);   // normb_e = 2norm of bcomp_e
	normb_t = bcomp_t.norm2(); // normb_t = 2norm of bcomp_t
	b.Norm2(&normb_exact);    // normb_exact = 2norm of b
	
	if(verbose) // we only need to print this out once, because norms are a global op
		cout << "\n2-norm of computed RHS (Epetra)                              = " << normb_e << endl
			<< "2-norm of computed RHS (Tpetra)                              = " << normb_t << endl
			<< "2-norm of exact RHS                                          = " << normb_exact << endl
			<< "2-norm of difference between computed and exact RHS (Epetra) = " << residual_e << endl
			<< "2-norm of difference between computed and exact RHS (Tpetra) = " << residual_t << endl;
	
}

//=========================================================================================
// helper function to handle outputing the test results (but not the residuals)
//=========================================================================================
void outputResults(bool const verbose, int niters, std::string const name, int neq, int nnz,
				   double epetraInsertTime, double epetraFillCompleteTime, 
				   double epetraMatvecTime, double epetraNumFlops, 
				   double tpetraInsertTime, double tpetraFillCompleteTime,
				   double tpetraMatvecTime, double tpetraNumFlops) {
#ifdef TPETRA_MPI
	Tpetra::MpiComm<int, double> commV(MPI_COMM_WORLD);
	Tpetra::MpiComm<int, int> commE(MPI_COMM_WORLD);
#else
	Tpetra::SerialComm<int, double> commV;
	Tpetra::SerialComm<int, int> commE;
#endif
	int const numProcs = commE.getNumImages();
	
	// vectors to hold all values on the root node
	Teuchos::Array<int> niters_g(numProcs); // same for Epetra and Tpetra
	
	Teuchos::Array<double> epetraInsertTime_g(numProcs);
	Teuchos::Array<double> epetraFillCompleteTime_g(numProcs);
	Teuchos::Array<double> epetraMatvecTime_g(numProcs);
	Teuchos::Array<double> epetraNumFlops_g(numProcs);
	
	Teuchos::Array<double> tpetraInsertTime_g(numProcs);
	Teuchos::Array<double> tpetraFillCompleteTime_g(numProcs);
	Teuchos::Array<double> tpetraMatvecTime_g(numProcs);
	Teuchos::Array<double> tpetraNumFlops_g(numProcs);
	
	// do the gathers
	commE.gatherAll(&niters, &niters_g[0], 1);
	commV.gatherAll(&epetraInsertTime, &epetraInsertTime_g[0], 1);
	commV.gatherAll(&epetraFillCompleteTime, &epetraFillCompleteTime_g[0], 1);
	commV.gatherAll(&epetraMatvecTime, &epetraMatvecTime_g[0], 1);
	commV.gatherAll(&epetraNumFlops, &epetraNumFlops_g[0], 1);
	
	commV.gatherAll(&tpetraInsertTime, &tpetraInsertTime_g[0], 1);
	commV.gatherAll(&tpetraFillCompleteTime, &tpetraFillCompleteTime_g[0], 1);
	commV.gatherAll(&tpetraMatvecTime, &tpetraMatvecTime_g[0], 1);
	commV.gatherAll(&tpetraNumFlops, &tpetraNumFlops_g[0], 1);
	
	double const nnzPerRow(static_cast<double>(nnz) / static_cast<double>(neq));
	
	if(verbose) {
		cout << "\n*************************************************************************************************" << endl;
		cout << "Package name, Matrix Name, NEQ, NNZ, NNZ/Row, PID, Insert Time, FillComplete Time, Matvec Time, MFLOP Rate" << endl;
		cout << "*************************************************************************************************" << endl;
		for(int i = 0; i < numProcs; i++) {
			cout << "Epetra\t" << name << "\t" << neq << "\t" << nnz << "\t" << nnzPerRow << "\t" << i << "\t" 
				 << epetraInsertTime_g[i] << "\t" << epetraFillCompleteTime_g[i] << "\t"
				 << epetraMatvecTime_g[i] << "\t" << (epetraNumFlops_g[i] / epetraMatvecTime_g[i] / 1000000.0) 
				 << endl;
			cout << "Tpetra\t" << name << "\t" << neq << "\t" << nnz << "\t" << nnzPerRow << "\t" << i << "\t" 
				 << tpetraInsertTime_g[i] << "\t" << tpetraFillCompleteTime_g[i] << "\t"
				 << tpetraMatvecTime_g[i] << "\t" << (tpetraNumFlops_g[i] / tpetraMatvecTime_g[i] / 1000000.0) 
				 << endl;
		}
		cout << "\n# Matvecs (niters) = " << niters_g << endl;
		cout << "# Flops (Epetra) = " << epetraNumFlops_g << endl;
		cout << "# Flops (Tpetra) = " << tpetraNumFlops_g << endl;
	}
	
}
