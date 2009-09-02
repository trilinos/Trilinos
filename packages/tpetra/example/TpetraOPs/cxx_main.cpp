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
#ifdef HAVE_TPETRA_TBB
#include "Tpetra_TBB_TaskScheduler.hpp"
#include "Tpetra_TBB_Vec_Kernels.hpp"
#endif

#ifdef TPETRA_MPI

#ifdef HAVE_TPETRA_TBB
#include "Tpetra_TBB_MpiPlatform.hpp"
#endif

#include "Tpetra_MpiPlatform.hpp"

#else
#include "Tpetra_SerialPlatform.hpp"
#endif

// Teuchos includes
#include "Teuchos_Array.hpp"

void test(Epetra_Comm& comm, Epetra_Map& map, Epetra_CrsMatrix& A, Epetra_Vector& xexact,
      Epetra_Vector& b, int dim, bool verbose, const std::string& name);

void outputResults(bool verbose, const std::string& name, int neq, int nnz,
                   std::vector<double>& results);

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
    cout << "== Tpetra OPerations" << endl;
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

  int dim = map.NumGlobalElements();
  
  // ------------------------------------------------------------------
  // start of performance testing
  // ------------------------------------------------------------------
  
  // convert dim and nnz from global values to local values?
  dim = map.NumMyElements();
  
  test(comm, map, A, xexact, b, dim, verbose, name);
  
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

enum {
  INSERT_TIME = 0,
  FILLCOMPLETE_TIME = 1,
  NUM_MATVECS = 2,
  MATVEC_TIME = 3,
  MATVEC_MFLOPS = 4,
  NUM_AXPYS = 5,
  AXPY_TIME = 6,
  AXPY_MFLOPS = 7,
  NUM_DOTS = 8,
  DOT_TIME = 9,
  DOT_MFLOPS = 10,
  NUM_VALUES = 11
};

//=========================================================================================
// main testing function: does performance testing on Tpetra matrix/vector objects
//=========================================================================================
void test(Epetra_Comm& comm, Epetra_Map& map,
          Epetra_CrsMatrix& A, Epetra_Vector& xexact, 
          Epetra_Vector& b, int dim, bool verbose,
          const std::string& name) {
  // ------------------------------------------------------------------
  // create Tpetra versions of map, xexact, and b
  // ------------------------------------------------------------------
  
  // create Tpetra VectorSpace<int, double> , named vectorspace
  // should be compatible with map.
  if(!map.LinearMap())
    cerr << "*** Epetra_Map is not contiguous, can't create VectorSpace (yet). ***" << endl;
#ifdef TPETRA_MPI

  int num_threads = 2;

#ifdef HAVE_TPETRA_TBB
  Tpetra::TBB_MpiPlatform<int, double> platformV(MPI_COMM_WORLD, num_threads);
  Tpetra::TBB_MpiPlatform<int, int> platformE(MPI_COMM_WORLD, num_threads);
#else
  Tpetra::MpiPlatform<int, double> platformV(MPI_COMM_WORLD);
  Tpetra::MpiPlatform<int, int> platformE(MPI_COMM_WORLD);
#endif

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
  
  std::vector<double> results(NUM_VALUES, 0.0);

  // ------------------------------------------------------------------
  // measure time to do creation and insertions
  // ------------------------------------------------------------------
  
  double startTime = timer.ElapsedTime();
  
  Tpetra::CisMatrix<int, double> At(vectorspace);
  for(int i = 0; i < dim; i++) {
    int GIDi = A.GRID(i);
    A.ExtractGlobalRowView(GIDi, numEntries, values, indices);
    At.submitEntries(Tpetra::Insert, GIDi, numEntries, values, indices);
  }
  double tpetraInsertTime = timer.ElapsedTime() - startTime;
  results[INSERT_TIME] = tpetraInsertTime;

  // ------------------------------------------------------------------
  // measure time to do fillComplete
  // ------------------------------------------------------------------
  
  At.fillComplete();
 
  results[FILLCOMPLETE_TIME] = timer.ElapsedTime() - tpetraInsertTime;

  // ------------------------------------------------------------------
  // measure time to do multiply/apply
  // ------------------------------------------------------------------
  
  // First we call A.FillComplete() so that we can call A.NumGlobalNonzeros().
  // This will give us the nnz value we use.
  A.FillComplete();
  int neq = A.NumGlobalRows();
  int nnz = A.NumGlobalNonzeros();
 
  // Next, compute how many times we should call the Multiply method, 
  // assuming a rate of 500 MFLOPS and a desired time of 1 second total.
  int nmatvecs = static_cast<int>(500000000.0 / (200.0* static_cast<double>(nnz)));  
 
  results[NUM_MATVECS] = nmatvecs;

  Tpetra::Vector<int, double> bcomp_t(vectorspace);
  Teuchos::Flops flops_t;
  At.setFlopCounter(flops_t);

  startTime = timer.ElapsedTime();

  for(int i = 0; i < nmatvecs; i++) 
    At.apply(xexact_t, bcomp_t); // At * xexact_t = bcomp_t

  results[MATVEC_TIME] = timer.ElapsedTime() - startTime;
  results[MATVEC_MFLOPS] = At.getFlops()/1.e+6; // Total number of Tpetra FLOPS in Multiplies
  results[MATVEC_MFLOPS] /= results[MATVEC_TIME];
 
  // ------------------------------------------------------------------
  // calculate & output residuals
  // ------------------------------------------------------------------
  
  // make level 2 deep copy, Tpetra::Vector cpy ctr would only make level 1 deep copy
  Tpetra::Vector<int, double> resid_t(bcomp_t.scalarPointer(), bcomp_t.getNumMyEntries(), bcomp_t.vectorSpace());
  
  // resid = b - bcomp
  resid_t.update(1.0, b_t, -1.0, bcomp_t, 0.0);
  double residual_t;
  residual_t = resid_t.norm2(); // residual_t = 2norm of resid_t
  double normb_t, normb_exact;
  normb_t = bcomp_t.norm2(); // normb_t = 2norm of bcomp_t
  b.Norm2(&normb_exact);    // normb_exact = 2norm of b
  
  if(verbose) // we only need to print this out once, because norms are a global op
    cout << "2-norm of computed RHS (Tpetra)                              = " << normb_t << endl
      << "2-norm of exact RHS                                          = " << normb_exact << endl
      << "2-norm of difference between computed and exact RHS (Tpetra) = " << residual_t << endl;
  
  int naxpys = static_cast<int>(500000000.0 / (100.0*static_cast<double>(neq)));

  results[NUM_AXPYS] = naxpys;

  startTime = timer.ElapsedTime();

  for(int i=0; i<naxpys; ++i) {
    resid_t.update(1.0, b_t, 0.5);
  }

  results[AXPY_TIME] = timer.ElapsedTime() - startTime;

  results[AXPY_MFLOPS] = 2.e-6*neq;
  results[AXPY_MFLOPS] *= naxpys;
  results[AXPY_MFLOPS] /= results[AXPY_TIME];

  int ndots = static_cast<int>(500000000.0 / (100.0*static_cast<double>(neq)));

  results[NUM_DOTS] = ndots;

  startTime = timer.ElapsedTime();

  double dot = 0.0;
  for(int i=0; i<ndots; ++i) {
    dot = resid_t.dotProduct(b_t);
  }

  results[DOT_TIME] = timer.ElapsedTime() - startTime;

  results[DOT_MFLOPS] = 2.e-6*neq;
  results[DOT_MFLOPS] *= ndots;
  results[DOT_MFLOPS] /= results[DOT_TIME];

  // ------------------------------------------------------------------
  // output results
  // ------------------------------------------------------------------
  
  outputResults(verbose, name, neq, nnz, results);
}

//=========================================================================================
// helper function to handle outputing the test results (but not the residuals)
//=========================================================================================
void outputResults(bool verbose, const std::string& name, int neq, int nnz,
                   std::vector<double>& results) {
#ifdef TPETRA_MPI
  Tpetra::MpiComm<int, double> commV(MPI_COMM_WORLD);
#else
  Tpetra::SerialComm<int, double> commV;
#endif
 
  // vectors to hold min/max values on the root node
  std::vector<double> results_min(NUM_VALUES);
  std::vector<double> results_max(NUM_VALUES);
 
  // do the gathers
 
  commV.minAll(&results[0], &results_min[0], results.size());
  commV.maxAll(&results[0], &results_max[0], results.size());
 
  unsigned name_len = name.size();
  if (name_len < 14) name_len = 14;

  if(verbose) {
    cout << "\n**************************************************************************************" << endl;
    cout.width(name_len);
    cout << "Matrix Name";
    cout << "   NEQ     NNZ" <<endl;
    cout << "****************************************************************************************" << endl;
    cout.width(name_len);
    cout << name;
    cout.width(6); cout << neq;
    cout.width(8); cout << nnz;
    cout << "\n" << endl;
    cout << "****************************************************************************************" << endl;
    cout << "Min. Insert Time  Max. Insert Time   Min. FillComplete  Max. FillComplete" << endl;
    cout << "****************************************************************************************" << endl;
    cout.width(14); cout << results_min[INSERT_TIME];
    cout.width(14); cout << results_max[INSERT_TIME];
    cout.width(20); cout << results_min[FILLCOMPLETE_TIME];
    cout.width(14); cout << results_max[FILLCOMPLETE_TIME];
    cout << "\n" << endl;
    cout << "\n**************************************************************************************" << endl;
    cout << "  Min. Matvec Time  Max. Matvec Time    Min. MFLOPs    Max. MFLOPs" << endl;
    cout << "\n**************************************************************************************" << endl;
    cout.width(14); cout << results_min[MATVEC_TIME];
    cout.width(14); cout << results_max[MATVEC_TIME];
    cout.width(20); cout << results_min[MATVEC_MFLOPS];
    cout.width(14); cout << results_max[MATVEC_MFLOPS];
    cout << "\n" << endl;
    cout << "\n**************************************************************************************" << endl;
    cout << "  Min. AXPY Time  Max. AXPY Time    Min. MFLOPs    Max. MFLOPs" << endl;
    cout << "\n**************************************************************************************" << endl;
    cout.width(14); cout << results_min[AXPY_TIME];
    cout.width(14); cout << results_max[AXPY_TIME];
    cout.width(20); cout << results_min[AXPY_MFLOPS];
    cout.width(14); cout << results_max[AXPY_MFLOPS];
    cout << "\n" << endl;
    cout << "\n**************************************************************************************" << endl;
    cout << "  Min. DOT Time  Max. DOT Time    Min. MFLOPs    Max. MFLOPs" << endl;
    cout << "\n**************************************************************************************" << endl;
    cout.width(14); cout << results_min[DOT_TIME];
    cout.width(14); cout << results_max[DOT_TIME];
    cout.width(20); cout << results_min[DOT_MFLOPS];
    cout.width(14); cout << results_max[DOT_MFLOPS];
    cout << "\n" << endl;
    cout << "# Matvecs = " << static_cast<int>(results_min[NUM_MATVECS]) << endl;
    cout << "# AXPYs = " << static_cast<int>(results_min[NUM_AXPYS]) << endl;
    cout << "# DOTs = " << static_cast<int>(results_min[NUM_DOTS]) << endl;
  }
}

