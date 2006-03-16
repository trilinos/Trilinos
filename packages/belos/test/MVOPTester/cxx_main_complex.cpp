//@HEADER
// ************************************************************************
// 
//
//                 Belos: Block Linear Solvers Package
//                 Copyright (2004) Sandia Corporation
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
//
//  This test instantiates the Belos classes using a complex scalar type
//  and checks functionality.
//

#include "BelosConfigDefs.hpp"
#include "BelosMVOPTester.hpp"
#include "Teuchos_CommandLineProcessor.hpp"

#ifdef HAVE_MPI
#include <mpi.h>
#endif

// I/O for Harwell-Boeing files
#ifdef HAVE_BELOS_TRIUTILS
#include "iohb.h"
#endif

#include "MyMultiVec.hpp"
#include "MyOperator.hpp"
#include "MyBetterOperator.hpp"

using namespace Teuchos;

int main(int argc, char *argv[])
{
  int ierr, gerr;
  gerr = 0;

#ifdef HAVE_MPI
  // Initialize MPI and setup an Epetra communicator
  MPI_Init(&argc,&argv);
#endif

  // PID info
  int MyPID = 0;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &MyPID);
#endif
  bool verbose = false;
  std::string filename("mhd1280b.cua");

  // number of global elements
  int blockSize = 5;

  CommandLineProcessor cmdp(false,true);
  cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
  cmdp.setOption("filename",&filename,"Filename for Harwell-Boeing test matrix.");
  if (cmdp.parse(argc,argv) != CommandLineProcessor::PARSE_SUCCESSFUL) {
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return -1;
  }

#ifdef HAVE_COMPLEX
  typedef std::complex<double> ST;
#elif HAVE_COMPLEX_H
  typedef ::complex<double> ST;
#else
  typedef double ST;
  // no complex. quit with failure.
  if (verbose && MyPID==0) {
    cout << "Not compiled with complex support." << endl;
    if (verbose && MyPID==0) {
      cout << "End Result: TEST FAILED" << endl;
    }
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return -1;
  }
#endif

  // Issue several useful typedefs;
  typedef Belos::MultiVec<ST> MV;
  typedef Belos::Operator<ST> OP;
  typedef Belos::MultiVecTraits<ST,MV> MVT;
  typedef Belos::OperatorTraits<ST,MV,OP> OPT;

  // Create an output manager to handle the I/O from the solver
  RefCountPtr<Belos::OutputManager<ST> > MyOM 
    = rcp( new Belos::OutputManager<ST>( MyPID ) );
  if (verbose) {
  }


#ifndef HAVE_BELOS_TRIUTILS
  cout << "This test requires Triutils. Please configure with --enable-triutils." << endl;
#ifdef HAVE_MPI
  MPI_Finalize() ;
#endif
  if (verbose && MyPID==0) {
    cout << "End Result: TEST FAILED" << endl;	
  }
  return -1;
#endif

  // Get the data from the HB file
  int info;
  int dim,dim2,nnz;
  double *dvals;
  int *colptr,*rowind;
  ST *cvals;
  nnz = -1;
  info = readHB_newmat_double(filename.c_str(),&dim,&dim2,&nnz,&colptr,&rowind,&dvals);
  if (info == 0 || nnz < 0) {
    if (MyOM->isVerbosityAndPrint( Belos::Errors )) {
      cout << "Error reading '" << filename << "'" << endl;
      cout << "End Result: TEST FAILED" << endl;
    }
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return -1;
  }
  // Convert interleaved doubles to complex values
  cvals = new ST[nnz];
  for (int ii=0; ii<nnz; ii++) {
    cvals[ii] = ST(dvals[ii*2],dvals[ii*2+1]);
  }

  // Build the problem matrix
  RefCountPtr< MyBetterOperator<ST> > A1
    = rcp( new MyBetterOperator<ST>(dim,colptr,nnz,rowind,cvals) );

  // Create a MyMultiVec for cloning
  std::vector<ScalarTraits<ST>::magnitudeType> v(blockSize);
  RefCountPtr< MyMultiVec<ST> > ivec = rcp( new MyMultiVec<ST>(dim,blockSize) );
  MVT::MvNorm(*ivec,&v);

  // Create a MyOperator for testing against
  RefCountPtr<MyOperator<ST> > A2 = rcp( new MyOperator<ST>(dim) );

  // test the multivector and its adapter
  ierr = Belos::TestMultiVecTraits<ST,MV>(MyOM,ivec);
  gerr |= ierr;
  switch (ierr) {
  case Belos::Ok:
    if ( verbose && MyPID==0 ) {
      cout << "*** MyMultiVec<complex> PASSED TestMultiVecTraits()" << endl;
    }
    break;
  case Belos::Error:
    if ( verbose && MyPID==0 ) {
      cout << "*** MyMultiVec<complex> FAILED TestMultiVecTraits() ***" 
           << endl << endl;
    }
    break;
  }

  // test the operator and its adapter
  ierr = Belos::TestOperatorTraits<ST,MV,OP>(MyOM,ivec,A2);
  gerr |= ierr;
  switch (ierr) {
  case Belos::Ok:
    if ( verbose && MyPID==0 ) {
      cout << "*** MyOperator<complex> PASSED TestOperatorTraits()" << endl;
    }
    break;
  case Belos::Error:
    if ( verbose && MyPID==0 ) {
      cout << "*** MyOperator<complex> FAILED TestOperatorTraits() ***" 
           << endl << endl;
    }
    break;
  }

  // test the operator and its adapter
  ierr = Belos::TestOperatorTraits<ST,MV,OP>(MyOM,ivec,A1);
  gerr |= ierr;
  switch (ierr) {
  case Belos::Ok:
    if ( verbose && MyPID==0 ) {
      cout << "*** MyBetterOperator<complex> PASSED TestOperatorTraits()" << endl;
    }
    break;
  case Belos::Error:
    if ( verbose && MyPID==0 ) {
      cout << "*** MyBetterOperator<complex> FAILED TestOperatorTraits() ***" 
           << endl << endl;
    }
    break;
  }

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  if (gerr) {
    if (verbose && MyPID==0)
      cout << "End Result: TEST FAILED" << endl;
    return -1;
  }
  //
  // Default return value
  //
  if (verbose && MyPID==0)
    cout << "End Result: TEST PASSED" << endl;
  return 0;

}
