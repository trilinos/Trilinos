//@HEADER
// ************************************************************************
// 
//
//                 Anasazi: Block Eigensolvers Package
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
//  This test instantiates the Anasazi classes using a complex scalar type
//  and checks functionality.
//

#include "AnasaziConfigDefs.hpp"
#include "AnasaziBasicSort.hpp"
#include "AnasaziEpetraAdapter.hpp"
#include "AnasaziMVOPTester.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "AnasaziBasicOutputManager.hpp"

#ifdef HAVE_MPI
#include <mpi.h>
#endif

// I/O for Harwell-Boeing files
#ifdef HAVE_ANASAZI_TRIUTILS
#include "iohb.h"
#endif

#include "MyMultiVec.hpp"
#include "MyOperator.hpp"
#include "MyBetterOperator.hpp"

using namespace Teuchos;

int main(int argc, char *argv[])
{
  bool ierr, gerr;
  gerr = true;

#ifdef HAVE_MPI
  // Initialize MPI and setup an Epetra communicator
  MPI_Init(&argc,&argv);
#endif

  int MyPID;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &MyPID);
#else 
  MyPID = 0;
#endif
  bool verbose = false;
  std::string filename("mhd1280b.cua");

  // number of global elements
  int blockSize = 5;

  CommandLineProcessor cmdp(false,true);
  cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
  cmdp.setOption("debug","quiet",&verbose,"Print messages and results.");
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
  //typedef Anasazi::MultiVec<ST> MV;
  typedef Anasazi::MultiVec<ST> MV;
  typedef Anasazi::Operator<ST> OP;
  typedef Anasazi::MultiVecTraits<ST,MV> MVT;
  typedef Anasazi::OperatorTraits<ST,MV,OP> OPT;

  // Create an output manager to handle the I/O from the solver
  RefCountPtr<Anasazi::OutputManager<ST> > MyOM 
    = rcp( new Anasazi::BasicOutputManager<ST>() );
  if (verbose) {
    MyOM->setVerbosity( Anasazi::Warnings );
  }


#ifndef HAVE_ANASAZI_TRIUTILS
  cout << "This test requires Triutils. Please configure with --enable-triutils." << endl;
#ifdef EPETRA_MPI
  MPI_Finalize() ;
#endif
  MyOM->print(Anasazi::Warnings,"End Result: TEST FAILED\n");
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
    MyOM->stream(Anasazi::Warnings) 
      << "Warning reading '" << filename << "'" << endl
      << "End Result: TEST FAILED" << endl;
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
  ierr = Anasazi::TestMultiVecTraits<ST,MV>(MyOM,ivec);
  gerr &= ierr;
  if (ierr) {
    MyOM->print(Anasazi::Warnings, "*** MyMultiVec<complex> PASSED TestMultiVecTraits()\n");
  }
  else {
    MyOM->print(Anasazi::Warnings, "*** MyMultiVec<complex> FAILED TestMultiVecTraits() ***\n\n");
  }

  // test the operator and its adapter
  ierr = Anasazi::TestOperatorTraits<ST,MV,OP>(MyOM,ivec,A2);
  gerr &= ierr;
  if (ierr) {
    MyOM->print(Anasazi::Warnings,"*** MyOperator<complex> PASSED TestOperatorTraits()\n");
  }
  else {
    MyOM->print(Anasazi::Warnings,"*** MyOperator<complex> FAILED TestOperatorTraits() ***\n\n");
  }

  // test the operator and its adapter
  ierr = Anasazi::TestOperatorTraits<ST,MV,OP>(MyOM,ivec,A1);
  gerr &= ierr;
  if (ierr) {
    MyOM->print(Anasazi::Warnings,"*** MyBetterOperator<complex> PASSED TestOperatorTraits()\n");
  }
  else {
    MyOM->print(Anasazi::Warnings,"*** MyBetterOperator<complex> FAILED TestOperatorTraits() ***\n\n");
  }

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  if (gerr == false) {
    MyOM->print(Anasazi::Warnings,"End Result: TEST FAILED\n");
    return -1;
  }
  //
  // Default return value
  //
  MyOM->print(Anasazi::Warnings,"End Result: TEST PASSED\n");
  return 0;

}
