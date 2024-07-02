// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
//
//  This test instantiates the Anasazi classes using a complex scalar type
//  and checks functionality.
//

#include "AnasaziConfigDefs.hpp"
#include "AnasaziBasicSort.hpp"
#include "AnasaziMVOPTester.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "AnasaziBasicOutputManager.hpp"

#ifdef HAVE_MPI
#include <mpi.h>
#endif

// I/O for Harwell-Boeing files
#ifdef HAVE_ANASAZI_TRIUTILS
#include "Trilinos_Util_iohb.h"
#endif

#include "MyMultiVec.hpp"
#include "MyOperator.hpp"
#include "MyBetterOperator.hpp"

using namespace Teuchos;

int main(int argc, char *argv[])
{
  using std::cout;
  using std::endl;
  bool ierr, gerr;
  gerr = true;

#ifdef HAVE_MPI
  // Initialize MPI
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

  // Issue several useful typedefs;
  typedef std::complex<double> ST;
  typedef Anasazi::MultiVec<ST> MV;
  typedef Anasazi::Operator<ST> OP;
  typedef Anasazi::MultiVecTraits<ST,MV> MVT;
  //typedef Anasazi::OperatorTraits<ST,MV,OP> OPT; // unused

  // Create an output manager to handle the I/O from the solver
  RCP<Anasazi::OutputManager<ST> > MyOM
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
  std::vector<ST> cvals(nnz);
  for (int ii=0; ii<nnz; ii++) {
    cvals[ii] = ST(dvals[ii*2],dvals[ii*2+1]);
  }
  // Build the problem matrix
  RCP< const MyBetterOperator<ST> > A1
    = rcp( new MyBetterOperator<ST>(dim,colptr,nnz,rowind,&cvals[0]) );


  // Create a MyMultiVec for cloning
  std::vector<ScalarTraits<ST>::magnitudeType> v(blockSize);
  RCP< MyMultiVec<ST> > ivec = rcp( new MyMultiVec<ST>(dim,blockSize) );
  MVT::MvNorm(*ivec,v);

  // Create a MyOperator for testing against
  RCP<MyOperator<ST> > A2 = rcp( new MyOperator<ST>(dim) );

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

  // Clean up.
  std::free( dvals );
  std::free( colptr );
  std::free( rowind );

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
