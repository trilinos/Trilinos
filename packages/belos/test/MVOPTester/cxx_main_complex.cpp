//@HEADER
// ************************************************************************
//
//                 Belos: Block Linear Solvers Package
//                  Copyright 2004 Sandia Corporation
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
//
//  This test instantiates the Belos classes using a std::complex scalar type
//  and checks functionality.
//

#include "BelosConfigDefs.hpp"
#include "BelosMVOPTester.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "BelosOutputManager.hpp"

#ifdef HAVE_MPI
#include <mpi.h>
#endif

// I/O for Harwell-Boeing files
#ifdef HAVE_BELOS_TRIUTILS
#include "Trilinos_Util_iohb.h"
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

  bool success = false;
  bool verbose = false;
  try {
    int MyPID;
#ifdef HAVE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &MyPID);
#else
    MyPID = 0;
#endif
    (void) MyPID; // forestall "set but not used" warnings

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
    typedef std::complex<double> ST;
#else
    typedef double ST;
    // no std::complex. quit with failure.
    if (verbose && MyPID==0) {
      std::cout << "Not compiled with std::complex support." << std::endl;
      if (verbose && MyPID==0) {
        std::cout << "End Result: TEST FAILED" << std::endl;
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
    //typedef Belos::OperatorTraits<ST,MV,OP> OPT;

    // Create an output manager to handle the I/O from the solver
    RCP<Belos::OutputManager<ST> > MyOM
      = rcp( new Belos::OutputManager<ST>() );
    if (verbose) {
      MyOM->setVerbosity( Belos::Warnings );
    }

#ifndef HAVE_BELOS_TRIUTILS
    std::cout << "This test requires Triutils. Please configure with --enable-triutils." << std::endl;
#ifdef EPETRA_MPI
    MPI_Finalize() ;
#endif
    MyOM->print(Belos::Warnings,"End Result: TEST FAILED\n");
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
      MyOM->stream(Belos::Warnings)
        << "Warning reading '" << filename << "'" << std::endl
        << "End Result: TEST FAILED" << std::endl;
#ifdef HAVE_MPI
      MPI_Finalize();
#endif
      return -1;
    }
    // Convert interleaved doubles to std::complex values
    std::vector<ST> cvals(nnz);
    for (int ii=0; ii<nnz; ii++) {
      cvals[ii] = ST(dvals[ii*2],dvals[ii*2+1]);
    }
    // Build the problem matrix
    RCP< MyBetterOperator<ST> > A1
      = rcp( new MyBetterOperator<ST>(dim,colptr,nnz,rowind,&cvals[0]) );


    // Create a MyMultiVec for cloning
    std::vector<ScalarTraits<ST>::magnitudeType> v(blockSize);
    RCP< MyMultiVec<ST> > ivec = rcp( new MyMultiVec<ST>(dim,blockSize) );
    MVT::MvNorm(*ivec,v);

    // Create a MyOperator for testing against
    RCP<MyOperator<ST> > A2 = rcp( new MyOperator<ST>(dim) );

    // test the multivector and its adapter
    ierr = Belos::TestMultiVecTraits<ST,MV>(MyOM,ivec);
    gerr &= ierr;
    if (ierr) {
      MyOM->print(Belos::Warnings, "*** MyMultiVec<std::complex> PASSED TestMultiVecTraits()\n");
    }
    else {
      MyOM->print(Belos::Warnings, "*** MyMultiVec<std::complex> FAILED TestMultiVecTraits() ***\n\n");
    }

    // test the operator and its adapter
    ierr = Belos::TestOperatorTraits<ST,MV,OP>(MyOM,ivec,A2);
    gerr &= ierr;
    if (ierr) {
      MyOM->print(Belos::Warnings,"*** MyOperator<std::complex> PASSED TestOperatorTraits()\n");
    }
    else {
      MyOM->print(Belos::Warnings,"*** MyOperator<std::complex> FAILED TestOperatorTraits() ***\n\n");
    }

    // test the operator and its adapter
    ierr = Belos::TestOperatorTraits<ST,MV,OP>(MyOM,ivec,A1);
    gerr &= ierr;
    if (ierr) {
      MyOM->print(Belos::Warnings,"*** MyBetterOperator<std::complex> PASSED TestOperatorTraits()\n");
    }
    else {
      MyOM->print(Belos::Warnings,"*** MyBetterOperator<std::complex> FAILED TestOperatorTraits() ***\n\n");
    }

    // Clean up.
    free( dvals );
    free( colptr );
    free( rowind );

    success = gerr;
    if (success) {
      MyOM->print(Belos::Warnings,"End Result: TEST PASSED\n");
    } else {
      MyOM->print(Belos::Warnings,"End Result: TEST FAILED\n");
    }
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

#ifdef HAVE_MPI
    MPI_Finalize();
#endif

  return ( success ? EXIT_SUCCESS : EXIT_FAILURE );
}
