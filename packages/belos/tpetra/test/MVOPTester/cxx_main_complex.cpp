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
//  This test instantiates the Belos classes using a std::complex scalar type
//  and checks functionality.
//

#include <Teuchos_UnitTestHarness.hpp>
#include <Tpetra_DefaultPlatform.hpp>
#include <Tpetra_CrsMatrix.hpp>

#include "BelosConfigDefs.hpp"
#include "BelosMVOPTester.hpp"
#include "BelosTpetraAdapter.hpp"
#include "BelosOutputManager.hpp"

// I/O for Harwell-Boeing files
#ifdef HAVE_BELOS_TRIUTILS
#include "iohb.h"
#endif

namespace {

  using Teuchos::as;
  using Teuchos::RCP;
  using Teuchos::ArrayRCP;
  using Teuchos::rcp;
  using Tpetra::Map;
  using Tpetra::DefaultPlatform;
  using Tpetra::Platform;
  using std::vector;
  using std::sort;
  using Teuchos::arrayViewFromVector;
  using Teuchos::broadcast;
  using Teuchos::OrdinalTraits;
  using Teuchos::ScalarTraits;
  using Teuchos::Comm;
  using Tpetra::MultiVector;
  using Tpetra::CrsMatrix;
  using std::endl;
  using Teuchos::Array;
  using Teuchos::ArrayView;
  using Teuchos::NO_TRANS;
  using Teuchos::TRANS;
  using Teuchos::CONJ_TRANS;
  using Belos::OutputManager;
  using Belos::Warnings;

  bool testMpi = true;
  double errorTolSlack = 1e+1;

#define PRINT_VECTOR(v) \
   { \
     out << #v << ": "; \
     copy(v.begin(), v.end(), ostream_iterator<Ordinal>(out," ")); \
     out << endl; \
   }

  TEUCHOS_STATIC_SETUP()
  {
    Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
    clp.addOutputSetupOptions(true);
    clp.setOption(
        "test-mpi", "test-serial", &testMpi,
        "Test MPI (if available) or force test of serial.  In a serial build,"
        " this option is ignord and a serial comm is always used." );
    clp.setOption(
        "error-tol-slack", &errorTolSlack,
        "Slack off of machine epsilon used to check test results" );
  }

  template<class Ordinal>
  RCP<const Platform<Ordinal> > getDefaultPlatform()
  {
    if (testMpi) {
      return DefaultPlatform<Ordinal>::getPlatform();
    }
    return rcp(new Tpetra::SerialPlatform<Ordinal>());
  }

  template<class Ordinal, class Scalar> 
  RCP<CrsMatrix<Ordinal,Scalar> > constructTriDiagMatrix(const Map<Ordinal> &map) 
  {
    RCP<CrsMatrix<Ordinal,Scalar> > op = rcp( new CrsMatrix<Ordinal,Scalar>(map) );
    op->fillComplete();
    return op;
  }

  //
  // UNIT TESTS
  // 

  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MultiVector, MVTestDist, Ordinal, Scalar )
  {
    typedef Tpetra::MultiVector<Ordinal,Scalar> MV;
    const Ordinal dim = 500;
    const Ordinal numVecs = 5;
    const Ordinal ZERO = OrdinalTraits<Ordinal>::zero();
    // Create an output manager to handle the I/O from the solver
    RCP<OutputManager<Scalar> > MyOM = rcp( new OutputManager<Scalar>(Warnings,rcp(&out,false)) );
    // create a platform  
    const Platform<Ordinal> & platform = *(getDefaultPlatform<Ordinal>());
    // create a comm  
    RCP<Comm<Ordinal> > comm = platform.createComm();
    // create a uniform contiguous map
    Map<Ordinal> map(dim,ZERO,platform);
    RCP<MV> mvec = rcp( new MV(map,numVecs,true) );
    bool res = Belos::TestMultiVecTraits<Scalar,MV>(MyOM,mvec);
    TEST_EQUALITY_CONST(res,true);
    // All procs fail if any proc fails
    int globalSuccess_int = -1;
    reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, &globalSuccess_int );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }

  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MultiVector, MVTestLocal, Ordinal, Scalar )
  {
    typedef Tpetra::MultiVector<Ordinal,Scalar> MV;
    const Ordinal dim = 500;
    const Ordinal numVecs = 5;
    const Ordinal ZERO = OrdinalTraits<Ordinal>::zero();
    // Create an output manager to handle the I/O from the solver
    RCP<OutputManager<Scalar> > MyOM = rcp( new OutputManager<Scalar>(Warnings,rcp(&out,false)) );
    // create a platform  
    const Platform<Ordinal> & platform = *(getDefaultPlatform<Ordinal>());
    // create a comm  
    RCP<Comm<Ordinal> > comm = platform.createComm();
    // create a uniform contiguous map
    Map<Ordinal> map(dim,ZERO,platform,true);
    RCP<MV> mvec = rcp( new MV(map,numVecs,true) );
    bool res = Belos::TestMultiVecTraits<Scalar,MV>(MyOM,mvec);
    TEST_EQUALITY_CONST(res,true);
    // All procs fail if any proc fails
    int globalSuccess_int = -1;
    reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, &globalSuccess_int );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }

  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MultiVector, OPTestLocal, Ordinal, Scalar )
  {
    typedef Tpetra::MultiVector<Ordinal,Scalar> MV;
    typedef Tpetra::Operator<Ordinal,Scalar>    OP;
    // const Ordinal dim = 500;
    const Ordinal dim = 10;
    const Ordinal numVecs = 5;
    const Ordinal ZERO = OrdinalTraits<Ordinal>::zero();
    // Create an output manager to handle the I/O from the solver
    RCP<OutputManager<Scalar> > MyOM = rcp( new OutputManager<Scalar>(Warnings,rcp(&out,false)) );
    // create a platform  
    const Platform<Ordinal> & platform = *(getDefaultPlatform<Ordinal>());
    // create a comm  
    RCP<Comm<Ordinal> > comm = platform.createComm();
    // create a uniform contiguous map (local)
    Map<Ordinal> map(dim,ZERO,platform,true);
    // create a CrsMatrix
    RCP<OP> op = constructTriDiagMatrix<Ordinal,Scalar>(map);
    // create a multivector
    RCP<MV> mvec = rcp( new MV(map,numVecs,true) );
    bool res = Belos::TestOperatorTraits<Scalar,MV,OP>(MyOM,mvec,op);
    TEST_EQUALITY_CONST(res,true);
    // All procs fail if any proc fails
    int globalSuccess_int = -1;
    reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, &globalSuccess_int );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }

  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MultiVector, OPTestDist, Ordinal, Scalar )
  {
    typedef Tpetra::MultiVector<Ordinal,Scalar> MV;
    typedef Tpetra::Operator<Ordinal,Scalar>    OP;
    // const Ordinal dim = 500;
    const Ordinal dim = 10;
    const Ordinal numVecs = 5;
    const Ordinal ZERO = OrdinalTraits<Ordinal>::zero();
    // Create an output manager to handle the I/O from the solver
    RCP<OutputManager<Scalar> > MyOM = rcp( new OutputManager<Scalar>(Warnings,rcp(&out,false)) );
    // create a platform  
    const Platform<Ordinal> & platform = *(getDefaultPlatform<Ordinal>());
    // create a comm  
    RCP<Comm<Ordinal> > comm = platform.createComm();
    // create a uniform contiguous map
    Map<Ordinal> map(dim,ZERO,platform);
    // create a CrsMatrix
    RCP<OP> op = constructTriDiagMatrix<Ordinal,Scalar>(map);
    // create a multivector
    RCP<MV> mvec = rcp( new MV(map,numVecs,true) );
    bool res = Belos::TestOperatorTraits<Scalar,MV,OP>(MyOM,mvec,op);
    TEST_EQUALITY_CONST(res,true);
    // All procs fail if any proc fails
    int globalSuccess_int = -1;
    reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, &globalSuccess_int );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }

  //
  // INSTANTIATIONS
  //

#ifdef HAVE_TEUCHOS_COMPLEX
#  define UNIT_TEST_GROUP_ORDINAL_COMPLEX_FLOAT(ORDINAL)\
     typedef std::complex<float> ComplexFloat; \
     UNIT_TEST_GROUP_ORDINAL_SCALAR(ORDINAL, ComplexFloat)
#  define UNIT_TEST_GROUP_ORDINAL_COMPLEX_DOUBLE(ORDINAL)\
     typedef std::complex<double> ComplexDouble; \
     UNIT_TEST_GROUP_ORDINAL_SCALAR(ORDINAL, ComplexDouble)
#else
#  define UNIT_TEST_GROUP_ORDINAL_COMPLEX_FLOAT(ORDINAL)
#  define UNIT_TEST_GROUP_ORDINAL_COMPLEX_DOUBLE(ORDINAL)
#endif

  // Uncomment this for really fast development cycles but make sure to comment
  // it back again before checking in so that we can test all the types.
  #define FAST_DEVELOPMENT_UNIT_TEST_BUILD

#define UNIT_TEST_GROUP_ORDINAL_SCALAR( ORDINAL, SCALAR ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MultiVector, MVTestDist, ORDINAL, SCALAR ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MultiVector, MVTestLocal, ORDINAL, SCALAR ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MultiVector, OPTestDist, ORDINAL, SCALAR ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MultiVector, OPTestLocal, ORDINAL, SCALAR )

# ifdef FAST_DEVELOPMENT_UNIT_TEST_BUILD
#    define UNIT_TEST_GROUP_ORDINAL( ORDINAL ) \
         /*UNIT_TEST_GROUP_ORDINAL_COMPLEX_FLOAT(ORDINAL)*/ \
         UNIT_TEST_GROUP_ORDINAL_SCALAR(ORDINAL, double)
     UNIT_TEST_GROUP_ORDINAL(int)
# else // not FAST_DEVELOPMENT_UNIT_TEST_BUILD

#    define UNIT_TEST_GROUP_ORDINAL( ORDINAL ) \
         UNIT_TEST_GROUP_ORDINAL_SCALAR(ORDINAL, char) \
         UNIT_TEST_GROUP_ORDINAL_SCALAR(ORDINAL, int) \
         UNIT_TEST_GROUP_ORDINAL_SCALAR(ORDINAL, float) \
         UNIT_TEST_GROUP_ORDINAL_SCALAR(ORDINAL, double) \
         UNIT_TEST_GROUP_ORDINAL_COMPLEX_FLOAT(ORDINAL) \
         UNIT_TEST_GROUP_ORDINAL_COMPLEX_DOUBLE(ORDINAL)
     UNIT_TEST_GROUP_ORDINAL(int)

     typedef short int ShortInt;
     UNIT_TEST_GROUP_ORDINAL(ShortInt)
     typedef long int LongInt;
     UNIT_TEST_GROUP_ORDINAL(LongInt)
#    ifdef HAVE_TEUCHOS_LONG_LONG_INT
        typedef long long int LongLongInt;
        UNIT_TEST_GROUP_ORDINAL(LongLongInt)
#    endif

# endif // FAST_DEVELOPMENT_UNIT_TEST_BUILD

}
