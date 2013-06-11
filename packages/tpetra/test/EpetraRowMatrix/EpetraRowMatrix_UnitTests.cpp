/*
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
*/

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_OrdinalTraits.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_VerboseObject.hpp>
#include <Teuchos_oblackholestream.hpp>
#include <Teuchos_FancyOStream.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_as.hpp>

#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_DefaultPlatform.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_EpetraRowMatrix.hpp"

#ifdef HAVE_MPI
#  include "Epetra_MpiComm.h"
#else
#  include "Epetra_SerialComm.h"
#endif


namespace Teuchos {
  template <>
    ScalarTraits<int>::magnitudeType
    relErr( const int &s1, const int &s2 )
    {
      typedef ScalarTraits<int> ST;
      return ST::magnitude(s1-s2);
    }

  template <>
    ScalarTraits<char>::magnitudeType
    relErr( const char &s1, const char &s2 )
    {
      typedef ScalarTraits<char> ST;
      return ST::magnitude(s1-s2);
    }
}

namespace {

  using Teuchos::RCP;
  using Teuchos::ArrayRCP;
  using Teuchos::rcp;
  using Teuchos::arcpClone;
  using Tpetra::Map;
  using Tpetra::OptimizeOption;
  using Tpetra::DoOptimizeStorage;
  using Tpetra::DoNotOptimizeStorage;
  using Tpetra::DefaultPlatform;
  using Tpetra::global_size_t;
  using std::sort;
  using Teuchos::arrayView;
  using Teuchos::broadcast;
  using Teuchos::OrdinalTraits;
  using Teuchos::ScalarTraits;
  using Teuchos::Comm;
  using Tpetra::MultiVector;
  using Tpetra::Vector;
  using std::endl;
  using std::swap;
  using Teuchos::Array;
  using Teuchos::ArrayView;
  using Tpetra::Operator;
  using Tpetra::CrsMatrix;
  using Tpetra::CrsGraph;
  using Tpetra::RowMatrix;
  using Tpetra::INSERT;
  using Tpetra::Import;
  using std::string;
  using Teuchos::tuple;
  using Teuchos::VERB_NONE;
  using Teuchos::VERB_LOW;
  using Teuchos::VERB_MEDIUM;
  using Teuchos::VERB_HIGH;
  using Teuchos::VERB_EXTREME;
  using Tpetra::ProfileType;
  using Tpetra::StaticProfile;
  using Tpetra::DynamicProfile;
  using Teuchos::ETransp;
  using Teuchos::NO_TRANS;
  using Teuchos::TRANS;
  using Teuchos::CONJ_TRANS;
  using Teuchos::EDiag;
  using Teuchos::UNIT_DIAG;
  using Teuchos::NON_UNIT_DIAG;
  using Teuchos::EUplo;
  using Teuchos::UPPER_TRI;
  using Teuchos::LOWER_TRI;
  using Tpetra::LocallyReplicated;
  using Tpetra::GloballyDistributed;

  typedef DefaultPlatform::DefaultPlatformType::NodeType Node;

  bool testMpi = true;
  double errorTolSlack = 1e+1;
  string filedir;


  TEUCHOS_STATIC_SETUP()
  {
    Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
    clp.setOption(
        "filedir",&filedir,"Directory of expected matrix files.");
    clp.addOutputSetupOptions(true);
    clp.setOption(
        "test-mpi", "test-serial", &testMpi,
        "Test MPI (if available) or force test of serial.  In a serial build,"
        " this option is ignored and a serial comm is always used." );
    clp.setOption(
        "error-tol-slack", &errorTolSlack,
        "Slack off of machine epsilon used to check test results" );
  }

  RCP<const Comm<int> > getDefaultComm()
  {
    RCP<const Comm<int> > ret;
    if (testMpi) {
      ret = DefaultPlatform::getDefaultPlatform().getComm();
    }
    else {
      ret = rcp(new Teuchos::SerialComm<int>());
    }
    return ret;
  }

  //
  // UNIT TESTS
  // 

  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( EpetraRowMatrix, BasicFunctionality, LO, GO, Scalar )
  {
    // generate a tridiagonal matrix
    typedef ScalarTraits<Scalar> ST;
    typedef CrsMatrix<Scalar,LO,GO,Node> MAT;
    typedef MultiVector<Scalar,LO,GO,Node> MV;
    typedef typename ST::magnitudeType Mag;
    typedef ScalarTraits<Mag> MT;
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    // get a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    const int numImages = size(*comm);
    // create a Map
    const size_t numLocal = 10;
    const size_t numVecs = 5;
    RCP<const Map<LO,GO> > map = Tpetra::createContigMap<LO,GO>(INVALID,numLocal,comm);
    // create a matrix, modeled closely on Chris' CrsMatrix unit-tests.
    RCP<MAT> matrix(new MAT(map, 3));
    for (GO r=map->getMinGlobalIndex(); r <= map->getMaxGlobalIndex(); ++r) {
      if (r == map->getMinAllGlobalIndex()) {
        matrix->insertGlobalValues(r, tuple(r,r+1), tuple(ST::one(),ST::one()) );
      }
      else if (r == map->getMaxAllGlobalIndex()) {
        matrix->insertGlobalValues(r, tuple(r-1,r), tuple(ST::one(),ST::one()) );
      }
      else {
        matrix->insertGlobalValues(r, tuple(r-1,r,r+1), tuple(ST::one(),ST::one(),ST::one()) );
      }
    }
    for (GO r=map->getMinGlobalIndex(); r <= map->getMaxGlobalIndex(); ++r) {
      // increment the diagonals
      matrix->sumIntoGlobalValues(r, tuple(r), tuple(ST::one()) );
    }
    matrix->fillComplete();
    TEST_EQUALITY( matrix->getNodeNumDiags(), numLocal );
    TEST_EQUALITY( matrix->getGlobalNumDiags(), numImages*numLocal );
    TEST_EQUALITY( matrix->getGlobalNumEntries(), 3*numImages*numLocal - 2 );

#ifdef HAVE_MPI
    Epetra_MpiComm ecomm(MPI_COMM_WORLD);
#else
    Epetra_SerialComm ecomm;
#endif

    //now create the EpetraRowMatrix class, which is the real target of this
    //unit-test:
    Tpetra::EpetraRowMatrix<MAT> erowmat(matrix, ecomm);

    int myRow = numLocal/2;
    int erowmat_NumMyRowEntries = -1;
    TEST_EQUALITY_CONST( erowmat.NumMyRowEntries( myRow, erowmat_NumMyRowEntries ) == 0, true );

    size_t numEntriesLocalRow = matrix->getNumEntriesInLocalRow(myRow);

    TEST_EQUALITY_CONST( erowmat_NumMyRowEntries == (int)numEntriesLocalRow, true );

    //test the matrix-vector product by comparing the result of
    //CrsMatrix::apply with the result of Epetra_RowMatrix::Multiply:

    MV tmv1(map,numVecs,true), tmv2(map,numVecs,true);
    tmv1.randomize();
    tmv2.randomize();
    tmv1.putScalar(1.0);
    tmv2.putScalar(0.0);
    matrix->apply(tmv1,tmv2);

    Epetra_BlockMap emap(int(numImages*numLocal), 1, 0, ecomm);
    Epetra_MultiVector emv1(emap,numVecs), emv2(emap,numVecs);
    emv1.PutScalar(1.0);
    emv2.PutScalar(0.0);
    erowmat.Multiply(false, emv1, emv2);

    ArrayView<Scalar> vals(emv2.Values(),numLocal*numVecs);
    MV tmvres(map,vals,numLocal,numVecs);
    tmvres.update(-ST::one(),tmv2,ST::one());
    Array<Mag> norms(numVecs), zeros(numVecs,MT::zero());
    tmvres.norm2(norms());
    TEST_COMPARE_FLOATING_ARRAYS(norms,zeros,MT::zero());
  }


  // 
  // INSTANTIATIONS
  //

  // Uncomment this for really fast development cycles but make sure to comment
  // it back again before checking in so that we can test all the types.
  #define FAST_DEVELOPMENT_UNIT_TEST_BUILD

#define UNIT_TEST_GROUP_ORDINAL_SCALAR( LO, GO, SCALAR ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( EpetraRowMatrix, BasicFunctionality     , LO, GO, SCALAR )

#define UNIT_TEST_GROUP_ORDINAL( ORDINAL ) \
    UNIT_TEST_GROUP_ORDINAL_ORDINAL( ORDINAL, ORDINAL )

# ifdef FAST_DEVELOPMENT_UNIT_TEST_BUILD
#    define UNIT_TEST_GROUP_ORDINAL_ORDINAL( LO, GO ) \
         UNIT_TEST_GROUP_ORDINAL_SCALAR( LO, GO, double)
     UNIT_TEST_GROUP_ORDINAL(int)

# else // not FAST_DEVELOPMENT_UNIT_TEST_BUILD

#    define UNIT_TEST_GROUP_ORDINAL_ORDINAL( LO, GO ) \
         UNIT_TEST_GROUP_ORDINAL_SCALAR(LO, GO, float)  \
         UNIT_TEST_GROUP_ORDINAL_SCALAR(LO, GO, double)

     UNIT_TEST_GROUP_ORDINAL(int)

     typedef long int LongInt;
     UNIT_TEST_GROUP_ORDINAL_ORDINAL( int, LongInt )
#    ifdef HAVE_TEUCHOS_LONG_LONG_INT
        typedef long long int LongLongInt;
        UNIT_TEST_GROUP_ORDINAL_ORDINAL( int, LongLongInt )
#    endif

# endif // FAST_DEVELOPMENT_UNIT_TEST_BUILD

}
