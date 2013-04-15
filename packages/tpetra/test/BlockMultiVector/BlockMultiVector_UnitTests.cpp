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

#include <Tpetra_TestingUtilities.hpp>

#include <Tpetra_BlockMultiVector.hpp>
#include <Tpetra_Vector.hpp>
#include <iterator>

// FINISH: add test for BlockMultiVector with a node containing zero local entries
// FINISH: add tests for local BlockMultiVectors

namespace Teuchos {
  template <>
  ScalarTraits<int>::magnitudeType
  relErr( const int &s1, const int &s2 ) {
    typedef ScalarTraits<int> ST;
    return ST::magnitude(s1-s2);
  }

  template <>
  ScalarTraits<char>::magnitudeType
  relErr( const char &s1, const char &s2 ) {
    typedef ScalarTraits<char> ST;
    return ST::magnitude(s1-s2);
  }
}

namespace {

  using Tpetra::TestingUtilities::getNode;
  using Tpetra::TestingUtilities::getDefaultComm;

  using std::endl;
  using std::copy;
  using std::ostream_iterator;
  using std::string;

  using Teuchos::TypeTraits::is_same;
  using Teuchos::RCP;
  using Teuchos::ArrayRCP;
  using Teuchos::rcp;
  using Teuchos::null;
  using Teuchos::Array;
  using Teuchos::ArrayView;
  using Teuchos::Comm;
  using Teuchos::Range1D;
  using Teuchos::Tuple;
  using Teuchos::as;
  using Teuchos::OrdinalTraits;
  using Teuchos::ScalarTraits;
  using Teuchos::arrayView;
  using Teuchos::tuple;
  using Teuchos::NO_TRANS;
  using Teuchos::TRANS;
  using Teuchos::CONJ_TRANS;
  using Teuchos::VERB_DEFAULT;
  using Teuchos::VERB_NONE;
  using Teuchos::VERB_LOW;
  using Teuchos::VERB_MEDIUM;
  using Teuchos::VERB_HIGH;
  using Teuchos::VERB_EXTREME;

  using Tpetra::BlockMap;
  using Tpetra::BlockMultiVector;
  using Tpetra::global_size_t;
  using Tpetra::DefaultPlatform;
  using Tpetra::GloballyDistributed;

  using Tpetra::createLocalMapWithNode;

  double errorTolSlack = 1e+1;

  TEUCHOS_STATIC_SETUP()
  {
    Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
    clp.addOutputSetupOptions(true);
    clp.setOption(
        "test-mpi", "test-serial", &Tpetra::TestingUtilities::testMpi,
        "Test MPI (if available) or force test of serial.  In a serial build,"
        " this option is ignored and a serial comm is always used." );
    clp.setOption(
        "error-tol-slack", &errorTolSlack,
        "Slack off of machine epsilon used to check test results" );
  }


  //
  // UNIT TESTS
  //

  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( BlockMultiVector, basic, LO, GO, Scalar , Node )
  {
    RCP<Node> node = getNode<Node>();
    typedef Tpetra::BlockMultiVector<Scalar,LO,GO,Node> BMV;
    typedef typename ScalarTraits<Scalar>::magnitudeType Magnitude;
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    // get a comm and node
    RCP<const Comm<int> > comm = getDefaultComm();
    const int numImages = comm->getSize();
    // create a BlockMap
    const size_t numLocal = 12;
    const size_t numLocalBlocks = 4;
    const LO blockSize = 3;
    const size_t numVecs  = 7;
    RCP<BlockMap<LO,GO,Node> > blkmap = rcp( new BlockMap<LO,GO,Node>(INVALID,numLocalBlocks,blockSize,0,comm,node) );
    BMV bmvec(blkmap,numVecs,true);
    TEST_EQUALITY( bmvec.getNumVectors(), numVecs );
    TEST_EQUALITY( bmvec.getLocalLength(), numLocal );
    TEST_EQUALITY( bmvec.getGlobalLength(), numImages*numLocal );
    // we zeroed it out in the constructor; all norms should be zero
    Array<Magnitude> norms(numVecs), zeros(numVecs);
    std::fill(zeros.begin(),zeros.end(),ScalarTraits<Magnitude>::zero());
    bmvec.norm2(norms);
    TEST_COMPARE_FLOATING_ARRAYS(norms,zeros,ScalarTraits<Magnitude>::zero());
    bmvec.norm1(norms);
    TEST_COMPARE_FLOATING_ARRAYS(norms,zeros,ScalarTraits<Magnitude>::zero());
    bmvec.normInf(norms);
    TEST_COMPARE_FLOATING_ARRAYS(norms,zeros,ScalarTraits<Magnitude>::zero());
    // print it
    out << bmvec << endl;
  }

  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( BlockMultiVector, OrthoDot, LO, GO, Scalar , Node )
  {
    RCP<Node> node = getNode<Node>();
    typedef typename ScalarTraits<Scalar>::magnitudeType Mag;
    typedef Tpetra::BlockMultiVector<Scalar,LO,GO,Node> BMV;
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    const Scalar S0 = ScalarTraits<Scalar>::zero();
    const Mag M0 = ScalarTraits<Mag>::zero();
    // get a comm and node
    RCP<const Comm<int> > comm = getDefaultComm();
    const int numImages = comm->getSize();
    // create a Map
    const LO indexBase = 0;
    const size_t numLocalBlocks = 4;
    const LO blockSize = 3;
    const size_t numVectors = 3;
    RCP<BlockMap<LO,GO,Node> > blkmap = rcp( new BlockMap<LO,GO,Node>(INVALID,numLocalBlocks,blockSize,indexBase,comm,node) );
    const bool zeroOut = true;
    BMV bmvec1(blkmap,numVectors,zeroOut),
        bmvec2(blkmap,numVectors,zeroOut);
    Array<Scalar> dots1(numVectors), dots2(numVectors), zeros(numVectors);
    Array<Mag>    norms1(numVectors), norms2(numVectors), ans(numVectors);
    std::fill(zeros.begin(),zeros.end(),ScalarTraits<Scalar>::zero());
    // these should be numerically orthogonal even in finite arithmetic, because both are zero. 1-norms are zero.
    bmvec1.dot(bmvec2,dots1());
    bmvec2.dot(bmvec1,dots2());
    TEST_COMPARE_FLOATING_ARRAYS(dots2,zeros,M0);
    TEST_COMPARE_FLOATING_ARRAYS(dots1,zeros,M0);
    TEST_EQUALITY_CONST( bmvec1.getVector(0)->dot(*bmvec2.getVector(0)), S0);
    bmvec1.norm1(norms1());
    bmvec2.norm1(norms2());
    std::fill(ans.begin(), ans.end(), M0);
    TEST_COMPARE_FLOATING_ARRAYS(norms1,ans,M0);
    TEST_COMPARE_FLOATING_ARRAYS(norms1,ans,M0);
    // replace local entries s.t.
    // mvec1 = [1 1]  and  mvec2 = [0 0]
    //         [0 0]               [1 1]
    // still numerically orthogonal even in finite arithmetic. norms are numImages.
    for (size_t j=0; j < numVectors; ++j) {
      bmvec1.replaceLocalValue(0, 1, j,ScalarTraits<Scalar>::one());
      bmvec2.replaceGlobalValue(blkmap->getGlobalBlockID(1), 1, j,ScalarTraits<Scalar>::one());
    }
    bmvec1.dot(bmvec2,dots1());
    bmvec2.dot(bmvec1,dots2());
    TEST_COMPARE_FLOATING_ARRAYS(dots2,zeros,M0);
    TEST_COMPARE_FLOATING_ARRAYS(dots1,zeros,M0);
    TEST_EQUALITY_CONST( bmvec1.getVector(0)->dot(*bmvec2.getVector(0)), S0);
    bmvec1.norm1(norms1());
    bmvec2.norm1(norms2());
    std::fill(ans.begin(), ans.end(), as<Mag>(numImages));
    TEST_COMPARE_FLOATING_ARRAYS(norms1,ans,M0);
    TEST_COMPARE_FLOATING_ARRAYS(norms2,ans,M0);
    // sum into local entries s.t.
    // mvec1 = [1 1]  and  mvec2 = [-1 -1]
    //         [1 1]               [ 1  1]
    // still numerically orthogonal even in finite arithmetic. norms are 2*numImages.
    for (size_t j=0; j < numVectors; ++j) {
      bmvec1.sumIntoLocalValue(1, 1, j,ScalarTraits<Scalar>::one());
      bmvec2.sumIntoGlobalValue(blkmap->getGlobalBlockID(0), 1, j,-ScalarTraits<Scalar>::one());
    }
    bmvec1.dot(bmvec2,dots1());
    bmvec2.dot(bmvec1,dots2());
    TEST_COMPARE_FLOATING_ARRAYS(dots2,zeros,M0);
    TEST_COMPARE_FLOATING_ARRAYS(dots1,zeros,M0);
    TEST_EQUALITY_CONST( bmvec1.getVector(0)->dot(*bmvec2.getVector(0)), S0);
    bmvec1.norm1(norms1());
    bmvec2.norm1(norms2());
    std::fill(ans.begin(), ans.end(), as<Mag>(2*numImages));
    TEST_COMPARE_FLOATING_ARRAYS(norms1,ans,M0);
    TEST_COMPARE_FLOATING_ARRAYS(norms2,ans,M0);
  }


//
// INSTANTIATIONS
//

#define UNIT_TEST_GROUP( SCALAR, LO, GO, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( BlockMultiVector, basic             , LO, GO, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( BlockMultiVector, OrthoDot          , LO, GO, SCALAR, NODE )

  TPETRA_ETI_MANGLING_TYPEDEFS()

  TPETRA_INSTANTIATE_TESTMV_NOGPU( UNIT_TEST_GROUP )

}
