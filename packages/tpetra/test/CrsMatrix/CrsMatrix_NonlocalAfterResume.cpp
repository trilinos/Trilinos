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
#include <Tpetra_CrsMatrix.hpp>

namespace {
  using Tpetra::TestingUtilities::getNode;
  using Tpetra::TestingUtilities::getDefaultComm;

  using Teuchos::Array;
  using Teuchos::ArrayRCP;
  using Teuchos::ArrayView;
  using Teuchos::as;
  using Teuchos::FancyOStream;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::outArg;
  using Teuchos::broadcast;
  using Teuchos::OrdinalTraits;
  using Teuchos::ScalarTraits;
  using Teuchos::Comm;
  using Teuchos::tuple;
  using Teuchos::null;
  using Teuchos::ParameterList;
  using Teuchos::parameterList;

  using Tpetra::Map;
  using Tpetra::CrsMatrix;
  using Tpetra::Import;
  using Tpetra::global_size_t;
  using Tpetra::createNonContigMapWithNode;
  using Tpetra::createUniformContigMapWithNode;
  using Tpetra::createContigMapWithNode;
  using Tpetra::createLocalMapWithNode;
  using Tpetra::createVector;
  using Tpetra::createCrsMatrix;
  using Tpetra::DefaultPlatform;
  using Tpetra::ProfileType;
  using Tpetra::StaticProfile;
  using Tpetra::DynamicProfile;
  using Tpetra::OptimizeOption;
  using Tpetra::DoOptimizeStorage;
  using Tpetra::DoNotOptimizeStorage;
  using Tpetra::GloballyDistributed;
  using Tpetra::INSERT;

  TEUCHOS_STATIC_SETUP()
  {
    Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
    clp.addOutputSetupOptions(true);
    clp.setOption(
        "test-mpi", "test-serial", &Tpetra::TestingUtilities::testMpi,
        "Test MPI (if available) or force test of serial.  In a serial build,"
        " this option is ignored and a serial comm is always used." );
  }

//
// UNIT TEST(S)
//

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( CrsMatrix, NonlocalAfterResume, LO, GO, Scalar, Node )
{
  RCP<Node> node = getNode<Node>();
  // test that an exception is thrown when we exceed statically allocated memory
  typedef ScalarTraits<Scalar> ST;
  const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
  // get a comm
  RCP<const Comm<int> > comm = getDefaultComm();
  const size_t numImages = size(*comm);
  const size_t myImageID = rank(*comm);
  // create a row Map, 5 rows per processor
  const GO numLocal = 5;
  RCP<const Map<LO,GO,Node> > rmap = createContigMapWithNode<LO,GO,Node>(INVALID,numLocal,comm,node);
  RCP<const Map<LO,GO,Node> > cmap;
  // create a column Map, with super- and sub-diagonal blocks
  {
    Array<GO> cols;
    for (GO c=rmap->getMinGlobalIndex(); c <= rmap->getMaxGlobalIndex(); ++c) cols.push_back(c);
    if (rmap->getMinGlobalIndex() >= rmap->getMinAllGlobalIndex() + numLocal) {
      for (GO c = rmap->getMinGlobalIndex()-numLocal; c < rmap->getMinGlobalIndex(); ++c) {
        cols.push_back(c);
      }
    }
    if (rmap->getMaxGlobalIndex()+numLocal <= rmap->getMaxAllGlobalIndex()) {
      for (GO c = rmap->getMaxGlobalIndex()+1; c <= rmap->getMaxGlobalIndex()+numLocal; ++c) {
        cols.push_back(c);
      }
    }
    cmap = createNonContigMapWithNode<LO,GO,Node>(cols(), comm, node);
  }
  {
    //----------------------------------------------------------------------
    // put in diagonal, locally
    //----------------------------------------------------------------------
    Tpetra::CrsMatrix<Scalar,LO,GO,Node> matrix(rmap,cmap,3,DynamicProfile);
    for (GO r=rmap->getMinGlobalIndex(); r <= rmap->getMaxGlobalIndex(); ++r) {
      matrix.insertGlobalValues(r,tuple(r),tuple(ST::one()));
    }
    // fill, but do not pack, because we will add new entries below
    RCP<ParameterList> params = parameterList(); 
    params->set("Optimize Storage",false);
    TEST_NOTHROW       ( matrix.fillComplete( params ) );
    TEST_EQUALITY_CONST( matrix.isFillComplete(),      true );
    TEST_EQUALITY_CONST( matrix.isStorageOptimized(), false );
    TEST_EQUALITY      ( matrix.getGlobalNumEntries(), numLocal*numImages );
    TEST_EQUALITY      ( matrix.getNodeNumEntries(),   (size_t)numLocal   );

    //----------------------------------------------------------------------
    // add super-diagonal, non-locally
    //----------------------------------------------------------------------
    // because fillComplete() was called above, we must call resumeFill() before adding new entries
    matrix.resumeFill();
    if (rmap->getMinGlobalIndex()+numLocal < rmap->getMaxAllGlobalIndex()) {
      for (GO r=rmap->getMinGlobalIndex(); r <= rmap->getMaxGlobalIndex(); ++r) {
        matrix.insertGlobalValues(r+numLocal,tuple(r),tuple(ST::one()));
      }
    }
    // fill, but do not pack, because we will add new entries below
    params->set("Optimize Storage",false);
    TEST_NOTHROW       ( matrix.fillComplete( params ) );
    TEST_EQUALITY_CONST( matrix.isFillComplete(),      true );
    TEST_EQUALITY_CONST( matrix.isStorageOptimized(), false );
    TEST_EQUALITY      ( matrix.getGlobalNumEntries(), 2*numLocal*numImages-numLocal );
    {
      size_t expected = numLocal;
      if (myImageID > 0) expected += numLocal; // super-diagonal
      TEST_EQUALITY( matrix.getNodeNumEntries(), expected );
    }

    //----------------------------------------------------------------------
    // add sub-diagonal block, non-locally
    //----------------------------------------------------------------------
    // because fillComplete() was called above, we must call resumeFill() before adding new entries
    matrix.resumeFill();
    if (rmap->getMinGlobalIndex() >= rmap->getMinAllGlobalIndex()+numLocal) {
      for (GO r=rmap->getMinGlobalIndex(); r <= rmap->getMaxGlobalIndex(); ++r) {
        matrix.insertGlobalValues(r-numLocal,tuple(r),tuple(ST::one()));
      }
    }
    // fill; it is okay to pack now
    params->set("Optimize Storage",true);
    TEST_NOTHROW       ( matrix.fillComplete( params ) );
    TEST_EQUALITY_CONST( matrix.isFillComplete(), true );
    TEST_EQUALITY_CONST( matrix.isStorageOptimized(), true );
    TEST_EQUALITY      ( matrix.getGlobalNumEntries(), 3*numLocal*numImages-2*numLocal );
    {
      size_t expected = numLocal;
      if (myImageID > 0)           expected += numLocal; // super-diagonal
      if (myImageID < numImages-1) expected += numLocal; // sub-diagonal
      TEST_EQUALITY( matrix.getNodeNumEntries(), expected );
    }
  }
  // All procs fail if any node fails
  int globalSuccess_int = -1;
  reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, outArg(globalSuccess_int) );
  TEST_EQUALITY_CONST( globalSuccess_int, 0 );
}

//
// INSTANTIATIONS
//

#define UNIT_TEST_GROUP( SCALAR, LO, GO, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( CrsMatrix, NonlocalAfterResume, LO, GO, SCALAR, NODE )

  TPETRA_ETI_MANGLING_TYPEDEFS()

  TPETRA_INSTANTIATE_SLGN( UNIT_TEST_GROUP )
}
