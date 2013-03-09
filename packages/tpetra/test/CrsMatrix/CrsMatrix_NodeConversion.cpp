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
  using std::endl;

  using Teuchos::TypeTraits::is_same;
  using Teuchos::as;
  using Teuchos::FancyOStream;
  using Teuchos::RCP;
  using Teuchos::ArrayRCP;
  using Teuchos::rcp;
  using Teuchos::arcp;
  using Teuchos::outArg;
  using Teuchos::arcpClone;
  using Teuchos::arrayView;
  using Teuchos::broadcast;
  using Teuchos::OrdinalTraits;
  using Teuchos::ScalarTraits;
  using Teuchos::Comm;
  using Teuchos::Array;
  using Teuchos::ArrayView;
  using Teuchos::tuple;
  using Teuchos::null;
  using Teuchos::ParameterList;
  using Teuchos::parameterList;

  using Tpetra::Map;
  using Tpetra::Vector;
  using Tpetra::CrsMatrix;
  using Tpetra::CrsGraph;
  using Tpetra::RowMatrix;
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
  // UNIT TESTS
  //

  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( CrsMatrix, NodeConversion, SCALAR, LO, GO, N2 )
  {
    typedef typename Kokkos::DefaultNode::DefaultNodeType N1;
    typedef Map<LO,GO,N1>              Map1;
    typedef CrsMatrix<SCALAR,LO,GO,N1> Mat1;
    typedef Map<LO,GO,N2>              Map2;
    typedef CrsMatrix<SCALAR,LO,GO,N2> Mat2;
    // create a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    const int numImages = comm->getSize();
    //const int myImageID = comm->getRank();
    const size_t        numLocal  = 10;
    const global_size_t numGlobal = numImages*numLocal;

    RCP<N1> n1 = getNode<N1>();
    RCP<N2> n2 = getNode<N2>();

    // create a contiguous uniform distributed map with numLocal entries per node
    RCP<const Map1> map1 = createUniformContigMapWithNode<LO,GO>(numGlobal,comm,n1);
    RCP<Mat1>       A1 = createCrsMatrix<SCALAR>(map1,3);

    // empty source, not filled
    {
      RCP<ParameterList> plClone = parameterList();
      // default: plClone->set("fillComplete clone",true);
      RCP<Mat2> A2 = A1->template clone<N2>(n2,plClone);
      TEST_EQUALITY_CONST( A2->isFillComplete(), true );
      TEST_EQUALITY_CONST( A2->isStorageOptimized(), true );
      TEST_EQUALITY_CONST( A2->getNodeNumEntries(), (size_t)0 );
      TEST_EQUALITY_CONST( A2->getCrsGraph()->getNodeAllocationSize(), (size_t)0 );
    }

    // one entry per row
    for (GO grow =map1->getMinGlobalIndex();
            grow<=map1->getMaxGlobalIndex();
            ++grow)
    {
      if (grow == map1->getMinGlobalIndex())      A1->insertGlobalValues(grow, tuple<GO>(grow,grow+1), tuple<SCALAR>(1.0,-1.0));
      else if (grow == map1->getMaxGlobalIndex()) A1->insertGlobalValues(grow, tuple<GO>(grow-1,grow), tuple<SCALAR>(-1.0,1.0));
      else                                        A1->insertGlobalValues(grow, tuple<GO>(grow-1,grow,grow+1), tuple<SCALAR>(-1.0,1.0,-1.0));
    }
    // source has global indices, not filled, dynamic profile
    {
      RCP<ParameterList> plClone = parameterList();
      plClone->set("fillComplete clone",false);
      plClone->set("Static profile clone",false);
      // default: plClone->set("Locally indexed clone",false);
      RCP<Mat2> A2 = A1->template clone<N2>(n2,plClone);
      TEST_EQUALITY_CONST( A2->hasColMap(), false );
      TEST_EQUALITY_CONST( A2->isFillComplete(), false );
      TEST_EQUALITY_CONST( A2->isGloballyIndexed(), true );
      TEST_EQUALITY_CONST( A2->getCrsGraph()->getNodeAllocationSize(), (size_t)(numLocal*3-2) );
      TEST_EQUALITY( A2->getNodeNumEntries(), A1->getNodeNumEntries() );
      TEST_NOTHROW( A2->insertGlobalValues(map1->getMaxLocalIndex(), tuple<GO>(map1->getMinLocalIndex()), tuple<SCALAR>(1.0)) );
      TEST_NOTHROW( A2->insertGlobalValues(map1->getMinLocalIndex(), tuple<GO>(map1->getMaxLocalIndex()), tuple<SCALAR>(1.0)) );
      TEST_NOTHROW( A2->fillComplete() );
      TEST_EQUALITY_CONST( A2->getNodeNumEntries(), A1->getNodeNumEntries()+2 );
    }

    // source has local indices
    A1->fillComplete();

    {
      RCP<ParameterList> plClone = parameterList();
      plClone->set("Static profile clone", false);
      RCP<ParameterList> plCloneFill = sublist(plClone,"fillComplete");
      plCloneFill->set("Optimize Storage",false);
      RCP<Mat2> A2 = A1->template clone<N2>(n2,plClone);
      TEST_EQUALITY_CONST( A2->isFillComplete(), true );
      TEST_EQUALITY_CONST( A2->isStorageOptimized(), false );
      A2->resumeFill();
      for (LO lrow = map1->getMinLocalIndex();
              lrow < map1->getMaxLocalIndex()-1;
              ++lrow)
      {
        TEST_NOTHROW( A2->insertLocalValues(lrow, tuple<LO>(lrow+2), tuple<SCALAR>(1.0)) );
      }
      A2->fillComplete();
      TEST_EQUALITY_CONST( A2->isFillComplete(), true );
      TEST_EQUALITY_CONST( A2->isStorageOptimized(), true );
      TEST_EQUALITY_CONST( A2->getNodeNumEntries(), A1->getNodeNumEntries()+numLocal-2 );
    }

  }

//
// INSTANTIATIONS
//

#define NC_TESTS(D,LO,GO,N2) \
    TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( CrsMatrix, NodeConversion, double, int, int, N2 )

  TPETRA_ETI_MANGLING_TYPEDEFS()

  TPETRA_INSTANTIATE_DOUBLE_INT_INT_N(NC_TESTS)
}
