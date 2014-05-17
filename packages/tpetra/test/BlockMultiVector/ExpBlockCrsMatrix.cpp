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

// Some Macro Magic to ensure that if CUDA and KokkosCompat is enabled
// only the .cu version of this file is actually compiled
#include <Tpetra_config.h>
#ifdef HAVE_TPETRA_KOKKOSCOMPAT
#include <KokkosCore_config.h>
#ifdef KOKKOS_USE_CUDA_BUILD
  #define DO_COMPILATION
#else
  #ifndef KOKKOS_HAVE_CUDA
    #define DO_COMPILATION
  #endif
#endif
#else
  #define DO_COMPILATION
#endif

#ifdef DO_COMPILATION

#include <Tpetra_TestingUtilities.hpp>
#include <Tpetra_Experimental_BlockCrsMatrix.hpp>

namespace {

  //
  // UNIT TESTS
  //

  // Test BlockCrsMatrix's constructors.
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( ExpBlockCrsMatrix, ctor, Scalar, LO, GO, Node )
  {
    using Tpetra::TestingUtilities::getNode;
    using Tpetra::TestingUtilities::getDefaultComm;
    using Teuchos::Comm;
    using Teuchos::RCP;
    using Teuchos::rcp;
    using std::endl;
    typedef Tpetra::Experimental::BlockMultiVector<Scalar, LO, GO, Node> BMV;
    typedef Tpetra::Experimental::BlockCrsMatrix<Scalar, LO, GO, Node> BCM;
    typedef Tpetra::CrsGraph<LO, GO, Node> graph_type;
    typedef Tpetra::Map<LO, GO, Node> map_type;
    typedef Tpetra::global_size_t GST;

    out << "Testing Tpetra::Experimental::BlockCrsMatrix" << endl;
    Teuchos::OSTab tab0 (out);

    RCP<const Comm<int> > comm = getDefaultComm ();
    RCP<Node> node = getNode<Node> ();
    const GST INVALID = Teuchos::OrdinalTraits<GST>::invalid ();

    out << "Creating mesh row Map" << endl;

    const size_t numLocalMeshPoints = 12;
    const GO indexBase = 0;
    // mfh 16 May 2014: Tpetra::CrsGraph still needs the row Map as an
    // RCP.  Later interface changes will let us pass in the Map by
    // const reference and assume view semantics.
    RCP<const map_type> meshRowMapPtr =
      rcp (new map_type (INVALID, numLocalMeshPoints, indexBase, comm, node));

    const LO blockSize = 4;
    //const LO numVecs = 3;

    // mfh 16 May 2014: Make a graph.  This graph is empty; that's OK.
    // We just want to check that the constructors of BlockCrsMatrix
    // work.
    out << "Creating mesh graph" << endl;
    graph_type graph (meshRowMapPtr, 0, Tpetra::StaticProfile);
    graph.fillComplete ();

    // Test the default constructor.  We can't wrap this in a
    // TEST_NOTHROW, because the variable need to be in scope for
    // tests below.
    out << "Testing default constructor" << endl;
    BCM blockMat;

    // Test the two-argument constructor.
    out << "Testing two-argument constructor" << endl;
    TEST_NOTHROW( blockMat = BCM (graph, blockSize) );

    // Test the four-argument constructor.
    // Doing this requires the domain and range point Maps.
    out << "Testing four-argument constructor" << endl;
    map_type pointDomainMap = BMV::makePointMap (* (graph.getDomainMap ()), blockSize);
    TEST_ASSERT( ! blockMat.getDomainMap ().is_null () &&
                 pointDomainMap.isSameAs (* (blockMat.getDomainMap ())) );
    map_type pointRangeMap = BMV::makePointMap (* (graph.getRangeMap ()), blockSize);
    TEST_ASSERT( ! blockMat.getRangeMap ().is_null () &&
                 pointRangeMap.isSameAs (* (blockMat.getRangeMap ())) );
    TEST_NOTHROW( blockMat = BCM (graph, pointDomainMap, pointRangeMap, blockSize ) );
  }

//
// INSTANTIATIONS
//

#define UNIT_TEST_GROUP( SCALAR, LO, GO, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( ExpBlockCrsMatrix, ctor, SCALAR, LO, GO, NODE )

  TPETRA_ETI_MANGLING_TYPEDEFS()

  TPETRA_INSTANTIATE_TESTMV_NOGPU( UNIT_TEST_GROUP )

} // namespace (anonymous)

#endif  //DO_COMPILATION

