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

#include <Tpetra_Experimental_BlockMultiVector.hpp>

namespace {

  //
  // UNIT TESTS
  //

  // Make sure that BlockMultiVector's constructor actually compiles.
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( ExpBlockMultiVector, ctor, Scalar, LO, GO, Node )
  {
    using Tpetra::TestingUtilities::getNode;
    using Tpetra::TestingUtilities::getDefaultComm;
    using Teuchos::Comm;
    using Teuchos::RCP;
    typedef Tpetra::Experimental::BlockMultiVector<Scalar, LO, GO, Node> BMV;
    typedef Tpetra::Map<LO, GO, Node> map_type;
    typedef Tpetra::global_size_t GST;

    RCP<const Comm<int> > comm = getDefaultComm ();
    RCP<Node> node = getNode<Node> ();
    const GST INVALID = Teuchos::OrdinalTraits<GST>::invalid ();

    const size_t numLocalMeshPoints = 12;
    const GO indexBase = 0;
    map_type meshMap (INVALID, numLocalMeshPoints, indexBase, comm, node);
    //RCP<const map_type> mapPtr = Teuchos::rcpFromRef (map); // nonowning RCP

    const LO blockSize = 4;
    const LO numVecs = 3;

    BMV X (meshMap, blockSize, numVecs);
    TEST_EQUALITY( X.getBlockSize (), blockSize );
    TEST_EQUALITY( X.getNumVectors (), numVecs );
    TEST_EQUALITY_CONST( X.getMap ().is_null (), false );
    TEST_EQUALITY_CONST( ! X.getMap ().is_null () &&
                         X.getMap ()->isSameAs (meshMap), true );

    map_type pointMap = BMV::makePointMap (meshMap, blockSize);
    TEST_EQUALITY_CONST( pointMap.isSameAs (X.getPointMap ()), true );

    BMV Y (meshMap, pointMap, blockSize, numVecs);
    TEST_EQUALITY( Y.getBlockSize (), blockSize );
    TEST_EQUALITY( Y.getNumVectors (), numVecs );
    TEST_EQUALITY_CONST( Y.getMap ().is_null (), false );
    TEST_EQUALITY_CONST( ! Y.getMap ().is_null () &&
                         Y.getMap ()->isSameAs (meshMap), true );

    TEST_EQUALITY_CONST( Y.getPointMap ().isSameAs (X.getPointMap ()), true );

    BMV Z; // Exercise the default constructor.

    TEST_EQUALITY_CONST( Z.getMap ().is_null (), true );
    TEST_EQUALITY_CONST( Z.getBlockSize (), static_cast<LO> (0) );
    TEST_EQUALITY_CONST( Z.getNumVectors (), static_cast<LO> (0) );
  }

//
// INSTANTIATIONS
//

#define UNIT_TEST_GROUP( SCALAR, LO, GO, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( ExpBlockMultiVector, ctor, SCALAR, LO, GO, NODE )

  TPETRA_ETI_MANGLING_TYPEDEFS()

  TPETRA_INSTANTIATE_TESTMV_NOGPU( UNIT_TEST_GROUP )

} // namespace (anonymous)

#endif  //DO_COMPILATION

