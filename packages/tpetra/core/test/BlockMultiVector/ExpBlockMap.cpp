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

#include <Tpetra_ConfigDefs.hpp>
#include <Tpetra_TestingUtilities.hpp>
#include <Tpetra_Experimental_BlockMultiVector.hpp>

namespace {

  //
  // UNIT TESTS
  //

  // Test makePointMap.  In particular, if a GID g_mesh is in the
  // input mesh Map on a process, and if the block size is b, then the
  // GID g_point := g_mesh * b must be in the output point Map on that
  // process.
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( ExpBlockMap, HilbertsHotel, Scalar, LO, GO, Node )
  {
    using Tpetra::TestingUtilities::getNode;
    using Tpetra::TestingUtilities::getDefaultComm;
    using Teuchos::Array;
    using Teuchos::ArrayView;
    using Teuchos::Comm;
    using Teuchos::outArg;
    using Teuchos::REDUCE_MIN;
    using Teuchos::reduceAll;
    using Teuchos::RCP;
    typedef Tpetra::Experimental::BlockMultiVector<Scalar, LO, GO, Node> BMV;
    typedef Tpetra::Map<LO, GO, Node> map_type;
    typedef Tpetra::global_size_t GST;
    typedef typename Array<GO>::size_type size_type;

    RCP<const Comm<int> > comm = getDefaultComm ();
    RCP<Node> node = getNode<Node> ();
    const GST INVALID = Teuchos::OrdinalTraits<GST>::invalid ();

    // Create a nonoverlapping, contiguous mesh Map.
    const size_t numLocalMeshPoints = 2;
    const LO blockSize = 5;
    const GO indexBase = 0;
    map_type meshMap (INVALID, numLocalMeshPoints, indexBase, comm, node);

    TEST_ASSERT( meshMap.isContiguous () );
    TEST_ASSERT( meshMap.isOneToOne () );

    // Create a nonoverlapping point Map from the contiguous mesh Map.
    map_type pointMap = BMV::makePointMap (meshMap, blockSize);

    // The point Map is contiguous and one to one, so the mesh Map
    // must also be.  Otherwise, the assignment of "work" to processes
    // wouldn't match between the mesh Map and the point Map.
    TEST_ASSERT( pointMap.isContiguous () );
    TEST_ASSERT( pointMap.isOneToOne () );

    // Check that the point Map has the right numbers of global and
    // local indices.
    TEST_EQUALITY( pointMap.getGlobalNumElements (),
                   meshMap.getGlobalNumElements () * static_cast<GST> (blockSize) );
    TEST_EQUALITY( pointMap.getNodeNumElements (),
                   meshMap.getNodeNumElements () * static_cast<GST> (blockSize) );

    // If a GID g_mesh is in the input mesh Map on a process, and if
    // the block size is b, then the GID g_point := g_mesh * b must be
    // in the output point Map on that process.  Furthermore, the
    // order must be the same as in the original Map.
    {
      ArrayView<const GO> gblMeshInds = meshMap.getNodeElementList ();
      TEST_EQUALITY( static_cast<size_t> (gblMeshInds.size ()),
                     static_cast<size_t> (meshMap.getNodeNumElements ()) );
      Array<GO> gblPointIndsIdeal (gblMeshInds.size () * static_cast<size_type> (blockSize));
      for (size_type m = 0; m < gblMeshInds.size (); ++m) {
        const size_type offset = m * static_cast<size_type> (blockSize);
        const GO pointGidStart = gblMeshInds[m] * static_cast<GO> (blockSize);
        for (LO k = 0; k < blockSize; ++k) {
          gblPointIndsIdeal[offset + static_cast<size_type> (k)] =
            pointGidStart + static_cast<GO> (k);
        }
      }

      TEST_EQUALITY( static_cast<size_t> (gblPointIndsIdeal.size ()),
                     static_cast<size_t> (pointMap.getNodeNumElements ()) );
      ArrayView<const GO> gblPointIndsActual = pointMap.getNodeElementList ();
      TEST_EQUALITY( gblPointIndsIdeal.size (), gblPointIndsActual.size () );
      if (static_cast<size_t> (gblPointIndsIdeal.size ()) ==
          static_cast<size_t> (pointMap.getNodeNumElements ())) {
        TEST_COMPARE_ARRAYS( gblPointIndsIdeal, gblPointIndsActual );
      }
    }

    // Create an overlapping mesh Map from the above nonoverlapping
    // mesh Map.  In the overlapping mesh Map, every process gets
    // exactly one overlapping GID.  In the (nonoverlapping) mesh Map
    // on the calling process, my min GID overlaps with exactly one
    // process ((myRank + 1) % numProcs) in the ovrlpng mesh Map.
    Teuchos::Array<GO> ovrlpngGIDs (numLocalMeshPoints + 1);
    for (LO lid = 0; lid < static_cast<LO> (numLocalMeshPoints); ++lid) {
      ovrlpngGIDs[lid] = meshMap.getGlobalElement (lid);
    }
    ovrlpngGIDs[numLocalMeshPoints] = ovrlpngGIDs[numLocalMeshPoints-1] %
      static_cast<GO> (meshMap.getGlobalNumElements ());
    map_type ovrlpngMeshMap (INVALID, ovrlpngGIDs (), indexBase, comm, node);

    // If the communicator has more than one process, then the
    // overlapping mesh Map must overlap, that is, not be one to one.
    TEST_ASSERT( comm->getSize () == 1 || ! ovrlpngMeshMap.isOneToOne () );

    // Create an overlapping point Map from the overlapping mesh Map.
    map_type ovrlpngPointMap = BMV::makePointMap (ovrlpngMeshMap, blockSize);

    // The overlapping point Map must not be one to one if the mesh
    // Map is not one to one.
    TEST_ASSERT( comm->getSize () == 1 || ! ovrlpngPointMap.isOneToOne () );

    // Check that the point Map has the right numbers of global and
    // local indices.
    TEST_EQUALITY( static_cast<GST> (ovrlpngPointMap.getGlobalNumElements ()),
                   static_cast<GST> (ovrlpngMeshMap.getGlobalNumElements ()) * static_cast<GST> (blockSize) );
    TEST_EQUALITY( static_cast<size_t> (ovrlpngPointMap.getNodeNumElements ()),
                   static_cast<size_t> (ovrlpngMeshMap.getNodeNumElements ()) * static_cast<size_t> (blockSize) );

    // If a GID g_mesh is in the input mesh Map on a process, and if
    // the block size is b, then the GID g_point := g_mesh * b must be
    // in the output point Map on that process.  Furthermore, the
    // order must be the same as in the original Map.
    {
      ArrayView<const GO> gblMeshInds = ovrlpngMeshMap.getNodeElementList ();
      TEST_EQUALITY( static_cast<size_t> (gblMeshInds.size ()),
                     static_cast<size_t> (ovrlpngMeshMap.getNodeNumElements ()) );
      Array<GO> gblPointIndsIdeal (gblMeshInds.size () * static_cast<size_type> (blockSize));
      for (size_type m = 0; m < gblMeshInds.size (); ++m) {
        const size_type offset = m * static_cast<size_type> (blockSize);
        const GO pointGidStart = gblMeshInds[m] * static_cast<GO> (blockSize);
        for (LO k = 0; k < blockSize; ++k) {
          gblPointIndsIdeal[offset + static_cast<size_type> (k)] =
            pointGidStart + static_cast<GO> (k);
        }
      }

      TEST_EQUALITY( static_cast<size_t> (gblPointIndsIdeal.size ()),
                     static_cast<size_t> (ovrlpngPointMap.getNodeNumElements ()) );
      ArrayView<const GO> gblPointIndsActual = ovrlpngPointMap.getNodeElementList ();
      TEST_EQUALITY( gblPointIndsIdeal.size (), gblPointIndsActual.size () );
      if (static_cast<size_t> (gblPointIndsIdeal.size ()) ==
          static_cast<size_t> (ovrlpngPointMap.getNodeNumElements ())) {
        TEST_COMPARE_ARRAYS( gblPointIndsIdeal, gblPointIndsActual );
      }
    }
  }

//
// INSTANTIATIONS
//

#define UNIT_TEST_GROUP( SCALAR, LO, GO, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( ExpBlockMap, HilbertsHotel, SCALAR, LO, GO, NODE )

  TPETRA_ETI_MANGLING_TYPEDEFS()

  TPETRA_INSTANTIATE_TESTMV_NOGPU( UNIT_TEST_GROUP )

} // namespace (anonymous)


