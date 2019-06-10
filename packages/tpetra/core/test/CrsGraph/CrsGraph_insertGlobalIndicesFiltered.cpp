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

#include "Tpetra_CrsGraph.hpp"
#include "Tpetra_Export.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_TestingUtilities.hpp"
#include <type_traits>

namespace { // (anonymous)
  using Teuchos::outArg;
  using Teuchos::RCP;
  using Teuchos::REDUCE_MIN;
  using Teuchos::reduceAll;
  using std::endl;
  typedef Tpetra::global_size_t GST;

  //
  // UNIT TESTS
  //

  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( CrsGraph, insertGlobalIndicesFiltered, LO, GO, NODE_TYPE )
  {
    static_assert (std::is_integral<LO>::value,
                   "LO must be a built-in integer type.");
    static_assert (std::is_integral<GO>::value,
                   "GO must be a built-in integer type.");
    // It's harder to test NODE_TYPE; without a "this is a valid node
    // type" compile-time query function, that would call for "is
    // Kokkos::Compat::KokkosDeviceWrapperNode<E, M> for some E which
    // is a valid Kokkos execution space, and some M which is a valid
    // Kokkos memory space."

    typedef Tpetra::CrsGraph<LO, GO, NODE_TYPE> crs_graph_type;
    typedef Tpetra::Map<LO, GO, NODE_TYPE> map_type;
    typedef Tpetra::Export<LO, GO, NODE_TYPE> export_type;
    int lclSuccess = 1; // to set below
    int gblSuccess = 0; // to set below
    const GST INVALID = Teuchos::OrdinalTraits<GST>::invalid ();

    auto comm = Tpetra::TestingUtilities::getDefaultComm ();
    const int numProcs = comm->getSize ();
    TEST_EQUALITY( numProcs, 2 );
    if (numProcs != 2) {
      out << "This test requires exactly two MPI processes!" << endl;
    }
    const int myRank = comm->getRank ();

    const GO indexBase = 0;
    const GO gblNumRows = 2;
    constexpr LO lclNumRows = 1;
    const GO gblRows_proc0_src[lclNumRows] = {0};
    const GO gblRows_proc1_src[lclNumRows] = {1};
    const GO gblRows_proc0_tgt[lclNumRows] = {1};
    const GO gblRows_proc1_tgt[lclNumRows] = {0};

    // Source and target graphs need different row Maps, so that
    // communication is nontrivial.  We do this by giving Process
    // (0,1) row (1,0) -- that is, by flipping row ownership.
    RCP<const map_type> rowMap_src (new map_type (gblNumRows,
                                                  (myRank == 0) ? gblRows_proc0_src : gblRows_proc1_src,
                                                  lclNumRows, indexBase, comm));
    RCP<const map_type> rowMap_tgt (new map_type (gblNumRows,
                                                  (myRank == 0) ? gblRows_proc0_tgt : gblRows_proc1_tgt,
                                                  lclNumRows, indexBase, comm));
    export_type theExport (rowMap_src, rowMap_tgt);

    // Construct the target column Map so that the input sequence of
    // column indices would have noncontiguous chunks filtered out.
    // It's OK if all local rows do the same thing.
    constexpr LO lclNumColMapInds_tgt = 6;
    const GO gblColMapInds_tgt[lclNumColMapInds_tgt] =
      {1, 2, 4, 5, 6, 8};
    RCP<const map_type> colMap_tgt (new map_type (INVALID, gblColMapInds_tgt,
                                                  lclNumColMapInds_tgt,
                                                  indexBase, comm));
    constexpr LO lclNumColMapInds_src = 10;
    const GO gblColMapInds_src[lclNumColMapInds_src] =
      {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
    RCP<const map_type> colMap_src (new map_type (INVALID, gblColMapInds_src,
                                                  lclNumColMapInds_src,
                                                  indexBase, comm));

    // The graph in either case is the graph of a 2 x 10 matrix.
    RCP<const map_type> ranMap (new map_type (gblNumRows, lclNumRows, indexBase, comm));
    const GO gblNumCols = 10;
    const GO lclNumCols = 5; // works since we know numProcs == 2
    RCP<const map_type> domMap (new map_type (gblNumCols, lclNumCols, indexBase, comm));

    // This must be enough for any row in the source or target matrix.
    const size_t maxNumEntPerRow = static_cast<size_t> (lclNumColMapInds_src);

    // Buffer for storing output of getGlobalRowCopy.
    Teuchos::Array<GO> gblColIndsBuf (maxNumEntPerRow);

#ifdef TPETRA_ENABLE_DEPRECATED_CODE
    const Tpetra::ProfileType profileTypes[2] =
      {Tpetra::DynamicProfile, Tpetra::StaticProfile};
#else
    const Tpetra::ProfileType profileTypes[1] = {Tpetra::StaticProfile};
#endif // TPETRA_ENABLE_DEPRECATED_CODE
    for (auto profileType_src : profileTypes) {
      crs_graph_type graph_src (rowMap_src, maxNumEntPerRow, profileType_src);
      for (LO lclRow = 0; lclRow < lclNumRows; ++lclRow) {
        const GO gblRow = rowMap_src->getGlobalElement (lclRow);
        // Every row of the source graph gets the same column indices.
        graph_src.insertGlobalIndices (gblRow, lclNumColMapInds_src, gblColMapInds_src);
      }
      graph_src.fillComplete (domMap, ranMap);

      for (auto profileType_tgt : profileTypes) {
        // Filtering only happens if the target graph has a column
        // Map, so we must give it a column Map.
        crs_graph_type graph_tgt (rowMap_tgt, colMap_tgt, maxNumEntPerRow, profileType_tgt);

        // mfh 20 Jul 2017: If we clone this test in order to test
        // CrsMatrix column index filtering, then we should include
        // both the separate doExport,fillComplete case, and the
        // exportAndFillComplete combined case.
        graph_tgt.doExport (graph_src, theExport, Tpetra::INSERT);
        graph_tgt.fillComplete (domMap, ranMap);

        for (LO lclRow = 0; lclRow < lclNumRows; ++lclRow) {
          const GO gblRow = rowMap_tgt->getGlobalElement (lclRow);
          size_t numEnt = 0; // output argument
          graph_tgt.getGlobalRowCopy (gblRow, gblColIndsBuf (), numEnt);
          TEST_EQUALITY( numEnt, static_cast<size_t> (lclNumColMapInds_tgt) );
          if (numEnt == static_cast<size_t> (lclNumColMapInds_tgt)) {
            for (LO k = 0; k < static_cast<LO> (numEnt); ++k) {
              TEST_EQUALITY( gblColIndsBuf[k], gblColMapInds_tgt[k] );
            }
          }
        }
      } // for each target graph profile type
    } // for each source graph profile type

    // Make sure that the test succeeded on all processes.
    lclSuccess = success ? 1 : 0;
    gblSuccess = 0; // output argument
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY_CONST( gblSuccess, 1 );
  }

//
// INSTANTIATIONS
//

// Tests to build and run.  We will instantiate them over all enabled
// local ordinal (LO), global ordinal (GO), and Kokkos Node (NT) types.
#define UNIT_TEST_GROUP( LO, GO, NT ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( CrsGraph, insertGlobalIndicesFiltered, LO, GO, NT )

  TPETRA_ETI_MANGLING_TYPEDEFS()

  TPETRA_INSTANTIATE_LGN( UNIT_TEST_GROUP )

} // namespace (anonymous)
