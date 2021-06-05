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
// ************************************************************************
// @HEADER
*/

#include "Tpetra_TestingUtilities.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_Import.hpp"
#include "Tpetra_Map.hpp"
#include "Kokkos_Core.hpp"

namespace { // (anonymous)

  using Tpetra::TestingUtilities::getDefaultComm;
  using Teuchos::ArrayView;
  using Teuchos::Comm;
  using Teuchos::outArg;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::REDUCE_MIN;
  using Teuchos::reduceAll;
  using Teuchos::tuple;
  using std::endl;
  using GST = Tpetra::global_size_t;

  // Both source and target matrices have one row on each process.
  //
  // Target matrix global column indices:
  // Proc 0: Global row index 0: [0, 1, 2, 3, 4, 5]
  // Proc 1: Global row index 1: [0, 1, 2, 3, 4, 5]
  //
  // Source matrix global column indices:
  // Proc 0: Global row index 1: []
  // Proc 1: Global row index 0: [3, 4, 5, 6, 7, 8, 9]
  //
  // After Import, target should look like this:
  // Proc 0: Global row index 0: [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
  // Proc 1: Global row index 1: [0, 1, 2, 3, 4, 5]

  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( CrsMatrix, UnpackMerge1, Scalar, Node )
  {
    using LO = Tpetra::Map<>::local_ordinal_type;
    using GO = Tpetra::Map<>::global_ordinal_type;
    using crs_matrix_type = Tpetra::CrsMatrix<Scalar, LO, GO, Node>;
    using import_type = Tpetra::Import<LO, GO, Node>;
    using map_type = Tpetra::Map<LO, GO, Node>;

    RCP<const Comm<int> > comm = getDefaultComm();
    const int myRank = comm->getRank();
    const int numProcs = comm->getSize();

    out << "Test that Tpetra::CrsMatrix::unpackAndCombine into a "
      "target matrix with a non-static graph merges column indices"
      << endl;
    Teuchos::OSTab tab1(out);

    TEST_ASSERT( numProcs == 2 );
    if (numProcs != 2) {
      out << "This test requires exactly two MPI processes, but you "
        "ran it with " << numProcs << " process(es)." << endl;
      return;
    }

    const GO gblNumRows (2);
    const GO indexBase (0);
    std::vector<GO> srcRowMapInds;
    std::vector<GO> tgtRowMapInds;
    if (myRank == 0) {
      srcRowMapInds.push_back(1);
      tgtRowMapInds.push_back(0);
    }
    else if (myRank == 1) {
      srcRowMapInds.push_back(0);
      tgtRowMapInds.push_back(1);
    }
    const LO srcLclNumRows (srcRowMapInds.size());
    const LO tgtLclNumRows (tgtRowMapInds.size());

    RCP<const map_type> srcRowMap =
      rcp(new map_type(static_cast<GST>(gblNumRows),
                       srcRowMapInds.data(), srcLclNumRows,
                       indexBase, comm));
    RCP<const map_type> tgtRowMap =
      rcp(new map_type(static_cast<GST>(gblNumRows),
                       tgtRowMapInds.data(), tgtLclNumRows,
                       indexBase, comm));

    const GO gblNumCols = 10;
    RCP<const map_type> colMap =
      rcp(new map_type(static_cast<GST>(gblNumCols),
                       indexBase, comm,
                       Tpetra::LocallyReplicated));
    RCP<const map_type> domMap =
      rcp(new map_type(static_cast<GST>(gblNumCols),
                       indexBase, comm,
                       Tpetra::GloballyDistributed));
    RCP<const map_type> ranMap = srcRowMap;
    import_type importer(srcRowMap, tgtRowMap);

    std::vector<GO> srcGblColInds;
    if (myRank == 1) {
      srcGblColInds = std::vector<GO>{{3, 4, 5, 6, 7, 8, 9}};
    }
    std::vector<GO> tgtGblColInds{{0, 1, 2, 3, 4, 5}};
    std::vector<Scalar> srcVals(srcGblColInds.size(), Scalar(1.0));
    std::vector<Scalar> tgtVals(tgtGblColInds.size(), Scalar(1.0));

    std::vector<Scalar> expectedTgtVals;
    if (myRank == 0) {
      expectedTgtVals.resize(10);
      for(LO k = 0; k < LO(3); ++k) {
        expectedTgtVals[k] = 1.0;
      }
      for(LO k = LO(3); k < LO(6); ++k) {
        expectedTgtVals[k] = 2.0;
      }
      for(LO k = LO(6); k < LO(10); ++k) {
        expectedTgtVals[k] = 1.0;
      }
    }
    else if (myRank == 1) {
      expectedTgtVals.resize(6);
      for(LO k = 0; k < LO(6); ++k) {
        expectedTgtVals[k] = 1.0;
      }
    }

    for (const bool A_src_is_fill_complete : {false, true}) {
      out << "A_src will" << (A_src_is_fill_complete ? "" : " NOT")
          << " be fill complete." << endl;
      crs_matrix_type A_src(srcRowMap, colMap, srcGblColInds.size());
      crs_matrix_type A_tgt(tgtRowMap, colMap, tgtGblColInds.size());

      for (LO lclRow = 0; lclRow < srcLclNumRows; ++lclRow) {
        const GO gblRow = srcRowMap->getGlobalElement(lclRow);
        A_tgt.insertGlobalValues(gblRow,
          Teuchos::ArrayView<const GO>(tgtGblColInds),
          Teuchos::ArrayView<const Scalar>(tgtVals));
        A_src.insertGlobalValues(gblRow,
          Teuchos::ArrayView<const GO>(srcGblColInds),
          Teuchos::ArrayView<const Scalar>(srcVals));
      }
      if (A_src_is_fill_complete) {
        A_src.fillComplete(domMap, ranMap);
      }

      out << "Finished A_src.fillComplete(domMap, ranMap)" << endl;

      TEST_ASSERT( ! A_tgt.isStaticGraph() );

      A_tgt.doImport(A_src, importer, Tpetra::INSERT);
      A_tgt.fillComplete(domMap, ranMap);

      Kokkos::fence(); // since we're accessing data on host now

      typename crs_matrix_type::local_inds_host_view_type lclColInds;
      typename crs_matrix_type::values_host_view_type vals;
      const LO lclRowToTest (0);
      A_tgt.getLocalRowView(lclRowToTest, lclColInds, vals);

      const LO expectedNumEnt = myRank == 0 ? LO(10) : LO(6);
      TEST_EQUALITY( LO(lclColInds.size()), expectedNumEnt );
      TEST_EQUALITY( LO(vals.size()), expectedNumEnt );

      if (success && myRank == 0) {
        for (LO k = 0; k < expectedNumEnt; ++k) {
          TEST_EQUALITY( lclColInds[k], LO(k) );
          const Scalar expectedVal = expectedTgtVals[k];
          TEST_EQUALITY( vals[k], expectedVal );
        }
      }

      // Test whether all processes passed the test.
      int lclSuccess = success ? 1 : 0;
      int gblSuccess = 0;
      reduceAll(*comm, REDUCE_MIN, lclSuccess, outArg(gblSuccess));
      TEST_EQUALITY_CONST( gblSuccess, 1 );
    }
  }

  //TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( CrsMatrix, UnpackMerge2, Scalar, Node )
  TEUCHOS_UNIT_TEST( CrsMatrix, UnpackMerge2 )
  {
    using Scalar = Tpetra::CrsMatrix<>::scalar_type;
    using LO = Tpetra::Map<>::local_ordinal_type;
    using GO = Tpetra::Map<>::global_ordinal_type;
    using Node = Tpetra::Map<>::node_type;
    using crs_matrix_type = Tpetra::CrsMatrix<Scalar, LO, GO, Node>;
    using import_type = Tpetra::Import<LO, GO, Node>;
    using map_type = Tpetra::Map<LO, GO, Node>;
    using gids_type = typename crs_matrix_type::nonconst_global_inds_host_view_type;
    using vals_type = typename crs_matrix_type::nonconst_values_host_view_type;
    int lclSuccess = 1;
    int gblSuccess = 0;

    RCP<const Comm<int> > comm = getDefaultComm();
    const int myRank = comm->getRank();
    const int numProcs = comm->getSize();

    out << "Regression test with a real-life example" << endl;
    Teuchos::OSTab tab1(out);

    TEST_ASSERT( numProcs == 2 );
    if (numProcs != 2) {
      out << "This test requires exactly two MPI processes, but you "
        "ran it with " << numProcs << " process(es)." << endl;
      return;
    }

    const GO gblNumRows (2);
    const GO indexBase (0);
    std::vector<GO> srcRowMapInds;
    std::vector<GO> tgtRowMapInds;
    if (myRank == 0) {
      srcRowMapInds.push_back(1);
      tgtRowMapInds.push_back(0);
    }
    else if (myRank == 1) {
      srcRowMapInds.push_back(0);
      tgtRowMapInds.push_back(1);
    }
    const LO srcLclNumRows (srcRowMapInds.size());
    const LO tgtLclNumRows (tgtRowMapInds.size());

    RCP<const map_type> srcRowMap =
      rcp(new map_type(static_cast<GST>(gblNumRows),
                       srcRowMapInds.data(), srcLclNumRows,
                       indexBase, comm));
    RCP<const map_type> tgtRowMap =
      rcp(new map_type(static_cast<GST>(gblNumRows),
                       tgtRowMapInds.data(), tgtLclNumRows,
                       indexBase, comm));
    // [0, ... 199,999]
    const GO gblNumCols = 200000;
    RCP<const map_type> colMap =
      rcp(new map_type(static_cast<GST>(gblNumCols),
                       indexBase, comm,
                       Tpetra::LocallyReplicated));
    RCP<const map_type> domMap =
      rcp(new map_type(static_cast<GST>(gblNumCols),
                       indexBase, comm,
                       Tpetra::GloballyDistributed));
    RCP<const map_type> ranMap = srcRowMap;
    import_type importer(srcRowMap, tgtRowMap);

    // Input to insert: 96 entries, sent from Proc 1.
    std::vector<GO> srcGblColInds {{
      142944, 142945, 142946, 142947, 142948, 142949, 142950, 142951,
      142952, 142953, 142954, 142955, 142959, 142960, 142961, 142965,
      142966, 142967, 142968, 142969, 142970, 143142, 143143, 143144,
      198279, 198280, 198281, 198282, 198283, 198284, 198291, 198292,
      198293, 198303, 198304, 198305, 198309, 198310, 198311, 198333,
      198334, 198335, 198336, 198337, 198338, 198339, 198340, 198341,
      198342, 198343, 198344, 198345, 198346, 198347, 198348, 198349,
      198350, 198351, 198352, 198353, 198354, 198355, 198356, 198699,
      198700, 198701, 198702, 198703, 198704, 198705, 198706, 198707,
      198708, 198709, 198710, 198711, 198712, 198713, 198729, 198730,
      198731, 198732, 198733, 198734, 198735, 198736, 198737, 198738,
      198739, 198740, 198741, 198742, 198743, 198744, 198745, 198746
    }};

    // Current contents of Proc 0 row: 96 entries.
    std::vector<GO> tgtGblColInds {{
      166215, 166216, 166217, 166218, 166219, 166220, 166221, 166222,
      166223, 166224, 166225, 166226, 166227, 166228, 166229, 166230,
      166231, 166232, 166233, 166234, 166235, 166236, 166237, 166238,
      166239, 166240, 166241, 166242, 166243, 166244, 166245, 166246,
      166247, 198279, 198280, 198281, 198282, 198283, 198284, 198285,
      198286, 198287, 198288, 198289, 198290, 198291, 198292, 198293,
      198294, 198295, 198296, 198297, 198298, 198299, 198300, 198301,
      198302, 198303, 198304, 198305, 198306, 198307, 198308, 198309,
      198310, 198311, 198312, 198313, 198314, 198315, 198316, 198317,
      198333, 198334, 198335, 198336, 198337, 198338, 198339, 198340,
      198341, 198342, 198343, 198344, 198345, 198346, 198347, 198348,
      198349, 198350, 198351, 198352, 198353, 198354, 198355, 198356
    }};

    TEST_EQUALITY( srcGblColInds.size(), size_t(96) );
    TEST_EQUALITY( tgtGblColInds.size(), size_t(96) );

    std::vector<GO> srcCpy (srcGblColInds);
    auto srcBeg = srcCpy.begin();
    auto srcEnd = srcCpy.end();
    std::sort(srcBeg, srcEnd);
    srcEnd = std::unique(srcBeg, srcEnd);

    std::vector<GO> tgtCpy (tgtGblColInds);
    auto tgtBeg = tgtCpy.begin();
    auto tgtEnd = tgtCpy.end();
    std::sort(tgtBeg, tgtEnd);
    tgtEnd = std::unique(tgtBeg, tgtEnd);

    std::vector<GO> unionGblColInds(srcGblColInds.size() +
                                    tgtGblColInds.size());
    auto unionEnd = std::set_union(srcBeg, srcEnd, tgtBeg, tgtEnd,
                                    unionGblColInds.begin());
    unionGblColInds.resize(unionEnd - unionGblColInds.begin());
    const size_t unionSize = unionGblColInds.size();

    out << "Number of elements in set union of column indices: "
        << unionSize << endl;

    std::vector<Scalar> srcVals(srcGblColInds.size(), Scalar(1.0));
    std::vector<Scalar> tgtVals(tgtGblColInds.size(), Scalar(1.0));

    std::vector<Scalar> expectedTgtVals;
    if (myRank == 0) {
      std::vector<GO> intersectionGblColInds
        (std::min(srcGblColInds.size(), tgtGblColInds.size()));
      auto intersectionEnd =
        std::set_intersection(srcBeg, srcEnd, tgtBeg, tgtEnd,
                              intersectionGblColInds.begin());
      intersectionGblColInds.resize
        (intersectionEnd - intersectionGblColInds.begin());

      expectedTgtVals.resize(unionSize);
      for(size_t k = 0; k < unionSize; ++k) {
        const GO resultGid = unionGblColInds[k];
        auto it = std::lower_bound(intersectionGblColInds.begin(),
                                   intersectionGblColInds.end(),
                                   resultGid);
        if (it != intersectionGblColInds.end() &&
            *it == resultGid) {
          expectedTgtVals[k] = 2.0;
        }
        else {
          expectedTgtVals[k] = 1.0;
        }
      }
    }

    crs_matrix_type A_src(srcRowMap, colMap, srcGblColInds.size());
    crs_matrix_type A_tgt(tgtRowMap, colMap, tgtGblColInds.size());

    for (LO lclRow = 0; lclRow < srcLclNumRows; ++lclRow) {
      const GO gblRow = srcRowMap->getGlobalElement(lclRow);
      A_tgt.insertGlobalValues(gblRow,
        Teuchos::ArrayView<const GO>(tgtGblColInds),
        Teuchos::ArrayView<const Scalar>(tgtVals));
      A_src.insertGlobalValues(gblRow,
        Teuchos::ArrayView<const GO>(srcGblColInds),
        Teuchos::ArrayView<const Scalar>(srcVals));
    }
    A_src.fillComplete(domMap, ranMap);

    TEST_ASSERT( ! A_tgt.isStaticGraph() );

    A_tgt.doImport(A_src, importer, Tpetra::INSERT);
    A_tgt.fillComplete(domMap, ranMap);

    Kokkos::fence(); // since we're accessing data on host now

    if (myRank == 0) {
      const GO gblRowToTest = tgtRowMap->getMinGlobalIndex();
      size_t numEnt = A_tgt.getNumEntriesInGlobalRow(gblRowToTest);
      gids_type gblColInds("gids",numEnt);
      vals_type vals("vals",numEnt);
      A_tgt.getGlobalRowCopy(gblRowToTest, gblColInds,
                             vals, numEnt);

      const LO expectedNumEnt(expectedTgtVals.size());
      TEST_EQUALITY( size_t(numEnt), size_t(expectedNumEnt) );
      TEST_EQUALITY( size_t(gblColInds.extent(0)),
                     size_t(expectedNumEnt) );
      TEST_EQUALITY( size_t(vals.extent(0)), size_t(expectedNumEnt) );

      if (success) {
        for (LO k = 0; k < expectedNumEnt; ++k) {
          TEST_EQUALITY( gblColInds[k], unionGblColInds[k] );
          const Scalar expectedVal = expectedTgtVals[k];
          TEST_EQUALITY( vals[k], expectedVal );
        }
      }
    }

    lclSuccess = success ? 1 : 0;
    gblSuccess = 0;
    reduceAll(*comm, REDUCE_MIN, lclSuccess, outArg(gblSuccess));
    TEST_EQUALITY_CONST( gblSuccess, 1 );
  }

//
// INSTANTIATIONS
//

#define UNIT_TEST_GROUP( SCALAR, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( CrsMatrix, UnpackMerge1, SCALAR, NODE )
  //TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( CrsMatrix, UnpackMerge2, SCALAR, NODE )

  TPETRA_ETI_MANGLING_TYPEDEFS()

  TPETRA_INSTANTIATE_SN( UNIT_TEST_GROUP )

} // namespace (anonymous)
