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

  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( CrsMatrix, UnpackMerge, Scalar, Node )
  {
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
    using LO = Tpetra::Map<>::local_ordinal_type;
    using GO = Tpetra::Map<>::global_ordinal_type;
    using crs_matrix_type = Tpetra::CrsMatrix<Scalar, LO, GO, Node>;
    using import_type = Tpetra::Import<LO, GO, Node>;
    using map_type = Tpetra::Map<LO, GO, Node>;
    using GST = Tpetra::global_size_t;
    using STS = Teuchos::ScalarTraits<Scalar>;

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
      srcRowMapInds = std::vector<GO>{{1}};
      tgtRowMapInds = std::vector<GO>{{0}};
    }
    else if (myRank == 1) {
      srcRowMapInds = std::vector<GO>{{0}};
      tgtRowMapInds = std::vector<GO>{{1}};
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

      Teuchos::ArrayView<const LO> lclColInds;
      Teuchos::ArrayView<const Scalar> vals;
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

//
// INSTANTIATIONS
//

#define UNIT_TEST_GROUP( SCALAR, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( CrsMatrix, UnpackMerge, SCALAR, NODE )

  TPETRA_ETI_MANGLING_TYPEDEFS()

  TPETRA_INSTANTIATE_SN( UNIT_TEST_GROUP )

} // namespace (anonymous)


