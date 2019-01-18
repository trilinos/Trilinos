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

#include "Teuchos_UnitTestHarness.hpp"
#include "Tpetra_TestingUtilities.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_Import.hpp"
#include "Tpetra_Vector.hpp"

#include "Tpetra_Distributor.hpp"
#include "Tpetra_BlockCrsMatrix.hpp"
#include "Tpetra_BlockCrsMatrix_Helpers.hpp"
#include "Tpetra_CrsGraph.hpp"
#include "Tpetra_Import_Util.hpp"

namespace {
  using Teuchos::Comm;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::outArg;
  using std::endl;

  // bool testMpi = true;

  TEUCHOS_STATIC_SETUP()
  {
    Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
    clp.addOutputSetupOptions(true);
  }

  //
  // UNIT TESTS
  //
  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(ImportExport,ImportConstructExpert,LO,GO,NT) {
    //
    using Teuchos::RCP;
    using Teuchos::rcp;
    using std::endl;
    typedef Teuchos::Array<int>::size_type size_type;
    typedef Tpetra::BlockCrsMatrix<>::scalar_type Scalar;
    typedef Tpetra::BlockCrsMatrix<Scalar,LO,GO,NT> matrix_type;
    typedef Tpetra::Map<LO,GO,NT> map_type;
    typedef Tpetra::CrsGraph<LO,GO,NT>  graph_type;
    typedef Tpetra::global_size_t GST;
    typedef typename matrix_type::device_type device_type;
    typedef typename Kokkos::View<Scalar**, Kokkos::LayoutRight, device_type>::HostMirror block_type;

    Teuchos::OSTab tab0 (out);
    Teuchos::OSTab tab1 (out);

    RCP<const Comm<int> > comm = Tpetra::getDefaultComm();

    //    int rank = comm->getRank();

    const LO NumRows = 32;
    const GST gblNumRows = static_cast<GST> ( NumRows * comm->getSize ());
    const GO indexBase = 0;
    const size_t numEntPerRow = 11;

    RCP<const map_type> rowMap =
      rcp (new map_type (gblNumRows, static_cast<size_t> (NumRows),
                         indexBase, comm));
    const GO gblNumCols = static_cast<GO> (rowMap->getGlobalNumElements ());

    RCP<graph_type> G =
      rcp (new graph_type (rowMap, numEntPerRow,
                           Tpetra::StaticProfile));

    Teuchos::Array<GO> gblColInds (numEntPerRow);
    for (LO lclRow = 0; lclRow < NumRows; ++lclRow) {
      const GO gblInd = rowMap->getGlobalElement (lclRow);
      // Put some entries in the graph.
      for (LO k = 0; k < static_cast<LO> (numEntPerRow); ++k) {
        const GO curColInd = (gblInd + static_cast<GO> (3*k)) % gblNumCols;
        gblColInds[k] = curColInd;
      }
      G->insertGlobalIndices (gblInd, gblColInds ());
    }
    // Make the graph ready for use by BlockCrsMatrix.
    G->fillComplete ();
    const auto& meshRowMap = * (G->getRowMap ());
    // Contrary to expectations, asking for the graph's number of
    // columns, or asking the column Map for the number of entries,
    // won't give the correct number of columns in the graph.
    // const GO gblNumCols = graph->getDomainMap ()->getGlobalNumElements ();
    const LO lclNumRows = meshRowMap.getNodeNumElements ();
    const LO blkSize = 16;

    RCP<matrix_type> A = rcp (new matrix_type (*G, blkSize));
    // Create a values to use when filling the sparse matrix. Use primes because it's cute.
    Scalar prime[]={
      2,    3,   5,   7,  11,   13,  17,  19,  23,  29,
      31,  37,  41,  43,  47,   53,  59,  61,  67,  71,
      73,  79,  83,  89,  97,  101, 103, 107, 109, 113,
      127, 131, 137, 139, 149, 151, 157, 163, 167, 173,
      179, 181, 191, 193, 197, 199, 211, 223, 227, 229,
      233, 239, 241, 251, 257, 263, 269, 271};

    int idx=0;
    block_type curBlk ("curBlk", blkSize, blkSize);

    for (LO j = 0; j < blkSize; ++j) {
      for (LO i = 0; i < blkSize; ++i) {
        curBlk(i,j) = 1.0* prime[idx++];;
        if (idx >=58) idx = 0;
      }
    }

    // Fill in the block sparse matrix.
    for (LO lclRow = 0; lclRow < lclNumRows; ++lclRow) { // for each of my rows
      Teuchos::ArrayView<const LO> lclColInds;
      G->getLocalRowView (lclRow, lclColInds);

      // Put some entries in the matrix.
      for (LO k = 0; k < static_cast<LO> (lclColInds.size ()); ++k) {
        const LO lclColInd = lclColInds[k];
        const LO err =
          A->replaceLocalValues (lclRow, &lclColInd, curBlk.data (), 1);
        TEUCHOS_TEST_FOR_EXCEPTION(err != 1, std::logic_error, "Bug");
      }
    }


    auto olddist = G->getImporter()->getDistributor();

    Teuchos::RCP<const map_type> source = G->getImporter()->getSourceMap ();
    Teuchos::RCP<const map_type> target = G->getImporter()->getTargetMap ();

    // Still need to get remote GID's, and export LID's.

    Teuchos::Array<LO> saveremoteLIDs = G->getImporter()->getRemoteLIDs();

    Teuchos::Array<LO> remoteLIDs = G->getImporter()->getRemoteLIDs();

    Teuchos::Array<GO> remoteGIDs(remoteLIDs.size());
    Teuchos::Array<int> userRemotePIDs(remoteLIDs.size());
    for(size_t i=0; i< (size_t)remoteLIDs.size(); i++)
      remoteGIDs[i] = target->getGlobalElement(G->getImporter()->getRemoteLIDs()[i]);
    Tpetra::Import_Util::getRemotePIDs(*G->getImporter(),userRemotePIDs);

    const Teuchos::ArrayView<const LO> exportLIDs      = G->getImporter()->getExportLIDs();
    const Teuchos::ArrayView<const int> userExportPIDs = G->getImporter()->getExportPIDs();

    Tpetra::Import<LO,GO,NT> newimport(source,
                                       target,
                                       userRemotePIDs,
                                       exportLIDs,
                                       userExportPIDs ,
                                       Teuchos::null,
                                       Teuchos::null ); // plist == null

    Teuchos::RCP<const map_type> newsource = newimport.getSourceMap ();
    Teuchos::RCP<const map_type> newtarget = newimport.getTargetMap ();


    Teuchos::Array<LO> newexportLIDs = newimport.getExportLIDs();

    const int rank = newsource->getComm()->getRank();
    if(newexportLIDs.size()!=exportLIDs.size())
      {

        out <<"Rank "<<rank<<" newexportLIDs.size does not match exportLIDs.size()"<<endl;
        out <<"Rank "<<rank<<" oldExportLIDs "<<exportLIDs<<endl;
        out <<"Rank "<<rank<<" newExportLIDs "<<newexportLIDs<<endl;
        success = false;
      }
    else
      for(size_type i=0;i<exportLIDs.size();++i)
        if(exportLIDs[i]!=newexportLIDs[i]) {
          out <<"Rank "<<rank<<" exportLIDs["<<i<<"] ="<<exportLIDs[i]<<" != newexportLIDs[i] = "<<newexportLIDs[i]<<endl;
          success = false;
          break;
        }

    Teuchos::Array<LO> newremoteLIDs = newimport.getRemoteLIDs();
    if(newremoteLIDs.size()!=saveremoteLIDs.size())
      {
        out <<"Rank "<<rank<<" newremoteLIDs.size does not match remoteLIDs.size()"<<endl;
        out <<"Rank "<<rank<<" oldRemoteLIDs "<<saveremoteLIDs<<endl;
        out <<"Rank "<<rank<<" newRemoteLIDs "<<newremoteLIDs<<endl;
        success = false;
      }
    else
      for(size_type i=0;i<saveremoteLIDs.size();++i)
        if(saveremoteLIDs[i]!=newremoteLIDs[i]) {
          out <<"Rank "<<rank<<" remoteLIDs["<<i<<"] ="<<remoteLIDs[i]<<" != newremoteLIDs[i] = "<<newremoteLIDs[i]<<endl;
          success = false;
          break;
        }

    int globalSuccess_int = -1;
    Teuchos::reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, outArg(globalSuccess_int) );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );

    out<<" GlobalSuccess = "<<globalSuccess_int<<std::endl;

  }

  //
  // INSTANTIATIONS
  //

#define UNIT_TEST_3( LO, GO, NT ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( ImportExport,ImportConstructExpert,LO,GO,NT)

  TPETRA_ETI_MANGLING_TYPEDEFS()

  TPETRA_INSTANTIATE_LGN( UNIT_TEST_3 )

} // namespace (anonymous)


