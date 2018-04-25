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

#include <Tpetra_ConfigDefs.hpp>
#include <Tpetra_TestingUtilities.hpp>
#include <Teuchos_UnitTestHarness.hpp>

#include <map>
#include <Teuchos_OrdinalTraits.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_VerboseObject.hpp>
#include <Teuchos_as.hpp>
#include <Teuchos_Tuple.hpp>
#include "Tpetra_CrsGraph.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_DefaultPlatform.hpp"
#include "Tpetra_Distributor.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_Util.hpp"
#include "Tpetra_Import.hpp"
#include "Tpetra_Export.hpp"
#include "Tpetra_FEMultiVector.hpp"
#include "Tpetra_Details_gathervPrint.hpp"

namespace {

  // UNIT TESTS
  //

  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( FEMultiVector, doImport, LO, GO, Scalar )
  {
    using Teuchos::VERB_EXTREME;
    using Teuchos::VERB_NONE;
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::Array;
    using Teuchos::Comm;
    typedef Tpetra::global_size_t GST;

    const GST INVALID = Teuchos::OrdinalTraits<GST>::invalid ();
    const RCP<const Comm<int> > comm = Tpetra::getDefaultComm();
    //    const int numProcs = comm->getSize();


    // Prepare for verbose output, if applicable.
    //    Teuchos::EVerbosityLevel verbLevel = verbose ? VERB_EXTREME : VERB_NONE;
    Teuchos::EVerbosityLevel verbLevel = VERB_EXTREME;
    const bool doPrint = includesVerbLevel (verbLevel, VERB_EXTREME, true);
    if (doPrint) {
      out << "FEMultiVector unit test" << endl;
    }
    OSTab tab1 (out); // Add one tab level

    std::ostringstream err;
    int lclErr = 0;

    try {
      OSTab tab2 (out);
      const int num_local_elements = 3;

      // create Map
      RCP<const Map<LO, GO> > map =
        createContigMap<LO, GO> (INVALID, num_local_elements, comm);

      // create CrsGraph object
      RCP<CrsGraph<LO, GO> > graph =
        rcp (new CrsGraph<LO, GO> (map, 3, DynamicProfile));

      // Create a simple tridiagonal source graph.
      if (doPrint) {
        out << "Filling source CrsGraph" << endl;
      }
      Array<GO> entry(1);
      LO row = 0;
      for (size_t i = 0; i < map->getNodeNumElements (); ++i, ++row) {
        const GO globalrow = map->getGlobalElement (row);
        entry[0] = globalrow;
        graph->insertGlobalIndices (globalrow, entry());
        if (globalrow != map->getMinGlobalIndex()) {
          entry[0] = globalrow-1;
          graph->insertGlobalIndices (globalrow, entry());
        }
        if (globalrow != map->getMaxGlobalIndex()) {
          entry[0] = globalrow+1;
          graph->insertGlobalIndices (globalrow+!, entry());
        }       
      }
      graph->fillComplete();


      Scalar ZERO = Teuchos::ScalarTraits<Scalar>::zero();
      RCP<const Map<LO,GO> > domainMap = graph->getDomainMap();
      RCP<const Import<LO,GO> > importer = graph->getImporter();
      MultiVector<Scalar,LO,GO> Vdomain(domainMap,1), Vcolumn(graph->getColMap(),1);
      FEMultiVector<Scalar,LO,GO> Vfe(domainMap,importer,1);

      // 1) Test domain -> column
      size_t Ndomain = domainMap->getNodeNumElements();

      for(size_t i=0; i<Ndomain; i++)
        Vdomain->getDataNonConst()[i] = domainMap->getGlobalElement(i);
      Vcolumn->PutScalar(zero);
      Vcolumn->doImport(*Vdomain,i*mporter,Tpetra::ADD);



    } catch (std::exception& e) {
      err << "Proc " << myRank << ": " << e.what () << endl;
      lclErr = 1;
    }

    int gblErr = 0;
    reduceAll<int, int> (*comm, REDUCE_MAX, lclErr, outArg (gblErr));
    TEST_EQUALITY_CONST( gblErr, 0 );
    if (gblErr != 0) {
      Tpetra::Details::gathervPrint (out, err.str (), *comm);
      out << "Above test failed; aborting further tests" << std::endl;
      return;
    }
  }




// ===============================================================================


#define UNIT_TEST_GROUP_SC_LO_GO( SC, LO, GO )                   \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( FEMultiVector, doImport, LO, GO, SC )

  TPETRA_ETI_MANGLING_TYPEDEFS()

  // Test CrsMatrix for all Scalar, LO, GO template parameter
  // combinations, and the default Node type.
  TPETRA_INSTANTIATE_SLG_NO_ORDINAL_SCALAR( UNIT_TEST_GROUP_SC_LO_GO )

}
