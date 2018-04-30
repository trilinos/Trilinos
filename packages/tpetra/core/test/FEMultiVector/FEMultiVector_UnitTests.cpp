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

  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( FEMultiVector, doImport, LO, GO, Scalar, NO )
  {
    using Teuchos::VERB_EXTREME;
    using Teuchos::VERB_NONE;
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::Array;
    using Teuchos::Comm;
    using Teuchos::OSTab;
    using Teuchos::outArg;
    using Teuchos::REDUCE_MAX;
    using Teuchos::reduceAll;
    typedef Tpetra::global_size_t GST;

    const GST INVALID = Teuchos::OrdinalTraits<GST>::invalid ();
    const RCP<const Comm<int> > comm = Tpetra::getDefaultComm();
    const int myRank = comm->getRank();
    const int numProcs = comm->getSize();
    if(numProcs==1) return;


    // Prepare for verbose output, if applicable.
    //    Teuchos::EVerbosityLevel verbLevel = verbose ? VERB_EXTREME : VERB_NONE;
    Teuchos::EVerbosityLevel verbLevel = VERB_EXTREME;
    const bool doPrint = includesVerbLevel (verbLevel, VERB_EXTREME, true);
    if (doPrint) {
      out << "FEMultiVector unit test" << std::endl;
    }
    OSTab tab1 (out); // Add one tab level

    std::ostringstream err;
    int lclErr = 0;

    try {
      OSTab tab2 (out);
      const int num_local_elements = 3;

      // create Map
      RCP<const Tpetra::Map<LO, GO, NO> > map =
        rcp( new Tpetra::Map<LO,GO,NO>(INVALID, num_local_elements, 0, comm));

      // create CrsGraph object
      RCP<Tpetra::CrsGraph<LO, GO, NO> > graph =
             rcp (new Tpetra::CrsGraph<LO, GO, NO> (map, 3, Tpetra::DynamicProfile));

      // Create a simple tridiagonal source graph.
      Array<GO> entry(1);
      for (size_t i = 0; i < map->getNodeNumElements (); i++) {
        const GO globalrow = map->getGlobalElement (i);
        entry[0] = globalrow;
        graph->insertGlobalIndices (globalrow, entry());
        if (myRank!=0) {
          entry[0] = globalrow-1;
          graph->insertGlobalIndices (globalrow, entry());
        }
        if (myRank!=numProcs-1) {
          entry[0] = globalrow+1;
          graph->insertGlobalIndices (globalrow, entry());
        }       
      }
      graph->fillComplete();


      Scalar ZERO = Teuchos::ScalarTraits<Scalar>::zero();
      RCP<const Tpetra::Map<LO,GO,NO> > domainMap = graph->getDomainMap();
      RCP<const Tpetra::Import<LO,GO,NO> > importer = graph->getImporter();
      Tpetra::MultiVector<Scalar,LO,GO,NO> Vdomain(domainMap,1), Vcolumn(graph->getColMap(),1);
      Tpetra::FEMultiVector<Scalar,LO,GO,NO> Vfe(domainMap,importer,1);
      if(importer.is_null()) throw std::runtime_error("No valid importer");

      // 1) Test domain -> column
      size_t Ndomain = domainMap->getNodeNumElements();
      for(size_t i=0; i<Ndomain; i++)
        Vdomain.getDataNonConst(0)[i] = domainMap->getGlobalElement(i);
      Vcolumn.putScalar(ZERO);
      Vcolumn.doImport(Vdomain,*importer,Tpetra::ADD);

      Vfe.putScalar(ZERO);
      //      for(size_t i=0; i<Ndomain; i++)
        //Vdomain->getDataNonConst(0)[i] = domainMap->getGlobalElement(i);
      


    } catch (std::exception& e) {
      err << "Proc " << myRank << ": " << e.what () << std::endl;
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


#define UNIT_TEST_GROUP( SC, LO, GO, NO  )                         \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( FEMultiVector, doImport, LO, GO, SC, NO )

  TPETRA_ETI_MANGLING_TYPEDEFS()

  TPETRA_INSTANTIATE_SLGN( UNIT_TEST_GROUP )

}
