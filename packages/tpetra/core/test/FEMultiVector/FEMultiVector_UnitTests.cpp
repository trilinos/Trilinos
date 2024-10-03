// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Tpetra_TestingUtilities.hpp"
#include "Tpetra_CrsGraph.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_Core.hpp"
#include "Tpetra_Distributor.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_Util.hpp"
#include "Tpetra_Import.hpp"
#include "Tpetra_FEMultiVector.hpp"
#include "Tpetra_Assembly_Helpers.hpp"
#include "Tpetra_Details_gathervPrint.hpp"
#include "Teuchos_OrdinalTraits.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_VerboseObject.hpp"
#include <map>

namespace {

  template<class T1, class T2>
  void vector_check(size_t N, T1 & v1, T2 & v2) {
    const int myRank = v1.getMap()->getComm()->getRank();
    for(size_t i=0; i<N; i++)
      if(v1.getDataNonConst(0)[i] != v2.getDataNonConst(0)[i]) {
        std::stringstream oss;
        oss<<"["<<myRank<<"] Error: Mismatch on unknown "<<i<<" "<<v1.getDataNonConst(0)[i]<<" != "<<v2.getDataNonConst(0)[i]<<std::endl;
        throw std::runtime_error(oss.str());
      }
  }

  using std::endl;

  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::Array;
  using Teuchos::Comm;
  using Teuchos::OrdinalTraits;
  using Teuchos::ScalarTraits;
  using Teuchos::VERB_DEFAULT;
  using Teuchos::VERB_NONE;
  using Teuchos::VERB_LOW;
  using Teuchos::VERB_MEDIUM;
  using Teuchos::VERB_HIGH;
  using Teuchos::VERB_EXTREME;

  using Tpetra::createContigMapWithNode;
  using GST = Tpetra::global_size_t;

  //
  // UNIT TESTS
  //

  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( FEMultiVector, doImport, LO, GO, Scalar, NO )
  {
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
    Teuchos::OSTab tab1 (out); // Add one tab level

    std::ostringstream err;
    int lclErr = 0;

    try {
      Teuchos::OSTab tab2 (out);
      const int num_local_elements = 3;

      // create Map
      RCP<const Tpetra::Map<LO, GO, NO> > map =
        rcp( new Tpetra::Map<LO,GO,NO>(INVALID, num_local_elements, 0, comm));

      // create CrsGraph object
      RCP<Tpetra::CrsGraph<LO, GO, NO> > graph =
             rcp (new Tpetra::CrsGraph<LO, GO, NO> (map, 3));

      // Create a simple tridiagonal source graph.
      Array<GO> entry(1);
      for (size_t i = 0; i < map->getLocalNumElements (); i++) {
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
      Scalar ONE = Teuchos::ScalarTraits<Scalar>::one();
      RCP<const Tpetra::Map<LO,GO,NO> > domainMap = graph->getDomainMap();
      RCP<const Tpetra::Map<LO,GO,NO> > columnMap = graph->getColMap();
      RCP<const Tpetra::Import<LO,GO,NO> > importer = graph->getImporter();
      Tpetra::MultiVector<Scalar,LO,GO,NO> Vdomain(domainMap,1), Vcolumn(columnMap,1);
      Tpetra::FEMultiVector<Scalar,LO,GO,NO> Vfe(domainMap,importer,1);
      size_t Ndomain = domainMap->getLocalNumElements();
      size_t Ncolumn = domainMap->getLocalNumElements();

      if(importer.is_null()) throw std::runtime_error("No valid importer");

      // 1) Test column -> domain (no off-proc addition)
      Vcolumn.putScalar(ZERO);
      for(size_t i=0; i<Ndomain; i++)
        Vcolumn.getDataNonConst(0)[i] = domainMap->getGlobalElement(i);
      Vdomain.doExport(Vcolumn,*importer,Tpetra::ADD);


      Vfe.beginAssembly();
      Vfe.putScalar(ZERO);
      for(size_t i=0; i<Ndomain; i++)
        Vfe.getDataNonConst(0)[i] = domainMap->getGlobalElement(i);
      Vfe.endAssembly();
      vector_check(Ndomain,Vfe,Vdomain);

      // 2) Test column -> domain (with off-proc addition)
      Vdomain.putScalar(ZERO);
      Vcolumn.putScalar(ONE);
      Vdomain.doExport(Vcolumn,*importer,Tpetra::ADD);

      Vfe.putScalar(ZERO);
      Vfe.beginAssembly();
      Vfe.putScalar(ONE);
      Vfe.endAssembly();
      vector_check(Ncolumn,Vfe,Vdomain);
    } catch (std::exception& e) {
      err << "Proc " << myRank << ": " << e.what () << std::endl;
      lclErr = 1;
    }


    int gblErr = 0;
    Teuchos::reduceAll<int, int> (*comm, Teuchos::REDUCE_MAX, lclErr, Teuchos::outArg (gblErr));
    TEST_EQUALITY_CONST( gblErr, 0 );
    if (gblErr != 0) {
      Tpetra::Details::gathervPrint (out, err.str (), *comm);
      out << "Above test failed; aborting further tests" << std::endl;
      return;
    }
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( FEMultiVector, AssemblyHelpers, LO , GO , Scalar , Node )
  {
    const GST INVALID = Teuchos::OrdinalTraits<GST>::invalid ();
    const RCP<const Comm<int> > comm = Tpetra::getDefaultComm();
    int rank = comm->getRank();
    int size = comm->getSize();

    // create maps
    const size_t numLocal = 10;
    const size_t numOverlap = numLocal + (rank!=0) + (rank!=size-1);
    RCP<const Tpetra::Map<LO,GO,Node> > map = createContigMapWithNode<LO,GO,Node>(INVALID,numLocal,comm);

    Array<GO> overlapList(numOverlap);
    for(size_t i = 0; i< numLocal; i++) {
      overlapList[i] = map->getGlobalElement(i);
    }
    size_t iii = numLocal;
    if(rank!=0)      {overlapList[iii] = overlapList[0]-1; iii++;}
    if(rank!=size-1) {overlapList[iii] = overlapList[numLocal-1]+1;}

    RCP<const Tpetra::Map<LO,GO,Node> > overlapMap    = rcp(new Tpetra::Map<LO,GO,Node>(INVALID,overlapList(),0,comm));
    RCP<const Tpetra::Import<LO,GO,Node> > importer   = rcp(new Tpetra::Import<LO,GO,Node>(map,overlapMap));

    Tpetra::FEMultiVector<Scalar,LO,GO,Node> v1(map,importer,1);
    Tpetra::FEMultiVector<Scalar,LO,GO,Node> v2(map,importer,1);
    Tpetra::FEMultiVector<Scalar,LO,GO,Node> v3(map,importer,1);

    // Just check to make sure beginAssembly() / endAssembly() compile
    Tpetra::beginAssembly(v1,v2,v3);

    Tpetra::endAssembly(v1,v2,v3);
  }

#define UNIT_TEST_GROUP( SC, LO, GO, NO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( FEMultiVector, doImport, LO, GO, SC, NO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( FEMultiVector, AssemblyHelpers, LO, GO, SC, NO )

  TPETRA_ETI_MANGLING_TYPEDEFS()

  TPETRA_INSTANTIATE_TESTMV( UNIT_TEST_GROUP )
}
