// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_DefaultMpiComm.hpp>

#include <string>
#include <iostream>
#include <vector>
#include <set>

#include "Thyra_SpmdVectorSpaceBase.hpp"

#include "Panzer_DOFManager.hpp"
#include "Panzer_Filtered_GlobalIndexer.hpp"
#include "Panzer_IntrepidFieldPattern.hpp"
#include "Panzer_PauseToAttach.hpp"
#include "Panzer_BlockedTpetraLinearObjFactory.hpp"
#include "Panzer_TpetraLinearObjFactory.hpp"
#include "Panzer_LOCPair_GlobalEvaluationData.hpp"

#include "UnitTest_ConnManager.hpp"

// include some intrepid basis functions
// 2D basis
#include "Intrepid2_HGRAD_QUAD_C1_FEM.hpp"

#include "Kokkos_DynRankView.hpp"

using Teuchos::rcp;
using Teuchos::RCP;
using Teuchos::rcp_dynamic_cast;
using Teuchos::rcpFromRef;

typedef Kokkos::DynRankView<double, PHX::Device> FieldContainer;

namespace panzer
{
  typedef double ScalarT;
  using LocalOrdinalT = panzer::LocalOrdinal;
  using GlobalOrdinalT = panzer::GlobalOrdinal;
  using TpetraVector = Tpetra::Vector<ScalarT, LocalOrdinalT, GlobalOrdinalT>;

  typedef panzer::TpetraLinearObjFactory<panzer::Traits, ScalarT, LocalOrdinalT, GlobalOrdinalT> TpetraLinObjFactoryType;
  typedef panzer::TpetraLinearObjContainer<ScalarT, LocalOrdinalT, GlobalOrdinalT> TpetraLinObjContainerType;
  typedef panzer::BlockedTpetraLinearObjFactory<panzer::Traits, ScalarT, LocalOrdinalT, GlobalOrdinalT> BlockedTpetraLinObjFactoryType;
  typedef panzer::BlockedTpetraLinearObjContainer<ScalarT, LocalOrdinalT, GlobalOrdinalT> BlockedTpetraLinObjContainerType;

  template <typename Intrepid2Type>
  Teuchos::RCP<const panzer::FieldPattern> buildFieldPattern()
  {
    // build a geometric pattern from a single basis
    Teuchos::RCP<Intrepid2::Basis<PHX::exec_space, double, double>> basis = rcp(new Intrepid2Type);
    Teuchos::RCP<const panzer::FieldPattern> pattern = rcp(new panzer::Intrepid2FieldPattern(basis));
    return pattern;
  }

  // this just excercises a bunch of functions
  TEUCHOS_UNIT_TEST(Tpetra_LOF_FilteredUGI, tpetra_lof)
  {
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::rcp_dynamic_cast;

    typedef Thyra::SpmdVectorSpaceBase<double> SpmdSpace;

    // build global (or serial communicator)
#ifdef HAVE_MPI
    Teuchos::RCP<const Teuchos::MpiComm<int>> tComm = Teuchos::rcp(new Teuchos::MpiComm<int>(MPI_COMM_WORLD));
#else
    Teuchos::RCP<const Teuchos::SerialComm<int>> tComm = Teuchos::rcp(new Teuchos::SerialComm<int>(MPI_COMM_WORLD));
#endif

    int myRank = tComm->getRank();
    int numProc = tComm->getSize();

    RCP<ConnManager> connManager = rcp(new unit_test::ConnManager(myRank, numProc));
    RCP<DOFManager> dofManager = rcp(new DOFManager);
    dofManager->setConnManager(connManager, MPI_COMM_WORLD);

    RCP<const panzer::FieldPattern> patternC1 = buildFieldPattern<Intrepid2::Basis_HGRAD_QUAD_C1_FEM<PHX::exec_space, double, double>>();

    dofManager->addField("T", patternC1); // add it to all three blocks
    dofManager->addField("block_0", "Ux", patternC1);
    dofManager->addField("block_0", "Uy", patternC1);
    dofManager->addField("block_0", "P", patternC1);
    dofManager->addField("block_2", "rho", patternC1);

    dofManager->buildGlobalUnknowns();

    // get GIDs on the "bottom" side of local cell 0, use those as the filtering
    // GIDs
    std::vector<panzer::GlobalOrdinal> filtered;
    {
      std::pair<std::vector<int>, std::vector<int>> fieldOffsets = dofManager->getGIDFieldOffsets_closure("block_0", dofManager->getFieldNum("Ux"), 1, 0);

      std::vector<panzer::GlobalOrdinal> gids;
      dofManager->getElementGIDs(0, gids);

      filtered.resize(fieldOffsets.first.size());
      for (std::size_t i = 0; i < fieldOffsets.first.size(); i++)
        filtered[i] = gids[fieldOffsets.first[i]];
    }

    TEST_EQUALITY(filtered.size(), 2);

    RCP<Filtered_GlobalIndexer> filtered_ugi = rcp(new Filtered_GlobalIndexer);
    filtered_ugi->initialize(dofManager, filtered);

    out << "check out ownsership" << std::endl;
    {
      std::vector<bool> isOwned;
      dofManager->ownedIndices(filtered, isOwned);

      TEST_EQUALITY(isOwned.size(), 2);
      out << "IS OWNED: " << isOwned[0] << " " << isOwned[1] << std::endl;

      filtered_ugi->ownedIndices(filtered, isOwned);

      TEST_EQUALITY(isOwned.size(), 2);
      TEST_EQUALITY(isOwned[0], false);
      TEST_EQUALITY(isOwned[1], false);
    }

    out << "test out sizes construction by LOF" << std::endl;
    {
      std::vector<panzer::GlobalOrdinal> indices_f;
      filtered_ugi->getOwnedIndices(indices_f);

      RCP<TpetraLinObjFactoryType> lof = rcp(new TpetraLinObjFactoryType(tComm, filtered_ugi));
      TEST_EQUALITY(rcp_dynamic_cast<const SpmdSpace>(lof->getThyraDomainSpace(), true)->localSubDim(), Teuchos::as<int>(indices_f.size()));
      TEST_EQUALITY(rcp_dynamic_cast<const SpmdSpace>(lof->getThyraRangeSpace(), true)->localSubDim(), Teuchos::as<int>(indices_f.size()));

      RCP<Thyra::LinearOpBase<double>> A = lof->getThyraMatrix();
      TEST_EQUALITY(rcp_dynamic_cast<const SpmdSpace>(A->range(), true)->localSubDim(), Teuchos::as<int>(indices_f.size()));
      TEST_EQUALITY(rcp_dynamic_cast<const SpmdSpace>(A->domain(), true)->localSubDim(), Teuchos::as<int>(indices_f.size()));

      // This next chunk of code tests to ensure parallel communication works as
      // expected, in particular that a filtered owned vector can be used with
      // an unfiltered ghosted vector and that the entries in the ghosted vector
      // will be preserved.
      ////////////////////////////////////////////////////////////////////////////////////

      // check that parallel communication works as expected
      auto importer = lof->getGhostedImport();
      auto ghostedMap = lof->getGhostedMap();
      auto ownedMap = lof->getMap();

      TEST_ASSERT(importer != Teuchos::null);
      TEST_ASSERT(ghostedMap != Teuchos::null);
      TEST_ASSERT(ownedMap != Teuchos::null);

      auto locPair = Teuchos::rcp(new LOCPair_GlobalEvaluationData(lof, panzer::LinearObjContainer::X));

      auto ghosted_t_loc = rcp_dynamic_cast<TpetraLinObjContainerType>(locPair->getGhostedLOC());
      auto ghosted_x_vec = ghosted_t_loc->get_x();
      ghosted_x_vec->putScalar(0.0);

      auto global_t_loc = rcp_dynamic_cast<TpetraLinObjContainerType>(locPair->getGlobalLOC());
      auto global_x_vec = global_t_loc->get_x();
      global_x_vec->putScalar(1.0);

      ghosted_x_vec->doImport(*global_x_vec, *importer, Tpetra::INSERT);

      // check all the filtered indices remain 0
      int count = 0;
      auto xGhostedData = ghosted_x_vec->getData(0);
      for (std::size_t i = 0; i < filtered.size(); i++)
      {
        auto lid = ghostedMap->getLocalElement(filtered[i]);
        if (lid >= 0)
        {
          count++;
          TEST_EQUALITY(xGhostedData[lid], 0.0);
        }
      }

      std::size_t sum = 0;
      for (std::size_t i = 0; i < ghostedMap->getLocalNumElements(); i++)
      {
        sum += xGhostedData[i];
      }

      // check that ones sum up to the number of ids
      // that were not filtered
      TEST_EQUALITY(sum, ghosted_x_vec->getLocalLength() - count);

      // do a lazy test to ensure you can construct an owned matrix
      RCP<Thyra::LinearOpBase<double>> ownedMatrix = lof->getThyraMatrix();
      TEST_ASSERT(ownedMatrix != Teuchos::null);
    }
  }

}
