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
#include <Teuchos_CommHelpers.hpp>

#include "Kokkos_Core.hpp"
#include "Kokkos_DynRankView.hpp"

#include "Intrepid2_config.h"
#include "Intrepid2_HGRAD_HEX_C1_FEM.hpp"
#include "Intrepid2_HGRAD_HEX_C2_FEM.hpp"
#include "Intrepid2_HDIV_HEX_I1_FEM.hpp"
#include "Intrepid2_HCURL_HEX_I1_FEM.hpp"

#include "PanzerCore_config.hpp"

#include "Panzer_ConnManager.hpp"
#include "Panzer_IntrepidFieldPattern.hpp"
#include "Panzer_DOFManager.hpp"
#include "Panzer_GlobalIndexer_Utilities.hpp"

#include "CartesianConnManager.hpp"

using Teuchos::rcp;
using Teuchos::rcp_dynamic_cast;
using Teuchos::RCP;
using Teuchos::rcpFromRef;

namespace panzer::unit_test {

using Triplet = CartesianConnManager::Triplet<panzer::GlobalOrdinal>;

RCP<const panzer::FieldPattern> buildFieldPattern(
  RCP<Intrepid2::Basis<PHX::Device, double, double>> basis
) {
  // build a geometric pattern from a single basis
  return Teuchos::make_rcp<panzer::Intrepid2FieldPattern>(basis);
}

std::string getElementBlock(
  const Triplet & element,
  const CartesianConnManager & connManager
) {
  int localElmtId = connManager.computeLocalBrickElementIndex(element);
  return connManager.getBlockId(localElmtId);
}

TEUCHOS_UNIT_TEST(tCartesianDOFMgr_DynRankView, threed)
{
  using CCM = CartesianConnManager;
  using DOFManager = panzer::DOFManager;

  // build global (or serial communicator)
  #ifdef HAVE_MPI
    Teuchos::MpiComm<int> comm(MPI_COMM_WORLD);
  #else
    THIS_REALLY_DOES_NOT_WORK
  #endif

  int np = comm.getSize(); // number of processors

  // mesh description
  panzer::GlobalOrdinal nx = 10, ny = 7, nz = 4;
  int px = np, py = 1, pz = 1;
  int bx =  1, by = 2, bz = 1;

  // build velocity, temperature and pressure fields
  RCP<const panzer::FieldPattern> pattern_U = buildFieldPattern(Teuchos::make_rcp<Intrepid2::Basis_HGRAD_HEX_C2_FEM<PHX::Device, double, double>>());
  RCP<const panzer::FieldPattern> pattern_P = buildFieldPattern(Teuchos::make_rcp<Intrepid2::Basis_HGRAD_HEX_C1_FEM<PHX::Device, double, double>>());
  RCP<const panzer::FieldPattern> pattern_T = buildFieldPattern(Teuchos::make_rcp<Intrepid2::Basis_HGRAD_HEX_C1_FEM<PHX::Device, double, double>>());
  RCP<const panzer::FieldPattern> pattern_B = buildFieldPattern(Teuchos::make_rcp<Intrepid2::Basis_HDIV_HEX_I1_FEM< PHX::Device, double, double>>());
  RCP<const panzer::FieldPattern> pattern_E = buildFieldPattern(Teuchos::make_rcp<Intrepid2::Basis_HCURL_HEX_I1_FEM<PHX::Device, double, double>>());

  // build the topology
  const auto connManager = Teuchos::make_rcp<CCM>();
  connManager->initialize(comm,nx,ny,nz,px,py,pz,bx,by,bz);

  // build the dof manager, and assocaite with the topology
  const auto dofManager = Teuchos::make_rcp<DOFManager>();
  dofManager->setConnManager(connManager,*comm.getRawMpiComm());

  // add TEMPERATURE field to all element blocks (MHD and solid)
  dofManager->addField("TEMPERATURE",pattern_T);

  // add velocity (U) and PRESSURE fields to the MHD element block
  dofManager->addField("eblock-0_0_0","UX",pattern_U);
  dofManager->addField("eblock-0_0_0","UY",pattern_U);
  dofManager->addField("eblock-0_0_0","UZ",pattern_U);
  dofManager->addField("eblock-0_0_0","PRESSURE",pattern_P);
  dofManager->addField("eblock-0_0_0","B",pattern_B);
  dofManager->addField("eblock-0_0_0","E",pattern_E);

  // add velocity (U) fields to the solid element block
  dofManager->addField("eblock-0_1_0","UX",pattern_U);
  dofManager->addField("eblock-0_1_0","UY",pattern_U);
  dofManager->addField("eblock-0_1_0","UZ",pattern_U);

  // build global unknowns (useful comment!)
  dofManager->buildGlobalUnknowns();

  // print out some diagnostic information
  ///////////////////////////////////////////////////////////

  dofManager->printFieldInformation(out);

  out << std::endl << "Load balancing: " << printUGILoadBalancingInformation(*dofManager) << std::endl;

  out << std::endl << "Mesh Topology: " << std::endl;
  printMeshTopology(out,*dofManager);

  auto myOffset   = connManager->getMyBrickOffsetTriplet();
  auto myElements = connManager->getMyBrickElementsTriplet();

  out << "My Offset   = " << myOffset.x << " " << myOffset.y << " " << myOffset.z << std::endl;
  out << "My myElements = " << myElements.x << " " << myElements.y << " " << myElements.z << std::endl;

  // Test the kokkos version of field offsets
  {
    const int fieldNumber = dofManager->getFieldNum("B");
    std::vector<std::string> elementBlockNames;
    dofManager->getElementBlockIds(elementBlockNames);
    TEST_ASSERT(elementBlockNames.size() > 0);
    TEST_ASSERT(fieldNumber >= 0);
    const auto& hostOffsetsStdVector = dofManager->getGIDFieldOffsets(elementBlockNames[0],fieldNumber);
    TEST_EQUALITY(hostOffsetsStdVector.size(),6);
    const auto kokkosOffsets = dofManager->getGIDFieldOffsetsKokkos(elementBlockNames[0],fieldNumber);
    const auto hostKokkosOffsets = Kokkos::create_mirror_view(kokkosOffsets);
    Kokkos::deep_copy(hostKokkosOffsets,kokkosOffsets);
    typename PHX::Device().fence();

    TEST_EQUALITY(hostOffsetsStdVector.size(),hostKokkosOffsets.size());
    for (size_t i=0; i < hostOffsetsStdVector.size(); ++i) {
      TEST_EQUALITY(hostOffsetsStdVector[i],hostKokkosOffsets(i));
    }
  }

}

} // namespace panzer::unit test
