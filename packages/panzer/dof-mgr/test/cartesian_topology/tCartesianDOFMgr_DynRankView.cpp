// @HEADER
// ***********************************************************************
//
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//                 Copyright (2011) Sandia Corporation
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
// Questions? Contact Roger P. Pawlowski (rppawlo@sandia.gov) and
// Eric C. Cyr (eccyr@sandia.gov)
// ***********************************************************************
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
  int rank = comm.getRank(); // processor rank

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

  // int temp_num = dofManager->getFieldNum("TEMPERATURE");
  int ux_num   = dofManager->getFieldNum("UX");
  int uy_num   = dofManager->getFieldNum("UY");
  // int uz_num   = dofManager->getFieldNum("UZ");
  // int p_num    = dofManager->getFieldNum("PRESSURE");
  // int b_num    = dofManager->getFieldNum("B");
  // int e_num    = dofManager->getFieldNum("E");

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
