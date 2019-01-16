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

using Teuchos::RCP;
using Teuchos::rcp;

#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include "Panzer_STK_Version.hpp"
#include "PanzerAdaptersSTK_config.hpp"

#include "Panzer_WorksetContainer.hpp"
#include "Panzer_IntrepidBasisFactory.hpp"
#include "Panzer_DOFManager.hpp"
#include "Panzer_OrientationsInterface.hpp"

#include "Panzer_STK_Interface.hpp"
#include "Panzer_STK_CubeHexMeshFactory.hpp"
#include "Panzer_STKConnManager.hpp"
#include "Panzer_STK_WorksetFactory.hpp"

namespace panzer {

  TEUCHOS_UNIT_TEST(workset_container, identifiers)
  {
    using Teuchos::RCP;
    using Teuchos::rcp;

    typedef Intrepid2::Basis<PHX::Device::execution_space,double,double> IntrepidBasis;

    std::string element_block = "eblock-0_0_0";
    std::string sideset = "left";
    int workset_size = 10;

    RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
    pl->set("X Elements",8);
    pl->set("Y Elements",8);
    pl->set("Z Elements",8);

    panzer_stk::CubeHexMeshFactory factory;
    factory.setParameterList(pl);
    RCP<panzer_stk::STK_Interface> mesh = factory.buildMesh(MPI_COMM_WORLD);
    RCP<const Teuchos::MpiComm<int> > tComm = Teuchos::rcp(new Teuchos::MpiComm<int>(MPI_COMM_WORLD));

    int myRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

    // build DOF Manager (with a single HDiv basis)
    /////////////////////////////////////////////////////////////
 
    // build the connection manager
    const auto conn_manager = rcp(new panzer_stk::STKConnManager(mesh));
    auto dof_manager = rcp(new panzer::DOFManager(conn_manager,MPI_COMM_WORLD));

    // build an intrepid basis and a related field pattern for seeding the DOFManager
    {
       RCP<IntrepidBasis> intrepid_basis
           = panzer::createIntrepid2Basis<PHX::Device::execution_space,double,double>("HGrad",1,
                                                                                      *mesh->getCellTopology(element_block));
      RCP<panzer::Intrepid2FieldPattern> field_pattern = rcp(new panzer::Intrepid2FieldPattern(intrepid_basis));

      dof_manager->addField(element_block, "T", field_pattern);
    }
    dof_manager->buildGlobalUnknowns();

    // build WorksetContainer
    //////////////////////////////////////////////////////////////

    auto wkstFactory = Teuchos::rcp(new panzer_stk::WorksetFactory(mesh));
    wkstFactory->setOrientationsInterface(Teuchos::rcp(new OrientationsInterface(dof_manager)));
    auto wkstContainer = Teuchos::rcp(new panzer::WorksetContainer(wkstFactory));

    // Test Volume Worksets
    //////////////////////////////////////////////////////////////

    // Double check nothing is generated for a bad input
    TEST_ASSERT(wkstContainer->getWorksets(WorksetDescriptor("element block",workset_size))->size() == 0)

    // Check worksets
    {
      // Grab the worksets for the given element block
      auto worksets_ptr = wkstContainer->getWorksets(WorksetDescriptor(element_block,workset_size));
      const std::vector<Workset> & worksets = *worksets_ptr;

      // Make sure we actually got something back
      TEST_ASSERT(worksets.size()>0)

      // Check for unique identifiers in worksets - no 0 identifiers allowed
      std::set<std::size_t> identifiers;
      for(const auto & wkst : worksets) {
        TEST_ASSERT(wkst.getIdentifier() != 0u)
        TEST_ASSERT(identifiers.find(wkst.getIdentifier()) == identifiers.end());
        identifiers.insert(wkst.getIdentifier());
      }
    }

    // Test Sideset (volume) Worksets
    //////////////////////////////////////////////////////////////

    // Double check nothing is generated for a bad input
    TEST_ASSERT(wkstContainer->getWorksets(WorksetDescriptor("element block", "sideset"))->size() == 0)

    // Check worksets
    {
      auto worksets_ptr = wkstContainer->getWorksets(WorksetDescriptor(element_block, sideset));
      const std::vector<Workset> & worksets = *worksets_ptr;

      // Make sure we actually got something back - but only for process 0 (process 1 should not touch the 'left' sideset)
      if(myRank == 0){
        TEST_ASSERT(worksets.size()>0)
      }

      // Check for unique identifiers in worksets - no 0 identifiers allowed
      std::set<std::size_t> identifiers;
      for(const auto & wkst : worksets) {
        TEST_ASSERT(wkst.getIdentifier() != 0u)
        TEST_ASSERT(identifiers.find(wkst.getIdentifier()) == identifiers.end());
        identifiers.insert(wkst.getIdentifier());
      }

    }

    // Test Sideset (surface) Worksets
    //////////////////////////////////////////////////////////////

    // Double check nothing is generated for a bad input
    TEST_ASSERT(wkstContainer->getWorksets(WorksetDescriptor("element block", "sideset", true))->size() == 0)

    // Check worksets
    {
      auto worksets_ptr = wkstContainer->getWorksets(WorksetDescriptor(element_block, sideset, true));
      const std::vector<Workset> & worksets = *worksets_ptr;

      // Make sure we actually got something back - but only for process 0 (process 1 should not touch the 'left' sideset)
      if(myRank == 0){
        TEST_ASSERT(worksets.size()>0)
      }

      // Check for unique identifiers in worksets - no 0 identifiers allowed
      std::set<std::size_t> identifiers;
      for(const auto & wkst : worksets) {
        TEST_ASSERT(wkst.getIdentifier() != 0u)
        TEST_ASSERT(identifiers.find(wkst.getIdentifier()) == identifiers.end());
        identifiers.insert(wkst.getIdentifier());
      }
    }

    // Test BC Sideset Worksets
    //////////////////////////////////////////////////////////////

    // Our mesh is too simple to test for it working, but we can test it not working just fine

    // If a block is found, but the sidesets are mismatched, throw error (sidesets must have same name - this is a current limitation in the mesh reader)
    TEST_THROW(wkstContainer->getWorksets(WorksetDescriptor("element block", "other block", "sideset", "other set")), std::logic_error)
    TEST_THROW(wkstContainer->getWorksets(WorksetDescriptor(element_block, "other block", "sideset", "other set")), std::logic_error)

    // If the blocks are not joined then don't return any worksets
    TEST_EQUALITY(wkstContainer->getWorksets(WorksetDescriptor(element_block, "other block", "sideset", "sideset"))->size(), 0)
    TEST_EQUALITY(wkstContainer->getWorksets(WorksetDescriptor(element_block, "other block", sideset, sideset))->size(), 0)
    TEST_EQUALITY(wkstContainer->getWorksets(WorksetDescriptor("other block", element_block, sideset, sideset))->size(), 0)

    //////////////////////////////////////////////////////////////
  }
}
