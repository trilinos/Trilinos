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
#include "Panzer_WorksetDescriptor.hpp"
#include "Panzer_WorksetNeeds.hpp"
#include "Panzer_IntrepidBasisFactory.hpp"
#include "Panzer_IntegrationDescriptor.hpp"
#include "Panzer_DOFManager.hpp"

#include "Panzer_STK_Interface.hpp"
#include "Panzer_STK_CubeHexMeshFactory.hpp"
#include "Panzer_STKConnManager.hpp"
#include "Panzer_STK_WorksetFactory.hpp"

namespace panzer {

  TEUCHOS_UNIT_TEST(check_rotation_matrices, basic)
  {
    using Teuchos::RCP;
    using Teuchos::rcp;

    typedef Intrepid2::Basis<PHX::Device::execution_space,double,double> IntrepidBasis;

    std::string element_block = "eblock-0_0_0";
    std::string sideset = "left";

    RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
    pl->set("X Elements",1);
    pl->set("Y Elements",1);
    pl->set("Z Elements",1);
    pl->set("X0",-0.001);
    pl->set("Xf", 0.001);
    pl->set("Y0",-0.0005);
    pl->set("Yf", 0.0005);
    pl->set("Z0",-0.0003333333333333333);
    pl->set("Zf", 0.0003333333333333333);

    panzer_stk::CubeHexMeshFactory factory;
    factory.setParameterList(pl);
    RCP<panzer_stk::STK_Interface> mesh = factory.buildMesh(MPI_COMM_WORLD);
    RCP<const Teuchos::MpiComm<int> > tComm = Teuchos::rcp(new Teuchos::MpiComm<int>(MPI_COMM_WORLD));

    // build DOF Manager (with a single HDiv basis)
    /////////////////////////////////////////////////////////////
 
    // build the connection manager 
    const RCP<panzer::ConnManager<int,panzer::Ordinal64> > 
      conn_manager = rcp(new panzer_stk::STKConnManager<panzer::Ordinal64>(mesh));

    RCP<panzer::DOFManager<int,panzer::Ordinal64> > dof_manager
        = rcp(new panzer::DOFManager<int,panzer::Ordinal64>(conn_manager,MPI_COMM_WORLD));

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
    
    panzer::IntegrationDescriptor sid(2*2, panzer::IntegrationDescriptor::SURFACE);
    std::map<std::string, panzer::WorksetNeeds> wkstRequirements;
    wkstRequirements[element_block].addIntegrator(sid);

    RCP<panzer_stk::WorksetFactory> wkstFactory 
       = rcp(new panzer_stk::WorksetFactory(mesh)); // build STK workset factory
    RCP<panzer::WorksetContainer> wkstContainer     // attach it to a workset container (uses lazy evaluation)
       = rcp(new panzer::WorksetContainer(wkstFactory,wkstRequirements));

    wkstContainer->setGlobalIndexer(dof_manager);

    panzer::WorksetDescriptor workset_descriptor(element_block, panzer::WorksetSizeType::ALL_ELEMENTS, true,false);

    auto worksets = wkstContainer->getWorksets(workset_descriptor);

    TEST_ASSERT(worksets->size()==1);
    
    auto rot_matrices = (*worksets)[0].getIntegrationValues(sid).surface_rotation_matrices;
    auto normals = (*worksets)[0].getIntegrationValues(sid).surface_normals;

    TEST_ASSERT(rot_matrices.rank()==4);
    TEST_ASSERT(rot_matrices.extent_int(0)==7); // 7 cells (6 virtual, one owned)
    TEST_ASSERT(rot_matrices.extent_int(2)==3);
    TEST_ASSERT(rot_matrices.extent_int(3)==3);
    out << "SIZES = " 
        << rot_matrices.extent(0) << " " 
        << rot_matrices.extent(1) << " " 
        << rot_matrices.extent(2) << " " 
        << rot_matrices.extent(3) << std::endl;

    out << std::endl;

    for(int c=0;c<rot_matrices.extent_int(0);c++) {
      for(int p=0;p<rot_matrices.extent_int(1);p++) {
        out << "Cell,Point = " << c << "," << p << std::endl;
        out << "   " << rot_matrices(c,p,0,0) << " " << rot_matrices(c,p,0,1) << " " << rot_matrices(c,p,0,2) << std::endl;
        out << "   " << rot_matrices(c,p,1,0) << " " << rot_matrices(c,p,1,1) << " " << rot_matrices(c,p,1,2) << std::endl;
        out << "   " << rot_matrices(c,p,2,0) << " " << rot_matrices(c,p,2,1) << " " << rot_matrices(c,p,2,2) << std::endl;
        out << std::endl;

        // mutually orthongonal
        TEST_ASSERT(std::fabs(rot_matrices(c,p,0,0) * rot_matrices(c,p,1,0) +
                              rot_matrices(c,p,0,1) * rot_matrices(c,p,1,1) +
                              rot_matrices(c,p,0,2) * rot_matrices(c,p,1,2))<=1e-14);
        TEST_ASSERT(std::fabs(rot_matrices(c,p,0,0) * rot_matrices(c,p,2,0) +
                              rot_matrices(c,p,0,1) * rot_matrices(c,p,2,1) +
                              rot_matrices(c,p,0,2) * rot_matrices(c,p,2,2))<=1e-14);
        TEST_ASSERT(std::fabs(rot_matrices(c,p,1,0) * rot_matrices(c,p,2,0) +
                              rot_matrices(c,p,1,1) * rot_matrices(c,p,2,1) +
                              rot_matrices(c,p,1,2) * rot_matrices(c,p,2,2))<=1e-14);

        // normalized
        TEST_FLOATING_EQUALITY(std::sqrt(rot_matrices(c,p,0,0) * rot_matrices(c,p,0,0) +
                                         rot_matrices(c,p,0,1) * rot_matrices(c,p,0,1) +
                                         rot_matrices(c,p,0,2) * rot_matrices(c,p,0,2)),1.0,1e-14);
        TEST_FLOATING_EQUALITY(std::sqrt(rot_matrices(c,p,1,0) * rot_matrices(c,p,1,0) +
                                         rot_matrices(c,p,1,1) * rot_matrices(c,p,1,1) +
                                         rot_matrices(c,p,1,2) * rot_matrices(c,p,1,2)),1.0,1e-14);
        TEST_FLOATING_EQUALITY(std::sqrt(rot_matrices(c,p,2,0) * rot_matrices(c,p,2,0) +
                                         rot_matrices(c,p,2,1) * rot_matrices(c,p,2,1) +
                                         rot_matrices(c,p,2,2) * rot_matrices(c,p,2,2)),1.0,1e-14);
      }
    }
  }
}
