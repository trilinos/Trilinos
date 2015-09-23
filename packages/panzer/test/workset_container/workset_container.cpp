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

#include "Panzer_WorksetContainer.hpp"

namespace panzer {

  class unit_test_Factory : public panzer::WorksetFactoryBase {
  public:
     unit_test_Factory() {}
  
     /** Build sets of volume worksets
       */ 
     virtual
     Teuchos::RCP<std::vector<panzer::Workset> >
     getVolumeWorksets(const std::string & eBlock,
                       const panzer::PhysicsBlock & pb,
                       std::size_t worksetSize) const;
  
     /** Build sets of boundary condition worksets
       */
     virtual
     Teuchos::RCP<std::map<unsigned,panzer::Workset> > 
     getSideWorksets(const panzer::BC & bc,
		     const panzer::PhysicsBlock& pb) const; 
  
  };


  BC getBC(const std::string & ss,const std::string & eblk)
  {
     return BC(0,BCT_Dirichlet,ss,eblk,"NONE","NONE");
  }


  bool calledVol, calledSide;

  Teuchos::RCP<std::vector<panzer::Workset> > unit_test_Factory::
  getVolumeWorksets(const std::string & eBlock,
                    const panzer::PhysicsBlock & pb,
                    std::size_t worksetSize) const
  {
     calledVol = true;
     return Teuchos::rcp(new std::vector<panzer::Workset>);
  }
  
  Teuchos::RCP<std::map<unsigned,panzer::Workset> > unit_test_Factory::
  getSideWorksets(const panzer::BC & bc,
                 const panzer::PhysicsBlock & pb) const
  {
     calledSide = true;
     return Teuchos::rcp(new std::map<unsigned,panzer::Workset>);
  }

/*
  // Refactor made this test obsolete, though it should (have) been updated
  TEUCHOS_UNIT_TEST(workset_container, basic)
  {
     std::map<std::string,panzer::InputPhysicsBlock> ebToIpb;
     ebToIpb["block_0"].physics_block_id = "picture";
     ebToIpb["block_1"].physics_block_id = "picture";
     ebToIpb["block_2"].physics_block_id = "swan";

     panzer::buildPhysicsBlocks(block_ids_to_physics_ids,
                                block_ids_to_cell_topo,
                                physics_id_to_input_physics_blocks,
                                2, 1,
                                eqset_factory,
                                false,
                                physicsBlocks);

     Teuchos::RCP<const WorksetFactoryBase> factory 
        = Teuchos::rcp(new unit_test_Factory);
     WorksetContainer wkstCont;
     wkstCont.setFactory(factory);
     wkstCont.setInputPhysicsBlockMap(ebToIpb);

     TEST_ASSERT(wkstCont.getFactory()==factory);
     TEST_ASSERT(wkstCont.getInputPhysicsBlockMap().size()==ebToIpb.size());

     calledVol = false; calledSide = false;     
     wkstCont.getVolumeWorksets("block_0");
     TEST_ASSERT(calledVol); TEST_ASSERT(!calledSide);

     calledVol = false; calledSide = false;     
     wkstCont.getVolumeWorksets("block_0");
     TEST_ASSERT(!calledVol); TEST_ASSERT(!calledSide);

     calledVol = false; calledSide = false;     
     wkstCont.getVolumeWorksets("block_1");
     TEST_ASSERT(calledVol); TEST_ASSERT(!calledSide);

     ///////////////////////

     calledVol = false; calledSide = false;     
     wkstCont.getSideWorksets(getBC("left","block_0"));
     TEST_ASSERT(!calledVol); TEST_ASSERT(calledSide);

     calledVol = false; calledSide = false;     
     wkstCont.getSideWorksets(getBC("left","block_0"));
     TEST_ASSERT(!calledVol); TEST_ASSERT(!calledSide);

     calledVol = false; calledSide = false;     
     wkstCont.getSideWorksets(getBC("left","block_1"));
     TEST_ASSERT(!calledVol); TEST_ASSERT(calledSide);
  }
*/

}
