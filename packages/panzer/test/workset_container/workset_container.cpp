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
                    const panzer::InputPhysicsBlock & pb) const; 
  
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
                 const panzer::InputPhysicsBlock & pb) const
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
