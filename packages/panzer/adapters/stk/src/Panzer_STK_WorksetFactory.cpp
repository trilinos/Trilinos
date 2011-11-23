#include "Panzer_STK_WorksetFactory.hpp"

#include "Panzer_WorksetFactoryBase.hpp"
#include "Panzer_STK_SetupUtilities.hpp"

namespace panzer_stk {

/** Build sets of volume worksets
  */ 
Teuchos::RCP<std::vector<panzer::Workset> > WorksetFactory::
getVolumeWorksets(const std::string & eBlock,
                  const panzer::PhysicsBlock & pb,
                  std::size_t worksetSize) const
{
   return panzer_stk::buildWorksets(*mesh_, pb, worksetSize);
}

/** Build sets of boundary condition worksets
  */
Teuchos::RCP<std::map<unsigned,panzer::Workset> > WorksetFactory::
getSideWorksets(const panzer::BC & bc,
              const panzer::InputPhysicsBlock & pb) const
{
   return panzer_stk::buildBCWorksets(*mesh_,pb,bc);
}

}
