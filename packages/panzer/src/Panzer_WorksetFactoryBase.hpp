#ifndef __Panzer_WorksetFactoryBase_hpp__
#define __Panzer_WorksetFactoryBase_hpp__

#include <string>
#include <vector>
#include <map>

#include "Panzer_Workset.hpp"
#include "Panzer_BC.hpp"
#include "Panzer_PhysicsBlock.hpp"

namespace panzer {

/** Pure virtual base class used to construct 
  * worksets on volumes and side sets.
  */
class WorksetFactoryBase {
public:
   virtual ~WorksetFactoryBase() {}

   /** Build sets of volume worksets
     */ 
   virtual
   Teuchos::RCP<std::vector<panzer::Workset> >
   getVolumeWorksets(const std::string & eBlock,
                     const panzer::PhysicsBlock & pb,
                     std::size_t worksetSize) const = 0;

   /** Build sets of boundary condition worksets
     */
   virtual
   Teuchos::RCP<std::map<unsigned,panzer::Workset> > 
   getSideWorksets(const panzer::BC & bc,
                 const panzer::PhysicsBlock & pb) const = 0;
};

}

#endif
