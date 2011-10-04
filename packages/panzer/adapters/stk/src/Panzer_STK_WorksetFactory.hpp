#ifndef __Panzer_STK_WorksetFactory_hpp__
#define __Panzer_STK_WorksetFactory_hpp__

#include "Panzer_WorksetFactoryBase.hpp"

#include "Panzer_STK_Interface.hpp"

namespace panzer_stk {

/** Pure virtual base class used to construct 
  * worksets on volumes and side sets.
  */
class WorksetFactory : public panzer::WorksetFactoryBase {
public:
   WorksetFactory(Teuchos::RCP<const STK_Interface> & mesh) : mesh_(mesh) {}
   virtual ~WorksetFactory() {}

   /** Build sets of volume worksets
     */ 
   virtual
   Teuchos::RCP<std::vector<panzer::Workset> >
   getVolumeWorksets(const std::string & eBlock,
                     const panzer::InputPhysicsBlock & pb,
                     std::size_t worksetSize) const;

   /** Build sets of boundary condition worksets
     */
   virtual
   Teuchos::RCP<std::map<unsigned,panzer::Workset> > 
   getSideWorksets(const panzer::BC & bc,
                 const panzer::InputPhysicsBlock & pb) const;

private:

   Teuchos::RCP<const STK_Interface> mesh_;

};

}

#endif
