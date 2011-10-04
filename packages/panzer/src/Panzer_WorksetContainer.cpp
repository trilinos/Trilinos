#include "Panzer_WorksetContainer.hpp"

namespace panzer {

//! Default contructor, starts with no workset factory objects
WorksetContainer::WorksetContainer()
   : worksetSize_(1)
{}

/** Instantiate a workset object with a specified factory and input physics block
  * map.
  *
  * \param[in] factory Factory to be used for constructing worksets
  * \param[in] eb_to_ipb Map from "element blocks" to "input physics block" objects
  */ 
WorksetContainer::WorksetContainer(const Teuchos::RCP<const WorksetFactoryBase> & factory,
                                   const std::map<std::string,InputPhysicsBlock> & ebToIpb,
                                   std::size_t wkstSz)
   : wkstFactory_(factory), ebToIpb_(ebToIpb), worksetSize_(wkstSz)
{
}

/** Copies the workset factory, the InputPhysicsBlock map, and the workset size, but not constructed
  * worksets.
  */
WorksetContainer::WorksetContainer(const WorksetContainer & wc)
   : wkstFactory_(wc.wkstFactory_), ebToIpb_(wc.ebToIpb_), worksetSize_(wc.worksetSize_)
{
}

/** Clear all allocated worksets, maintain the workset factory and element to physics
  * block map.
  */ 
void WorksetContainer::clear()
{
   volWorksets_.clear();
   sideWorksets_.clear();
}

//! Look up an input physics block, throws an exception if it can be found.
const InputPhysicsBlock & WorksetContainer::lookupInputPhysicsBlock(const std::string & eBlock) const
{
   std::map<std::string,InputPhysicsBlock>::const_iterator itr = ebToIpb_.find(eBlock);
 
   TEST_FOR_EXCEPTION(itr==ebToIpb_.end(),std::logic_error, 
                      "WorksetContainer::lookupInputPhysicsBlock no InputPhysicsBlock object is associated "
                      "with the element block \""+eBlock+"\".");

   return itr->second;
}

//! Access, and construction of volume worksets
std::vector<Workset> & 
WorksetContainer::getVolumeWorksets(const std::string & eBlock)
{
   Teuchos::RCP<std::vector<Workset> > worksetVector;
   VolumeMap::iterator itr = volWorksets_.find(eBlock);
   if(itr==volWorksets_.end()) {
      // couldn't find workset, build it!
      const InputPhysicsBlock & ipb = lookupInputPhysicsBlock(eBlock);
      worksetVector = wkstFactory_->getVolumeWorksets(eBlock,ipb,worksetSize_);

      // store vector for reuse in the future
      volWorksets_[eBlock] = worksetVector;
   }
   else 
      worksetVector = itr->second;

   return *worksetVector;
}
 
//! Access, and construction of side worksets
std::map<unsigned,Workset> & 
WorksetContainer::getSideWorksets(const BC & bc)
{
   Teuchos::RCP<std::map<unsigned,Workset> > worksetMap;
   BCMap::iterator itr = sideWorksets_.find(bc);
   if(itr==sideWorksets_.end()) {
      // couldn't find workset, build it!
      const std::string & eBlock = bc.elementBlockID();
      const InputPhysicsBlock & ipb = lookupInputPhysicsBlock(eBlock);
      worksetMap = wkstFactory_->getSideWorksets(bc,ipb);

      // store map for reuse in the future
      sideWorksets_[bc] = worksetMap;
   }
   else 
      worksetMap = itr->second;

   return *worksetMap;
}

void WorksetContainer::allocateVolumeWorksets(const std::vector<std::string> & eBlocks)
{
   for(std::size_t i=0;i<eBlocks.size();i++) {
      // couldn't find workset, build it!
      const std::string & eBlock = eBlocks[i];
      const InputPhysicsBlock & ipb = lookupInputPhysicsBlock(eBlock);

      // store vector for reuse in the future
      volWorksets_[eBlock] = wkstFactory_->getVolumeWorksets(eBlock,ipb,worksetSize_);
   }
}

void WorksetContainer::allocateSideWorksets(const std::vector<BC> & bcs)
{
   for(std::size_t i=0;i<bcs.size();i++) {
      // couldn't find workset, build it!
      const BC & bc = bcs[i];
      const std::string & eBlock = bc.elementBlockID();
      const InputPhysicsBlock & ipb = lookupInputPhysicsBlock(eBlock);

      // store map for reuse in the future
      sideWorksets_[bc] = wkstFactory_->getSideWorksets(bc,ipb);
   }
}

}
