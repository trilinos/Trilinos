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

#include "Panzer_WorksetContainer.hpp"

namespace panzer {

//! Default contructor, starts with no workset factory objects
WorksetContainer::WorksetContainer()
   : worksetSize_(1)
{}

WorksetContainer::WorksetContainer(const Teuchos::RCP<const WorksetFactoryBase> & factory,
                                   const std::vector<Teuchos::RCP<PhysicsBlock> > & physicsBlocks,
                                   std::size_t wkstSz)
   : wkstFactory_(factory), worksetSize_(wkstSz)
{
   setPhysicsBlockVector(physicsBlocks);
}

/** Copies the workset factory and the workset size, but not constructed
  * worksets.
  */
WorksetContainer::WorksetContainer(const WorksetContainer & wc)
   : wkstFactory_(wc.wkstFactory_), ebToPb_(wc.ebToPb_), worksetSize_(wc.worksetSize_)
{
}

void WorksetContainer::setPhysicsBlockVector(const std::vector<Teuchos::RCP<PhysicsBlock> > & physicsBlocks)
{
   for(std::size_t i=0;i<physicsBlocks.size();i++) {
      ebToPb_[physicsBlocks[i]->elementBlockID()] = physicsBlocks[i];
   }
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
const PhysicsBlock & WorksetContainer::lookupPhysicsBlock(const std::string & eBlock) const
{
   std::map<std::string,Teuchos::RCP<PhysicsBlock> >::const_iterator itr = ebToPb_.find(eBlock);
 
   TEUCHOS_TEST_FOR_EXCEPTION(itr==ebToPb_.end(),std::logic_error, 
                      "WorksetContainer::lookupPhysicsBlock no PhysicsBlock object is associated "
                      "with the element block \""+eBlock+"\".");

   return *itr->second;
}

//! Access, and construction of volume worksets
Teuchos::RCP<std::vector<Workset> >  
WorksetContainer::getVolumeWorksets(const std::string & eBlock)
{
   Teuchos::RCP<std::vector<Workset> > worksetVector;
   VolumeMap::iterator itr = volWorksets_.find(eBlock);
   if(itr==volWorksets_.end()) {
      // couldn't find workset, build it!
      const PhysicsBlock & pb = lookupPhysicsBlock(eBlock);
      worksetVector = wkstFactory_->getVolumeWorksets(eBlock,pb);

      // store vector for reuse in the future
      volWorksets_[eBlock] = worksetVector;
   }
   else 
      worksetVector = itr->second;

   return worksetVector;
}

Teuchos::RCP<std::vector<Teuchos::RCP<std::vector<Workset> > > > 
WorksetContainer::getVolumeWorksets() const
{
   Teuchos::RCP<std::vector<Teuchos::RCP<std::vector<Workset> > > > worksets =
      Teuchos::rcp(new std::vector<Teuchos::RCP<std::vector<Workset> > >);

   // fill vector with RCP pointers
   for(VolumeMap::const_iterator itr=volWorksets_.begin();itr!=volWorksets_.end();itr++)
      worksets->push_back(itr->second);
 
   return worksets;
}
 
//! Access, and construction of side worksets
Teuchos::RCP<std::map<unsigned,Workset> > 
WorksetContainer::getSideWorksets(const BC & bc)
{
   Teuchos::RCP<std::map<unsigned,Workset> > worksetMap;
   SideId side(bc);
   SideMap::iterator itr = sideWorksets_.find(side);
   if(itr==sideWorksets_.end()) {
      // couldn't find workset, build it!
      const std::string & eBlock = side.eblk_id;
      const PhysicsBlock & pb = lookupPhysicsBlock(eBlock);
      worksetMap = wkstFactory_->getSideWorksets(bc,pb);

      // store map for reuse in the future
      sideWorksets_[side] = worksetMap;
   }
   else { 
      worksetMap = itr->second;
   }

   return worksetMap;
}

void WorksetContainer::allocateVolumeWorksets(const std::vector<std::string> & eBlocks)
{
   for(std::size_t i=0;i<eBlocks.size();i++) {
      // couldn't find workset, build it!
      const std::string & eBlock = eBlocks[i];
      const PhysicsBlock & pb = lookupPhysicsBlock(eBlock);

      // store vector for reuse in the future
      volWorksets_[eBlock] = wkstFactory_->getVolumeWorksets(eBlock,pb);
   }
}

void WorksetContainer::allocateSideWorksets(const std::vector<BC> & bcs)
{
   for(std::size_t i=0;i<bcs.size();i++) {
      // couldn't find workset, build it!
      const BC & bc = bcs[i];
      SideId side(bc);
      const std::string & eBlock = bc.elementBlockID();
      const PhysicsBlock & pb = lookupPhysicsBlock(eBlock);

      // store map for reuse in the future
      sideWorksets_[side] = wkstFactory_->getSideWorksets(bc,pb);
   }
}

void getVolumeWorksetsFromContainer(WorksetContainer & wc,
                                    const std::vector<std::string> & elementBlockNames,
                                    std::map<std::string,Teuchos::RCP<std::vector<Workset> > > & volumeWksts) 
{
   for(std::size_t i=0;i<elementBlockNames.size();i++)
      volumeWksts[elementBlockNames[i]] = wc.getVolumeWorksets(elementBlockNames[i]);
}

void getSideWorksetsFromContainer(WorksetContainer & wc,
                                  const std::vector<BC> & bcs,
                                  std::map<BC,Teuchos::RCP<std::map<unsigned,Workset> >,LessBC> & sideWksts)
{
   for(std::size_t i=0;i<bcs.size();i++) {
      Teuchos::RCP<std::map<unsigned,Workset> > wksts = wc.getSideWorksets(bcs[i]);
      if(wksts!=Teuchos::null)
         sideWksts[bcs[i]] = wksts;
   }
}

}
