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

#include "Panzer_Workset_Utilities.hpp"
#include "Panzer_CommonArrayFactories.hpp"
#include "Panzer_OrientationContainer.hpp"
#include "Panzer_Dimension.hpp"

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
   const WorksetDescriptor wd = blockDescriptor(eBlock);

   return getWorksets(wd);
}

Teuchos::RCP<std::vector<Workset> >  
WorksetContainer::getWorksets(const WorksetDescriptor & wd)
{
   Teuchos::RCP<std::vector<Workset> > worksetVector;
   VolumeMap::iterator itr = volWorksets_.find(wd);
   if(itr==volWorksets_.end()) {
      // couldn't find workset, build it!
      const PhysicsBlock & pb = lookupPhysicsBlock(wd.getElementBlock());
      worksetVector = wkstFactory_->getWorksets(wd,pb);

      // apply orientations to the just constructed worksets
      if(worksetVector!=Teuchos::null)
        applyOrientations(wd.getElementBlock(),*worksetVector);

      // store vector for reuse in the future
      volWorksets_[wd] = worksetVector;
   }
   else 
      worksetVector = itr->second;

   return worksetVector;
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

      // apply orientations to the worksets for this side
      if(worksetMap!=Teuchos::null)
        applyOrientations(side,*worksetMap);

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
      const WorksetDescriptor wd = blockDescriptor(eBlock);
      volWorksets_[eBlock] = wkstFactory_->getWorksets(wd,pb);

      // apply orientations to the worksets for this side
      if(volWorksets_[eBlock]!=Teuchos::null)
        applyOrientations(eBlock,*volWorksets_[eBlock]);
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

      // apply orientations to the worksets for this side
      if(sideWorksets_[side]!=Teuchos::null)
        applyOrientations(side,*sideWorksets_[side]);
   }
}

void WorksetContainer::
setGlobalIndexer(const Teuchos::RCP<const panzer::UniqueGlobalIndexerBase> & ugi)
{
  // apply the orientations for stored worksets
  applyOrientations(ugi);
}

void WorksetContainer::
applyOrientations(const Teuchos::RCP<const panzer::UniqueGlobalIndexerBase> & ugi)
{
  // this gurantees orientations won't accidently be applied twice.
  TEUCHOS_ASSERT(globalIndexer_==Teuchos::null);

  globalIndexer_ = ugi;

  // loop over volume worksets, apply orientations to each
  for(VolumeMap::iterator itr=volWorksets_.begin();
      itr!=volWorksets_.end();++itr) {
    std::string eBlock = itr->first.getElementBlock();
   
    applyOrientations(eBlock,*itr->second);
  }

  // loop over side worksets, apply orientations to each
  for(SideMap::iterator itr=sideWorksets_.begin();
      itr!=sideWorksets_.end();itr++) {
    SideId sideId = itr->first;

    applyOrientations(sideId,*itr->second);
  }
}

void WorksetContainer::
applyOrientations(const std::string & eBlock,std::vector<Workset> & worksets) const
{
  using Teuchos::RCP;

  typedef double Scalar;                          // orientation container scalar type
  typedef Intrepid::FieldContainer<Scalar> Array; // orientation container array type
  typedef std::pair<std::string,Teuchos::RCP<const PureBasis> > StrConstBasisPair;

  /////////////////////////////////
  // this is for volume worksets //
  /////////////////////////////////

  // short circuit if no global indexer exists
  if(globalIndexer_==Teuchos::null) { 
    Teuchos::FancyOStream fout(Teuchos::rcpFromRef(std::cout));
    fout.setOutputToRootOnly(0);
 
    fout << "Panzer Warning: No global indexer assigned to a workset container. "
         << "Orientation of the basis for edge basis functions cannot be applied, "
         << "if those basis functions are used, there will be problems!" << std::endl;
    return;
  }

  // extract a map from basis name to string/basis pairs that require orientations
  //////////////////////////////////////////////////////////////////////////////////
  std::map<std::string, StrConstBasisPair> orientedBasisFields;
  {
    const PhysicsBlock & pb = lookupPhysicsBlock(eBlock);

    // get all field/basis pairs
    std::vector<StrConstBasisPair> fieldBasisPairs;
    pb.getFieldLibraryBase()->basisPairs(fieldBasisPairs);
  
    // now use a map to extract only unique basis functions that
    // require orientations
    for(std::size_t i=0;i<fieldBasisPairs.size();i++) {
      StrConstBasisPair pair = fieldBasisPairs[i];
  
      if(pair.second->requiresOrientations()) {
        // name is guranteed to be unique for basis. 
        // just removes duplicate basis objects
        orientedBasisFields[pair.second->name()] = pair;
      }
    }
  }

  // loop over each basis requiring orientations, then apply them
  //////////////////////////////////////////////////////////////////////////////////

  // Note: It may be faster to loop over the basis pairs on the inside (not really sure)
  
  for(std::map<std::string, StrConstBasisPair>::const_iterator itr=orientedBasisFields.begin();
      itr!=orientedBasisFields.end();++itr) {
    // extract field name and basis name 
    const std::string & fieldName = itr->second.first;
    const PureBasis & basis = *itr->second.second;

    // build accessors for orientation fields
    RCP<const OrientationContainerBase<Scalar,Array> > orientationContainer 
        = buildOrientationContainer<Scalar,Array>(globalIndexer_,fieldName); 

    IntrepidFieldContainerFactory fc_factory;
 
    // loop over worksets compute and apply orientations
    for(std::size_t i=0;i<worksets.size();i++) {
      for(std::size_t j=0;j<worksets[i].details.size();j++) {

        // break out of the workset loop
        if(worksets[i].num_cells<=0) continue;

        int array0_sz = worksets[i].num_cells;
        int array1_sz = basis.functional->dimension(1);
        Array orientations = fc_factory.buildArray<double,panzer::Cell,panzer::BASIS>("orientations",array0_sz,array1_sz);

        WorksetDetails & details = *worksets[i].details[j];

        // compute orientations using the orientation container (and global indexer eventually)
        orientationContainer->getOrientations(eBlock,details.cell_local_ids,orientations);

        for(std::size_t basis_index=0;basis_index<details.bases.size();basis_index++) {
          Teuchos::RCP<const BasisIRLayout> layout = details.bases[basis_index]->basis_layout;
          TEUCHOS_ASSERT(layout!=Teuchos::null);
          TEUCHOS_ASSERT(layout->getBasis()!=Teuchos::null);
          if(layout->getBasis()->name()==basis.name()) {
            // apply orientations for this basis
            details.bases[basis_index]->applyOrientations(orientations);
          }
        }
      }
    }
  }
}

void WorksetContainer::
applyOrientations(const SideId & sideId,std::map<unsigned,Workset> & worksets) const
{
  using Teuchos::RCP;

  typedef double Scalar;                          // orientation container scalar type
  typedef Intrepid::FieldContainer<Scalar> Array; // orientation container array type
  typedef std::pair<std::string,Teuchos::RCP<const PureBasis> > StrConstBasisPair;

  /////////////////////////////////
  // this is for side worksets //
  /////////////////////////////////

  // short circuit if no global indexer exists
  if(globalIndexer_==Teuchos::null) { 
    Teuchos::FancyOStream fout(Teuchos::rcpFromRef(std::cout));
    fout.setOutputToRootOnly(0);
 
    fout << "Panzer Warning: No global indexer assigned to a workset container. "
         << "Orientation of the basis for edge basis functions cannot be applied, "
         << "if those basis functions are used, there will be problems!";
    return;
  }

  // extract a map from basis name to string/basis pairs that require orientations
  //////////////////////////////////////////////////////////////////////////////////
  std::map<std::string, StrConstBasisPair> orientedBasisFields;
  {
    const PhysicsBlock & pb = lookupPhysicsBlock(sideId.eblk_id);

    // get all field/basis pairs
    std::vector<StrConstBasisPair> fieldBasisPairs;
    pb.getFieldLibraryBase()->basisPairs(fieldBasisPairs);
  
    // now use a map to extract only unique basis functions that
    // require orientations
    for(std::size_t i=0;i<fieldBasisPairs.size();i++) {
      StrConstBasisPair pair = fieldBasisPairs[i];
  
      if(pair.second->requiresOrientations()) {
        // name is guranteed to be unique for basis. 
        // just removes duplicate basis objects
        orientedBasisFields[pair.second->name()] = pair;
      }
    }
  }

  // loop over each basis requiring orientations, then apply them
  //////////////////////////////////////////////////////////////////////////////////

  // Note: It may be faster to loop over the basis pairs on the inside (not really sure)
  
  for(std::map<std::string, StrConstBasisPair>::const_iterator itr=orientedBasisFields.begin();
      itr!=orientedBasisFields.end();++itr) {
    // extract field name and basis name 
    const std::string & fieldName = itr->second.first;
    const PureBasis & basis = *itr->second.second;

    // build accessors for orientation fields
    RCP<const OrientationContainerBase<Scalar,Array> > orientationContainer 
        = buildOrientationContainer<Scalar,Array>(globalIndexer_,fieldName); 
 
    // loop over worksets compute and apply orientations
    for(std::map<unsigned,Workset>::iterator itr=worksets.begin();
        itr!=worksets.end();++itr) {

      // break out of the workset loop
      if(itr->second.num_cells<=0) continue;

      int array0_sz = itr->second.num_cells;
      int array1_sz = basis.functional->dimension(1);
      IntrepidFieldContainerFactory fc_factory;
      Array orientations = fc_factory.buildArray<double,panzer::Cell,panzer::BASIS>("orientations",array0_sz,array1_sz);

      for(std::size_t j=0;j<itr->second.details.size();j++) {

        WorksetDetails & details = *itr->second.details[j];

        // compute orientations using the orientation container (and global indexer eventually)
        orientationContainer->getOrientations(sideId.eblk_id,details.cell_local_ids,orientations);

        for(std::size_t basis_index=0;basis_index<details.bases.size();basis_index++) {
          Teuchos::RCP<const BasisIRLayout> layout = details.bases[basis_index]->basis_layout;
          TEUCHOS_ASSERT(layout!=Teuchos::null);
          TEUCHOS_ASSERT(layout->getBasis()!=Teuchos::null);
          if(layout->getBasis()->name()==basis.name()) {
            // apply orientations for this basis
            details.bases[basis_index]->applyOrientations(orientations);
          }
        }
      }
    }
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
