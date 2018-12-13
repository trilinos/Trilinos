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

#include "Panzer_IntrepidOrientation.hpp"

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
                                   const std::map<std::string,WorksetNeeds> & needs)
   : wkstFactory_(factory), worksetSize_(-1)
{
  // thats all!
  ebToNeeds_ = needs;
}

/** Copies the workset factory and the workset size, but not constructed
  * worksets.
  */
WorksetContainer::WorksetContainer(const WorksetContainer & wc)
   : wkstFactory_(wc.wkstFactory_)
   , worksetSize_(wc.worksetSize_)
{
}

/** Clear all allocated worksets, maintain the workset factory and element to physics
  * block map.
  */
void WorksetContainer::clear()
{
   worksets_.clear();
   sideWorksets_.clear();
}

void WorksetContainer::
setNeeds(const std::string & eBlock,const WorksetNeeds & needs)
{
  clear(); // clear out old worksets
  ebToNeeds_[eBlock] = needs;
}

//! Look up an input physics block, throws an exception if it can be found.
const WorksetNeeds & WorksetContainer::lookupNeeds(const std::string & eBlock) const
{
   std::map<std::string,WorksetNeeds>::const_iterator itr = ebToNeeds_.find(eBlock);

   TEUCHOS_TEST_FOR_EXCEPTION(itr==ebToNeeds_.end(),std::logic_error,
                      "WorksetContainer::lookupNeeds no WorksetNeeds object is associated "
                      "with the element block \""+eBlock+"\".");

   return itr->second;
}

Teuchos::RCP<std::vector<Workset> >
WorksetContainer::getWorksets(const WorksetDescriptor & wd)
{
   Teuchos::RCP<std::vector<Workset> > worksetVector;
   WorksetMap::iterator itr = worksets_.find(wd);
   if(itr==worksets_.end()) {
      // couldn't find workset, build it!
      const WorksetNeeds & needs = lookupNeeds(wd.getElementBlock());
      worksetVector = wkstFactory_->getWorksets(wd,needs);

      // apply orientations to the just constructed worksets
      if(worksetVector!=Teuchos::null && wd.applyOrientations()) {
        applyOrientations(wd.getElementBlock(),*worksetVector);
      }

      if(worksetVector!=Teuchos::null)
        setIdentifiers(wd,*worksetVector);

      // store vector for reuse in the future
      worksets_[wd] = worksetVector;
   }
   else
      worksetVector = itr->second;

   return worksetVector;
}

//! Access, and construction of side worksets
Teuchos::RCP<std::map<unsigned,Workset> >
WorksetContainer::getSideWorksets(const WorksetDescriptor & desc)
{
   Teuchos::RCP<std::map<unsigned,Workset> > worksetMap;

   // this is the key for the workset map
   SideMap::iterator itr = sideWorksets_.find(desc);

   if(itr==sideWorksets_.end()) {
      // couldn't find workset, build it!
      if (desc.connectsElementBlocks()) {
        worksetMap = wkstFactory_->getSideWorksets(desc, lookupNeeds(desc.getElementBlock(0)),
                                                         lookupNeeds(desc.getElementBlock(1)));
      }
      else {
        worksetMap = wkstFactory_->getSideWorksets(desc,lookupNeeds(desc.getElementBlock(0)));
      }

      // apply orientations to the worksets for this side
      if(worksetMap!=Teuchos::null)
        applyOrientations(desc,*worksetMap);

      if(worksetMap!=Teuchos::null)
        setIdentifiers(desc,*worksetMap);

      // store map for reuse in the future
      sideWorksets_[desc] = worksetMap;
   }
   else {
      worksetMap = itr->second;
   }

   return worksetMap;
}


void WorksetContainer::
setGlobalIndexer(const Teuchos::RCP<const panzer::UniqueGlobalIndexerBase> & ugi)
{
  // apply the orientations for stored worksets
  applyOrientations(ugi);
}

void WorksetContainer::
addBasis(const std::string & type,int order,const std::string & rep_field)
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  for(auto itr=ebToNeeds_.begin();itr!=ebToNeeds_.end();++itr) {
    WorksetNeeds & needs = itr->second;
    RCP<PureBasis> basis = rcp(new PureBasis(type,order,needs.cellData));

    // add in the new basis
    needs.bases.push_back(basis);
    needs.rep_field_name.push_back(rep_field);
  }

  // clear all arrays, lazy evaluation means it will be rebuilt
  clear();
}

void WorksetContainer::
applyOrientations(const Teuchos::RCP<const panzer::UniqueGlobalIndexerBase> & ugi)
{
  // this gurantees orientations won't accidently be applied twice.
  TEUCHOS_ASSERT(globalIndexer_==Teuchos::null);

  globalIndexer_ = ugi;

  // this should be created once and stored in an appropriate place
  TEUCHOS_TEST_FOR_EXCEPTION(globalIndexer_ == Teuchos::null, std::logic_error,
                             "global indexer is not set yet");
  orientations_ = buildIntrepidOrientation(globalIndexer_);

  // loop over volume worksets, apply orientations to each
  for(WorksetMap::iterator itr=worksets_.begin();
      itr!=worksets_.end();++itr) {
    std::string eBlock = itr->first.getElementBlock();

    applyOrientations(eBlock,*itr->second);
  }

  // loop over side worksets, apply orientations to each
  for(SideMap::iterator itr=sideWorksets_.begin();
      itr!=sideWorksets_.end();itr++) {

    applyOrientations(itr->first,*itr->second);
  }
}

void WorksetContainer::
setIdentifiers(const WorksetDescriptor & wd,std::vector<Workset> & worksets)
{
  std::size_t hash = std::hash<WorksetDescriptor>()(wd); // this is really ugly, is this really a C++ standard?
  for(std::size_t i=0;i<worksets.size();i++)
    worksets[i].setIdentifier(hash+i);
}

void WorksetContainer::
setIdentifiers(const WorksetDescriptor & wd,std::map<unsigned,Workset> & workset_map)
{
  std::size_t hash = std::hash<WorksetDescriptor>()(wd); // this is really ugly, is this really a C++ standard?
  std::size_t offset = 0;
  for(auto itr : workset_map) {
    // itr.second.setIdentifier(hash+offset);
    workset_map[itr.first].setIdentifier(hash+offset);

    offset++;
  }
}

void WorksetContainer::
applyOrientations(const std::string & eBlock, std::vector<Workset> & worksets) const
{
  using Teuchos::RCP;

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

  // this should be matched to global indexer size (not sure how to retrive it)
  TEUCHOS_TEST_FOR_EXCEPTION(orientations_ == Teuchos::null, std::logic_error,
                             "intrepid2 orientation is not constructed");

  // loop over each basis requiring orientations, then apply them
  //////////////////////////////////////////////////////////////////////////////////
  
  // Note: It may be faster to loop over the basis pairs on the inside (not really sure)

  const WorksetNeeds & needs = lookupNeeds(eBlock);

  if(needs.bases.size()>0) {
    // sanity check that we aren't missing something (the old and new "needs" should not be used together)
    TEUCHOS_ASSERT(needs.getBases().size()==0);

    for(std::size_t w=0;w<needs.bases.size();w++) {
      const PureBasis & basis = *needs.bases[w];
  
      // no need for this if orientations are not required!
      if(!basis.requiresOrientations())
        continue;
  
      // build accessors for orientation fields
      std::vector<Intrepid2::Orientation> ortsPerBlock;
  
      // loop over worksets compute and apply orientations
      for(std::size_t i=0;i<worksets.size();i++) {
        // break out of the workset loop
        if(worksets[i].num_cells<=0) continue;
  
        for(std::size_t j=0;j<worksets[i].numDetails();j++) {
          WorksetDetails & details = worksets[i](j);
  
          ortsPerBlock.clear();
          for (int k=0;k<worksets[i].num_cells;++k) {
            ortsPerBlock.push_back((*orientations_)[details.cell_local_ids[k]]);
          }
  
          for(std::size_t basis_index=0;basis_index<details.bases.size();basis_index++) {
            Teuchos::RCP<const BasisIRLayout> layout = details.bases[basis_index]->basis_layout;
  
            // only apply orientations if its relevant to the current needs
            if(layout->getBasis()->name()!=basis.name())
              continue;
  
            TEUCHOS_ASSERT(layout!=Teuchos::null);
            TEUCHOS_ASSERT(layout->getBasis()!=Teuchos::null);
            if(layout->getBasis()->requiresOrientations()) {
              // apply orientations for this basis
              details.bases[basis_index]->applyOrientations(ortsPerBlock,(int) worksets[i].num_cells);
            }
          }
        }
      }
    } // end for w
  }
  else if(needs.getBases().size()>0) {
    // sanity check that we aren't missing something (the old and new "needs" should not be used together)
    TEUCHOS_ASSERT(needs.bases.size()==0);

    // This is for forwards compatibility, the needs now use "getBasis" calls as opposed
    // to director accessors.
    for(const auto & bd : needs.getBases()) {
  
      // build accessors for orientation fields
      std::vector<Intrepid2::Orientation> ortsPerBlock;
  
      // loop over worksets compute and apply orientations
      for(std::size_t i=0;i<worksets.size();i++) {
        // break out of the workset loop
        if(worksets[i].num_cells<=0) continue;
  
        for(std::size_t j=0;j<worksets[i].numDetails();j++) {
          WorksetDetails & details = worksets[i](j);
  
          ortsPerBlock.clear();
          // for (int k=0;k<worksets[i].num_cells;++k) {
          for (int k=0;k<details.numOwnedCells();++k) {
            ortsPerBlock.push_back((*orientations_)[details.cell_local_ids[k]]);
          }
  
          for(const auto & id : needs.getIntegrators()) {
            // apply orientations for this basis
            details.getBasisValues(bd,id).applyOrientations(ortsPerBlock,(int) worksets[i].num_cells);
          }
        }
      }
    } // end for w
  }
}
  
  
void WorksetContainer::
applyOrientations(const WorksetDescriptor & desc,std::map<unsigned,Workset> & worksets) const
{
  using Teuchos::RCP;
  
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
  
  // loop over each basis requiring orientations, then apply them
  //////////////////////////////////////////////////////////////////////////////////
  
  // Note: It may be faster to loop over the basis pairs on the inside (not really sure)
  const WorksetNeeds & needs = lookupNeeds(desc.getElementBlock());

  if(needs.bases.size()>0) {
    // sanity check that we aren't missing something (the old and new "needs" should not be used together)
    TEUCHOS_ASSERT(needs.getBases().size()==0);
    for(std::size_t i=0;i<needs.bases.size();i++) {
      const PureBasis & basis = *needs.bases[i];
      
      // no need for this if orientations are not required!
      if(!basis.requiresOrientations()) continue;
      
      // build accessors for orientation fields
      std::vector<Intrepid2::Orientation> ortsPerBlock;  
      
      // loop over worksets compute and apply orientations
      for(std::map<unsigned,Workset>::iterator itr=worksets.begin(); 
          itr!=worksets.end();++itr) { 
        
        // break out of the workset loop
        if(itr->second.num_cells<=0) continue;
        
        for(std::size_t j=0;j<itr->second.numDetails();j++) {  
          WorksetDetails & details = itr->second(j);
  
          ortsPerBlock.clear();
          for (int k=0;k<itr->second.num_cells;++k) {
            ortsPerBlock.push_back((*orientations_)[details.cell_local_ids[k]]);
          }
  
          for(std::size_t basis_index=0;basis_index<details.bases.size();basis_index++) {
            Teuchos::RCP<const BasisIRLayout> layout = details.bases[basis_index]->basis_layout;
  
            // only apply orientations if its relevant to the current needs
            if(layout->getBasis()->name()!=basis.name())
              continue;
  
            TEUCHOS_ASSERT(layout!=Teuchos::null);
            TEUCHOS_ASSERT(layout->getBasis()!=Teuchos::null);
            if(layout->getBasis()->requiresOrientations()) {
              // apply orientations for this basis
              details.bases[basis_index]->applyOrientations(ortsPerBlock,(int) itr->second.num_cells);
            }
          }
        }
      }
    } // end for i
  }
  else if(needs.getBases().size()>0) {
    // sanity check that we aren't missing something (the old and new "needs" should not be used together)
    TEUCHOS_ASSERT(needs.bases.size()==0);

    // This is for forwards compatibility, the needs now use "getBasis" calls as opposed
    // to director accessors.
    for(const auto & bd : needs.getBases()) {
  
      // build accessors for orientation fields
      std::vector<Intrepid2::Orientation> ortsPerBlock;
  
      // loop over worksets compute and apply orientations
      for(std::map<unsigned,Workset>::iterator itr=worksets.begin(); 
          itr!=worksets.end();++itr) { 
        
        // break out of the workset loop
        if(itr->second.num_cells<=0) continue;
        
        for(std::size_t j=0;j<itr->second.numDetails();j++) {  
          WorksetDetails & details = itr->second(j);
  
          ortsPerBlock.clear();
          for (int k=0;k<itr->second.num_cells;++k) {
            ortsPerBlock.push_back((*orientations_)[details.cell_local_ids[k]]);
          }
  
          for(const auto & id : needs.getIntegrators()) {
            // apply orientations for this basis
            details.getBasisValues(bd,id).applyOrientations(ortsPerBlock,(int) itr->second.num_cells);
          }
        }
      }
    } // end for w
  }
}

void getVolumeWorksetsFromContainer(WorksetContainer & wc,
                                    const std::vector<std::string> & elementBlockNames,
                                    std::map<std::string,Teuchos::RCP<std::vector<Workset> > > & volumeWksts)
{
   for(std::size_t i=0;i<elementBlockNames.size();i++) {
      WorksetDescriptor wd = blockDescriptor(elementBlockNames[i]);
      volumeWksts[elementBlockNames[i]] = wc.getWorksets(wd);
   }
}

void getSideWorksetsFromContainer(WorksetContainer & wc,
                                  const std::vector<BC> & bcs,
                                  std::map<BC,Teuchos::RCP<std::map<unsigned,Workset> >,LessBC> & sideWksts)
{
   for(std::size_t i=0;i<bcs.size();i++) {
      WorksetDescriptor wd(bcs[i].elementBlockID(),bcs[i].sidesetID());
      Teuchos::RCP<std::map<unsigned,Workset> > wksts = wc.getSideWorksets(wd);
      if(wksts!=Teuchos::null)
         sideWksts[bcs[i]] = wksts;
   }
}

}
