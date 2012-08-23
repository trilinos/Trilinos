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

#ifndef __Panzer_SingleBlockDOFManager_impl_hpp__
#define __Panzer_SingleBlockDOFManager_impl_hpp__

#include <map>

#include "mpi.h"

#include "Panzer_config.hpp"
#include "Panzer_FieldPattern.hpp"
#include "Panzer_FieldAggPattern.hpp"
#include "Panzer_GeometricAggFieldPattern.hpp"
#include "Panzer_ConnManager.hpp"
#include "Panzer_UniqueGlobalIndexer.hpp"
#include "Panzer_UniqueGlobalIndexer_Utilities.hpp"

#include "Panzer_SingleBlockDOFManager_decl.hpp"

#include "Teuchos_RCP.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_ArrayView.hpp"

#include "Tpetra_DefaultPlatform.hpp"
#include "Tpetra_Map.hpp"

namespace panzer {

using Teuchos::RCP;

template <typename LO, typename GO>
SingleBlockDOFManager<LO,GO>::SingleBlockDOFManager()
  : numFields_(0),buildConnectivityRun_(false),requireOrientations_(false)
{ }

template <typename LO, typename GO>
void SingleBlockDOFManager<LO,GO>::setConnManager(const Teuchos::RCP<ConnManager<LO,GO> > & connMngr, MPI_Comm mpiComm)
{
  TEUCHOS_TEST_FOR_EXCEPTION(buildConnectivityRun_,std::logic_error,
                      "SingleBlockDOFManager::setConnManager: setConnManager cannot be called after "
                      "buildGlobalUnknowns has been called"); 
  connMngr_ = connMngr;
  communicator_ = Teuchos::rcp(new Teuchos::MpiComm<int>(Teuchos::opaqueWrapper(mpiComm)));
}


  //Adds a field to be used in creating the Global Numbering
  //Returns the index for the field pattern
template <typename LO, typename GO>
int SingleBlockDOFManager<LO,GO>::addField(const std::string & str, const Teuchos::RCP<const FieldPattern> & pattern)
{

  TEUCHOS_TEST_FOR_EXCEPTION(buildConnectivityRun_,std::logic_error,
                      "SingleBlockDOFManager::addField: addField cannot be called after "
                      "buildGlobalUnknowns has been called"); 
  if(fieldPatterns_.size()>0)
    TEUCHOS_TEST_FOR_EXCEPTION(!(pattern->equals(*fieldPatterns_[0])),std::logic_error,"Patterns must be equal!");
  fieldPatterns_.push_back(pattern);
  fieldStringToIndex_.insert(std::map<std::string,int>::value_type(str, numFields_));
  //The default values for IDs are the sequential order they are added in.
  fieldStringToID_.insert(std::map<std::string,int>::value_type(str,numFields_));
  fieldOrder_.push_back(str);
  ++numFields_;
  return numFields_-1;
}


  //TO-WRITE: A method that differentiates between elementblocks
  //TO-WRITE: Field Ordering Method...I'm not sure why this is that important.

  //Returns the fieldpattern of the given name
  //This could also be done using the number you'd get from getFieldNum which
  //isn't yet included.
template <typename LO, typename GO>
Teuchos::RCP<const FieldPattern> SingleBlockDOFManager<LO,GO>::getFieldPattern(const std::string & name)
{
  return fieldPatterns_[fieldStringToIndex_[name]];
}

template <typename LO, typename GO>
void SingleBlockDOFManager<LO,GO>::getOwnedIndices(std::vector<GO> & indicies) const{
  indicies.resize(owned_.size());
  for (size_t i = 0; i < owned_.size(); ++i) {
    indicies[i]=owned_[i];
  }
}

template <typename LO, typename GO>
void SingleBlockDOFManager<LO,GO>::getOwnedAndSharedIndices(std::vector<GO> & indicies) const{
  indicies.resize(owned_and_ghosted_.size());
  for (size_t i = 0; i < owned_and_ghosted_.size(); ++i) {
    indicies[i]=owned_and_ghosted_[i];
  }
}
  
  //gets the number of fields
template <typename LO, typename GO>
int SingleBlockDOFManager<LO,GO>::getNumFields() const
{
  return numFields_;
}

template <typename LO, typename GO>
const std::vector<int> & SingleBlockDOFManager<LO,GO>::getGIDFieldOffsets(const std::string & blockID, int fieldNum) const{
  return fa_fp->localOffsets(fieldNum);
}

template <typename LO, typename GO>
void SingleBlockDOFManager<LO,GO>::getElementGIDs(LO localElementID, std::vector<GO> & gids, const std::string & blockIdHint) const {
  std::vector<GO> which = elementGIDs_[localElementID];
  gids.resize(which.size());
  for (size_t i = 0; i < which.size(); ++i) {
    gids[i]=which[i];
  }
}

//builds the global unknowns array
template <typename LO, typename GO>
void SingleBlockDOFManager<LO,GO>::buildGlobalUnknowns()
{
  TEUCHOS_TEST_FOR_EXCEPTION(buildConnectivityRun_,std::logic_error,
                      "SingleBlockDOFManager::buildGlobalUnknowns: buildGlobalUnknowns cannot be called again "
                      "after buildGlobalUnknowns has been called"); 
  //Some stuff for the Map dealios.
  typedef Tpetra::DefaultPlatform::DefaultPlatformType Platform;
  typedef Tpetra::DefaultPlatform::DefaultPlatformType::NodeType Node;
  typedef Tpetra::Map<LO, GO, Node> Map;
  //I'm not sure if this should go here, or if this will work,
  //but I'm going to need this information to make the maps.
  Tpetra::DefaultPlatform::DefaultPlatformType &platform = Tpetra::DefaultPlatform::getDefaultPlatform();
  RCP<const Teuchos::Comm<int> > comm = platform.getComm();
  RCP<Node> node = platform.getNode();

  std::vector<RCP<const panzer::FieldPattern> > field_patterns;
  //I needed this before I changed how I'm storring things...
  std::vector< std::pair< int, RCP<const panzer::FieldPattern > > > paired_patterns(numFields_);


//  int i=0;
//  for (std::map<std::string,int>::const_iterator itr=fieldStringToIndex_.begin(); itr!=fieldStringToIndex_.end(); ++i,itr++) {
//    int ind=itr->find(fieldOrder_[i])->second;
//    field_patterns.push_back(fieldPatterns_[ind]);
//    paired_patterns[i]=std::make_pair(i,fieldPatterns_[ind]);
//  }
  for(int i=0;i<numFields_;++i){
    int ind=fieldStringToIndex_.find(fieldOrder_[i])->second;
    field_patterns.push_back(fieldPatterns_[ind]);
    paired_patterns[i]=std::make_pair(i,fieldPatterns_[ind]);
    //fieldStringToID_.insert(std::map<std::string,int>::value_type(fieldOrder_[i],i));
    //It has been correctly set above, we just need to reassign in case there has been a new.
    //fieldStringToID_[fieldOrder_[i]]=i;
    //None of this is necessary. This map is updated everytime field order is set.
  }
  ga_fp = Teuchos::rcp(new GeometricAggFieldPattern(field_patterns));
  fa_fp = Teuchos::rcp(new FieldAggPattern(paired_patterns));
  connMngr_->buildConnectivity(*ga_fp);

  std::set<GO> concerned_set;
  
  const std::vector<LO> & myElements=connMngr_->getElementBlock("eblock-0_0");
  for (std::size_t e = 0; e<myElements.size(); ++e) {
    LO connSize=connMngr_->getConnectivitySize(myElements[e]);
    const GO* elementConn = connMngr_->getConnectivity(myElements[e]);
    for (LO i = 0; i < connSize; ++i) {
      concerned_set.insert(elementConn[i]);
    }
  }

  for(typename std::set<GO>::const_iterator itr=concerned_set.begin(); itr!=concerned_set.end(); ++itr){
    owned_and_ghosted_.push_back(*itr);
  }


  RCP<const Map> all_nodes = Tpetra::createNonContigMapWithNode<LO,GO>(owned_and_ghosted_,comm,node);
  RCP<const Map> assigned_nodes = Tpetra::createOneToOne<LO,GO,Node>(all_nodes);
  
  owned_=createVector(assigned_nodes->getNodeElementList());

  //Here we generate all of the content for elementGIDs_.

  //Using myElements from above.
  //I make an assumption here that the Local IDs given by myElements
  //will always be contigious.
  for (std::size_t i = 0; i < myElements.size(); ++i) {
    LO connSize=connMngr_->getConnectivitySize(myElements[i]);
    const GO* elementConn = connMngr_->getConnectivity(myElements[i]);
    std::vector<GO> gids(connSize*numFields_);
    for (LO j = 0; j < connSize; ++j) {
      for (int k = 0; k < numFields_; ++k) {
        //this needs to be edited to take into account field ordering.
        int location=(j*numFields_)+k;
        int value=(elementConn[j]*numFields_)+k;
        gids[location]=value;
      }
    }
    //gids is now filled and ready to go into  elementGIDs_.
    //Again: because we are going through contigiously numbered
    //local elements in order, I assume pushing and then
    //referencing by order will be fine.
    elementGIDs_.push_back(gids);
  }

  buildConnectivityRun_=true;
  if(requireOrientations_)
    buildUnknownsOrientation();
}


template <typename LO, typename GO>
int SingleBlockDOFManager<LO,GO>::getFieldNum(const std::string & string) const{
  std::map<std::string,int>::const_iterator itr = fieldStringToID_.find(string);
  if(itr==fieldStringToID_.end()){
    TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,"INVALID:FIELD NAME"); 
  }
  //NOTE THIS DOESN"T ERROR CHECK AT ALL! BADD!
  return itr->second;
}

  
template <typename LO, typename GO>
bool SingleBlockDOFManager<LO,GO>::fieldInBlock(const std::string & field, const std::string & block) const{
  TEUCHOS_TEST_FOR_EXCEPTION(block!="eblock-0_0",std::logic_error,"SingleBlockDOFManager::fieldInBlock: Invalid block name.");
  return fieldStringToIndex_.find(field) != fieldStringToIndex_.end();
}

template <typename LO, typename GO>
const std::vector<int> & SingleBlockDOFManager<LO,GO>::getBlockFieldNumbers(const std::string & blockId) const{
  TEUCHOS_TEST_FOR_EXCEPTION(blockId!="eblock-0_0",std::logic_error,"SingleBlockDOFManager::getBlockFieldNumbers: Invalid block name.");
  //TODO:CHECK BLOCK NAME.
  //{ return fieldAggPattern_.find(block)->second->fieldIds(); }
//  std::vector<int> to_return;
//  for(std::map<std::string,int>::const_iterator itr=fieldStringToIndex_.begin(); itr!=fieldStringToIndex_.end();itr++){
//    to_return.push_back(itr->second);
//  }

  return fa_fp->fieldIds();

}

template <typename LO, typename GO>
void SingleBlockDOFManager<LO,GO>::ownedIndices(const std::vector<GO> & indices,std::vector<bool> & isOwned) const{
  //Resizes the isOwned array.
  if(indices.size()!=isOwned.size())
    isOwned.resize(indices.size(),false);
  typename std::vector<GO>::const_iterator endOf = owned_.end();
  for (std::size_t i = 0; i < indices.size(); ++i) {
    isOwned[i] = ( std::find(owned_.begin(), owned_.end(), indices[i])!=endOf );
  }
}

template <typename LO, typename GO>
void SingleBlockDOFManager<LO,GO>::setFieldOrder(const std::vector<std::string> & fieldOrder){
  TEUCHOS_TEST_FOR_EXCEPTION(buildConnectivityRun_,std::logic_error,
                      "SingleBlockDOFManager::setFieldOrder: setFieldOrder cannot be called after "
                      "buildGlobalUnknowns has been called"); 
  if(validFieldOrder(fieldOrder)){
    fieldOrder_=fieldOrder;
    //We also need to reassign the IDs to this ordering.
    for (size_t i = 0; i < fieldOrder_.size(); ++i) {
      fieldStringToID_[fieldOrder_[i]]=i;
    }
  }
  else
    TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,"Invalid Field Ordering!");
    
}


template <typename LO, typename GO>
bool SingleBlockDOFManager<LO,GO>::validFieldOrder(const std::vector<std::string> & proposed_fieldOrder){
  if(fieldOrder_.size()!=proposed_fieldOrder.size())
    return false;
  //I'm using a not very efficient way of doing this, but there should never
  //be that many fields, so it isn't that big of a deal.
  for (size_t i = 0; i < fieldOrder_.size(); ++i) {
    bool found=false;
    for (size_t j = 0; j < proposed_fieldOrder.size(); ++j) {
      if(fieldOrder_[i]==proposed_fieldOrder[j])
        found=true;
    }
    if(!found)
      return false;
  }
  return true;
}

template <typename LO, typename GO>
void SingleBlockDOFManager<LO,GO>::buildUnknownsOrientation(){
  orientation_.clear();
  //unlike the dof manager, this simple SingleBlockDOFManager only has one block.
  std::size_t myElementCount = connMngr_->getElementBlock("eblock-0_0").size();

  orientation_.resize(myElementCount);
  
  //Here the DOF code has a search for a certain fieldAgg pattern,
  //but I only have one, so it shouldn't matter.
  std::vector<std::pair<int,int> > topEdgeIndices;
  orientation_helpers::computePatternEdgeIndices(*ga_fp,topEdgeIndices);
  //The element GID count for a specific block is just the GID count
  //for the proc. Remember, one element block.
  //Also, this will always be run after buildGlobalUnknowns has been run.
  std::size_t numGIDs = owned_.size();
  const std::vector<LO> & elmts = getElementBlock("eblock-0_0");
  for (std::size_t e = 0; e < elmts.size(); ++e) {
    std::vector<char> & eOrientation = orientation_[elmts[e]];
    eOrientation.resize(numGIDs);
    for (std::size_t s = 0; s < eOrientation.size(); ++s) {
      eOrientation[s]=1;//1 by default
    }
    LO connSz = connMngr_->getConnectivitySize(elmts[e]);
    const GO * connPtr = connMngr_->getConnectivity(elmts[e]);
    const std::vector<GO> connectivity(connPtr, connPtr+connSz);

    orientation_helpers::computeCellEdgeOrientations(topEdgeIndices,
                            connectivity, *fa_fp, eOrientation);

  }


  
}

template <typename LO, typename GO>
const std::string & SingleBlockDOFManager<LO,GO>::getFieldString(int num) const{
  //I need the string associated with this int in fieldstringtoid.
  std::string toRet;
  for(std::map<std::string,int>::const_iterator itr=fieldStringToID_.begin(); itr!=fieldStringToID_.end();itr++){
    if(itr->second==num)
      return itr->first;
  }
  TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,"SingleBlockDOFManager::getFieldString: invalid number.");
}


} /*panzer*/

#endif
