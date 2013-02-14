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

#ifndef PANZER_BLOCKED_DOF_MANAGER_IMPL_HPP
#define PANZER_BLOCKED_DOF_MANAGER_IMPL_HPP

#include <map>

#include "Panzer_GeometricAggFieldPattern.hpp"

#include "Teuchos_DefaultMpiComm.hpp"

using Teuchos::RCP;

namespace panzer {

// ************************************************************
// class BlockedDOFManager
// ************************************************************

template <typename LocalOrdinalT,typename GlobalOrdinalT>
BlockedDOFManager<LocalOrdinalT,GlobalOrdinalT>::BlockedDOFManager()
   : fieldsRegistered_(false), maxSubFieldNum_(-1), requireOrientations_(false), useDOFManagerFEI_(true)
{ }

template <typename LocalOrdinalT,typename GlobalOrdinalT>
BlockedDOFManager<LocalOrdinalT,GlobalOrdinalT>::BlockedDOFManager(const Teuchos::RCP<ConnManager<LocalOrdinalT,GlobalOrdinalT> > & connMngr,MPI_Comm mpiComm)
   : fieldsRegistered_(false), maxSubFieldNum_(-1), requireOrientations_(false), useDOFManagerFEI_(true)
{
   setConnManager(connMngr,mpiComm);
}

////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////

template <typename LocalOrdinalT,typename GlobalOrdinalT>
int BlockedDOFManager<LocalOrdinalT,GlobalOrdinalT>::getFieldNum(const std::string & str) const
{
   std::map<std::string,int>::const_iterator itr = fieldStrToNum_.find(str);

   // return based on what was found
   if(itr==fieldStrToNum_.end()) {
      // incorrect field name
      TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,
                         "BlockedDOFManager::getFieldNum No field with the name \"" + str + "\" has been added");
   }
   else {
      return itr->second;
   }
}

template <typename LocalOrdinalT,typename GlobalOrdinalT>
const std::string & BlockedDOFManager<LocalOrdinalT,GlobalOrdinalT>::getFieldString(int number) const
{
   std::map<int,std::string>::const_iterator itr = fieldNumToStr_.find(number);

   // return based on what was found
   if(itr==fieldNumToStr_.end()) {
      std::stringstream ss; ss << number; // itoa() in c-language
      // incorrect field name
      TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,
                         "BlockedDOFManager::getFieldString No field with number \"" + ss.str() + "\" has been added");
   }
   else {
      return itr->second;
   }
}

template <typename LocalOrdinalT,typename GlobalOrdinalT>
bool BlockedDOFManager<LocalOrdinalT,GlobalOrdinalT>::fieldInBlock(const std::string & field,const std::string & block) const
{
   // try  to find element block
   std::map<std::string,std::set<std::string> >::const_iterator fieldsItr = blockIdToFieldStrings_.find(block);
   TEUCHOS_TEST_FOR_EXCEPTION(fieldsItr==blockIdToFieldStrings_.end(),std::logic_error,
                      "BlockedDOFManager::fieldInBlock could not find the element block \""+block+"\"");

   // find field in element block 
   const std::set<std::string> & fields = fieldsItr->second;
   std::set<std::string>::const_iterator itr = fields.find(field);
   return itr!=fields.end(); 
}

/** Get field numbers associated with a particular element block.
  */
template <typename LocalOrdinalT,typename GlobalOrdinalT>
const std::vector<int> & BlockedDOFManager<LocalOrdinalT,GlobalOrdinalT>::getBlockFieldNumbers(const std::string & block) const
{
   // try to find element block
   std::map<std::string,std::vector<int> >::const_iterator fieldsItr = blockIdToFieldNumbers_.find(block);
   TEUCHOS_TEST_FOR_EXCEPTION(fieldsItr==blockIdToFieldNumbers_.end(),std::logic_error,
                      "BlockedDOFManager::getBlockFieldNumbers cannot field elemenet block, has registerFields() been called?");

   return fieldsItr->second;
}

template <typename LocalOrdinalT,typename GlobalOrdinalT>
void BlockedDOFManager<LocalOrdinalT,GlobalOrdinalT>::getElementGIDs(LocalOrdinalT localElmtId,std::vector<GlobalOrdinal> & gids,const std::string & blockIdHint) const
{
   // WARNING: there is an assumed ordering being used here it
   // corresponds directly to the blockGIDOffset_ map and (as
   // a result) the getBlockGIDOffset function. However for 
   // the sake of speed this conversion is implicit. 
   //  
   // Any changes to the order should be reflected in the
   // blockGIDOffset_ map.

   gids.resize(0);

   // loop over field block manager and grab indices
   for(std::size_t fbm=0;fbm<fieldBlockManagers_.size();fbm++) {
      std::vector<GlobalOrdinalT> fieldBlockOwned;

      fieldBlockManagers_[fbm]->getElementGIDs(localElmtId,fieldBlockOwned,blockIdHint);

      for(std::size_t i=0;i<fieldBlockOwned.size();i++) 
         gids.push_back(std::make_pair(fbm,fieldBlockOwned[i]));
   }
}

template <typename LocalOrdinalT,typename GlobalOrdinalT>
void BlockedDOFManager<LocalOrdinalT,GlobalOrdinalT>::getElementOrientation(LocalOrdinalT localElmtId,std::vector<double> & gidsOrientation) const
{
   // WARNING: there is an assumed ordering being used here it
   // corresponds directly to the blockGIDOffset_ map and (as
   // a result) the getBlockGIDOffset function. However for 
   // the sake of speed this conversion is implicit. 
   //  
   // Any changes to the order should be reflected in the
   // blockGIDOffset_ map.

   gidsOrientation.resize(0);

   // loop over field block manager and grab indices
   for(std::size_t fbm=0;fbm<fieldBlockManagers_.size();fbm++) {
      std::vector<double> blkOrientation;

      fieldBlockManagers_[fbm]->getElementOrientation(localElmtId,blkOrientation);

      for(std::size_t i=0;i<blkOrientation.size();i++) 
         gidsOrientation.push_back(blkOrientation[i]);
   }
}

template <typename LocalOrdinalT,typename GlobalOrdinalT>
const std::vector<int> & BlockedDOFManager<LocalOrdinalT,GlobalOrdinalT>::getGIDFieldOffsets(const std::string & blockId,int fieldNum) const
{
   typedef std::map<std::string,std::map<int,std::vector<int> > > FieldOffsetsMap;
  
   FieldOffsetsMap::iterator blockItr = gidFieldOffsets_.find(blockId);
   if(blockItr!=gidFieldOffsets_.end()) {
      std::map<int,std::vector<int> > & fieldToVectorMap = blockItr->second;
      std::map<int,std::vector<int> >::const_iterator itr = fieldToVectorMap.find(fieldNum);

      // we have found the vector, return the precomputed one
      if(itr!=fieldToVectorMap.end())
         return itr->second;
   }
   else {
      std::vector<std::string> elementBlocks;
      getElementBlockIds(elementBlocks);
      TEUCHOS_TEST_FOR_EXCEPTION(std::find(elementBlocks.begin(),elementBlocks.end(),blockId)==elementBlocks.end(),std::logic_error,
                         "BlockedDOFManager::getGIDFieldOffsets: Block ID \""+blockId+"\" does not exist");

      gidFieldOffsets_[blockId] = std::map<int,std::vector<int> >();
      blockItr = gidFieldOffsets_.find(blockId);
   }

   // grab relevant map from iterator
   std::map<int,std::vector<int> > & fieldToVectorMap = blockItr->second;
  
   // we have not found the vector, now we need to build one
   ////////////////////////////////////////////////////////////////

   // first grab all pieces that are needed for extracting GIDs from sub system
   int fieldBlock = getFieldBlock(fieldNum);
   Teuchos::RCP<const UniqueGlobalIndexer<LocalOrdinalT,GlobalOrdinalT> > dofManager = fieldBlockManagers_[fieldBlock];

   // grab offsets for sub dof manager. Notice you must convert to field number used by sub manager!
   const std::vector<int> & subGIDOffsets 
         = dofManager->getGIDFieldOffsets(blockId,dofManager->getFieldNum(getFieldString(fieldNum)));

   // increment offsets to correspond with blocked system
   int gidOffset = getBlockGIDOffset(blockId,fieldBlock); 
   std::vector<int> & finalFieldOffsets = fieldToVectorMap[fieldNum];
   finalFieldOffsets.resize(subGIDOffsets.size());
   for(std::size_t i=0;i<finalFieldOffsets.size();i++)
      finalFieldOffsets[i] = gidOffset+subGIDOffsets[i];

   return finalFieldOffsets;
}

template <typename LocalOrdinalT,typename GlobalOrdinalT>
bool BlockedDOFManager<LocalOrdinalT,GlobalOrdinalT>::LessThan
::operator()(const Teuchos::Tuple<int,3> & a,const Teuchos::Tuple<int,3> & b) const 
{
   if(a[0] < b[0]) return true;
   if(a[0] > b[0]) return false;

   // a[0]==b[0]  
   if(a[1] < b[1]) return true;
   if(a[1] > b[1]) return false;

   // a[1]==b[1] && a[0]==b[0] 
   if(a[2] < b[2]) return true;
   if(a[2] > b[2]) return false;

   // a[2]==b[2] && a[1]==b[1] && a[0]==b[0]
   return false; // these are equal to, but not less than!
}

template <typename LocalOrdinalT,typename GlobalOrdinalT>
const std::pair<std::vector<int>,std::vector<int> > & 
BlockedDOFManager<LocalOrdinalT,GlobalOrdinalT>::getGIDFieldOffsets_closure(const std::string & blockId,int fieldNum,int subcellDim,int subcellId) const
{
   typename std::map<std::string,TupleToVectorPairMap>::iterator blockItr = gidFieldOffsets_closure_.find(blockId);
   if(blockItr!=gidFieldOffsets_closure_.end()) {
      TupleToVectorPairMap & fieldToTupleMap = blockItr->second;
      typename TupleToVectorPairMap::const_iterator itr =
            fieldToTupleMap.find(Teuchos::tuple(fieldNum,subcellDim,subcellId));

      // we have found the vector, return the precomputed one
      if(itr!=fieldToTupleMap.end())
         return itr->second;
   }
   else {
      std::vector<std::string> elementBlocks;
      getElementBlockIds(elementBlocks);
      TEUCHOS_TEST_FOR_EXCEPTION(std::find(elementBlocks.begin(),elementBlocks.end(),blockId)==elementBlocks.end(),std::logic_error,
                         "BlockedDOFManager::getGIDFieldOffsets: Block ID \""+blockId+"\" does not exist");

      gidFieldOffsets_closure_[blockId] = TupleToVectorPairMap();
      blockItr = gidFieldOffsets_closure_.find(blockId);
   }

   // grab relevant map from iterator
   TupleToVectorPairMap & fieldToTupleMap = blockItr->second;
  
   // we have not found the vector, now we need to build one
   ////////////////////////////////////////////////////////////////

   // first grab all pieces that are needed for extracting GIDs from sub system
   int fieldBlock = getFieldBlock(fieldNum);
   Teuchos::RCP<const UniqueGlobalIndexer<LocalOrdinalT,GlobalOrdinalT> > dofManager = fieldBlockManagers_[fieldBlock];

   // grab offsets for sub dof manager. Notice you must convert to field number used by sub manager!
   const std::pair<std::vector<int>,std::vector<int> > & subGIDOffsets_closure
         = dofManager->getGIDFieldOffsets_closure(blockId,dofManager->getFieldNum(getFieldString(fieldNum)),subcellDim,subcellId);

   // increment offsets to correspond with blocked system
   int gidOffset = getBlockGIDOffset(blockId,fieldBlock); 
   std::pair<std::vector<int>,std::vector<int> > & finalFieldOffsets = fieldToTupleMap[Teuchos::tuple(fieldNum,subcellDim,subcellId)];
   finalFieldOffsets.first.resize(subGIDOffsets_closure.first.size());
   finalFieldOffsets.second = subGIDOffsets_closure.second;
   for(std::size_t i=0;i<finalFieldOffsets.first.size();i++)
      finalFieldOffsets.first[i] = gidOffset+subGIDOffsets_closure.first[i];

   return finalFieldOffsets;
}

template <typename LocalOrdinalT,typename GlobalOrdinalT>
void BlockedDOFManager<LocalOrdinalT,GlobalOrdinalT>::getOwnedIndices(std::vector<GlobalOrdinal> & indices) const
{
   // loop over field block manager and grab indices
   for(std::size_t fbm=0;fbm<fieldBlockManagers_.size();fbm++) {
      std::vector<GlobalOrdinalT> fieldBlockOwned;

      fieldBlockManagers_[fbm]->getOwnedIndices(fieldBlockOwned);

      for(std::size_t i=0;i<fieldBlockOwned.size();i++) 
         indices.push_back(std::make_pair(fbm,fieldBlockOwned[i]));
   }
}

template <typename LocalOrdinalT,typename GlobalOrdinalT>
void BlockedDOFManager<LocalOrdinalT,GlobalOrdinalT>::getOwnedAndSharedIndices(std::vector<GlobalOrdinal> & indices) const
{
   // loop over field block manager and grab indices
   for(std::size_t fbm=0;fbm<fieldBlockManagers_.size();fbm++) {
      std::vector<GlobalOrdinalT> fieldBlockOwned;

      fieldBlockManagers_[fbm]->getOwnedAndSharedIndices(fieldBlockOwned);

      for(std::size_t i=0;i<fieldBlockOwned.size();i++) 
         indices.push_back(std::make_pair(fbm,fieldBlockOwned[i]));
   }
}

template <typename LocalOrdinalT,typename GlobalOrdinalT>
void BlockedDOFManager<LocalOrdinalT,GlobalOrdinalT>::ownedIndices(const std::vector<GlobalOrdinal> & indices,std::vector<bool> & isOwned) const
{
   isOwned.resize(0);

   std::vector<std::vector<GlobalOrdinalT> > blockIndices(fieldBlockManagers_.size());
   for(std::size_t i=0;i<indices.size();i++)
      blockIndices[indices[i].first].push_back(indices[i].second); 
 
   // build bool vector stating if each sub block is owned
   std::vector<std::vector<bool> > blockIsOwned(fieldBlockManagers_.size());
   std::vector<std::vector<bool>::const_iterator> blockItrs;
   for(std::size_t fbm=0;fbm<fieldBlockManagers_.size();fbm++) {
      fieldBlockManagers_[fbm]->ownedIndices(blockIndices[fbm],blockIsOwned[fbm]);

      // setup iterators to boolean vectors
      blockItrs.push_back(blockIsOwned[fbm].begin());
   }

   // loop over indices, consider their block and look it up
   // in iterator vector
   for(std::size_t i=0;i<indices.size();i++) {
      int block = indices[i].first;

      // set owned status from iterator of block
      bool owned = *blockItrs[block];
      isOwned.push_back(owned);

      // increment block iterator
      blockItrs[block]++;
   }

   // quick error sanity check
   TEUCHOS_ASSERT(isOwned.size()==indices.size());
   for(std::size_t fbm=0;fbm<fieldBlockManagers_.size();fbm++) {
      TEUCHOS_TEST_FOR_EXCEPTION(blockItrs[fbm]!=blockIsOwned[fbm].end(),std::logic_error,
                       "BlockedDOFManager::ownedIndices: Did not consume all sub block boolean entries as expected.");
   }
    
}


////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////


/** \brief Set the connection manager and MPI_Comm objects.
  *
  * Set the connection manager and MPI_Comm objects. If this method
  * is called more than once, the behavior is to reset the indices in
  * the DOF manager.  However, the fields will be the same (this assumes
  * that the element blocks are consistent with the fields). The indices
  * will need to be rebuilt by calling <code>buildGlobalUnknowns</code>.
  *
  * \param[in] connMngr Connection manager to use.
  * \param[in] mpiComm  Communicator to use.
  */
template <typename LocalOrdinalT,typename GlobalOrdinalT>
void BlockedDOFManager<LocalOrdinalT,GlobalOrdinalT>::setConnManager(const Teuchos::RCP<ConnManager<LocalOrdinalT,GlobalOrdinalT> > & connMngr,MPI_Comm mpiComm)
{
   communicator_ = Teuchos::rcp(new Teuchos::MpiComm<int>(Teuchos::opaqueWrapper(mpiComm)));

   // this kills any old connection manager as well as the old FEI objects
   resetIndices();

   connMngr_ = connMngr;

   mpiComm_ = *communicator_->getRawMpiComm();
}

/** \brief Reset the indicies for this DOF manager.
  *
  * This method resets the indices and wipes out internal state. This method
  * does preserve the fields and the patterns added to the object. Also the
  * old connection manager is returned.
  *
  * \returns Old connection manager.
  */
template <typename LocalOrdinalT,typename GlobalOrdinalT>
Teuchos::RCP<ConnManager<LocalOrdinalT,GlobalOrdinalT> > BlockedDOFManager<LocalOrdinalT,GlobalOrdinalT>::resetIndices()
{
   Teuchos::RCP<ConnManager<LocalOrdinalT,GlobalOrdinalT> > connMngr = connMngr_;

   connMngr_ = Teuchos::null;
   ownedGIDHashTable_.clear(); 
   blockGIDOffset_.clear();

   for(std::size_t fbm=0;fbm<fieldBlockManagers_.size();fbm++) {
     Teuchos::RCP<DOFManagerFEI<LocalOrdinalT,GlobalOrdinalT> > dofMngr 
         = Teuchos::rcp_dynamic_cast<DOFManagerFEI<LocalOrdinalT,GlobalOrdinalT> >(fieldBlockManagers_[fbm]);
     if(dofMngr!=Teuchos::null)
       dofMngr->resetIndices();
   }

   return connMngr;
}

template <typename LocalOrdinalT,typename GlobalOrdinalT>
void BlockedDOFManager<LocalOrdinalT,GlobalOrdinalT>::addField(const std::string & str,
                                                        const Teuchos::RCP<const FieldPattern> & pattern)
{
   std::vector<std::string> elementBlockIds;
   connMngr_->getElementBlockIds(elementBlockIds);

   // loop over blocks adding field pattern to each 
   for(std::size_t i=0;i<elementBlockIds.size();i++)
      addField(elementBlockIds[i],str,pattern);
}

template <typename LocalOrdinalT,typename GlobalOrdinalT>
void BlockedDOFManager<LocalOrdinalT,GlobalOrdinalT>::addField(const std::string & blockId,const std::string & str,
                                                        const Teuchos::RCP<const FieldPattern> & pattern)
{
   TEUCHOS_TEST_FOR_EXCEPTION(fieldsRegistered(),std::logic_error,
                      "BlockedDOFManager::addField: addField cannot be called after registerFields or"
                      "buildGlobalUnknowns has been called"); 

   fieldStringToPattern_[std::make_pair(blockId,str)] = pattern;
   blockIdToFieldStrings_[blockId].insert(str);
}

template <typename LocalOrdinalT,typename GlobalOrdinalT>
void BlockedDOFManager<LocalOrdinalT,GlobalOrdinalT>::registerFields() 
{
   fieldBlockManagers_.clear();
   fieldStrToNum_.clear();
   fieldNumToStr_.clear();
   fieldNumToFieldBlk_.clear();
   maxSubFieldNum_ = -1;
   
   fieldsRegistered_ = false;

   // test validity of the field order, build default if none is provided
   {
      // build a unique set of fields, so we can compare validate the ordered list
      std::set<std::string> fields;
      for(std::map<std::pair<std::string,std::string>,Teuchos::RCP<const FieldPattern> >::const_iterator
          fieldItr=fieldStringToPattern_.begin(); fieldItr!=fieldStringToPattern_.end();++fieldItr) {
         std::string fieldName = fieldItr->first.second;
         fields.insert(fieldName);
      }

      // construct default field order if neccessary
      if(fieldOrder_.size()==0) {
         std::set<std::string>::const_iterator itr;
         for(itr=fields.begin();itr!=fields.end();itr++) {
            std::vector<std::string> block;
            block.push_back(*itr);
            fieldOrder_.push_back(block);
         }
      }

      // check validity of field order: no repeats, and everything is accounted for
      bool validOrder = validFieldOrder(fieldOrder_,fields);
      if(!validOrder) {
         // for outputing
         std::stringstream ss;

         ss << "BlockedDOFManager::registerFields - Field order is invalid!\n";

         ss << "   fields = [ ";
         for(std::set<std::string>::const_iterator itr=fields.begin();
             itr!=fields.end();++itr)
            ss << "\"" << *itr << "\" ";
         ss << " ]\n";

         ss << "   fieldOrder = [ ";
         for(std::vector<std::vector<std::string> >::const_iterator bitr=fieldOrder_.begin();
             bitr!=fieldOrder_.end();++bitr) {
            ss << "[ ";
            for(std::vector<std::string>::const_iterator itr=bitr->begin();
                itr!=bitr->end();++itr) {
               ss << "\"" << *itr << "\" ";
            }
            ss << " ], ";
         }
         ss << " ]\n";

         TEUCHOS_TEST_FOR_EXCEPTION(!validOrder,std::logic_error,ss.str());
      }
   }

   // build sub DOFManagers for each field block
   for(std::size_t fldBlk=0;fldBlk<fieldOrder_.size();fldBlk++) {
      Teuchos::RCP<panzer::UniqueGlobalIndexer<LocalOrdinalT,GlobalOrdinalT> > dofManager = buildNewIndexer(getConnManager(),mpiComm_);

      // add in these fields to the new manager
      this->addFieldsToFieldBlockManager(fieldOrder_[fldBlk],*dofManager);

      fieldBlockManagers_.push_back(dofManager); 
   }
    
   ////////////////////////////////
   // build field numbers: two stage algorithm

   // 1st Stage: Extract field numbers used by each sub DOFManager.
   //            determine largest of these
   //
   // - note at this point since "validFieldOrder" has
   //   been called we are gurranteed to not have repeated fields
   maxSubFieldNum_ = -1;
   std::map<std::string,int> tempStrToNum;
   for(std::size_t fldBlk=0;fldBlk<fieldBlockManagers_.size();fldBlk++) {
      Teuchos::RCP<panzer::UniqueGlobalIndexer<LocalOrdinalT,GlobalOrdinalT> > dofManager =
         fieldBlockManagers_[fldBlk];
      const std::vector<std::string> & activeFields = fieldOrder_[fldBlk];

      int fieldNum = 0;
      for(std::size_t f=0;f<activeFields.size();f++) {
         fieldNum = dofManager->getFieldNum(activeFields[f]); 
         tempStrToNum[activeFields[f]] = fieldNum;

         maxSubFieldNum_ = (fieldNum>maxSubFieldNum_) ? fieldNum : maxSubFieldNum_;
      }
   }

   // 2nd Stage: Using field block index, field number and largest field number
   //            build a up fieldStrToNum_ map and fieldNumToFieldBlk_
   int numOffset = 0;
   for(std::size_t fldBlk=0;fldBlk<fieldBlockManagers_.size();fldBlk++) {
      const std::vector<std::string> & activeFields = fieldOrder_[fldBlk];
      for(std::size_t f=0;f<activeFields.size();f++) {
         // compute offset field number
         int fieldNum = tempStrToNum[activeFields[f]]+numOffset;

         // build up map data
         fieldStrToNum_[activeFields[f]] = fieldNum;
         fieldNumToStr_[fieldNum] = activeFields[f];
         fieldNumToFieldBlk_[fieldNum] = fldBlk; 
      }

      // increament field number offset based on largest sub field number
      numOffset += (maxSubFieldNum_+1);
   }

   // end build field numbers 
   ////////////////////////////////

   // build block to field numbers: this requires field numbers have been built
   // and that "getFieldNum" behaves correctly
   for(std::map<std::string,std::set<std::string> >::const_iterator itr=blockIdToFieldStrings_.begin();
       itr!=blockIdToFieldStrings_.end();++itr) {
      const std::set<std::string> & fields = itr->second;

      std::vector<int> & fieldNums = blockIdToFieldNumbers_[itr->first];
      for(std::set<std::string>::const_iterator fldItr=fields.begin();  
          fldItr!=fields.end();++fldItr) {
         fieldNums.push_back(getFieldNum(*fldItr));
      }
   }

   // everything completed, mark as fields registered
   fieldsRegistered_ = true;
}

template <typename LocalOrdinalT,typename GlobalOrdinalT>
Teuchos::RCP<UniqueGlobalIndexer<LocalOrdinalT,GlobalOrdinalT> > 
BlockedDOFManager<LocalOrdinalT,GlobalOrdinalT>::
buildNewIndexer(const Teuchos::RCP<ConnManager<LocalOrdinalT,GlobalOrdinalT> > & connManager,MPI_Comm mpiComm) const
{
  if(getUseDOFManagerFEI()) {
    Teuchos::RCP<panzer::DOFManagerFEI<LocalOrdinalT,GlobalOrdinalT> > dofManager = Teuchos::rcp(new panzer::DOFManagerFEI<LocalOrdinalT,GlobalOrdinalT>);
    dofManager->setConnManager(connManager,mpiComm);

    return dofManager;
  }
  else {
    Teuchos::RCP<panzer::DOFManager<LocalOrdinalT,GlobalOrdinalT> > dofManager = Teuchos::rcp(new panzer::DOFManager<LocalOrdinalT,GlobalOrdinalT>);
    dofManager->setConnManager(connManager,mpiComm);

    return dofManager;
  }

}

template <typename LocalOrdinalT,typename GlobalOrdinalT>
void BlockedDOFManager<LocalOrdinalT,GlobalOrdinalT>::
setOrientationsRequired(const Teuchos::RCP<UniqueGlobalIndexer<LocalOrdinalT,GlobalOrdinalT> > & indexer,bool required) const
{
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;

  // standard version
  {
    RCP<DOFManager<LocalOrdinalT,GlobalOrdinalT> > dofManager = rcp_dynamic_cast<DOFManager<LocalOrdinalT,GlobalOrdinalT> >(indexer);

    if(dofManager!=Teuchos::null) {
      dofManager->setOrientationsRequired(required);
      return;
    }
  }

  // now the FEI version
  {
    RCP<DOFManagerFEI<LocalOrdinalT,GlobalOrdinalT> > dofManager = rcp_dynamic_cast<DOFManagerFEI<LocalOrdinalT,GlobalOrdinalT> >(indexer);

    if(dofManager!=Teuchos::null) {
      dofManager->setOrientationsRequired(required);
      return;
    }
  }

  // you should never get here!
  TEUCHOS_ASSERT(false);
}

template <typename LocalOrdinalT,typename GlobalOrdinalT>
void BlockedDOFManager<LocalOrdinalT,GlobalOrdinalT>::
buildGlobalUnknowns(const Teuchos::RCP<UniqueGlobalIndexer<LocalOrdinalT,GlobalOrdinalT> > & indexer,const Teuchos::RCP<const FieldPattern> & geomPattern) const
{
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;

  // standard version
  {
    RCP<DOFManager<LocalOrdinalT,GlobalOrdinalT> > dofManager = rcp_dynamic_cast<DOFManager<LocalOrdinalT,GlobalOrdinalT> >(indexer);

    if(dofManager!=Teuchos::null) {
      dofManager->buildGlobalUnknowns(geomPattern);
      return;
    }
  }

  // now the FEI version
  {
    RCP<DOFManagerFEI<LocalOrdinalT,GlobalOrdinalT> > dofManager = rcp_dynamic_cast<DOFManagerFEI<LocalOrdinalT,GlobalOrdinalT> >(indexer);

    if(dofManager!=Teuchos::null) {
      dofManager->buildGlobalUnknowns(geomPattern);
      return;
    }
  }

  // you should never get here!
  TEUCHOS_ASSERT(false);
}

template <typename LocalOrdinalT,typename GlobalOrdinalT>
void BlockedDOFManager<LocalOrdinalT,GlobalOrdinalT>::
printFieldInformation(const Teuchos::RCP<UniqueGlobalIndexer<LocalOrdinalT,GlobalOrdinalT> > & indexer,std::ostream & os) const
{
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;

  // standard version
  {
    RCP<DOFManager<LocalOrdinalT,GlobalOrdinalT> > dofManager = rcp_dynamic_cast<DOFManager<LocalOrdinalT,GlobalOrdinalT> >(indexer);

    if(dofManager!=Teuchos::null) {
      dofManager->printFieldInformation(os);
      return;
    }
  }

  // now the FEI version
  {
    RCP<DOFManagerFEI<LocalOrdinalT,GlobalOrdinalT> > dofManager = rcp_dynamic_cast<DOFManagerFEI<LocalOrdinalT,GlobalOrdinalT> >(indexer);

    if(dofManager!=Teuchos::null) {
      dofManager->printFieldInformation(os);
      return;
    }
  }

  // you should never get here!
  TEUCHOS_ASSERT(false);
}

template <typename LocalOrdinalT,typename GlobalOrdinalT>
int BlockedDOFManager<LocalOrdinalT,GlobalOrdinalT>::
getElementBlockGIDCount(const Teuchos::RCP<UniqueGlobalIndexer<LocalOrdinalT,GlobalOrdinalT> > & indexer,const std::string & elementBlock) const
{
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;

  // standard version
  {
    RCP<DOFManager<LocalOrdinalT,GlobalOrdinalT> > dofManager = rcp_dynamic_cast<DOFManager<LocalOrdinalT,GlobalOrdinalT> >(indexer);

    if(dofManager!=Teuchos::null) 
      return dofManager->getElementBlockGIDCount(elementBlock);
  }

  // now the FEI version
  {
    RCP<DOFManagerFEI<LocalOrdinalT,GlobalOrdinalT> > dofManager = rcp_dynamic_cast<DOFManagerFEI<LocalOrdinalT,GlobalOrdinalT> >(indexer);

    if(dofManager!=Teuchos::null)
      return dofManager->getElementBlockGIDCount(elementBlock);
  }

  // you should never get here!
  TEUCHOS_ASSERT(false);

  return -1;
}

template <typename LocalOrdinalT,typename GlobalOrdinalT>
void BlockedDOFManager<LocalOrdinalT,GlobalOrdinalT>::
addFieldsToFieldBlockManager(const std::vector<std::string> & activeFields,
                             UniqueGlobalIndexer<LocalOrdinalT,GlobalOrdinalT> & fieldBlockManager) const
{
  using Teuchos::Ptr;
  using Teuchos::ptrFromRef;
  using Teuchos::ptr_dynamic_cast;

  Ptr<UniqueGlobalIndexer<LocalOrdinalT,GlobalOrdinalT> > ugi_ptr = ptrFromRef(fieldBlockManager);

  // standard version
  {
    Ptr<DOFManager<LocalOrdinalT,GlobalOrdinalT> > dofManager_ptr = ptr_dynamic_cast<DOFManager<LocalOrdinalT,GlobalOrdinalT> >(ugi_ptr);

    if(dofManager_ptr!=Teuchos::null) {
      addFieldsToFieldBlockManager(activeFields,*dofManager_ptr);
      return;
    }
  }

  // now the FEI version
  {
    Ptr<DOFManagerFEI<LocalOrdinalT,GlobalOrdinalT> > dofManager_ptr = ptr_dynamic_cast<DOFManagerFEI<LocalOrdinalT,GlobalOrdinalT> >(ugi_ptr);

    if(dofManager_ptr!=Teuchos::null) {
      addFieldsToFieldBlockManager(activeFields,*dofManager_ptr);
      return;
    }
  }

  // you should never get here!
  TEUCHOS_ASSERT(false);
}

template <typename LocalOrdinalT,typename GlobalOrdinalT>
void BlockedDOFManager<LocalOrdinalT,GlobalOrdinalT>::
addFieldsToFieldBlockManager(const std::vector<std::string> & activeFields,
                             DOFManagerFEI<LocalOrdinalT,GlobalOrdinalT> & fieldBlockManager) const
{
   std::vector<std::size_t> correctnessCheck(activeFields.size(),0);
   std::vector<std::string> elementBlocks;
   this->getElementBlockIds(elementBlocks);

   // loop over element blocks adding each field in this element block and this field block
   for(std::size_t eb=0;eb<elementBlocks.size();eb++) {
      std::string elementBlock = elementBlocks[eb];

      // loop over active fields extracting those that are associated with this element block
      for(std::size_t f=0;f<activeFields.size();f++) {
         std::string fieldName = activeFields[f];
         Teuchos::RCP<const FieldPattern> fp = this->getFieldPattern(elementBlock,fieldName);

         if(fp!=Teuchos::null) {
            fieldBlockManager.addField(elementBlock,fieldName,fp);
            correctnessCheck[f] = 1; // all active fields should be placed in DOFManager
         }
      }
   }

   // verify correctness check
   std::size_t correctFlag = std::accumulate(correctnessCheck.begin(),correctnessCheck.end(),0);
   TEUCHOS_TEST_FOR_EXCEPTION(correctFlag!=activeFields.size(),std::logic_error,
                      "BlockedDOFManager::addFieldsToFieldBlockManager detected inconsistincies in the active fields.");

   // set field order
   fieldBlockManager.setFieldOrder(activeFields);

   // register added fields
   fieldBlockManager.registerFields();
}

template <typename LocalOrdinalT,typename GlobalOrdinalT>
void BlockedDOFManager<LocalOrdinalT,GlobalOrdinalT>::
addFieldsToFieldBlockManager(const std::vector<std::string> & activeFields,
                             DOFManager<LocalOrdinalT,GlobalOrdinalT> & fieldBlockManager) const
{
   std::vector<std::size_t> correctnessCheck(activeFields.size(),0);
   std::vector<std::string> elementBlocks;
   this->getElementBlockIds(elementBlocks);

   // loop over element blocks adding each field in this element block and this field block
   for(std::size_t eb=0;eb<elementBlocks.size();eb++) {
      std::string elementBlock = elementBlocks[eb];

      // loop over active fields extracting those that are associated with this element block
      for(std::size_t f=0;f<activeFields.size();f++) {
         std::string fieldName = activeFields[f];
         Teuchos::RCP<const FieldPattern> fp = this->getFieldPattern(elementBlock,fieldName);

         if(fp!=Teuchos::null) {
            fieldBlockManager.addField(elementBlock,fieldName,fp);
            correctnessCheck[f] = 1; // all active fields should be placed in DOFManager
         }
      }
   }

   // verify correctness check
   std::size_t correctFlag = std::accumulate(correctnessCheck.begin(),correctnessCheck.end(),0);
   TEUCHOS_TEST_FOR_EXCEPTION(correctFlag!=activeFields.size(),std::logic_error,
                      "BlockedDOFManager::addFieldsToFieldBlockManager detected inconsistincies in the active fields.");

   // set field order
   fieldBlockManager.setFieldOrder(activeFields);
}

template <typename LocalOrdinalT,typename GlobalOrdinalT>
void BlockedDOFManager<LocalOrdinalT,GlobalOrdinalT>::setFieldOrder(const std::vector<std::vector<std::string> > & fieldOrder)
{
   fieldOrder_ = fieldOrder;
}

/** Get the field order used. Return the field strings.
  */
template <typename LocalOrdinalT,typename GlobalOrdinalT>
void BlockedDOFManager<LocalOrdinalT,GlobalOrdinalT>::getFieldOrder(std::vector<std::vector<std::string> > & fieldOrder) const
{
   fieldOrder = fieldOrder_;
}

template <typename LocalOrdinalT,typename GlobalOrdinalT>
int BlockedDOFManager<LocalOrdinalT,GlobalOrdinalT>::getNumFields() const
{
   if(fieldsRegistered())
      return fieldStrToNum_.size();

   // more work needs to be done if the fields have not yet been registered
   // pull it from the (block id x field name) ==> pattern map
   std::set<std::string> fields;
   std::map<std::pair<std::string,std::string>,Teuchos::RCP<const FieldPattern> >::const_iterator itr;
   for(itr=fieldStringToPattern_.begin();itr!=fieldStringToPattern_.end();++itr)
      fields.insert(itr->first.second);

   return fields.size();
}

// build the global unknown numberings
//   1. this builds the pattens
//   2. initializes the connectivity
//   3. calls initComplete
template <typename LocalOrdinalT,typename GlobalOrdinalT>
void BlockedDOFManager<LocalOrdinalT,GlobalOrdinalT>::buildGlobalUnknowns(const Teuchos::RCP<const FieldPattern> & geomPattern)
{
   if(!fieldsRegistered()) {
      // std::cout << "register fields in" << std::endl;
      registerFields();
   }

   // save the geometry pattern
   geomPattern_ = geomPattern;

   // build global unknowns for each field block
   for(std::size_t fb=0;fb<fieldBlockManagers_.size();fb++) {
      // std::cout << "building field block fb = " << fb << std::endl;
      // fieldBlockManagers_[fb]->setOrientationsRequired(getOrientationsRequired());
      // fieldBlockManagers_[fb]->buildGlobalUnknowns(geomPattern_);

      setOrientationsRequired(fieldBlockManagers_[fb],getOrientationsRequired());
      buildGlobalUnknowns(fieldBlockManagers_[fb],geomPattern_);
   }

   // build field block offsets: this helps fast construction
   // of GID offset vectors. GIDs are ordering by field block.
   std::vector<std::string> elementBlocks;
   getElementBlockIds(elementBlocks);
   for(std::size_t eb=0;eb<elementBlocks.size();eb++) {
      int offset = 0;
      for(std::size_t fb=0;fb<fieldBlockManagers_.size();fb++) {
         // int cnt = fieldBlockManagers_[fb]->getElementBlockGIDCount(elementBlocks[eb]);
         int cnt = getElementBlockGIDCount(fieldBlockManagers_[fb],elementBlocks[eb]);
         blockGIDOffset_[std::make_pair(elementBlocks[eb],fb)] = offset;
         offset += cnt;
      }
   }
}

// build the global unknown numberings
//   1. this builds the pattens
//   2. initializes the connectivity
//   3. calls initComplete
template <typename LocalOrdinalT,typename GlobalOrdinalT>
void BlockedDOFManager<LocalOrdinalT,GlobalOrdinalT>::buildGlobalUnknowns()
{
   if(!fieldsRegistered())
      registerFields();

   // build the pattern for the ID layout on the mesh
   std::vector<RCP<const FieldPattern> > patVector;
   RCP<GeometricAggFieldPattern> aggFieldPattern = Teuchos::rcp(new GeometricAggFieldPattern);;
   std::map<std::pair<std::string,std::string>,Teuchos::RCP<const FieldPattern> >::iterator f2p_itr;
   for(f2p_itr=fieldStringToPattern_.begin();f2p_itr!=fieldStringToPattern_.end();f2p_itr++)
      patVector.push_back(f2p_itr->second);
   aggFieldPattern->buildPattern(patVector);

   // setup connectivity mesh
   connMngr_->buildConnectivity(*aggFieldPattern);

   // using new geometric pattern, build global unknowns
   buildGlobalUnknowns(aggFieldPattern);
}

template <typename LocalOrdinalT,typename GlobalOrdinalT>
void BlockedDOFManager<LocalOrdinalT,GlobalOrdinalT>::printFieldInformation(std::ostream & os) const
{
   os << "BlockedDOFManager Field Information: " << std::endl;
   
   if(fieldsRegistered()) {
      // Print field block DOF managers
      for(std::size_t fbm=0;fbm<fieldBlockManagers_.size();fbm++) {
         os << "*************************************************\n";
         os << "Field Block Index = " << fbm << std::endl;
         printFieldInformation(fieldBlockManagers_[fbm],os);

         // print out mapping between sub field IDs and blocked field IDs
         os << "   Field String to Field Id (blocked/sub):\n";
         for(std::size_t i=0;i<fieldOrder_[fbm].size();i++) {
            std::string fieldString = fieldOrder_[fbm][i];
            int fieldNum = getFieldNum(fieldString);
            os << "      \"" << fieldString << "\" is field ID " << fieldNum 
               << "/" << fieldBlockManagers_[fbm]->getFieldNum(fieldString) << std::endl;
         }
         os << std::endl;
      }
   }
   else {
      // fields are not registered
      os << "Fields not yet registered! Unknowns not built (call registerFields or buildGlobalUnknowns)" << std::endl;
   }
}

template <typename LocalOrdinalT,typename GlobalOrdinalT>
Teuchos::RCP<const FieldPattern> 
BlockedDOFManager<LocalOrdinalT,GlobalOrdinalT>::getFieldPattern(const std::string & blockId, const std::string & fieldName) const
{
   std::map<std::pair<std::string,std::string>,Teuchos::RCP<const FieldPattern> >::const_iterator itr;
   itr = fieldStringToPattern_.find(std::make_pair(blockId,fieldName));

   if(itr==fieldStringToPattern_.end()) // not found
      return Teuchos::null;
   else // found
      return itr->second;
}

/** Check the validity of a field order. This is used internally
  * as a sanity check. Checks for no repeats, bogus fields, and all fields
  * being included.
  *
  * \param[in] fieldOrder_ut Field order vector under test (ut).
  *
  * \returns true if the vector is valid, false otherwise.
  */
template <typename LocalOrdinalT,typename GlobalOrdinalT>
bool BlockedDOFManager<LocalOrdinalT,GlobalOrdinalT>::validFieldOrder(const std::vector<std::vector<std::string> > & fieldOrder_ut,
                                                                      const std::set<std::string> & fields) const
{
   std::set<std::string> orderedFields;
   std::size_t numberInOrder = 0;

   for(std::size_t b=0;b<fieldOrder_ut.size();b++) {
      numberInOrder += fieldOrder_ut[b].size();
      orderedFields.insert(fieldOrder_ut[b].begin(),
                           fieldOrder_ut[b].end());
   }

   bool correctCount = (numberInOrder==fields.size());
   bool sameFields = (orderedFields==fields);
 
   return correctCount && sameFields;
}

template <typename LocalOrdinalT,typename GlobalOrdinalT>
int BlockedDOFManager<LocalOrdinalT,GlobalOrdinalT>::getNumFieldBlocks() const
{
   if(fieldOrder_.size()==0) 
      return 1; // only one field block
   return fieldOrder_.size();
}

///////////////////////////////////////////////////////////////////////////

}

#endif
