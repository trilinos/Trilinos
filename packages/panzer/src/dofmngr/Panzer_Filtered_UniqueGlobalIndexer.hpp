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

#ifndef __Panzer_Filtered_UniqueGlobalIndexer_hpp__
#define __Panzer_Filtered_UniqueGlobalIndexer_hpp__

#include "Panzer_UniqueGlobalIndexer.hpp"

namespace panzer {

/**
  */
template <typename LocalOrdinalT,typename GlobalOrdinalT>
class Filtered_UniqueGlobalIndexer : public UniqueGlobalIndexer<LocalOrdinalT,GlobalOrdinalT> {
public:
   // These functions are unique to this class (including constructors)
   
   /** Default constructor
     */
   Filtered_UniqueGlobalIndexer();

   /** Initialize method that allows use of the default constructor and
     * may help with further inheritence.
     *
     * \param[in] ugi The global indexer to filter the global IDs of
     * \param[in] filteredIndices Indices to be filtered out of the <code>ugi</code> argument
     *
     * \note Repeated or unused (not in <code>ugi.getOwnedIndices</code>)indices in 
     *       <code>filteredIndices</code> are ignored without detection or impact.
     */
   void initialize(const Teuchos::RCP<const UniqueGlobalIndexer<LocalOrdinalT,GlobalOrdinalT> > & ugi,
                   const std::vector<GlobalOrdinalT> & filteredIndices);

   // This functions are overriden, and the filtered indices removed
  
   virtual void getOwnedIndices(std::vector<GlobalOrdinalT> & indices) const;

   virtual void ownedIndices(const std::vector<GlobalOrdinalT> & indices,std::vector<bool> & isOwned) const;

   // The following functions are simply part of the decorator pattern and
   // are simple pass throughs
   
   virtual ~Filtered_UniqueGlobalIndexer() {}

   virtual Teuchos::RCP<Teuchos::Comm<int> > getComm() const 
   { return base_->getComm(); }

   virtual int getNumFields() const 
   { return base_->getNumFields(); }

   virtual const std::string & getFieldString(int fieldNum) const
   { return base_->getFieldString(fieldNum); }

   virtual int getFieldNum(const std::string & str) const 
   { return base_->getFieldNum(str); }

   virtual void getFieldOrder(std::vector<std::string> & fieldOrder) const 
   { base_->getFieldOrder(fieldOrder); }

   virtual void getElementBlockIds(std::vector<std::string> & elementBlockIds) const 
   { base_->getElementBlockIds(elementBlockIds); }

   virtual bool fieldInBlock(const std::string & field, const std::string & block) const 
   { return base_->fieldInBlock(field,block); }

   virtual const std::vector<int> & getBlockFieldNumbers(const std::string & blockId) const 
   { return base_->getBlockFieldNumbers(blockId); }

   virtual const std::vector<int> & getGIDFieldOffsets(const std::string & blockId,int fieldNum) const 
   { return base_->getGIDFieldOffsets(blockId,fieldNum); }

   virtual const std::pair<std::vector<int>,std::vector<int> > & 
   getGIDFieldOffsets_closure(const std::string & blockId, int fieldNum,
                              int subcellDim,int subcellId) const 
   { return base_->getGIDFieldOffsets_closure(blockId,fieldNum,subcellDim,subcellId); }

   virtual void getElementOrientation(LocalOrdinalT localElmtId,std::vector<double> & gidsOrientation) const 
   { base_->getElementOrientation(localElmtId,gidsOrientation); }

   virtual const std::vector<LocalOrdinalT> & getElementBlock(const std::string & blockId) const 
   { return base_->getElementBlock(blockId); }

   virtual void getElementGIDs(LocalOrdinalT localElmtId,std::vector<GlobalOrdinalT> & gids,const std::string & blockIdHint="") const 
   { base_->getElementGIDs(localElmtId,gids,blockIdHint); }

   virtual void getOwnedAndSharedIndices(std::vector<GlobalOrdinalT> & indices) const 
   { base_->getOwnedAndSharedIndices(indices); }

   virtual int getElementBlockGIDCount(const std::string & blockId) const 
   { return base_->getElementBlockGIDCount(blockId); }

   virtual int getElementBlockGIDCount(const std::size_t & blockIndex) const 
   { return base_->getElementBlockGIDCount(blockIndex); }

   virtual void getElementAndAssociatedLIDs(LocalOrdinalT localElmtId, std::vector<LocalOrdinalT>& lids) const
   { base_->getElementAndAssociatedLIDs(localElmtId,lids); }

   virtual void getElementAndAssociatedGIDs(LocalOrdinalT localElmtId, std::vector<GlobalOrdinalT>& gids) const
   { base_->getElementAndAssociatedGIDs(localElmtId,gids); }

private:

   Teuchos::RCP<const UniqueGlobalIndexer<LocalOrdinalT,GlobalOrdinalT> > base_;

   std::vector<GlobalOrdinalT> owned_;
};

}

#endif
