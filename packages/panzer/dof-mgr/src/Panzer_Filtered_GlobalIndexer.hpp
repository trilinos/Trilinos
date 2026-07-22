// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Panzer_Filtered_GlobalIndexer_hpp__
#define __Panzer_Filtered_GlobalIndexer_hpp__

#include "Panzer_GlobalIndexer.hpp"

namespace panzer {

/** This class wraps a DOFManager and removes certain DOFs from the
    linear system. This is used to filter out or remove boundary
    conditions from a DOFManager.
  */
class Filtered_GlobalIndexer : public GlobalIndexer {
public:
   // These functions are unique to this class (including constructors)
   
   /** Default constructor */
   Filtered_GlobalIndexer();

   /** Initialize method that allows use of the default constructor and
     * may help with further inheritence.
     *
     * \param[in] ugi The global indexer to filter the global IDs of
     * \param[in] filteredIndices Indices to be filtered out of the <code>ugi</code> argument
     *
     * \note Repeated or unused (not in <code>ugi.getOwnedIndices</code>)indices in 
     *       <code>filteredIndices</code> are ignored without detection or impact.
     */
   void initialize(const Teuchos::RCP<const GlobalIndexer> & ugi,
                   const std::vector<panzer::GlobalOrdinal> & filteredIndices);

   /** Get an indicator describing if a particular local GID has been filtered. This method
     * requires communication. 
     *
     * \param[out] indicator Vector the same length of output argument of
     *                       <code>getOwendAndGhostedIndices</code>. If a value is one it
     *                       is included (not filtered), if it is zero then the GID has
     *                       been filtered out.
     */
   void getOwnedAndGhostedNotFilteredIndicator(std::vector<int> & indicator) const;

   /** Get the set of filtered indices that are owned and ghosted. 
     *
     * \param[out] indices Set of filtered indices
     */
   void getFilteredOwnedAndGhostedIndices(std::vector<panzer::GlobalOrdinal> & indices) const ;

   // This functions are overriden, and the filtered indices removed
  
   /**
    *  \brief Get the set of indices owned by this processor.
    *
    *  \note This is the set of owned indices from the base
    *        `GlobalIndexer` with the filtered indices removed.
    *
    *  \param[out] indices A `vector` that will be filled with the indices
    *                      owned by this processor.
    */
   virtual void getOwnedIndices(std::vector<panzer::GlobalOrdinal>& indices) const
   {
     indices = owned_;
   }

   /**
    *  \brief Get the set of indices ghosted for this processor.
    *
    *  \note This is the set of owned indices from the base
    *        `GlobalIndexer` (UGI) that have been filtered out, combined
    *        with the ghosted indices from the base UGI.
    *
    *  \param[out] indices A `vector` that will be filled with the indices
    *                      ghosted for this processor.
    */
   virtual void getGhostedIndices(std::vector<panzer::GlobalOrdinal>& indices) const
   { 
     indices = ghosted_;
   }

   /**
    *  \brief Get the set of owned and ghosted indices for this processor.
    *
    *  \note This is the set of owned and ghosted indices from the base
    *        `GlobalIndexer`, regardless of filtering.
    *
    *  \param[out] indices A `vector` that will be filled with the owned and
    *                      ghosted indices for this processor.
    */
   virtual void
   getOwnedAndGhostedIndices(std::vector<panzer::GlobalOrdinal>& indices) const 
   { 
     using std::size_t;
     indices.resize(owned_.size() + ghosted_.size());
     for (size_t i(0); i < owned_.size(); ++i)
       indices[i] = owned_[i];
     for (size_t i(0); i < ghosted_.size(); ++i)
       indices[owned_.size() + i] = ghosted_[i];
   }

   // For backwards compatibility with Epetra. Will be deprecated.
   virtual void getElementGIDsAsInt(panzer::LocalOrdinal localElmtId,std::vector<int> & gids,const std::string & blockIdHint="") const 
   { base_->getElementGIDsAsInt(localElmtId,gids,blockIdHint); }

   // For backwards compatibility with Epetra. Will be deprecated.
   virtual void getOwnedIndicesAsInt(std::vector<int>& indices) const
   {
     indices.resize(owned_.size());
     for (std::size_t i=0; i < owned_.size(); ++i)
       indices[i] = owned_[i];
   }
 
   // For backwards compatibility with Epetra. Will be deprecated.
   virtual void getGhostedIndicesAsInt(std::vector<int>& indices) const
   { 
     indices.resize(ghosted_.size());
     for (std::size_t i=0; i < ghosted_.size(); ++i)
       indices[i] = ghosted_[i];
   }

   // For backwards compatibility with Epetra. Will be deprecated.
   virtual void
   getOwnedAndGhostedIndicesAsInt(std::vector<int>& indices) const 
   { 
     indices.resize(owned_.size() + ghosted_.size());
     for (std::size_t i=0; i < owned_.size(); ++i)
       indices[i] = owned_[i];
     for (std::size_t i=0; i < ghosted_.size(); ++i)
       indices[owned_.size() + i] = ghosted_[i];
   }


   /**
    *  \brief Get the number of indices owned by this processor.
    *
    *  \note This is the number of owned indices from the base
    *        `GlobalIndexer`, less the number of filtered indices.
    *
    *  \returns The number of indices owned by this processor.
    */
   virtual int getNumOwned() const
   { return owned_.size(); }

   /**
    *  \brief Get the number of indices ghosted for this processor.
    *
    *  \note This is the number of owned indices from the base
    *        `GlobalIndexer` (UGI) that have been filtered out, plus the
    *        number of ghosted indices from the base UGI.
    *
    *  \returns The number of indices ghosted for this processor.
    */
   virtual int getNumGhosted() const
   { return ghosted_.size(); }

   /**
    *  \brief Get the number of owned and ghosted indices for this processor.
    *
    *  \note This is the number of owned and ghosted indices from the base
    *        `GlobalIndexer`, regardless of filtering.
    *
    *  \returns The number of owned and ghosted indices for this processor.
    */
   virtual int getNumOwnedAndGhosted() const
   { return owned_.size() + ghosted_.size(); }

   virtual void ownedIndices(const std::vector<panzer::GlobalOrdinal> & indices,std::vector<bool> & isOwned) const;

   // The following functions are simply part of the decorator pattern and
   // are simple pass throughs
   
   virtual ~Filtered_GlobalIndexer() {}

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

   virtual void getElementOrientation(panzer::LocalOrdinal localElmtId,std::vector<double> & gidsOrientation) const 
   { base_->getElementOrientation(localElmtId,gidsOrientation); }

   virtual const std::vector<panzer::LocalOrdinal> & getElementBlock(const std::string & blockId) const 
   { return base_->getElementBlock(blockId); }

   virtual void getElementGIDs(panzer::LocalOrdinal localElmtId,std::vector<panzer::GlobalOrdinal> & gids,const std::string & blockIdHint="") const 
   { base_->getElementGIDs(localElmtId,gids,blockIdHint); }

   virtual int getElementBlockGIDCount(const std::string & blockId) const 
   { return base_->getElementBlockGIDCount(blockId); }

   virtual int getElementBlockGIDCount(const std::size_t & blockIndex) const 
   { return base_->getElementBlockGIDCount(blockIndex); }

   virtual Teuchos::RCP<const ConnManager> getConnManager() const
   { return base_->getConnManager(); }

private:

   Teuchos::RCP<const GlobalIndexer> base_;

   /**
    *  \brief The list of owned indices.
    *
    *  The list of the owned indices from the base `GlobalIndexer` with
    *  the filtered indices removed.
    */
   std::vector<panzer::GlobalOrdinal> owned_;

   /**
    *  \brief The list of ghosted indices.
    *
    *  The list of the owned indices from the base `GlobalIndexer` (UGI)
    *  that have been filtered out, combined with the ghosted indices from the
    *  base UGI.
    */
   std::vector<panzer::GlobalOrdinal> ghosted_;
};

}

#endif
