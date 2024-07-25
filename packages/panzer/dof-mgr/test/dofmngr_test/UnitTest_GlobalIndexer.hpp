// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __UnitTest_GlobalIndexer_hpp__
#define __UnitTest_GlobalIndexer_hpp__

#include <vector>
#include <string>

#include "Teuchos_RCP.hpp"
#include "Panzer_GlobalIndexer.hpp"

namespace panzer {
namespace unit_test {

/** This a dummy global indexer used for testing. It has two fields
  *    "U" and "T"
  * There is one element block called "block_0" with two elements. Each element
  * resides on a different processor.
  */
class GlobalIndexer : public virtual panzer::GlobalIndexer {
public:
   GlobalIndexer(int rank,int procCount);

   ~GlobalIndexer() {}

   /** \brief Get the number used for access to this
     *        field
     *
     * Get the number used for access to this
     * field. This is used as the input parameter
     * to the other functions that provide access
     * to the global unknowns.
     *
     * \param[in] str Human readable name of the field
     *
     * \returns A unique integer associated with the
     *          field if the field exisits. Otherwise
     *          a -1 is returned.
     */
   virtual int getFieldNum(const std::string & str) const;

   virtual int getNumFields() const { return 2; }
   virtual void getFieldOrder(std::vector<std::string> & order) const 
   { order.push_back("U"); order.push_back("T"); } 

   virtual const std::string & getFieldString(int fieldNum) const;

   /** Get the communicator 
     */
   virtual Teuchos::RCP<Teuchos::Comm<int> > getComm() const;

   /** What are the blockIds included in this connection manager?
     */
   virtual void getElementBlockIds(std::vector<std::string> & elementBlockIds) const; 

   /** Get field numbers associated with a particular element block.
     */
   virtual const std::vector<int> & getBlockFieldNumbers(const std::string & block) const;

   /** Is the specified field in the element block? 
     */
   virtual bool fieldInBlock(const std::string & field, const std::string & block) const;

   /** Get the local element IDs for a paricular element
     * block.
     *
     * \param[in] blockId Block ID
     *
     * \returns Vector of local element IDs.
     */
  virtual const std::vector<panzer::LocalOrdinal> & getElementBlock(const std::string & blockId) const;

   /** \brief Get the global IDs for a particular element. This function
     * overwrites the <code>gids</code> variable.
     */
   virtual void getElementGIDs(panzer::LocalOrdinal localElmtId,std::vector<panzer::GlobalOrdinal> & gids,const std::string & blockId="") const;

   virtual void getElementOrientation(panzer::LocalOrdinal /* localElmtId */, std::vector<double> & /* gidsOrientation */) const
   { TEUCHOS_ASSERT(false); }

   /** \brief Use the field pattern so that you can find a particular
     *        field in the GIDs array.
     */
   virtual const std::vector<int> & getGIDFieldOffsets(const std::string & blockId,int fieldNum) const;

   /** \brief Use the field pattern so that you can find a particular
     *        field in the GIDs array. This version lets you specify the sub
     *        cell you are interested in and gets the closure. Meaning all the
     *        IDs of equal or lesser sub cell dimension that are contained within
     *        the specified sub cell. For instance for an edge, this function would
     *        return offsets for the edge and the nodes on that edge.
     *
     * \param[in] blockId
     * \param[in] fieldNum
     * \param[in] subcellDim
     * \param[in] subcellId
     */
   // virtual const std::vector<int> & 
   virtual const std::pair<std::vector<int>,std::vector<int> > & 
   getGIDFieldOffsets_closure(const std::string & blockId, int fieldNum,
                              int subcellDim,int subcellId) const;

   /**
    *  \brief Get the set of indices owned by this processor.
    *
    *  \param[out] indices A `vector` that will be filled with the indices
    *                      owned by this processor.
    *
    *  \throws `std::logic_error` If the processor rank is neither 0 nor 1.
    */
   virtual void getOwnedIndices(std::vector<panzer::GlobalOrdinal>& indices) const;

   /**
    *  \brief Get the set of indices ghosted for this processor.
    *
    *  \param[out] indices A `vector` that will be filled with the indices
    *                      ghosted for this processor.
    *
    *  \throws `std::logic_error` If the processor rank is neither 0 nor 1.
    */
   virtual void getGhostedIndices(std::vector<panzer::GlobalOrdinal>& indices) const;

   /**
    *  \brief Get the set of owned and ghosted indices for this processor.
    *
    *  \param[out] indices A `vector` that will be filled with the owned and
    *                      ghosted indices for this processor.
    *
    *  \throws `std::logic_error` If the processor rank is neither 0 nor 1.
    */
   virtual void getOwnedAndGhostedIndices(std::vector<panzer::GlobalOrdinal>& indices) const;

   // Epetra accessors. Will be deprecated.
   virtual void getElementGIDsAsInt(panzer::LocalOrdinal localElmtId,std::vector<int> & gids,const std::string & blockId="") const;
   virtual void getOwnedIndicesAsInt(std::vector<int>& indices) const;
   virtual void getGhostedIndicesAsInt(std::vector<int>& indices) const;
   virtual void getOwnedAndGhostedIndicesAsInt(std::vector<int>& indices) const;

   /**
    *  \brief Get the number of indices owned by this processor.
    *
    *  \throws `std::logic_error` If the processor rank is neither 0 nor 1.
    *
    *  \returns The number of indices owned by this processor.
    */
   virtual int getNumOwned() const;

   /**
    *  \brief Get the number of indices ghosted for this processor.
    *
    *  \throws `std::logic_error` If the processor rank is neither 0 nor 1.
    *
    *  \returns The number of indices ghosted for this processor.
    */
   virtual int getNumGhosted() const;

   /**
    *  \brief Get the number of owned and ghosted indices for this processor.
    *
    *  \throws `std::logic_error` If the processor rank is neither 0 nor 1.
    *
    *  \returns The number of owned and ghosted indices for this processor.
    */
   virtual int getNumOwnedAndGhosted() const;

   /** Get a yes/no on ownership for each index in a vector
     */
   virtual void ownedIndices(const std::vector<panzer::GlobalOrdinal> & indices,std::vector<bool> & isOwned) const;

   int getElementBlockGIDCount(const std::string &) const;
   int getElementBlockGIDCount(const std::size_t &) const;

   /** Not needed, so return null.
     */
   virtual Teuchos::RCP<const ConnManager> getConnManager() const
   { return Teuchos::null; }

private:
   int procRank_;
   mutable Teuchos::RCP<std::vector<panzer::LocalOrdinal> > elements_b0_; // local element IDs
   mutable Teuchos::RCP<std::vector<panzer::LocalOrdinal> > elements_b1_; // local element IDs
   std::vector<panzer::LocalOrdinal> emptyElements_; // empty element IDs
   mutable Teuchos::RCP<std::vector<int> > field0Offset_b0_; // local element IDs
   mutable Teuchos::RCP<std::vector<int> > field1Offset_b0_; // local element IDs
   mutable Teuchos::RCP<std::vector<int> > field0Offset_b1_; // local element IDs
   mutable Teuchos::RCP<std::vector<int> > field1Offset_b1_; // local element IDs

   mutable Teuchos::RCP<std::vector<int> > block0Fields_; 
   mutable Teuchos::RCP<std::vector<int> > block1Fields_; 
};

} // end unit test
} // end panzer

#endif
