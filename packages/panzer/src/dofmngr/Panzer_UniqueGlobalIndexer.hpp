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

#ifndef __Panzer_UniqueGlobalIndexer_hpp__
#define __Panzer_UniqueGlobalIndexer_hpp__

#include <vector>
#include <string>

#include <boost/unordered_map.hpp> // a hash table for buildLocalIds()

#include "Teuchos_RCP.hpp"
#include "Teuchos_Comm.hpp"

namespace panzer {

class UniqueGlobalIndexerBase {
public:
   //! Pure virtual destructor: prevents warnings with inline empty implementation 
   virtual ~UniqueGlobalIndexerBase() = 0;

   /** Get communicator associated with this global indexer.
     */
   virtual Teuchos::RCP<Teuchos::Comm<int> > getComm() const = 0;

   /** Get the number of fields (total) stored by this DOF manager
     */
   virtual int getNumFields() const = 0;

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
   virtual int getFieldNum(const std::string & str) const = 0;

   /** Get the field order used by this global indexer.
     */
   virtual void getFieldOrder(std::vector<std::string> & fieldOrder) const = 0;

   /** \brief Reverse lookup of the field string from
     *        a field number.
     *
     * \param[in] num Field number. Assumed to be 
     *                a valid field number.  Computed
     *                from <code>getFieldNum</code>.
     *
     * \returns Field name. 
     */
   virtual const std::string & getFieldString(int num) const = 0;

   /** What are the blockIds included in this connection manager?
     */
   virtual void getElementBlockIds(std::vector<std::string> & elementBlockIds) const = 0; 

   /** Is the specified field in the element block? 
     */
   virtual bool fieldInBlock(const std::string & field, const std::string & block) const = 0;

   /** Get field numbers associated with a particular element block.
     */
   virtual const std::vector<int> & getBlockFieldNumbers(const std::string & blockId) const = 0;

   /** \brief Use the field pattern so that you can find a particular
     *        field in the GIDs array.
     */
   virtual const std::vector<int> & getGIDFieldOffsets(const std::string & blockId,int fieldNum) const = 0;

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
   virtual const std::pair<std::vector<int>,std::vector<int> > & 
   getGIDFieldOffsets_closure(const std::string & blockId, int fieldNum,
                              int subcellDim,int subcellId) const = 0;
};

template <typename LocalOrdinalT,typename GlobalOrdinalT>
class UniqueGlobalIndexer : public UniqueGlobalIndexerBase {
public:
   //! Pure virtual destructor: prevents warnings with inline empty implementation 
   virtual ~UniqueGlobalIndexer() = 0;

   /** Get communicator associated with this global indexer.
     */
   virtual Teuchos::RCP<Teuchos::Comm<int> > getComm() const = 0;

   /** Get the number of fields (total) stored by this DOF manager
     */
   virtual int getNumFields() const = 0;

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
   virtual int getFieldNum(const std::string & str) const = 0;

   /** Get the field order used by this global indexer.
     */
   virtual void getFieldOrder(std::vector<std::string> & fieldOrder) const = 0;

   /** What are the blockIds included in this connection manager?
     */
   virtual void getElementBlockIds(std::vector<std::string> & elementBlockIds) const = 0; 

   /** Is the specified field in the element block? 
     */
   virtual bool fieldInBlock(const std::string & field, const std::string & block) const = 0;

   /** Get field numbers associated with a particular element block.
     */
   virtual const std::vector<int> & getBlockFieldNumbers(const std::string & blockId) const = 0;

   /** \brief Use the field pattern so that you can find a particular
     *        field in the GIDs array.
     */
   virtual const std::vector<int> & getGIDFieldOffsets(const std::string & blockId,int fieldNum) const = 0;

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
   virtual const std::pair<std::vector<int>,std::vector<int> > & 
   getGIDFieldOffsets_closure(const std::string & blockId, int fieldNum,
                              int subcellDim,int subcellId) const = 0;

   // Methods requiring Local or Global OrdinalT
   ////////////////////////////////////////////////////////////////////////

   /** \brief Get a vector containg the orientation of the GIDs relative to the neighbors.
     */
   virtual void getElementOrientation(LocalOrdinalT localElmtId,std::vector<double> & gidsOrientation) const = 0;

   /** Get the local element IDs for a paricular element
     * block.
     *
     * \param[in] blockId Block ID
     *
     * \returns Vector of local element IDs.
     */
   virtual const std::vector<LocalOrdinalT> & getElementBlock(const std::string & blockId) const = 0;

   /** \brief Get the global IDs for a particular element. This function
     * overwrites the <code>gids</code> variable.
     */
   virtual void getElementGIDs(LocalOrdinalT localElmtId,std::vector<GlobalOrdinalT> & gids,const std::string & blockIdHint="") const = 0;

   /** Get set of indices owned by this processor
     */
   virtual void getOwnedIndices(std::vector<GlobalOrdinalT> & indices) const = 0;

   /** Get set of indices owned and shared by this processor.
     * This can be thought of as the ``ghosted'' indices.
     */
   virtual void getOwnedAndSharedIndices(std::vector<GlobalOrdinalT> & indices) const = 0;

   /** Get a yes/no on ownership for each index in a vector
     */
   virtual void ownedIndices(const std::vector<GlobalOrdinalT> & indices,std::vector<bool> & isOwned) const = 0;

   /** Access the local IDs for an element. The local ordering is according to
     * the <code>getOwnedAndSharedIndices</code> method.
     */
   const std::vector<LocalOrdinalT> & getElementLIDs(LocalOrdinalT localElmtId) const
   { return localIDs_[localElmtId]; }

protected:

   /** This method is used by derived classes to the construct the local IDs from 
     * the <code>getOwnedAndSharedIndices</code> method.
     */
   void buildLocalIds();

private:
   std::vector<std::vector<LocalOrdinalT> > localIDs_; 
};

// prevents a warning because a destructor does not exist
inline UniqueGlobalIndexerBase::~UniqueGlobalIndexerBase() {}

// prevents a warning because a destructor does not exist
template <typename LocalOrdinalT,typename GlobalOrdinalT>
inline UniqueGlobalIndexer<LocalOrdinalT,GlobalOrdinalT>::~UniqueGlobalIndexer() {}

template <typename LocalOrdinalT,typename GlobalOrdinalT>
inline void UniqueGlobalIndexer<LocalOrdinalT,GlobalOrdinalT>::buildLocalIds() 
{
  std::vector<GlobalOrdinalT> ownedAndShared;
  this->getOwnedAndSharedIndices(ownedAndShared);
   
  // build global to local hash map (temporary and used only once)
  boost::unordered_map<GlobalOrdinalT,LocalOrdinalT> hashMap;
  for(std::size_t i=0;i<ownedAndShared.size();i++)
    hashMap[ownedAndShared[i]] = i;

  std::vector<std::string> elementBlocks;
  this->getElementBlockIds(elementBlocks);
 
  // compute total number of elements
  std::size_t numElmts = 0;
  for(std::size_t eb=0;eb<elementBlocks.size();eb++)
    numElmts += this->getElementBlock(elementBlocks[eb]).size();
  localIDs_.resize(numElmts); // allocate local ids

  // perform computation of local ids
  for(std::size_t eb=0;eb<elementBlocks.size();eb++) {
    std::vector<GlobalOrdinalT> gids;
    const std::vector<LocalOrdinalT> & elmts = this->getElementBlock(elementBlocks[eb]);

    for(std::size_t e=0;e<elmts.size();e++) {
      this->getElementGIDs(elmts[e],gids,elementBlocks[eb]);
      std::vector<LocalOrdinalT> & lids = localIDs_[elmts[e]];
      lids.resize(gids.size());
 
      for(std::size_t g=0;g<gids.size();g++)
        lids[g] = hashMap[gids[g]];
    }
  } 
}

}

#endif
