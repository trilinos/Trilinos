#ifndef __UnitTest_UniqueGlobalIndexer_hpp__
#define __UnitTest_UniqueGlobalIndexer_hpp__

#include <vector>
#include <string>

#include "Teuchos_RCP.hpp"

#include "Panzer_UniqueGlobalIndexer.hpp"

namespace panzer {
namespace unit_test {

/** This a dummy global indexer used for testing. It has two fields
  *    "U" and "T"
  * There is one element block called "block_0" with two elements. Each element
  * resides on a different processor.
  */
class UniqueGlobalIndexer : public virtual panzer::UniqueGlobalIndexer<short,int> {
public:
   UniqueGlobalIndexer(int rank,int procCount);

   ~UniqueGlobalIndexer() {}

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

   /** What are the blockIds included in this connection manager?
     */
   virtual void getElementBlockIds(std::vector<std::string> & elementBlockIds) const; 

   /** Get the local element IDs for a paricular element
     * block.
     *
     * \param[in] blockId Block ID
     *
     * \returns Vector of local element IDs.
     */
   virtual const std::vector<short> & getElementBlock(const std::string & blockId) const;

   /** \brief Get the global IDs for a particular element. This function
     * overwrites the <code>gids</code> variable.
     */
   virtual void getElementGIDs(short localElmtId,std::vector<int> & gids) const;

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

   /** Get set of indices owned by this processor
     */
   virtual void getOwnedIndices(std::vector<int> & indices) const;

   /** Get set of indices owned and shared by this processor.
     * This can be thought of as the ``ghosted'' indices.
     */
   virtual void getOwnedAndSharedIndices(std::vector<int> & indices) const;

   /** Get a yes/no on ownership for each index in a vector
     */
   virtual void ownedIndices(const std::vector<int> & indices,std::vector<bool> & isOwned) const;

private:
   int procRank_;
   mutable Teuchos::RCP<std::vector<short> > elements_; // local element IDs
   mutable Teuchos::RCP<std::vector<int> > field0Offset_; // local element IDs
   mutable Teuchos::RCP<std::vector<int> > field1Offset_; // local element IDs
};

} // end unit test
} // end panzer

#endif
