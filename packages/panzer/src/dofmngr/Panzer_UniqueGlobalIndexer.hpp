#ifndef __Panzer_UniqueGlobalIndexer_hpp__
#define __Panzer_UniqueGlobalIndexer_hpp__

#include <vector>
#include <string>

#include "Teuchos_RCP.hpp"

#include "Panzer_ConnManager.hpp"

namespace panzer {

template <typename LocalOrdinalT,typename GlobalOrdinalT>
class UniqueGlobalIndexer {
public:
   //! Pure virtual destructor: prevents warnings with inline empty implementation 
   virtual ~UniqueGlobalIndexer() = 0;

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

   /** What are the blockIds included in this connection manager?
     */
   virtual void getElementBlockIds(std::vector<std::string> & elementBlockIds) const = 0; 

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
   virtual void getElementGIDs(LocalOrdinalT localElmtId,std::vector<GlobalOrdinalT> & gids) const = 0;

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
   // virtual const std::vector<int> & 
   virtual const std::pair<std::vector<int>,std::vector<int> > & 
   getGIDFieldOffsets_closure(const std::string & blockId, int fieldNum,
                                                               int subcellDim,int subcellId) const = 0;

   /** Get set of indices owned by this processor
     */
   virtual void getOwnedIndices(std::vector<GlobalOrdinalT> & indices) const = 0;

   /** Get set of indices owned and shared by this processor.
     * This can be thought of as the ``ghosted'' indices.
     */
   virtual void getOwnedAndSharedIndices(std::vector<GlobalOrdinalT> & indices) const = 0;

};

// prevents a warning because a destructor does not exist
template <typename LocalOrdinalT,typename GlobalOrdinalT>
inline UniqueGlobalIndexer<LocalOrdinalT,GlobalOrdinalT>::~UniqueGlobalIndexer() {}

}

#endif
