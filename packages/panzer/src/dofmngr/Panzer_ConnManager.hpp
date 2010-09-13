#ifndef __Panzer_ConnManager_hpp__
#define __Panzer_ConnManager_hpp__

#include <vector>

// Teuchos includes
#include "Teuchos_RCP.hpp"

namespace panzer {

class FieldPattern; // from DOFManager

/** Pure abstract base class templated on the
  * global and local ordinal types. It is assumed
  * that element blocks are number by a GlobalOrdinal
  * and local element IDs use the LocalOrdinal.
  */
template <typename LocalOrdinal,typename GlobalOrdinal>
class ConnManager {
public:
   virtual ~ConnManager() {}

   /** Tell the connection manager to build the connectivity assuming
     * a particular field pattern.
     *
     * \param[in] fp Field pattern to build connectivity for
     */
   virtual void buildConnectivity(const FieldPattern & fp) = 0;

   /** Get ID connectivity for a particular element
     *
     * \param[in] localElmtId Local element ID
     *
     * \returns Pointer to beginning of indices, with total size
     *          equal to <code>getConnectivitySize(localElmtId)</code>
     */
   virtual const GlobalOrdinal * getConnectivity(LocalOrdinal localElmtId) const = 0;

   /** How many mesh IDs are associated with this element?
     *
     * \param[in] localElmtId Local element ID
     *
     * \returns Number of mesh IDs that are associated with this element.
     */
   virtual int getConnectivitySize(LocalOrdinal localElmtId) const = 0;

   /** Get the block ID for a particular element.
     *
     * \param[in] localElmtId Local element ID
     */
   virtual GlobalOrdinal getBlockId(LocalOrdinal localElmtId) const = 0;

   /** How many element blocks in this mesh?
     */
   virtual GlobalOrdinal numElementBlocks() const = 0;

   /** Get the local element IDs for a paricular element
     * block.
     *
     * \param[in] blockID Block ID
     *
     * \returns Vector of local element IDs.
     */
   virtual const std::vector<LocalOrdinal> & getElementBlock(GlobalOrdinal blockID) const = 0;
};

};

#endif
