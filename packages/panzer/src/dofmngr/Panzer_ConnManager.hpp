#ifndef __Panzer_ConnManager_hpp__
#define __Panzer_ConnManager_hpp__

#include <vector>

// Teuchos includes
#include "Teuchos_RCP.hpp"

namespace panzer {

class FieldPattern; // from DOFManager

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
   virtual const int * getConnectivity(int localElmtId) const = 0;

   /** How many mesh IDs are associated with this element?
     *
     * \param[in] localElmtId Local element ID
     *
     * \returns Number of mesh IDs that are associated with this element.
     */
   virtual int getConnectivitySize(int localElmtId) const = 0;

   /** Get the block ID for a particular element.
     *
     * \param[in] localElmtId Local element ID
     */
   virtual int getBlockId(int localElmtId) const = 0;

   /** How many element blocks in this mesh?
     */
   virtual int numElementBlocks() const = 0;

   /** Get the local element IDs for a paricular element
     * block.
     *
     * \param[in] blockID Block ID
     *
     * \returns Vector of local element IDs.
     */
   virtual const std::vector<int> & getElementBlock(int blockID) const = 0;
};

};

#endif
