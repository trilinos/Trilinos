// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Panzer_ConnManager_hpp__
#define __Panzer_ConnManager_hpp__

#include <vector>

// Teuchos includes
#include "Teuchos_RCP.hpp"
#include "Shards_CellTopology.hpp"
#include "PanzerDofMgr_config.hpp"

namespace panzer {

class FieldPattern; // from DOFManager

  /// Pure virtual base class for supplying mesh connectivity information to the DOF Manager.
  class ConnManager {
  public:

    using GlobalOrdinal = panzer::GlobalOrdinal;
    using LocalOrdinal = int;

    virtual ~ConnManager() {}

    /** Tell the connection manager to build the connectivity assuming
     * a particular field pattern.
     *
     * \param[in] fp Field pattern to build connectivity for
     */
    virtual void buildConnectivity(const FieldPattern & fp) = 0;

    /** Build a clone of this connection manager, without any assumptions
     * about the required connectivity (i.e. <code>buildConnectivity</code>
     * has never been called).
     */
    virtual Teuchos::RCP<ConnManager> noConnectivityClone() const = 0;

    /** How many mesh IDs are associated with this element?
     *
     * \param[in] localElmtId Local element ID
     *
     * \returns Number of mesh IDs that are associated with this element.
     */
    virtual LocalOrdinal getConnectivitySize(LocalOrdinal localElmtId) const = 0;

    /** Get ID connectivity for a particular element
     *
     * \param[in] localElmtId Local element ID
     *
     * \returns Pointer to beginning of indices, with total size
     *          equal to <code>getConnectivitySize(localElmtId)</code>
     */
    virtual const GlobalOrdinal * getConnectivity(LocalOrdinal localElmtId) const = 0;

    /** Get the block ID for a particular element.
     *
     * \param[in] localElmtId Local element ID
     */
    virtual std::string getBlockId(LocalOrdinal localElmtId) const = 0;

    /** Returns the number of element blocks in this mesh */
    virtual std::size_t numElementBlocks() const = 0;

    /** What are the blockIds included in this connection manager */
    virtual void getElementBlockIds(std::vector<std::string> & elementBlockIds) const = 0;

    /** Returns the cellTopologies linked to element blocks in this connection manager */
    virtual void getElementBlockTopologies(std::vector<shards::CellTopology> & elementBlockTopologies) const = 0;

    /** Get the local element IDs for a paricular element
     * block.
     *
     * \param[in] blockID Block ID
     *
     * \returns Vector of local element IDs.
     */
    virtual const std::vector<LocalOrdinal> & getElementBlock(const std::string & blockID) const = 0;

    /** Get the local element IDs for all "neighbor" elements that reside in a paricular element
     * block (An element is a neighbor if it is in the one ring of owned elements).
     *
     * \param[in] blockID Block ID
     *
     * \returns Vector of local element IDs.
     */
    virtual const std::vector<LocalOrdinal> & getNeighborElementBlock(const std::string & blockID) const = 0;

    /** Get elements, if any, associated with <code>el</code>, excluding
     * <code>el</code> itself.
     */
    virtual const std::vector<LocalOrdinal>& getAssociatedNeighbors(const LocalOrdinal& el) const = 0;

    /** Return whether getAssociatedNeighbors will return true for at least one
     * input.
     */
    virtual bool hasAssociatedNeighbors() const = 0;
  };

}

#endif
