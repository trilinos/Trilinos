// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __UnitTest_ConnManager_hpp__
#define __UnitTest_ConnManager_hpp__

#include <vector>
#include <string>

#include "Teuchos_RCP.hpp"

#include "Panzer_ConnManager.hpp"
#include "Panzer_FieldPattern.hpp"

#ifdef HAVE_MPI
   #include "Teuchos_DefaultMpiComm.hpp"
   #include <mpi.h>
#else
   #include "Teuchos_DefaultSerialComm.hpp"
#endif

namespace panzer {
namespace unit_test {

// simple debugging call back
class ConnCallback { public: virtual void buildConnectivity(const FieldPattern & fp) = 0; };

/** This a dummy connection manager used for testing. It has three element blocks
  *    "block_0", "block_1", and "block_2".
  * There are 4 elements in blocks 0 and 1 and two in block 2.
  * This connection manager works on two processors and currently handles only nodal
  * field patterns.
  *
  * Here is a drawing showing the blocks, double lines distinguishes block boundaries
    where as single lines denote element boundaries.
 
       +------+------+
       |      |      |
       |   BLOCK 2   |
       |      |      |
       +======+======+-------+------+
       |      |      ||      |      |
       |      |      ||      |      |
       |      |      ||      |      |
       +---BLOCK 0---+----BLOCK 1---+
       |      |      ||      |      |
       |      |      ||      |      |
       |      |      ||      |      |
       +------+------+-------+------+

    This shows element global/local ids, and node global ids

       15-----16----17
       |      |      |
       |  8/2 |  9/3 |
       |      |      |
       10-----11-----12-----13----14
       |      |      |      |      |
       |  4/1 |  5/1 |  6/4 |  7/5 |
       |      |      |      |      |
       5------6------7------8------9
       |      |      |      |      |
       |  0/0 |  1/0 |  2/2 |  3/3 |
       |      |      |      |      |
       0------1------2------3------4

    This shows processor ownership

       +------+------+
       |      |      |
       |  P0  |  P0  |
       |      |      |
       +------+------+------+------+
       |      |      |      |      |
       |  P0  |  P1  |  P0  |  P0  |
       |      |      |      |      |
       +------+------+------+------+
       |      |      |      |      |
       |  P0  |  P1  |  P1  |  P1  |
       |      |      |      |      |
       +------+------+------+------+
  */
class ConnManager : public virtual panzer::ConnManager {
public:

   ConnManager(int rank,int procCount);

   ~ConnManager() {}

   /** Tell the connection manager to build the connectivity assuming
     * a particular field pattern.
     *
     * \param[in] fp Field pattern to build connectivity for
     */
   virtual void buildConnectivity(const FieldPattern & fp);

   /** Tell the connection manager to build the connectivity assuming
     * a particular field pattern.
     *
     * \param[in] fp Field pattern to build connectivity for
     */
   virtual Teuchos::RCP<panzer::ConnManager> noConnectivityClone() const
   { return Teuchos::rcp(new ConnManager(procRank_,2)); }

   /** Get ID connectivity for a particular element
     *
     * \param[in] localElmtId Local element ID
     *
     * \returns Pointer to beginning of indices, with total size
     *          equal to <code>getConnectivitySize(localElmtId)</code>
     */
   virtual const GlobalOrdinal * getConnectivity(LocalOrdinal localElmtId) const;

   /** How many mesh IDs are associated with this element?
     *
     * \param[in] localElmtId Local element ID
     *
     * \returns Number of mesh IDs that are associated with this element.
     */
   virtual LocalOrdinal getConnectivitySize(LocalOrdinal localElmtId) const;

   /** Get the block ID for a particular element.
     *
     * \param[in] localElmtId Local element ID
     */
   virtual std::string getBlockId(LocalOrdinal localElmtId) const;

   /** How many element blocks in this mesh?
     */
   virtual std::size_t numElementBlocks() const;

   /** What are the blockIds included in this connection manager?
     */
   virtual void getElementBlockIds(std::vector<std::string> & elementBlockIds) const;
   /** What are the cellTopologies linked to element blocks in this connection manager?
     */
   virtual void getElementBlockTopologies(std::vector<shards::CellTopology> & elementBlockTopologies) const;

   /** Get the local element IDs for a paricular element
     * block.
     *
     * \param[in] blockID Block ID
     *
     * \returns Vector of local element IDs.
     */
   virtual const std::vector<LocalOrdinal> & getElementBlock(const std::string & blockID) const;

   void setBuildConnectivityCallback(const Teuchos::RCP<ConnCallback> & callback)
   { callback_ = callback; }

   //! A dummy getNeighborElementBlock method included but not used
   virtual const std::vector<int> & getNeighborElementBlock(const std::string&) const
   { static std::vector<int> empty; return empty; }

   //! Not used.
   virtual const std::vector<int>& getAssociatedNeighbors(const LocalOrdinal& /* el */) const
   { static std::vector<int> empty; return empty; }
   //! Not used.
   virtual bool hasAssociatedNeighbors() const { return false; }

private:
   int procRank_;
   
   Teuchos::RCP<ConnCallback> callback_;
   std::map<std::string,std::vector<int> > elements_; // local element IDs
   std::vector<std::vector<GlobalOrdinal> > connectivity_;
};

} // end unit test
} // end panzer

#endif
