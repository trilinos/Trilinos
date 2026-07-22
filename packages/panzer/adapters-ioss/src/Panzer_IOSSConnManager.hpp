// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Panzer_IOSSConnManager_hpp__
#define __Panzer_IOSSConnManager_hpp__

#include <string>
#include <vector>
#include <array>

// Teuchos includes
#include "Teuchos_RCP.hpp"

// Kokkos includes
#include "Kokkos_DynRankView.hpp"

// Phalanx includes
#include "Phalanx_KokkosDeviceTypes.hpp"

// Panzer includes
#include "Panzer_ConnManager.hpp"
#include "Panzer_FieldPattern.hpp"
#include "Panzer_IntrepidFieldPattern.hpp"

// Ioss includes
#include "Ioss_CodeTypes.h"
#include "Ioss_ParallelUtils.h"
#include "Ioss_DBUsage.h"
#include "Ioss_PropertyManager.h"
#include "Ioss_Property.h"
#include "Ioss_IOFactory.h"
#include "Ioss_DatabaseIO.h"
#include "Ioss_GroupingEntity.h"
#include "Ioss_Region.h"
#include "Ioss_EntityBlock.h"
#include "Ioss_ElementTopology.h"
#include "Ioss_NodeBlock.h"
#include "Ioss_ElementBlock.h"

namespace panzer_ioss {

class IOSSConnManager : public panzer::ConnManager {

public:
   using LocalOrdinal = panzer::ConnManager::LocalOrdinal;
   using GlobalOrdinal = panzer::ConnManager::GlobalOrdinal;

   typedef typename std::vector<Ioss::NodeBlock*> NodeBlockContainer;
   typedef typename std::vector<Ioss::ElementBlock*> ElementBlockContainer;

   IOSSConnManager(Ioss::DatabaseIO * iossMeshDB);

   virtual ~IOSSConnManager() {}

   /** Tell the connection manager to build the connectivity assuming
     * a particular field pattern.
     *
     * \param[in] fp Field pattern to build connectivity for
     */
   virtual void buildConnectivity(const panzer::FieldPattern & fp);

   /** Build a clone of this connection manager, without any assumptions
     * about the required connectivity (e.g. <code>buildConnectivity</code>
     * has never been called).
     * This default version assumes an exodus-type database and a property
     * manager containing the single property DECOMPOSITION_METHOD=LINEAR
     */
   virtual Teuchos::RCP<panzer::ConnManager> noConnectivityClone() const {
     std::string type = "exodus";
     Ioss::PropertyManager properties;
     Ioss::Property decomp_prop("DECOMPOSITION_METHOD", "LINEAR");
     properties.add(decomp_prop);
     return noConnectivityClone(type, properties);
   };

   /** Build a clone of this connection manager, without any assumptions
        * about the required connectivity (e.g. <code>buildConnectivity</code>
        * has never been called).
        */
   virtual Teuchos::RCP<panzer::ConnManager> noConnectivityClone(std::string & type, Ioss::PropertyManager & properties) const;

   /** Get ID connectivity for a particular element
     *
     * \param[in] localElmtId Local element ID
     *
     * \returns Pointer to beginning of indices, with total size
     *          equal to <code>getConnectivitySize(localElmtId)</code>
     */
   virtual const GlobalOrdinal * getConnectivity(LocalOrdinal localElmtId) const
    { return &connectivity_[elmtLidToConn_[localElmtId]]; }

    /** Get ID connectivity for a particular element
      *
      * \param[in] localElmtId Local element ID
      *
      * \returns Pointer to beginning of indices, with total size
      *          equal to <code>getConnectivitySize(localElmtId)</code>
      */
    virtual GlobalOrdinal * getConnectivity(LocalOrdinal localElmtId)
    { return &connectivity_[elmtLidToConn_[localElmtId]]; }

    /** How many mesh IDs are associated with this element?
      *
      * \param[in] localElmtId Local element ID
      *
      * \returns Number of mesh IDs that are associated with this element.
      */
    virtual LocalOrdinal getConnectivitySize(LocalOrdinal localElmtId) const
    { return connSize_[localElmtId]; }

    /** Get the block ID for a particular element.
      *
      * \param[in] localElmtId Local element ID
      */
    virtual std::string getBlockId(LocalOrdinal localElmtId) const;

    /** How many element blocks in this mesh?
      */
    virtual std::size_t numElementBlocks() const {
       return iossElementBlocks_.size();
    };


    /** Get block IDs from IOSS mesh object
      */
    virtual void getElementBlockIds(std::vector<std::string> & elementBlockIds) const {
        elementBlockIds.clear();
        for (Ioss::ElementBlock * iossElementBlock : iossElementBlocks_)
           elementBlockIds.push_back(iossElementBlock->name());
    };


    /** What are the cellTopologies linked to element blocks in this connection manager?
     */
    virtual void getElementBlockTopologies(std::vector<shards::CellTopology> & elementBlockTopologies) const {
      elementBlockTopologies.clear();
      for (Ioss::ElementBlock * iossElementBlock : iossElementBlocks_) {
        elementBlockTopologies.push_back((iossToShardsTopology_.find(iossElementBlock->topology()->name()))->second);
      }
    };

    /** Get the local element IDs for a paricular element
      * block. These are only the owned element ids.
      *
      * \param[in] blockIndex Block Index
      *
      * \returns Vector of local element IDs.
      */
    virtual const std::vector<LocalOrdinal> & getElementBlock(const std::string & blockId) const {
        return *(elementBlocks_.find(blockId)->second);
    };

    virtual const std::vector<LocalOrdinal> & getNeighborElementBlock(const std::string & blockId) const
    { return *(neighborElementBlocks_.find(blockId)->second); }

    /** Get the coordinates (with local cell ids) for a specified element block and field pattern.
      *
      * \param[in] blockId Block containing the cells
      * \param[in] coordProvider Field pattern that builds the coordinates
      * \param[out] localCellIds Local cell Ids (indices)
      * \param[out] Resizable field container that contains the coordinates
      *             of the points on exit.
      */
    virtual void getDofCoords(const std::string & blockId,
                              const panzer::Intrepid2FieldPattern & coordProvider,
                              std::vector<LocalOrdinal> & localCellIds,
                              Kokkos::DynRankView<double,PHX::Device> & points) const;

     /** Get IOSS mesh database that this connection manager is built on.
       */
     Ioss::DatabaseIO *  getIossMeshDB() const
     { return iossMeshDB_; }

     /** How many elements are owned by this processor. Further,
       * the ordering of the local ids is such that the first
       * <code>getOwnedElementCount()</code> elements are owned
       * by this processor.
       */
     std::size_t getOwnedElementCount() const
     { return ownedElementCount_; }

     /** Before calling buildConnectivity, provide sideset IDs from which to
       * extract associated elements.
       */
     void associateElementsInSideset(const std::string sideset_id);

     /** After calling <code>buildConnectivity</code>, optionally check which
       * sidesets yielded no element associations in this communicator. This is a
       * parallel operation. In many applications, the outcome indicating
       * correctness is that the returned vector is empty.
       */
     std::vector<std::string> checkAssociateElementsInSidesets(const Teuchos::Comm<int>& /* comm */) const {
       std::vector<std::string> placeholder;
       return placeholder;
     };

     /** Get elements, if any, associated with <code>el</code>, excluding
       * <code>el</code> itself.
       */
     virtual const std::vector<LocalOrdinal>& getAssociatedNeighbors(const LocalOrdinal& /* el */) const {
       return placeholder_;
     };

     /** Return whether getAssociatedNeighbors will return true for at least one
       * input. Default implementation returns false.
       */
     virtual bool hasAssociatedNeighbors() const {
       return false;
     };


protected:

   /** Create a mapping from IOSS topology to Shards topology.
    */
   void createTopologyMapping();

   /** Get element info from IOSS database and build a local element mapping.
     */
   void buildLocalElementMapping();

   /** Erase element mapping info.
    */
   void clearLocalElementMapping();

   /** Build mapping from edge and face corner node numbers to unique global edge and face numbers.
    */
   void buildEdgeFaceCornerNodeMapping();

   void buildOffsetsAndIdCounts(const panzer::FieldPattern & fp,
                                  LocalOrdinal & nodeIdCnt, LocalOrdinal & edgeIdCnt,
                                  LocalOrdinal & faceIdCnt, LocalOrdinal & cellIdCnt,
                                  GlobalOrdinal & nodeOffset, GlobalOrdinal & edgeOffset,
                                  GlobalOrdinal & faceOffset, GlobalOrdinal & cellOffset) const;

   LocalOrdinal addSubcellConnectivities(const panzer::FieldPattern & fp, std::string & blockId,
                                         std::size_t elmtIdInBlock, std::size_t elmtLid,
                                         unsigned subcellRank, LocalOrdinal idCnt,GlobalOrdinal offset);

   /* Determine whether a FieldPattern object is compatible with the Ioss::ElementTopology
    * of every block in the mesh.
    *
    *  \param[in] fp The field pattern
    *
    *  \returns true if the field pattern is compatible with all element topologies in the mesh.
    *
    */
   bool compatibleTopology(const panzer::FieldPattern & fp) const;


   Ioss::DatabaseIO * iossMeshDB_;
   Teuchos::RCP<Ioss::Region> iossRegion_;
   NodeBlockContainer iossNodeBlocks_;
   ElementBlockContainer iossElementBlocks_;
   std::map<std::string,Teuchos::RCP<std::vector<GlobalOrdinal>>> iossConnectivity_;
   std::map<std::string,shards::CellTopology> iossToShardsTopology_;
   std::map<std::string, const Ioss::ElementTopology *> iossElementBlockTopologies_;

   // element block information
   std::map<std::string,Teuchos::RCP<std::vector<LocalOrdinal> > > elementBlocks_;
   std::map<std::string,Teuchos::RCP<std::vector<LocalOrdinal> > > neighborElementBlocks_;

   static const int MAX_SUBCELL_CORNER_NODES_ = 6;
   GlobalOrdinal numUniqueEdges_, numUniqueFaces_;
   std::map<std::array<GlobalOrdinal,MAX_SUBCELL_CORNER_NODES_>,GlobalOrdinal> edgeNodeToEdgeMap_;
   std::map<std::array<GlobalOrdinal,MAX_SUBCELL_CORNER_NODES_>,GlobalOrdinal> faceNodeToFaceMap_;
   std::map<std::string,Teuchos::RCP<std::vector<std::vector<std::array<GlobalOrdinal,MAX_SUBCELL_CORNER_NODES_>>>>> elementEdgeNodes_;
   std::map<std::string,Teuchos::RCP<std::vector<std::vector<std::array<GlobalOrdinal,MAX_SUBCELL_CORNER_NODES_>>>>> elementFaceNodes_;

   std::vector<GlobalOrdinal> elmtLidToGid_; // element LID to GID map.
   std::vector<LocalOrdinal> elmtLidToConn_; // element LID to starting index in connectivity_
   std::vector<LocalOrdinal> connSize_; // Number of mesh IDs that are associated with each element
   std::vector<GlobalOrdinal> connectivity_; // Connectivity

   std::size_t ownedElementCount_;

   std::vector<std::string> sidesetsToAssociate_;
   std::vector<bool> sidesetYieldedAssociations_;
   std::vector<std::vector<LocalOrdinal> > elmtToAssociatedElmts_;

   std::vector<LocalOrdinal> placeholder_;
};

}

#endif
