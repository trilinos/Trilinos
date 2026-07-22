// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __CartesianConnManager_hpp__
#define __CartesianConnManager_hpp__

#include <vector>

#include "Teuchos_DefaultMpiComm.hpp"

#include "Panzer_ConnManager.hpp"

#include "Shards_CellTopology.hpp"
#include "Shards_BasicTopologies.hpp"


namespace panzer::unit_test {

/** This class is used to construct the topology of a mesh that is cartesian. The coordinates
  * and ordinals follow the right hand rule in both 2D and 3D. Thus the indices look like
  *
  *      k(z)
  *      |
  *      |___ i(x)
  *      /
  *     /
  *    j(y)
  *
  * and in 2D
  *
  *      j(y)
  *      |
  *      |___ i(x)
  *
  * You can also construct multiple element blocks. In 3D they will be named "eblock-i_j_k"
  * and in 2D "eblock-i_j".
  *
  * In 3D, the mesh can either have Hexahedral (brick) elements, or tetrahedral elements,
  * where each brick element is split into 6 tetrahedra
  * In 2D, the mesh can either have Quadrilateral (brick) elements, or triangular elements,
  * where each brick element is split into 2 triangles
  */
class CartesianConnManager : public virtual panzer::ConnManager {
public:
   using connectivity_entries_host_view_type = Kokkos::View<GlobalOrdinal*, Kokkos::HostSpace>;

public:
   // A utility structure for storing triplet indices
   template <typename T>
   struct Triplet {
      Triplet() : x(-1),y(-1),z(-1) {}
      Triplet(T a,T b,T c) : x(a),y(b),z(c) {}
      T x, y, z;
   };

   CartesianConnManager() {isInitialized = false, areBoundarySidesBuilt = false;}

   ~CartesianConnManager() = default;

   /** Initialize a 2D topology.
     *
     * \param[in] comm Communicator
     * \param[in] nx Number of elements in the x direction per block
     * \param[in] ny Number of elements in the y direction per block
     * \param[in] px Number of processors in the x direction
     * \param[in] py Number of processors in the y direction
     * \param[in] bx Number of blocks in the x direction
     * \param[in] by Number of blocks in the y direction
     * \param[in] elemTopo Topology of the mesh element (either Quadrilateral<4> or Triangular<3>)
     */
   void initialize(const Teuchos::MpiComm<int> & comm,GlobalOrdinal nx, GlobalOrdinal ny,
                   int px, int py,
                   int bx, int by,
                   const shards::CellTopology elemTopo=shards::getCellTopologyData<shards::Quadrilateral<4>>());

   /** Initialize a 3D topology.
     *
     * \param[in] comm Communicator
     * \param[in] nx Number of elements in the x direction per block
     * \param[in] ny Number of elements in the y direction per block
     * \param[in] nz Number of elements in the z direction per block
     * \param[in] px Number of processors in the x direction
     * \param[in] py Number of processors in the y direction
     * \param[in] pz Number of processors in the z direction
     * \param[in] bx Number of blocks in the x direction
     * \param[in] by Number of blocks in the y direction
     * \param[in] bz Number of blocks in the z direction
     * \param[in] elemTopo Topology of the mesh element (either Hexahedron<8> or Tetrahedron<4>)
     */
   void initialize(const Teuchos::MpiComm<int> & comm,GlobalOrdinal nx, GlobalOrdinal ny, GlobalOrdinal nz,
                   int px, int py, int pz,
                   int bx, int by, int bz,
                   const shards::CellTopology elemTopo=shards::getCellTopologyData<shards::Hexahedron<8>>());

   /** Get the triplet corresponding to the number of brick elements (Hexahedra or Quadrilaterals) this processor owns.
     * Publically exposed for testing.
     */
   Triplet<GlobalOrdinal> getMyBrickElementsTriplet() const { return myBrickElements_; }

   /** Get the triplet corresponding to the offset of brick elements beginning with what this processor owns.
     * Publically exposed for testing.
     */
   Triplet<GlobalOrdinal> getMyBrickOffsetTriplet() const { return myBrickOffset_; }

   /** This computes the local index of a brick element, from the global triplet. If the index isn't
     * on this processor a -1 is returned.
     */
   LocalOrdinal computeLocalBrickElementIndex(const Triplet<GlobalOrdinal> & brickElement) const;

   /** Compute the global index of a brick element from a global triplet.
     */
   GlobalOrdinal computeGlobalBrickElementIndex(const Triplet<GlobalOrdinal> & brickElement) const;

   /** Utility function for computing the processor i,j,k ranking. Publically exposed for testing.
     */
   static Triplet<int> computeMyRankTriplet(int myRank,int dim,const Triplet<int> & procs);

   /** Utility function for computing a global triplet of a brick from a local brick element index
     * Publically exposed for testing.
     */
   static Triplet<GlobalOrdinal> computeLocalBrickElementGlobalTriplet(int index,const Triplet<GlobalOrdinal> & myBrickElements,
                                                                       const Triplet<GlobalOrdinal> & myBrickOffset);

   /** Compute the global brick index from a global brick triplet.
     */
   static GlobalOrdinal computeGlobalBrickElementIndex(const Triplet<GlobalOrdinal> & brickElement,
                                                  const Triplet<GlobalOrdinal> & brickShape);

   /** Compute the global brick index for a triplet.
     */
   static LocalOrdinal computeLocalBrickElementIndex(const Triplet<GlobalOrdinal> & brickElement,
                                                const Triplet<GlobalOrdinal> & myBrickElements,
                                                const Triplet<GlobalOrdinal> & myBrickOffset);

   /** Number of sub elements per brick (e.g. 6 tets per Hex, 2 triangles per quad)
     */
   int numSubElemsPerBrickElement(){return numSubElemsPerBrickElement_;}

   /** Local node Ids of a brick from the local Id of the sub element
     */
   int getLocalBrickNodeFromSubElemNode(int subElement, int node) {
     return subElemToBrickElementNodesMap_[subElement][node];
   }

   // The next functions are all pure virtual from ConnManager
   ////////////////////////////////////////////////////////////////////////////////

   /** Tell the connection manager to build the connectivity assuming
     * a particular field pattern.
     *
     * \param[in] fp Field pattern to build connectivity for
     */
   virtual void buildConnectivity(const FieldPattern & fp);

   /** Build a clone of this connection manager, without any assumptions
     * about the required connectivity (e.g. <code>buildConnectivity</code>
     * has never been called).
     */
   virtual Teuchos::RCP<panzer::ConnManager> noConnectivityClone() const;

   /** Get ID connectivity for a particular element
     *
     * \param[in] localElmtId Local element ID
     *
     * \returns Pointer to beginning of indices, with total size
     *          equal to <code>getConnectivitySize(localElmtId)</code>
     */
   virtual const GlobalOrdinal * getConnectivity(LocalOrdinal localElmtId) const { return connectivity_.data() + localElmtId*numIds_; }

   /** How many mesh IDs are associated with this element?
     *
     * \param[in] localElmtId Local element ID
     *
     * \returns Number of mesh IDs that are associated with this element.
     */
   virtual LocalOrdinal getConnectivitySize(LocalOrdinal localElmtId) const
   { return numIds_; }


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
   virtual const std::vector<LocalOrdinal> & getElementBlock(const std::string & blockId) const;

   virtual const std::vector<LocalOrdinal> & getNeighborElementBlock(const std::string & /* s */) const
   { return emptyVector_; }

   virtual const std::vector<LocalOrdinal> & getAssociatedNeighbors(const LocalOrdinal& /* el */) const
   { return emptyVector_; }

   virtual bool hasAssociatedNeighbors() const
   { return false; }

   /** For each of the 6 faces (4 edges in 2d) of the brick, builds a vector sides on that brick face (edge in 2d).
       Each side is described with a pair of integers: the cell containing the side, and side ordinal to the cell
   */
   void buildLocalBoundarySides();

   /** Returns the boundary sides associated to one of the 6 faces (4 edges in 2D) of the brick.
       Each side is described with pair of integers, the cell containing the side, and side ordinal to the cell   
       \param[in] brickSide, the face (or edge in 2D) of the brick, whose sides are returned
   */
   virtual const std::vector<std::pair<LocalOrdinal,int>> & getBoundarySides(int brickSide) const;

private:

   // For each element block allocate owned local elements in that block
   void buildLocalElements();

   // For const returns associated with neighbor access (not implemented)
   const std::vector<LocalOrdinal> emptyVector_;

   // Update the connectivity vector with the field pattern. The connectivity is specified
   // here using the i,j,k index of the local element. Also, this is where the ordering
   // of a cell is embedded into the system.
   void updateConnectivity_2d(const panzer::FieldPattern & fp,int subcellDim,int localElementId,
                              LocalOrdinal& offset) const;

   // Update the connectivity vector with the field pattern. The connectivity is specified
   // here using the i,j,k index of the local element. Also, this is where the ordering
   // of a cell is embedded into the system.
   void updateConnectivity_3d(const panzer::FieldPattern & fp,int subcellDim,int localElementId,
                              LocalOrdinal& offset) const;

   int numProc_;
   int myRank_;
   Triplet<int> myRankIndex_;

   int dim_;
   Triplet<GlobalOrdinal> totalBrickElements_; // over all blocks and processors how many elements
   Triplet<GlobalOrdinal> brickElements_; // per block element counts
   Triplet<int> processors_;         // full mesh
   Triplet<int> blocks_;             // number of element blocks

   Triplet<GlobalOrdinal> myBrickElements_;
   Triplet<GlobalOrdinal> myBrickOffset_;

   std::map<std::string,std::vector<LocalOrdinal> > localElements_;

   // element vector to connectivity
   connectivity_entries_host_view_type connectivity_;
   int numIds_;

   GlobalOrdinal totalNodes_;
   GlobalOrdinal totalEdges_;
   GlobalOrdinal totalFaces_;
   GlobalOrdinal totalElements_;
   int numSubElemsPerBrickElement_;
   std::vector<int> numSubSidesPerBrickSide_;

   // element to brick map

   //note: nodes of a sub element are a subset of nodes of a brick element
   std::vector<std::vector<int> > subElemToBrickElementNodesMap_;

   //note: the edges of a sub elements can include edges of the brick element,
   //      diagonal edges of brick faces, and in 3D, diagonal edges of the brick element
   std::vector<std::vector<int> > subElemToBrickElementEdgesMap_;

   //note: the faces of a sub elements can include faces of the brick element,
   //      and in 3D, faces internal to the brick element
   std::vector<std::vector<int> > subElemToBrickElementFacesMap_;

   //  subElemToBrickElementEdgesMap_ in 2D,  subElemToBrickElementFacesMap_ in 3D
   std::vector<std::vector<int> > subElemToBrickElementSidesMap_;

   std::vector<std::vector<std::pair<LocalOrdinal,int>>> boundarySides_;
   bool isInitialized, areBoundarySidesBuilt;

   shards::CellTopology elemTopology_;


};

} // namespace panzer::unit_test

#endif // __CartesianConnManager_hpp__
