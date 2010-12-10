#ifndef stk_mesh_FEMInterface_hpp
#define stk_mesh_FEMInterface_hpp

#include <stk_mesh/fem/CellTopology.hpp>
#include <stk_mesh/base/Types.hpp>

namespace stk {
namespace mesh {
namespace fem {

static const EntityRank NODE_RANK = 0;
static const EntityRank INVALID_RANK = stk::mesh::InvalidEntityRank;

///
/// FEMInterface defines the interface to acquire finite element data.
/// An object with this interface is attached to the MetaData via
/// attributes and can be obtained from the metadata using the
/// stk::mesh::fem::get_fem_interface(meta_data) function.
///
/// The FEMInterface maintains three type of information, the spatial
/// dimension, the cell topology to entity rank mapping, and the
/// association of cell topology to part.
///
/// Spatial Dimension and Entity Ranks:
///
///   Spatial dimension defines the dimensional 
///   Discretization entities: Node, Edge, Side, Element
///   Node - Node rank must be 0, because the relationss from rank element rank to node rank are indexes by node not mapped, Point
///   Edge - Node rank + 1, Tangent, only for 2d and 3d
///   Side - Element rank - 1, Normal
///   Element - Fundamental unit for discretization
///
/// Topological entities: Node, Edge, Face, Solid;  2d: Node, Edge, Face;  3d: Node, Edge, Face, Solid
/// 
/// Topological dimension is not the same as entity rank.  For example
/// a particle has topological dimension of 1, surprised?  Since a
/// particle has a node which is topological rank 0, and shards defines
/// a topological dimension as one greater than any contained cell.
/// Therefore, particle has topological rank 1.  However, particle is
/// an element in the finite element model, and therefore has.
///
///
///
/// Coordinate system may need to be included.
///

/** 
 * @brief <b>FEMInterface</b> defines the interface for a finite
 * element model to the STK mesh.
 * 
 */
class FEMInterface
{
public:
  virtual ~FEMInterface()
  {}


  /** 
   * @brief <b>set_spatial_dimension</b> sets the spatial dimension of
   * the finite element model.  The spatial dimension is utilized to
   * determine the entity ranks.
   * 
   * @param spatial_dimension 
   */
  virtual void set_spatial_dimension(size_t spatial_dimension) = 0;
  
  virtual void register_cell_topology(const CellTopology cell_topology, EntityRank entity_rank) = 0;

  virtual void set_cell_topology(Part &part, CellTopology cell_topology) = 0;

  /** 
   * @brief <b>get_spatial_dimension</b> 
   * 
   * 
   * @return 
   */
  virtual size_t get_spatial_dimension() const = 0;
  
  virtual CellTopology get_cell_topology(const Part &part) const = 0;

  
  virtual EntityRank get_entity_rank(const CellTopology cell_topology) const = 0;


  virtual Part &get_part(const CellTopology cell_topology) const = 0;
};


void set_fem_interface(MetaData &meta_data, FEMInterface *fem_interface);

// These are order get_fem_interface() O(log n), where n is the number of attribute defined on the mesh.
FEMInterface &get_fem_interface(const MetaData &meta_data);
FEMInterface &get_fem_interface(const BulkData &bulk_data);
FEMInterface &get_fem_interface(const Part &part);
FEMInterface &get_fem_interface(const Bucket &bucket);
FEMInterface &get_fem_interface(const Entity &entity);

void register_cell_topology(MetaData &meta_data, const CellTopology cell_topology, EntityRank entity_rank);

void set_cell_topology(Part & part, CellTopology cell_topology);

template< class TopologyTraits >
void set_cell_topology(Part & p) {
  return set_cell_topology(p, shards::getCellTopologyData<TopologyTraits>());
}


CellTopology get_cell_topology(const Part & part);
CellTopology get_cell_topology(const Bucket & bucket);
CellTopology get_cell_topology(const Entity & entity);


EntityRank get_entity_rank(const MetaData &meta_data, CellTopology cell_topology);
Part &get_part(const MetaData &meta_data, const CellTopology cell_topology);


/// Rank of nodes (always zero)      
inline EntityRank node_rank(size_t spatial_dimension) {
  return NODE_RANK;
}

/// Rank of edges (1 for 2D and 3D)
inline EntityRank edge_rank(size_t spatial_dimension) {
  return spatial_dimension > 1 ? 1 : InvalidEntityRank;
}

/// Rank of faces (2 for 3D)
inline EntityRank face_rank(size_t spatial_dimension) {
  return spatial_dimension > 2 ? 2 : InvalidEntityRank;
}

/// Rank of sides (1 for 2D, 2 for 3D)
inline EntityRank side_rank(size_t spatial_dimension) {
  return spatial_dimension - 1;
}

/// Rank of elements (spatial_dimension)
inline EntityRank element_rank(size_t spatial_dimension) {
  return spatial_dimension;
}

inline EntityRank node_rank(const FEMInterface &fem) {
  return node_rank(fem.get_spatial_dimension());
}

inline EntityRank edge_rank(const FEMInterface &fem) {
  return edge_rank(fem.get_spatial_dimension());
}

inline EntityRank face_rank(const FEMInterface &fem) {
  return face_rank(fem.get_spatial_dimension());
}

inline EntityRank side_rank(const FEMInterface &fem) {
  return side_rank(fem.get_spatial_dimension());
}

inline EntityRank element_rank(const FEMInterface &fem) {
  return element_rank(fem.get_spatial_dimension());
}

} // namespace fem

/** 
 * @brief <b>declare_part</b> 
 * 
 * @param meta_data 
 * @param name 
 * 
 * @return 
 */
Part &declare_part(MetaData &meta_data, const std::string &name);
Part &declare_part(MetaData &meta_data, const std::string &name, EntityRank entity_rank);
Part &declare_part(MetaData &meta_data, const PartVector &part_vector);

Part &declare_part(MetaData &meta_data, const std::string &name, fem::CellTopology cell_topology);

template< class Top >
Part &declare_part(MetaData &meta_data, const std::string &name) {
  return declare_part(meta_data, name, shards::getCellTopologyData<Top>());
}

} // namespace mesh
} // namespace stk


/** 11/8/10 Design meeting to decide on FEM layer interfaces
 * Attendance:  Carter, Alan, Greg, Jim, Dave, Todd, Dan
 * Issues to resolve:
 *   * EntityRanks defined by discretization names vs topological names
 *   * discretization names:  vertex, side, element
 *   * topological names:  node, edge, face, solid (3D)
 *   Element Rank depends on (equals) spatial dimension  [ essential complexity ]
 *   Entity Rank is not 1:1 with respect to topological dimension
 *
 * Three fundamental components of FEM: (agreed on by all in attendance)
 *   1.  Cell Topology
 *   2.  Entity Ranks
 *   3.  Spatial Dimension
 *
 * Where do coordinates go?  In Fields...
 * Where does axi-symmetric information go?
 *
 * Issue of attribute mechanism on MetaData giving const only access:
 * Carter said this is an arbitrary restriction and it should be removed now that we have a need.
 * All in agreement that attributes should be const/nonconst accessible.
 *
 * Issue of Pure Abstract Class design for FEM layer.
 *   * Initial idea was to allow applications store mapping of parts to cell
 *     topologies in different ways.  As a map on MetaData or as attributes on
 *     parts.  Alan suggested we store them in an std::vector<cell_topology*>
 *     with one index into vector for every Part (ordinal).  This would make
 *     lookup constant which is about as fast as you can make it. 
 *   * What about a template for this variation?
 *     Concern that this would move implementation into headers and require
 *     templatization of stk_io.  
 *   * How does refinement topology get added into the mix?  Is this an
 *     extension of FEM or a different plugin entirely?  Carter said a
 *     complaint in the framework has been "Why do I have to see refinement
 *     topology when I don't need it?"  This is an argument for a different
 *     plug-in.
 *   * Pro/Con for abstract base class:
 *     Pro:  specialization allowed easily (axi-symmetric information,
 *       refinement topology, alternative implementation of cell_topology
 *       lookups, etc.), user-defined ranks, Java style "interface" API
 *     Con:  unessential complexity added, "abstract", no other examples in
 *       stk_mesh, no need for it now, virtual symbol table, no inline calls,
 *       non-obvious performance concerns (e.g. fem.element_rank() looks
 *       light-weight, but may not be).
 *   * Favor composition over inheritance.  Where does this come into play here?
 *   * What is lesson learned from Framework use of inheritance?  Dave & Carter
 *     discussed this and it appeared that in the framework the base classes
 *     were not pure-virtual "interfaces" and the resulting extensive
 *     inheritance caused great problems.  The idea here is that a pure-virtual
 *     base class is safer.  There is also the idea that this class should only
 *     be inherited from once.
 *   * Voting on abstract vs. concrete implementation:
 *     Abstract:  Todd, Alan, Greg, Dan.
 *     Concrete:  Carter, Jim.
 *     Abstain:  Dave.
 *     Greg comment:  concrete classes that need specializations lead to work-arounds.
 *     Dan comment:  inner loop performance issues are easy to resolve.
 *   * We also talked about how to access the element_rank EntityRank
 *     information.  Currently Dave has some free functions that use MetaData
 *     to get the attribute for the FEM interface and it gets spatial_dimension
 *     off of that class and computes the element rank.  Alan suggested that we
 *     move the element_rank function to the FEM interface class instead and
 *     then we will use a block of code like:
 *       FEMInterface & fem = get_fem(part);
 *       fem.element_rank();
 *     In this way, the DefaultFEM class would set the internal element_rank
 *     value at construction with the spatial dimension and then the accessor
 *     would simply return the value.
 *
 *  Issue of constructor for MetaData:
 *    The main issue is that MetaData requires a list of entity rank names at
 *      construction and it uses this internally to check invariants related to
 *      the maximum entity rank and it is used in BulkData for both the maximum
 *      entity rank and the actual entity rank names in error messages.  When
 *      reading from an exodus file, you need to pass MetaData into stk_io but
 *      you don't know if the problem is 2D or 3D yet and this is necessary to
 *      create the correct list of entity rank names.
 *    Dave added a default MetaData constructor and an additional state for
 *      before/after setting entity rank names.  This allowed him to create a
 *      MetaData and pass it to stk_io to fill up with MetaData and entity rank
 *      names based on spatial dimension from the exodus file.
 *    Dan/Greg agree that stk_io needs better code for interacting with exodus,
 *      MetaData, and BulkData.  Currenlty, there is a use_case_mesh class that
 *      has a nice set of features but because its in a use case its not really
 *      part of the library, yet it is being used heavily clients of stk_mesh to
 *      read in exodus files.
 *    Alan suggested that we should have a factory in stk_io which creates
 *      MetaData and passes it back through a shared pointer.  
 *    We had a discussion about whether entity ranks must be contiguous or not.
 *      It appears that base is okay with non-contiguous entty ranks, but FEM
 *      requires them to be contiguous.
 *    We then had a discussion about who is a client of base and not FEM.  It
 *      appears that arbitrary polyhedra is building on top of FEM and
 *      peridynamics is also building on top of FEM.  Those were the two big
 *      examples that were supposed to build on top of base directly and not
 *      use the FEM layer.  So at this point we don't have even a single
 *      use-case for base all by itself.  This raised the question of whether
 *      we should just move the FEM interface calls into MetaData directly and
 *      access an internal FEM class to store all the mappings.  This part of
 *      the conversation occured at the end of 2.5 hours of discussion and
 *      trailed off with no conclusion.
 *    Dave is going to implement an additional state in MetaData for the
 *      before/after setting entity rank names and implement checks for this
 *      state in the appropriate member functions in MetaData.  This will allow
 *      his current implementation with stk_io to move forward easily.
 *    
 **/

#endif // stk_mesh_FEMInterface_hpp
