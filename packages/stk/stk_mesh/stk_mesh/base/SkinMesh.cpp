#include <map>
#include <set>
#include <algorithm>

#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/BulkModification.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/base/Relation.hpp>
#include <stk_mesh/base/BoundaryAnalysis.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>
#include <stk_mesh/base/SkinMesh.hpp>

#include <stk_util/parallel/ParallelComm.hpp>

namespace stk {
namespace mesh {

namespace {

typedef std::pair< const CellTopologyData *, EntityVector > SideKey;
typedef std::vector< EntitySideComponent >                  SideVector;
typedef std::map< SideKey, SideVector>                      BoundaryMap;

//Comparator class to sort EntitySideComponent
//first by entity identifier and then by side_ordinal
class EntitySideComponentLess {
  public:
    EntitySideComponentLess(const BulkData &mesh) : m_mesh(mesh) { }

    bool operator () (const EntitySideComponent & lhs, const EntitySideComponent &rhs) const
    {
      const EntityId lhs_elem_id = m_mesh.identifier(lhs.entity);
      const EntityId rhs_elem_id = m_mesh.identifier(rhs.entity);

      return (lhs_elem_id != rhs_elem_id) ?
        (lhs_elem_id < rhs_elem_id) :
        (lhs.side_ordinal < rhs.side_ordinal);
    }
  private:
    EntitySideComponentLess();
    const BulkData &m_mesh;
};

//Convience class to help with communication.
//sort first by element identifier and then by side_ordinal
class ElementIdSide {
  public:
    EntityId        elem_id;
    RelationIdentifier side_ordinal;

    ElementIdSide() :
      elem_id(0), side_ordinal(0) {}

    ElementIdSide( EntityId id, RelationIdentifier ordinal) :
      elem_id(id), side_ordinal(ordinal) {}

    ElementIdSide(const BulkData &mesh, const EntitySideComponent & esc) :
      elem_id(mesh.identifier(esc.entity)), side_ordinal(esc.side_ordinal) {}

    ElementIdSide & operator = ( const ElementIdSide & rhs) {
      elem_id      = rhs.elem_id;
      side_ordinal = rhs.side_ordinal;
      return *this;
    }

    //compare first by elem_id and then by side_ordinal
    bool operator < ( const ElementIdSide & rhs) const {
      const ElementIdSide & lhs = *this;

      return (lhs.elem_id != rhs.elem_id) ?
        (lhs.elem_id < rhs.elem_id) :
        (lhs.side_ordinal < rhs.side_ordinal);
    }
};

//a reverse map used in unpacking to determine what side needs
//to be created
typedef std::map< ElementIdSide, SideKey>   ReverseBoundaryMap;

//a convience class to help with packing the comm buffer
class SideCommHelper {
  public:
    unsigned       proc_to;
    ElementIdSide  creating_elem_id_side; //used as a look up in the reverse boundary map
    EntityId       generated_side_id;

    SideCommHelper() :
      proc_to(0), creating_elem_id_side(), generated_side_id(0) {};

    SideCommHelper( unsigned p, const ElementIdSide & creating_eid, const EntityId side_id) :
      proc_to(p), creating_elem_id_side(creating_eid), generated_side_id(side_id) {}

    SideCommHelper( const SideCommHelper & sch) :
      proc_to(sch.proc_to),
      creating_elem_id_side(sch.creating_elem_id_side),
      generated_side_id(sch.generated_side_id)
  {}

    SideCommHelper & operator = (const SideCommHelper & rhs) {
      proc_to                = rhs.proc_to;
      creating_elem_id_side  = rhs.creating_elem_id_side;
      generated_side_id      = rhs.generated_side_id;
      return *this;
    }

    //compare first by proc_to, then by creating elem_id
    bool operator < ( const SideCommHelper & rhs) const {
      const SideCommHelper & lhs = *this;

      return (lhs.proc_to != rhs.proc_to) ?
        (lhs.proc_to < rhs.proc_to) :
        (lhs.creating_elem_id_side < rhs.creating_elem_id_side);
    }
};

// Use permutation that starts with lowest entity id
void ensure_consistent_order(const stk::mesh::BulkData &mesh, EntityVector & side_entities)
{
  ThrowRequire( !side_entities.empty() );

  EntityId lowest_id = mesh.identifier(side_entities.front());
  unsigned idx_of_lowest_id = 0;

  for (unsigned idx = 1; idx < side_entities.size(); ++idx) {
    EntityId curr_id = mesh.identifier(side_entities[idx]);
    if (curr_id < lowest_id) {
      idx_of_lowest_id = idx;
      lowest_id = curr_id;
    }
  }

  if (idx_of_lowest_id != 0) {
    std::rotate(side_entities.begin(),
                side_entities.begin() + idx_of_lowest_id,
                side_entities.end());
  }
}

// populate the side_map with 'owned' sides that need to be created
//
// a side needs to be created if the outside is NULL and the element
// does not have a current relation to a side for the side_ordinal
void add_owned_sides_to_map(
    const BulkData & mesh,
    const EntityRank element_rank,
    const EntitySideVector & boundary,
    BoundaryMap & side_map)
{
  const EntityRank side_rank = mesh.mesh_meta_data().side_rank();
  for (stk::mesh::EntitySideVector::const_iterator itr = boundary.begin();
      itr != boundary.end(); ++itr) {
    const EntitySideComponent & inside = itr->inside;
    const EntitySideComponent & outside = itr->outside;
    const RelationIdentifier side_ordinal = inside.side_ordinal;
    const Entity inside_entity = inside.entity;

    if ( mesh.parallel_owner_rank(inside_entity) == mesh.parallel_rank() &&
         !mesh.is_valid(outside.entity) )
    {
      // search through existing sides
      const int num_sides = mesh.num_connectivity(inside_entity, side_rank);
      Entity const * existing_side = mesh.begin(inside_entity, side_rank);
      ConnectivityOrdinal const *existing_ords = mesh.begin_ordinals(inside_entity, side_rank);
      // See if we already have a side matching side_ordinal
      int i=0;
      for (; i < num_sides ; ++i)
      {
        if (mesh.is_valid(existing_side[i]) && existing_ords[i] == static_cast<ConnectivityOrdinal>(side_ordinal))
        {
          break;
        }
      }

      // a relation for the side was not found
      if (i == num_sides) {
        //create the side_key
        //side_key.first := CellTopologyData * of the side topology
        //side_key.second := EntityVector * of the side_nodes with the nodes in the correct
        //                                            permutation for the side starting with
        //                                            the node with the smallest identifier
        SideKey side_key;

        side_key.first = get_subcell_nodes(mesh,
            inside_entity,
            side_rank, // subcell rank
            side_ordinal,     // subcell identifier
            side_key.second  // subcell nodes
            );
        ensure_consistent_order(mesh, side_key.second);

        //add this side to the side_map
        side_map[side_key].push_back(inside);
      }
    }
  }
}

// populate the side_map with 'non-owned' sides that need to be created
// who's side_key is already present in the side_map.  This process
// may need to communicate with the process which own these elements
//
// a side needs to be created if the outside is NULL and the element
// does not have a current relation to a side for the side_ordinal
void add_non_owned_sides_to_map(
    const BulkData & mesh,
    const EntityRank element_rank,
    const EntitySideVector & boundary,
    BoundaryMap & side_map)
{
  const EntityRank side_rank = mesh.mesh_meta_data().side_rank();
  for (stk::mesh::EntitySideVector::const_iterator itr = boundary.begin();
      itr != boundary.end(); ++itr) {
    const EntitySideComponent & inside = itr->inside;
    const EntitySideComponent & outside = itr->outside;
    const RelationIdentifier side_ordinal = inside.side_ordinal;
    const Entity inside_entity = inside.entity;

    // If this process does NOT own the inside and the outside entity does not exist
    if ( mesh.parallel_owner_rank(inside_entity) != mesh.parallel_rank() &&
         !mesh.is_valid(outside.entity) )
    {
      // search through existing sides
      const int num_sides = mesh.num_connectivity(inside_entity, side_rank);
      Entity const * existing_entity = mesh.begin(inside_entity, side_rank);
      ConnectivityOrdinal const *existing_ords = mesh.begin_ordinals(inside_entity, side_rank);
      // See if we already have a side matching side_ordinal
      int i=0;
      for (; i < num_sides ; ++i)
      {
        if (mesh.is_valid(existing_entity[i]) && existing_ords[i] == static_cast<ConnectivityOrdinal>(side_ordinal))
        {
            break;
        }
      }

      if (i == num_sides) {
        // Get the nodes for the inside entity
        SideKey side_key;

        side_key.first = get_subcell_nodes(mesh,
            inside_entity,
            side_rank, // subcell rank
            side_ordinal,     // subcell identifier
            side_key.second  // subcell nodes
            );
        ensure_consistent_order(mesh, side_key.second);

        //only add the side if the side_key currently exist in the map
        if ( side_map.find(side_key) != side_map.end()) {
          side_map[side_key].push_back(inside);
        }
      }
    }
  }
}

//sort the SideVector for each side_key.
//the process who owns the element that appears
//first is responsible for generating the identifier
//for the side and communicating the identifier to the
//remaining process
//
//the ElementIdSide of the first side
//is enough to uniquely identify a side.
//the reverse_map is used to find the side_key
//from this ElementIdSide.
//
//return the number of sides this process
//is responsible for creating
size_t determine_creating_processes(
    const BulkData & mesh,
    BoundaryMap & side_map,
    ReverseBoundaryMap & reverse_side_map)
{
  EntitySideComponentLess esc_lesser(mesh);

  int num_sides_to_create = 0;
  for (BoundaryMap::iterator i = side_map.begin();
      i != side_map.end();
      ++i)
  {
    const SideKey & side_key = i->first;
    SideVector & side_vector = i->second;

    //sort the side vectors base on entity identifier and side ordinal
    std::sort( side_vector.begin(), side_vector.end(), esc_lesser );

    const EntitySideComponent & first_side = *side_vector.begin();

    //does this process create the first side
    //if so it needs to create a side_id
    if (mesh.parallel_owner_rank(first_side.entity) == mesh.parallel_rank()) {
      ++num_sides_to_create;
    }
    const ElementIdSide elem_id_side(mesh, first_side);

    reverse_side_map[elem_id_side] = side_key;
  }

  return num_sides_to_create;
}

} //end un-named namespace

void skin_mesh( BulkData & mesh, EntityRank element_rank, Part * skin_part) {
  ThrowErrorMsgIf( mesh.synchronized_state() == BulkData::MODIFIABLE,
                   "mesh is not SYNCHRONIZED" );

  EntityVector owned_elements;

  // select owned
  Selector owned = MetaData::get(mesh).locally_owned_part();
  get_selected_entities( owned,
                         mesh.buckets(element_rank),
                         owned_elements);

  reskin_mesh(mesh, element_rank, owned_elements, skin_part);
}

void reskin_mesh( BulkData & mesh, EntityRank element_rank, EntityVector & owned_elements, Part * skin_part) {
  ThrowErrorMsgIf( mesh.synchronized_state() == BulkData::MODIFIABLE,
                   "mesh is not SYNCHRONIZED" );

  const EntityRank side_rank = mesh.mesh_meta_data().side_rank();
  EntityVector elements_closure;

  // compute owned closure
  find_closure( mesh, owned_elements, elements_closure );

  // compute boundary
  EntitySideVector boundary;
  boundary_analysis( mesh, elements_closure, element_rank, boundary);

  BoundaryMap side_map;

  add_owned_sides_to_map(mesh, element_rank, boundary, side_map);

  add_non_owned_sides_to_map(mesh, element_rank, boundary, side_map);

  ReverseBoundaryMap reverse_side_map;

  size_t num_sides_to_create = determine_creating_processes(mesh, side_map, reverse_side_map);

  //begin modification
  mesh.modification_begin();

  // formulate request ids for the new sides
  std::vector<size_t> requests(MetaData::get(mesh).entity_rank_count(), 0);
  requests[side_rank] = num_sides_to_create;

  // create the new sides
  EntityVector requested_sides;
  mesh.generate_new_entities(requests, requested_sides);

  //set to aid comm packing
  std::set<SideCommHelper> side_comm_helper_set;


  size_t current_side = 0;
  for ( BoundaryMap::iterator map_itr = side_map.begin();
        map_itr!= side_map.end();
        ++map_itr)
  {
    const SideKey & side_key    = map_itr->first;
    SideVector    & side_vector = map_itr->second;

    // Only generated keys for sides in which this process
    // owns the first element in the side vector
    const EntitySideComponent & first_side = *(side_vector.begin());
    if ( mesh.parallel_owner_rank(first_side.entity) == mesh.parallel_rank()) {

      //to be used as the key in the reverse boundary map
      //a side can be identified in two ways
      //  1. vector of nodes in a correct permutation and a side-topology
      //  2. element-id side-ordinal for the side
      // the reverse boundary map is used to go from (2) to (1).
      const ElementIdSide elem_id_side(mesh, first_side);

      Entity side = requested_sides[current_side];
      EntityId side_id = mesh.identifier(side);


      PartVector add_parts ;
      {
        MetaData &fem_meta_data = MetaData::get(mesh);
        Part * topo_part = &fem_meta_data.get_cell_topology_root_part(side_key.first);
        add_parts.push_back( topo_part);
        if (skin_part) {
          add_parts.push_back(skin_part);
          CellTopology topo = fem_meta_data.get_cell_topology(*topo_part);
          const PartVector topo_parts =  skin_part->subsets();
          for (stk::mesh::PartVector::const_iterator i=topo_parts.begin(); i!=topo_parts.end(); ++i) {
            // Decode side to element topology. Specific to tri and quad on linear hex and wedges
            if (topo == fem_meta_data.get_cell_topology(**i)) {
              if (std::string(side_key.first->name) == "Quadrilateral_4") { // A enum would be nice
                // Quad could be the face of a hex or wedge.
                if (std::string("Hexahedron_8") == stk::mesh::get_cell_topology(mesh.bucket(side_vector.front().entity)).getName() &&
                                                   std::string::npos != (*i)->name().find("hex8")) { // Magic string
                  add_parts.push_back(*i);
                }
                else if (std::string("Wedge_6") == stk::mesh::get_cell_topology(mesh.bucket(side_vector.front().entity)).getName() &&
                                                   std::string::npos != (*i)->name().find("wedge6")) { // Magic string
                  add_parts.push_back(*i);
                }
              }
              else {
                add_parts.push_back(*i);
              }
            }
          }
        }
      }
      mesh.change_entity_parts(side, add_parts);


      //declare the side->node relations
      const EntityVector & nodes = side_key.second;

      for (size_t i = 0; i<nodes.size(); ++i) {
        Entity node = nodes[i];
        mesh.declare_relation( side, node, i);
      }

      //declare the elem->side relations
      for (SideVector::iterator side_itr = side_vector.begin();
          side_itr != side_vector.end();
          ++side_itr)
      {
        Entity elem = side_itr->entity;
        //only declare relations for owned elements
        if (mesh.parallel_owner_rank(elem) == mesh.parallel_rank()) {
          const RelationIdentifier elem_side_ordinal = side_itr->side_ordinal;
          mesh.declare_relation( elem, side, elem_side_ordinal);
        }
        else {
          //add this side to the communication set
          side_comm_helper_set.insert( SideCommHelper( mesh.parallel_owner_rank(elem), elem_id_side, side_id));
        }
      }
      ++current_side;
    }
  }

  CommAll comm( mesh.parallel() );

  //pack send buffers
  for (int allocation_pass=0; allocation_pass<2; ++allocation_pass) {
    if (allocation_pass==1) {
     comm.allocate_buffers( mesh.parallel_size() /4, 0);
    }

    for (std::set<SideCommHelper>::const_iterator i=side_comm_helper_set.begin();
         i != side_comm_helper_set.end();
         ++i)
    {
      const SideCommHelper & side_helper = *i;
      comm.send_buffer(side_helper.proc_to)
        .pack<EntityId>(side_helper.creating_elem_id_side.elem_id)
        .pack<RelationIdentifier>(side_helper.creating_elem_id_side.side_ordinal)
        .pack<EntityId>(side_helper.generated_side_id);
    }
  }

  comm.communicate();

  for ( int ip = 0 ; ip < mesh.parallel_size() ; ++ip ) {
    CommBuffer & buf = comm.recv_buffer( ip );
    while ( buf.remaining() ) {
      ElementIdSide creating_elem_id_side;
      EntityId      generated_side_id;

      buf.unpack<EntityId>(creating_elem_id_side.elem_id)
         .unpack<RelationIdentifier>(creating_elem_id_side.side_ordinal)
         .unpack<EntityId>(generated_side_id);

      //find the side in the reverse boundary map
      const SideKey & side_key = reverse_side_map[creating_elem_id_side];

      //get the SideVector for the corresponding side_key
      SideVector & side_vector = side_map[side_key];

      PartVector add_parts ;
      {
        MetaData &fem_meta_data =  MetaData::get(mesh);
        Part * topo_part = &fem_meta_data.get_cell_topology_root_part(side_key.first);
        add_parts.push_back( topo_part);
        if (skin_part) {
          add_parts.push_back(skin_part);
          CellTopology topo = fem_meta_data.get_cell_topology(*topo_part);
          PartVector topo_parts =  skin_part->subsets();
          for (stk::mesh::PartVector::const_iterator i=topo_parts.begin(); i!=topo_parts.end(); ++i) {
            // Decode side to element topology. Specific to tri and quad on linear hex and wedges
            if (topo == fem_meta_data.get_cell_topology(**i)) {
              if (std::string(side_key.first->name) == "Quadrilateral_4") { // A enum would be nice
                // Quad could be the face of a hex or wedge.
                if (std::string("Hexahedron_8") == stk::mesh::get_cell_topology(mesh.bucket(side_vector.front().entity)).getName() &&
                                                   std::string::npos != (*i)->name().find("hex8")) { // Magic string
                  add_parts.push_back(*i);
                }
                else if (std::string("Wedge_6") == stk::mesh::get_cell_topology(mesh.bucket(side_vector.front().entity)).getName() &&
                                                   std::string::npos != (*i)->name().find("wedge6")) { // Magic string
                  add_parts.push_back(*i);
                }
              }
              else {
                add_parts.push_back(*i);
              }
            }
          }
        }
      }
      Entity side = mesh.declare_entity(side_rank,generated_side_id,add_parts);

      //declare the side->node relations
      const EntityVector & nodes = side_key.second;

      for (size_t i = 0; i<nodes.size(); ++i) {
        Entity node = nodes[i];
        mesh.declare_relation( side, node, i);
      }

      //declare the elem->side relations
      for (SideVector::iterator side_itr = side_vector.begin();
           side_itr != side_vector.end();
           ++side_itr)
      {
        Entity elem = side_itr->entity;
        //only declare relations for owned elements
        if (mesh.parallel_owner_rank(elem) == mesh.parallel_rank()) {
          const RelationIdentifier elem_side_ordinal = side_itr->side_ordinal;
          mesh.declare_relation( elem, side, elem_side_ordinal);
        }
      }
    }
  }
  mesh.modification_end();
}

}
}
