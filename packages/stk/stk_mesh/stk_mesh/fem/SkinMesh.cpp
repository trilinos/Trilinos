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

#include <stk_mesh/fem/BoundaryAnalysis.hpp>
#include <stk_mesh/fem/TopologyHelpers.hpp>
#include <stk_mesh/fem/SkinMesh.hpp>

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
    bool operator () (const EntitySideComponent & lhs, const EntitySideComponent &rhs)
    {
      const EntityId lhs_elem_id = lhs.entity->identifier();
      const EntityId rhs_elem_id = rhs.entity->identifier();

      return (lhs_elem_id != rhs_elem_id) ?
        (lhs_elem_id < rhs_elem_id) :
        (lhs.side_ordinal < rhs.side_ordinal);
    }
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

    ElementIdSide( const EntitySideComponent & esc) :
      elem_id(esc.entity->identifier()), side_ordinal(esc.side_ordinal) {}

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
void ensure_consistent_order(EntityVector & side_entities)
{
  ThrowRequire( !side_entities.empty() );

  EntityId lowest_id = side_entities.front()->identifier();
  unsigned idx_of_lowest_id = 0;

  for (unsigned idx = 1; idx < side_entities.size(); ++idx) {
    EntityId curr_id = side_entities[idx]->identifier();
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
  for (stk::mesh::EntitySideVector::const_iterator itr = boundary.begin();
      itr != boundary.end(); ++itr) {
    const EntitySideComponent & inside = itr->inside;
    const EntitySideComponent & outside = itr->outside;
    const RelationIdentifier side_ordinal = inside.side_ordinal;
    const Entity& inside_entity = *(inside.entity);

    if ( inside_entity.owner_rank() == mesh.parallel_rank() &&
        outside.entity == NULL ) {
      // search through existing sides
      PairIterRelation existing_sides = inside_entity.relations(element_rank -1);
      for (; existing_sides.first != existing_sides.second &&
          existing_sides.first->identifier() != side_ordinal ;
          ++existing_sides.first);

      // a relation the side was not found
      if (existing_sides.first == existing_sides.second) {
        //create the side_key
        //side_key.first := CellTopologyData * of the side topology
        //side_key.second := EntityVector * of the side_nodes with the nodes in the correct
        //                                            permutation for the side starting with
        //                                            the node with the smallest identifier
        SideKey side_key;

        side_key.first = get_subcell_nodes(
            inside_entity,
            element_rank - 1, // subcell rank
            side_ordinal,     // subcell identifier
            side_key.second  // subcell nodes
            );
        ensure_consistent_order(side_key.second);

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
  for (stk::mesh::EntitySideVector::const_iterator itr = boundary.begin();
      itr != boundary.end(); ++itr) {
    const EntitySideComponent & inside = itr->inside;
    const EntitySideComponent & outside = itr->outside;
    const RelationIdentifier side_ordinal = inside.side_ordinal;
    const Entity& inside_entity = *(inside.entity);

    // If this process does NOT own the inside and the outside entity does not exist
    if ( inside_entity.owner_rank() != mesh.parallel_rank() &&
        outside.entity == NULL ) {
      // search through existing sides
      PairIterRelation existing_sides = inside_entity.relations(element_rank -1);
      for (; existing_sides.first != existing_sides.second &&
          existing_sides.first->identifier() != side_ordinal ;
          ++existing_sides.first);

      // a relation to the side was not found
      if (existing_sides.first == existing_sides.second) {
        // Get the nodes for the inside entity
        SideKey side_key;

        side_key.first = get_subcell_nodes(
            inside_entity,
            element_rank - 1, // subcell rank
            side_ordinal,     // subcell identifier
            side_key.second  // subcell nodes
            );
        ensure_consistent_order(side_key.second);

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
  int num_sides_to_create = 0;
  for (BoundaryMap::iterator i = side_map.begin();
      i != side_map.end();
      ++i)
  {
    const SideKey & side_key = i->first;
    SideVector & side_vector = i->second;

    //sort the side vectors base on entity identifier and side ordinal
    std::sort( side_vector.begin(), side_vector.end(), EntitySideComponentLess() );

    const EntitySideComponent & first_side = *side_vector.begin();

    //does this process create the first side
    //if so it needs to create a side_id
    if (first_side.entity->owner_rank() == mesh.parallel_rank()) {
      ++num_sides_to_create;
    }
    const ElementIdSide elem_id_side(first_side);

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
  requests[element_rank -1] = num_sides_to_create;

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
    if ( first_side.entity->owner_rank() == mesh.parallel_rank()) {

      //to be used as the key in the reverse boundary map
      //a side can be identified in two ways
      //  1. vector of nodes in a correct permutation and a side-topology
      //  2. element-id side-ordinal for the side
      // the reverse boundary map is used to go from (2) to (1).
      const ElementIdSide elem_id_side( first_side);

      Entity & side = *(requested_sides[current_side]);
      EntityId side_id = side.identifier();


      PartVector add_parts ;
      {
        Part * topo_part = & fem::get_part(MetaData::get(mesh), side_key.first);
        add_parts.push_back( topo_part);
        if (skin_part) {
          add_parts.push_back(skin_part);
        }
      }
      mesh.change_entity_parts(side, add_parts);


      //declare the side->node relations
      const EntityVector & nodes = side_key.second;

      for (size_t i = 0; i<nodes.size(); ++i) {
        Entity & node = *nodes[i];
        mesh.declare_relation( side, node, i);
      }

      //declare the elem->side relations
      for (SideVector::iterator side_itr = side_vector.begin();
          side_itr != side_vector.end();
          ++side_itr)
      {
        Entity & elem = *(side_itr->entity);
        //only declare relations for owned elements
        if (elem.owner_rank() == mesh.parallel_rank()) {
          const RelationIdentifier elem_side_ordinal = side_itr->side_ordinal;
          mesh.declare_relation( elem, side, elem_side_ordinal);
        }
        else {
          //add this side to the communication set
          side_comm_helper_set.insert( SideCommHelper( elem.owner_rank(), elem_id_side, side_id));
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

  for ( unsigned ip = 0 ; ip < mesh.parallel_size() ; ++ip ) {
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
        Part * topo_part = & fem::get_part(MetaData::get(mesh), side_key.first);
        add_parts.push_back( topo_part);
        if (skin_part) {
          add_parts.push_back(skin_part);
        }
      }
      Entity & side = mesh.declare_entity(element_rank-1,generated_side_id,add_parts);

      //declare the side->node relations
      const EntityVector & nodes = side_key.second;

      for (size_t i = 0; i<nodes.size(); ++i) {
        Entity & node = *nodes[i];
        mesh.declare_relation( side, node, i);
      }

      //declare the elem->side relations
      for (SideVector::iterator side_itr = side_vector.begin();
           side_itr != side_vector.end();
           ++side_itr)
      {
        Entity & elem = *(side_itr->entity);
        //only declare relations for owned elements
        if (elem.owner_rank() == mesh.parallel_rank()) {
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
