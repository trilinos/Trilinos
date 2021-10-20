// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef Akri_MeshHelpers_h
#define Akri_MeshHelpers_h

#include <Akri_FieldRef.hpp>
#include <Akri_Vec.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/EntityKey.hpp>
#include <stk_mesh/base/Part.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_topology/topology.hpp>
#include <stk_util/parallel/CommSparse.hpp>
#include <iostream>
#include <set>
#include <string>
#include <utility>
#include <vector>

namespace krino { class FieldRef; }

namespace krino {

typedef std::pair<stk::topology, stk::mesh::Part*>  TopologyPartPair;
typedef std::vector<TopologyPartPair> TopologyPartVector;

struct StkMeshEntities
{
    typedef stk::mesh::Entity value_type;
    const value_type *mBegin;
    const value_type *mEnd;
    const value_type * begin() const { return mBegin; }
    const value_type * end() const { return mEnd; }
    size_t size() const { return mEnd - mBegin; }
    bool empty() const { return mEnd == mBegin; }
    value_type operator[](int i) const { return *(mBegin + i); }
};

void fill_element_node_coordinates(const stk::mesh::BulkData & mesh, stk::mesh::Entity element, const FieldRef coordsField, std::vector<Vector3d> & elementNodeCoords);
void fill_procs_owning_or_sharing_or_ghosting_node(const stk::mesh::BulkData& bulkData, stk::mesh::Entity node, std::vector<int> & procsOwningSharingOrGhostingNode);
double compute_maximum_element_size(stk::mesh::BulkData& mesh);
void compute_element_quality(const stk::mesh::BulkData & mesh, double & minEdgeLength, double & maxEdgeLength, double & minVolume, double & maxVolume);
void delete_all_entities_using_nodes_with_nodal_volume_below_threshold(stk::mesh::BulkData & mesh, const stk::mesh::Selector & blockSelector, const double threshold);
std::vector<unsigned> get_side_permutation(stk::topology topology, stk::mesh::Permutation node_permutation);
const stk::mesh::Part & find_element_part(const stk::mesh::BulkData& mesh, stk::mesh::Entity elem);
bool check_induced_parts(const stk::mesh::BulkData & mesh);
void attach_sides_to_elements(stk::mesh::BulkData & mesh);
void attach_entity_to_elements(stk::mesh::BulkData & mesh, stk::mesh::Entity entity);
void unpack_entities_from_other_procs(const stk::mesh::BulkData & mesh, std::set<stk::mesh::Entity> & entities, stk::CommSparse &commSparse);
void update_node_activation(stk::mesh::BulkData & mesh, stk::mesh::Part & active_part);
void activate_all_entities(stk::mesh::BulkData & mesh, stk::mesh::Part & active_part);
void destroy_custom_ghostings(stk::mesh::BulkData & mesh);
void delete_mesh_entities(stk::mesh::BulkData & mesh, std::vector<stk::mesh::Entity> & child_elems);
void debug_print_selector_parts(const stk::mesh::Selector & selector);
stk::mesh::PartVector filter_non_io_parts(const stk::mesh::PartVector & all_parts);
void activate_selected_sides_touching_active_elements(stk::mesh::BulkData & mesh, const stk::mesh::Selector & side_selector, stk::mesh::Part & active_part);
void get_partially_and_fully_coincident_elements(const stk::mesh::BulkData & mesh, stk::mesh::Entity elem, std::vector<stk::mesh::Entity> & coincident_elems);
bool check_element_side_connectivity(const stk::mesh::BulkData & mesh, const stk::mesh::Part & exterior_boundary_part, const stk::mesh::Part & active_part);
bool check_coincident_elements(const stk::mesh::BulkData & mesh, const stk::mesh::Part & active_part);
bool fix_coincident_element_ownership(stk::mesh::BulkData & mesh);
bool fix_face_and_edge_ownership(stk::mesh::BulkData & mesh);
bool check_face_and_edge_ownership(const stk::mesh::BulkData & mesh);
bool check_face_and_edge_relations(const stk::mesh::BulkData & mesh);
bool check_shared_entity_nodes(const stk::mesh::BulkData & mesh, stk::mesh::EntityKey remote_entity_key, std::vector<stk::mesh::EntityId> & remote_entity_node_ids);
bool check_shared_entity_nodes(const stk::mesh::BulkData & mesh, std::vector<stk::mesh::Entity> & entities);
bool check_shared_entity_nodes(const stk::mesh::BulkData & mesh);
void disconnect_entity(stk::mesh::BulkData & mesh, stk::mesh::Entity entity);
bool disconnect_and_destroy_entity(stk::mesh::BulkData & mesh, stk::mesh::Entity entity);

stk::mesh::PartVector get_common_io_parts(const stk::mesh::BulkData & mesh, const std::vector<stk::mesh::Entity> entities);
stk::mesh::PartVector get_removable_parts(const stk::mesh::BulkData & mesh, const stk::mesh::Bucket & bucket);
stk::mesh::PartVector get_removable_parts(const stk::mesh::BulkData & mesh, const stk::mesh::Entity entity);

void
store_edge_node_parent_ids(const stk::mesh::BulkData & mesh,
    const FieldRef & parent_id_field,
    stk::mesh::Entity edge_node_entity,
    stk::mesh::EntityId parent0_id,
    stk::mesh::EntityId parent1_id);

std::array<stk::mesh::EntityId, 2>
get_edge_node_parent_ids(const stk::mesh::BulkData & mesh,
    const FieldRef & parent_id_field,
    const stk::mesh::Entity edge_node_entity);

void get_parent_nodes_from_child(const stk::mesh::BulkData & mesh,
    stk::mesh::Entity child, const FieldRef & parent_id_field,
    std::set<stk::mesh::Entity> & parent_nodes);

double compute_child_position(const stk::mesh::BulkData & mesh, stk::mesh::Entity child, stk::mesh::Entity parent0, stk::mesh::Entity parent1);

// topology helpers
// NOTE: These use static storage, but it does not depend on anything so it should be ok for nested or multithreaded usage.
const unsigned * get_side_node_ordinals(stk::topology topology, unsigned side_ordinal);
const unsigned * get_edge_node_ordinals(stk::topology topology, unsigned edge_ordinal);

std::string debug_entity(const stk::mesh::BulkData & mesh, stk::mesh::Entity entity);
std::string debug_entity(const stk::mesh::BulkData & mesh, stk::mesh::Entity entity, const bool includeFields);

struct ChildNodeRequest
{
    std::vector<stk::mesh::Entity*> parents;
    stk::mesh::Entity* child;

    ChildNodeRequest(const std::vector<stk::mesh::Entity*> & in_parents, stk::mesh::Entity *in_child)
    : parents(in_parents), child(in_child) {}
};
struct SideRequest
{
    stk::mesh::Entity element;
    unsigned element_side_ordinal;
    stk::mesh::PartVector side_parts;

    SideRequest(stk::mesh::Entity in_element, unsigned in_element_side_ordinal, const stk::mesh::PartVector & in_side_parts)
    : element(in_element), element_side_ordinal(in_element_side_ordinal), side_parts(in_side_parts) {}
};

void batch_create_child_nodes(stk::mesh::BulkData & mesh, const std::vector< ChildNodeRequest > & child_node_requests, const stk::mesh::PartVector & node_parts, bool assert_32bit_ids, bool make_64bit_ids);
void batch_create_sides(stk::mesh::BulkData & mesh, const std::vector< SideRequest > & side_requests);
void make_side_ids_consistent_with_stk_convention(stk::mesh::BulkData & mesh);

double compute_element_volume_to_edge_ratio(stk::mesh::BulkData & mesh, stk::mesh::Entity element, const stk::mesh::Field<double> * const coords_field);

bool is_refinement_child(const stk::mesh::BulkData & stk_bulk, stk::mesh::Entity entity);
bool has_refinement_children(const stk::mesh::BulkData& stk_bulk, stk::mesh::Entity parent);
void get_refinement_immediate_children(const stk::mesh::BulkData& stk_bulk, stk::mesh::Entity parent, std::vector<stk::mesh::Entity> & children);
void get_refinement_leaf_children(const stk::mesh::BulkData& stk_bulk, stk::mesh::Entity entity, std::vector<stk::mesh::Entity> & leaf_children);
void get_refinement_all_children(const stk::mesh::BulkData& stk_bulk, stk::mesh::Entity entity, std::vector<stk::mesh::Entity> & children);

// Temporary method for manually correcting the relation permutation
void set_relation_permutation(stk::mesh::BulkData & mesh, stk::mesh::Entity from, stk::mesh::Entity to, stk::mesh::ConnectivityOrdinal to_ord, stk::mesh::Permutation to_permutation);

std::pair<stk::mesh::ConnectivityOrdinal, stk::mesh::Permutation>
determine_shell_side_ordinal_and_permutation(const stk::mesh::BulkData & mesh, stk::mesh::Entity shell, stk::mesh::Entity side);

stk::mesh::Permutation
determine_permutation(const stk::mesh::BulkData & mesh, const stk::mesh::Entity entity, const stk::mesh::Entity relative, const stk::mesh::ConnectivityOrdinal ordinal);

std::pair<stk::mesh::ConnectivityOrdinal, stk::mesh::Permutation>
determine_ordinal_and_permutation(const stk::mesh::BulkData & mesh, const stk::mesh::Entity entity, const stk::mesh::Entity relative);

inline stk::mesh::Entity find_entity_by_ordinal(const stk::mesh::BulkData &mesh, stk::mesh::Entity entity, stk::mesh::EntityRank rank, const unsigned ordinal)
{
  stk::mesh::ConnectivityOrdinal const* relative_ordinals = mesh.begin_ordinals(entity, rank);
  stk::mesh::Entity const* relatives = mesh.begin(entity, rank);
  const int num_relatives = mesh.num_connectivity(entity, rank);
  for (int i = 0; i < num_relatives; ++i)
  {
    if (relative_ordinals[i] == ordinal) {
      return relatives[i];
    }
  }
  return stk::mesh::Entity();
}

template<class CONTAINER>
void pack_entities_for_owning_proc(const stk::mesh::BulkData & mesh,
    const CONTAINER & entities,
    stk::CommSparse &commSparse)
{
  stk::pack_and_communicate(commSparse,[&]()
  {
    for (auto entity : entities)
    {
      const int entityOwner = mesh.parallel_owner_rank(entity);
      if (commSparse.parallel_rank() != entityOwner)
        commSparse.send_buffer(entityOwner).pack(mesh.entity_key(entity));
    }
  });
}

} // namespace krino

#endif // Akri_MeshHelpers_h
