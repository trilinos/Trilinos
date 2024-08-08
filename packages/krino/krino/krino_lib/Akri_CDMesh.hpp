// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef Akri_CDMesh_h
#define Akri_CDMesh_h

#include <Akri_TypeDefs.hpp>

#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/HashEntityAndEntityKey.hpp>
#include <Akri_EntityIdPool.hpp>

#include <vector>
#include <map>
#include <memory>

#include <Akri_PhaseTag.hpp>
#include <Akri_CDFEM_Support.hpp>
#include <Akri_CDFEM_Snapper.hpp>
#include <Akri_InterfaceID.hpp>
#include <Akri_NodeToCapturedDomains.hpp>
#include <Akri_OrderedIdPair.hpp>
#include <Akri_ParentsToChildMapper.hpp>
#include <Akri_SearchTree.hpp>
#include "Akri_ProlongationData.hpp"

namespace krino {

class SubElement;
class ElementObj;
class Mesh_Element;
class SubElementNode;
class SubElementMeshNode;
class ProlongationNodeData;
class ProlongationPointData;
class ProlongationElementData;
class ProlongationFacet;
class InterpolationEdge;
class LevelSet;
class Phase_Support;
class RefinementSupport;
class AuxMetaData;
struct SideDescription;
class InterfaceGeometry;

typedef std::vector< const SubElementNode * > NodeVec;
typedef std::set< const SubElementNode * > NodeSet;
typedef std::vector< const ProlongationFacet * > ProlongFacetVec;
typedef std::unordered_map<stk::mesh::EntityId, ProlongationNodeData *> EntityProlongationNodeMap;
typedef std::unordered_map<stk::mesh::EntityId, ProlongationElementData *> EntityProlongationElementMap;
typedef std::map<std::vector<unsigned>,std::unique_ptr<SearchTree<const ProlongationFacet*>>> PhaseProlongTreeMap;

enum CDMeshStatus
{
  MESH_MODIFIED = 1,
  COORDINATES_MAY_BE_MODIFIED = 2
};

class CDMesh {
public:

  CDMesh(stk::mesh::BulkData & mesh);

  virtual ~CDMesh();

  static bool decomposition_needs_update(const InterfaceGeometry & interfaceGeometry,
      const std::vector<std::pair<stk::mesh::Entity, stk::mesh::Entity>> & periodic_node_pairs);
  static void handle_possible_failed_time_step( stk::mesh::BulkData & mesh, const int step_count );
  static int decompose_mesh( stk::mesh::BulkData & mesh,
      const InterfaceGeometry & interfaceGeometry,
      const int step_count = 0,
      const std::vector<std::pair<stk::mesh::Entity, stk::mesh::Entity>> & periodic_node_pairs = {} );
  static void reset_mesh_to_original_undecomposed_state(stk::mesh::BulkData & mesh);
  static void nonconformal_adaptivity(stk::mesh::BulkData & mesh, const FieldRef coordsField, const InterfaceGeometry & interfaceGeometry);
  static void mark_interface_elements_for_adaptivity(stk::mesh::BulkData & mesh, const FieldRef coordsField, const RefinementSupport & refinementSupport, const InterfaceGeometry & interfaceGeometry, const int num_refinements);
  static void fixup_adapted_element_parts(stk::mesh::BulkData & mesh);
  static void rebuild_from_restart_mesh(stk::mesh::BulkData & mesh);
  static void prepare_for_resnapping(const stk::mesh::BulkData & mesh, const InterfaceGeometry & interfaceGeometry);
  void rebuild_after_rebalance_or_failed_step();

  static CDMesh* get_new_mesh() { return the_new_mesh.get(); }

  void snap_and_update_fields_and_captured_domains(const InterfaceGeometry & interfaceGeometry,
    NodeToCapturedDomainsMap & nodesToCapturedDomains) const;

  void communicate_prolongation_facet_fields() const;
  ProlongationQuery find_prolongation_node(const SubElementNode & node) const;
  const SubElementNode * find_new_node_with_common_ancestry_as_existing_node_with_given_id(const stk::mesh::EntityId nodeId) const;
  const SubElementNode * find_new_node_with_common_ancestry_as_existing_child_node(const stk::mesh::Entity node) const;

  PhaseTag determine_entity_phase(stk::mesh::Entity obj) const;

  std::vector<stk::mesh::Entity> get_nonconformal_elements() const;
  void generate_nonconformal_elements();

  bool triangulate(const InterfaceGeometry & interfaceGeometry); //return value indicates if any changes were made
  void decompose(const InterfaceGeometry & interfaceGeometry);
  void cut_sharp_features(const InterfaceGeometry & interfaceGeometry);
  void snap_nearby_intersections_to_nodes(const InterfaceGeometry & interfaceGeometry, NodeToCapturedDomainsMap & domainsAtNodes);
  void set_phase_of_uncut_elements(const InterfaceGeometry & interfaceGeometry);

  int spatial_dim() const { return my_spatial_dim; }

  const Phase_Support & get_phase_support() const { return my_phase_support; }
  const CDFEM_Support & get_cdfem_support() const { return my_cdfem_support; }
  CDFEM_Support & get_cdfem_support() { return my_cdfem_support; }
  bool need_nodes_for_prolongation() const { return INTERPOLATION != get_prolongation_model() && was_mesh_previously_decomposed(); }
  bool need_facets_for_prolongation() const { return ALE_CLOSEST_POINT == get_prolongation_model() && was_mesh_previously_decomposed(); }
  Prolongation_Model get_prolongation_model() const { return my_cdfem_support.get_prolongation_model(); }
  Edge_Interpolation_Model get_edge_interpolation_model() const { return my_cdfem_support.get_edge_interpolation_model(); }
  const std::vector<InterfaceID> & all_interface_ids(const std::vector<Surface_Identifier> & surfaceIdentifiers) const;
  std::vector<InterfaceID> active_interface_ids(const std::vector<Surface_Identifier> & surfaceIdentifiers) const;
  // This should really only be used for unit test purposes.
  void add_interface_id(const InterfaceID id) { crossing_keys.push_back(id); }

  const FieldRef get_coords_field() const { return my_cdfem_support.get_coords_field(); }
  const FieldRef get_cdfem_displacements_field() const { return my_cdfem_support.get_cdfem_displacements_field(); }
  const FieldSet & get_ale_prolongation_fields() const { return my_cdfem_support.get_ale_prolongation_fields(); }
  const FieldSet & get_interpolation_fields() const { return my_cdfem_support.get_interpolation_fields(); }
  const FieldSet & get_edge_interpolation_fields() const { return my_cdfem_support.get_edge_interpolation_fields(); }
  const FieldSet & get_zeroed_fields() const { return my_cdfem_support.get_zeroed_fields(); }
  const FieldSet & get_element_fields() const { return my_cdfem_support.get_element_fields(); }
  const CDFEM_Snapper & get_snapper() const { return my_cdfem_support.get_snapper(); }

  stk::mesh::Part & get_parent_part() const { return my_cdfem_support.get_parent_part(); }
  stk::mesh::Part & get_child_part() const { return my_cdfem_support.get_child_part(); }
  stk::mesh::Part & get_internal_side_part() const { return my_cdfem_support.get_internal_side_part(); }
  stk::mesh::Part & get_active_part() const;
  stk::mesh::Part & get_locally_owned_part() const;
  stk::mesh::Part & get_globally_shared_part() const;
  stk::mesh::Part & get_block_boundary_part() const;

  // methods for finding or creating subelement nodes
  const SubElementNode * create_mesh_node( const Mesh_Element * owner,
      const int lnn,
      stk::mesh::Entity nodeEntity );
  const SubElementNode * create_steiner_node( const Mesh_Element * owner,
    const NodeVec & parents,
    const std::vector<double> & weights );
  const SubElementNode * create_child_internal_or_face_node( const Mesh_Element * owner,
    const NodeVec & parents,
    const std::vector<double> & weights );
  const SubElementNode * create_edge_node( const Mesh_Element * owner,
    const SubElementNode * parent1,  const SubElementNode * parent2,
    const double position);
  const SubElementNode * create_midside_node( const Mesh_Element * owner,
    const SubElementNode * parent1,  const SubElementNode * parent2, stk::mesh::Entity entity = stk::mesh::Entity());
  void create_node_obj( const SubElementNode * node, NodeVec & parents );

  void determine_node_signs(const InterfaceID & interface);
  void decompose_edges(const InterfaceID & interface);
  void determine_node_scores(const InterfaceID & interface);

  const AuxMetaData& aux_meta() const { return my_aux_meta; }
  AuxMetaData& aux_meta() { return my_aux_meta; }
  const stk::mesh::MetaData& stk_meta() const { return my_meta; }
  stk::mesh::MetaData& stk_meta() { return my_meta; }
  const stk::mesh::BulkData& stk_bulk() const { return my_meta.mesh_bulk_data(); }
  stk::mesh::BulkData& stk_bulk() { return my_meta.mesh_bulk_data(); }

  void add_periodic_node_pair(stk::mesh::Entity node1, stk::mesh::Entity node2);

  // just search for existing mesh node and return NULL if not found
  const SubElementMeshNode * get_mesh_node( stk::mesh::EntityId node_id ) const;
  // just search for existing mesh node and return NULL if not found, starting with subelement node (possibly from another mesh)
  const SubElementMeshNode * get_mesh_node(const SubElementNode * new_node) const;

  bool check_element_side_parts() const;

  const SubElement * find_child_element(stk::mesh::Entity elem_mesh_obj) const;
  static const Mesh_Element * find_mesh_element(stk::mesh::EntityId elemId, const std::vector<std::unique_ptr<Mesh_Element> > & searchElements);
  const Mesh_Element * find_mesh_element(stk::mesh::EntityId elemId) const { return find_mesh_element(elemId, elements); }
  Mesh_Element * find_mesh_element(stk::mesh::EntityId elemId) { return const_cast<Mesh_Element*>(find_mesh_element(elemId, elements)); } // Scott Meyers says this is the right way
  bool get_parent_child_coord_transformation(stk::mesh::Entity elem_mesh_obj, double * dParentdChild) const;
  void get_parent_nodes_and_weights(stk::mesh::Entity child, stk::mesh::Entity & parent0, stk::mesh::Entity & parent1, double & position) const;
  stk::mesh::Entity get_parent_element(stk::mesh::Entity elem_mesh_obj) const;

  double compute_cdfem_cfl(const Interface_CFL_Length_Scale lengthScaleType, const std::function<stk::math::Vector3d(stk::mesh::Entity)> & get_side_displacement) const;
  double compute_cdfem_displacement_cfl() const;
  double compute_non_rebased_cdfem_displacement_cfl() const;
  double compute_interface_velocity_cfl(const FieldRef velocityField, const double dt) const;

  int stash_step_count() const { return my_stash_step_count; }
  bool was_mesh_previously_decomposed() const { return my_stash_step_count >= 0; }
  const ProlongationNodeData * fetch_prolong_node(stk::mesh::EntityId node_id) const { EntityProlongationNodeMap::const_iterator it = my_prolong_node_map.find(node_id); return ( (it == my_prolong_node_map.end()) ? NULL : (it->second) ); }
  const ProlongationElementData * fetch_prolong_element(stk::mesh::EntityId elem_id) const { EntityProlongationElementMap::const_iterator it = my_prolong_element_map.find(elem_id); return ( (it == my_prolong_element_map.end()) ? NULL : (it->second) ); }
  const PartAndFieldCollections & get_prolong_part_and_field_collections() const { return myProlongPartAndFieldCollections; }

  void debug_output() const;
  void print_conformal_volumes_and_surface_areas() const;
  std::vector<std::unique_ptr<Mesh_Element>> & get_elements() { return elements; }
  const std::vector<std::unique_ptr<Mesh_Element>> & get_elements() const { return elements; }

  // FIXME: Make these private.
  SubElementNode * add_managed_node(std::unique_ptr<SubElementNode> node);
  // expose vector of nodes
  std::vector<std::unique_ptr<SubElementNode> > nodes;
  // expose vector of elements
  std::vector<std::unique_ptr<Mesh_Element> > elements;
  // expose vector of sorted child elements
  std::vector<const ElementObj *> child_elements;

public: // for unit testing
  void update_adaptivity_parent_entities();
  void determine_conformal_parts(stk::mesh::Entity entity, const PhaseTag & phase, stk::mesh::PartVector & add_parts, stk::mesh::PartVector & remove_parts) const;
  void clear();

private:
  //: Default constructor not allowed
  CDMesh();

  template <class MESH_FIXTURE, class LS_FIELD_POLICY, unsigned NUM_LS>
  friend class CompleteDecompositionFixture;

  template <class MESH_FIXTURE>
  friend class AnalyticDecompositionFixture;

  template <class MESH_FIXTURE, class LS_FIELD_POLICY, unsigned NUM_LS>
  friend class DecompositionFixture;

  void build_parallel_hanging_edge_nodes();
  void handle_hanging_children(const InterfaceID & interface);
  void parallel_sync_nodes_on_interface();

  void sync_node_signs_on_constrained_nodes();
  void parallel_sync_node_signs_on_shared_nodes();
  void sync_node_scores_on_constrained_nodes();
  void parallel_sync_node_scores_on_shared_nodes();

  bool decomposition_has_changed(const InterfaceGeometry & interfaceGeometry);
  bool elem_io_part_changed(const ElementObj & elem) const;
  void determine_nonconformal_parts(stk::mesh::Entity entity, stk::mesh::PartVector & add_parts, stk::mesh::PartVector & remove_parts) const;
  void determine_conformal_parts(const stk::mesh::PartVector & current_parts, const stk::mesh::EntityRank entity_rank, const PhaseTag & phase, stk::mesh::PartVector & add_parts, stk::mesh::PartVector & remove_parts) const;
  void determine_child_conformal_parts(stk::topology topology, const stk::mesh::PartVector & parent_parts, const PhaseTag & phase, stk::mesh::PartVector & child_parts) const;
  void determine_element_side_parts(const stk::mesh::Entity side, stk::mesh::PartVector & add_parts, stk::mesh::PartVector & remove_parts) const;
  bool element_side_should_be_active(const stk::mesh::Entity side) const;
  void stash_field_data(const int step_count) const;
  void stash_nodal_field_data() const;
  void stash_elemental_field_data() const;
  std::vector<std::vector<stk::mesh::Entity>> get_subelements_for_CDFEM_parents(const std::vector<stk::mesh::Entity> & sortedCdfemParentElems) const;

  void clear_prolongation_trees() const;
  void build_prolongation_trees() const;

  void clear_prolongation_data() const;

  Mesh_Element * create_mesh_element(stk::mesh::Entity mesh_obj);
  SubElementMeshNode * add_managed_mesh_node( std::unique_ptr<SubElementMeshNode> node );

  void store_child_node_parent_ids_and_weights_fields() const;
  bool modify_mesh();
  void handle_single_coincident_subelement(const Mesh_Element & elem, const SubElement * subelem, std::vector<SideDescription> & side_requests);
  void create_subelement_mesh_entities(const Mesh_Element & elem,
    const std::vector<const SubElement *> conformal_subelems);
  void attach_existing_and_identify_missing_subelement_sides(const Mesh_Element & elem,
    const std::vector<const SubElement *> conformal_subelems,
    std::vector<SideDescription> & side_requests);
  void update_uncut_element(const Mesh_Element & elem);

  std::vector<stk::mesh::Entity> get_owned_unused_old_child_elements_and_clear_child_elements();
  void set_entities_for_child_nodes_with_common_ancestry_as_existing_child_nodes();
  bool set_entities_for_existing_child_elements();
  void create_new_mesh_entities();
  void create_node_entities();
  void create_element_and_side_entities(std::vector<SideDescription> & side_requests);

  void determine_processor_prolongation_bounding_box(const bool guessAndCheckProcPadding, const double maxCFLGuess, BoundingBox & procBbox) const;
  void prolongation();
  void rebase_cdfem_displacements();
  double get_maximum_cdfem_displacement() const;

  void add_possible_interface_sides(std::vector<SideDescription> & sideRequests) const;
  bool check_element_side_parts(const std::vector<stk::mesh::Entity> & side_nodes) const;
  void update_element_side_parts();

  void parallel_communicate_elemental_death_fields() const;

  void generate_sorted_child_elements();
  void cache_node_ids();

  stk::mesh::Part & get_child_edge_node_part() const { return my_cdfem_support.get_child_node_part(); }
  FieldRef get_parent_node_ids_field() const { return my_cdfem_support.get_parent_node_ids_field(); }
  FieldRef get_parent_node_weights_field() const { return my_cdfem_support.get_parent_node_weights_field(); }

  void rebuild_child_part();
  void rebuild_parent_and_active_parts_using_nonconformal_and_child_parts();
  void restore_subelements();
  const SubElementNode * build_subelement_child_node(const stk::mesh::Entity node, const Mesh_Element & ownerMeshElem, std::map<stk::mesh::EntityId, const SubElementNode*> & idToSubElementNode);
  const SubElementNode * find_or_build_subelement_node_with_id(const stk::mesh::EntityId nodeId, const Mesh_Element & ownerMeshElem, std::map<stk::mesh::EntityId, const SubElementNode*> & idToSubElementNode);
  const SubElementNode * find_or_build_subelement_node(const stk::mesh::Entity node, const Mesh_Element & ownerMeshElem, std::map<stk::mesh::EntityId, const SubElementNode*> & idToSubElementNode);
  void find_or_build_midside_nodes(const stk::topology & elemTopo, const Mesh_Element & ownerMeshElem, const stk::mesh::Entity * elemNodes, const NodeVec & subelemNodes);

  stk::mesh::MetaData& my_meta;
  AuxMetaData& my_aux_meta;
  EntityIdPool my_entity_id_pool;

  const int my_spatial_dim;
  CDFEM_Support & my_cdfem_support;
  Phase_Support & my_phase_support;
  RefinementSupport & myRefinementSupport;
  typedef std::unordered_map<stk::mesh::EntityId, const SubElementMeshNode*> NodeMap;
  NodeMap mesh_node_map;

  static std::unique_ptr<CDMesh> the_new_mesh;

  std::unordered_map<stk::mesh::EntityId, std::vector<stk::mesh::EntityId> > my_periodic_node_id_map;

  mutable int my_stash_step_count;
  mutable PhaseProlongTreeMap my_phase_prolong_tree_map;

  mutable bool my_missing_remote_prolong_facets;
  mutable PartAndFieldCollections myProlongPartAndFieldCollections;
  mutable EntityProlongationNodeMap my_prolong_node_map;
  mutable EntityProlongationElementMap my_prolong_element_map;
  mutable ProlongFacetVec my_prolong_facets;

  mutable std::vector<InterfaceID> crossing_keys;

  mutable stk::diag::Timer my_timer_decompose;
  mutable stk::diag::Timer my_timer_decomposition_has_changed;
  mutable stk::diag::Timer my_timer_snap;
  mutable stk::diag::Timer my_timer_stash_field_data;
  mutable stk::diag::Timer my_timer_modify_mesh;
  mutable stk::diag::Timer my_timer_prolongation;
  mutable stk::diag::Timer my_timer_compute_CFL;
  stk::mesh::PartVector my_attribute_parts;

  mutable std::map<std::pair<const SubElementNode*,const SubElementNode*>,const SubElementNode*> my_midside_node_map;
};

} // namespace krino

#endif // Akri_CDMesh_h
