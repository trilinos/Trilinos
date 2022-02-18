// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef Akri_ProlongationData_h
#define Akri_ProlongationData_h

#include <vector>
#include <set>
#include <map>

#include <Akri_Facet.hpp>
#include <Akri_InterfaceID.hpp>
#include <Akri_CDMesh.hpp>
#include <Akri_Vec.hpp>

namespace krino {

class ProlongationData {
public:
  ProlongationData() {}
  ~ProlongationData() {}

  const double * get_field_data( const stk::mesh::FieldBase& state_field ) const { const int field_index = my_field_indices[state_field.mesh_meta_data_ordinal()]; return (field_index < 0) ? nullptr : &my_field_data[field_index]; }
  std::string missing_prolongation_fields_for_entity( const CDMesh & mesh, const stk::mesh::Entity dst ) const;
  void restore_fields(const stk::mesh::BulkData& stk_mesh, stk::mesh::Entity entity) const;

protected:
  void save_fields(const stk::mesh::BulkData& stk_mesh, stk::mesh::Entity entity);
  const std::vector<double> & get_field_data() const { return my_field_data; }
  std::vector<double> & get_field_data() { return my_field_data; }
  const std::vector<int> & get_field_indices() const { return my_field_indices; }
  std::vector<int> & get_field_indices() { return my_field_indices; }

protected:
  mutable std::vector<int> my_field_indices;
  mutable std::vector<double> my_field_data;
};

class ProlongationPointData : public ProlongationData {
public:
  ProlongationPointData(const Vector3d & coordinates) : my_coordinates(coordinates) {}
  ProlongationPointData(const CDMesh & mesh, const FacetDistanceQuery & facet_dist_query, const std::vector<const ProlongationNodeData *> & facet_nodes);
  ~ProlongationPointData() {}

  const Vector3d & get_coordinates() const { return my_coordinates; }
  Vector3d & get_coordinates() { return my_coordinates; }

protected:
  Vector3d my_coordinates;
};

class ProlongationElementData : public ProlongationData {
public:
  ProlongationElementData(const stk::mesh::BulkData& stk_mesh, stk::mesh::Entity element) : ProlongationData() { save_fields(stk_mesh, element); }
  ProlongationElementData(const stk::mesh::BulkData& stk_mesh, const std::vector<const ProlongationElementData *> & children_data, const std::vector< std::vector<double> > & children_intg_wts);
  ~ProlongationElementData() {}

private:
  //: copy constructor not allowed
  ProlongationElementData(const ProlongationElementData & copy);
};

class ProlongationNodeData : public ProlongationPointData {
public:
  ProlongationNodeData(const CDMesh & mesh, stk::mesh::Entity node, bool communicate_me_to_all_sharers);
  ProlongationNodeData(const stk::mesh::EntityId in_entityId, const Vector3d & coordinates, const std::vector<unsigned> & fields, const std::vector<unsigned> & ioparts)
    : ProlongationPointData(coordinates), my_entityId(in_entityId), myCommunicateMeToAllSharersFlag(false), my_fields(fields), my_ioparts(ioparts) {}

  ~ProlongationNodeData() {}

  static std::vector<unsigned> get_node_io_parts(const stk::mesh::BulkData& stk_mesh, stk::mesh::Entity entity);
  static std::vector<unsigned> get_fields_on_node(const stk::mesh::BulkData& stk_mesh, stk::mesh::Entity entity);
  static Vector3d get_node_coordinates(const CDMesh & mesh, stk::mesh::Entity node);

  static ProlongationNodeData * unpack_from_buffer( stk::CommBuffer & b, const stk::mesh::MetaData & stk_meta ); // static method that builds surface from data in buffer for off-processor communication
  void pack_into_buffer(stk::CommBuffer & b) const;
  bool communicate_me_to_all_sharers() const { return myCommunicateMeToAllSharersFlag; }

  const std::vector<unsigned> & get_fields() const { ThrowRequireMsg(!my_fields.empty(), "Fields not set for prolongation node."); return my_fields; }
  const std::vector<unsigned> & get_io_parts() const { ThrowRequireMsg(!my_ioparts.empty(), "IO Parts not set for prolongation node."); return my_ioparts; }
  stk::mesh::EntityId entityId() const { return my_entityId; }

protected:
  const stk::mesh::EntityId my_entityId;
  const bool myCommunicateMeToAllSharersFlag;
  mutable std::vector<unsigned> my_fields;
  mutable std::vector<unsigned> my_ioparts;

private:
  //: copy constructor not allowed
  ProlongationNodeData(const ProlongationNodeData & copy);
};

class ProlongationFacet {
public:
  ProlongationFacet(const CDMesh & mesh, const std::vector<const ProlongationNodeData *> & prolong_nodes, const std::vector<unsigned> & common_fields);
  ProlongationFacet(const CDMesh & mesh, stk::mesh::Entity side);

  static void communicate( const CDMesh & mesh, ProlongFacetVec & proc_prolong_facets, EntityProlongationNodeMap & proc_prolong_nodes, const std::vector<BoundingBox> & proc_target_bboxes );
  static std::set<const ProlongationNodeData *> get_facet_nodes_to_communicate( const ProlongFacetVec & proc_prolong_facets, const BoundingBox & proc_target_bbox );

  void pack_into_buffer(stk::CommBuffer & b) const;
  static ProlongationFacet * unpack_from_buffer( const CDMesh & mesh, stk::CommBuffer & b ); // static method that builds surface from data in buffer for off-processor communication
  void compute_common_fields();

  static const BoundingBox & get_bounding_box(const ProlongationFacet * prolong_facet) { return prolong_facet->get_facet()->bounding_box(); }

  Facet * get_facet() const { return my_facet.get(); }
  const std::vector<unsigned> & get_common_fields() const { return my_common_fields; }
  const ProlongationPointData * get_prolongation_point_data(const FacetDistanceQuery & dist_query) const { update_prolongation_point_data(dist_query); return my_prolongation_point_data.get(); }
  const std::vector<const ProlongationNodeData *> & get_prolongation_nodes() const { return my_prolong_nodes; }
  bool communicate_me(const BoundingBox & proc_target_bbox) const;

protected:
  void update_prolongation_point_data(const FacetDistanceQuery & dist_query) const;

protected:
  const CDMesh & my_mesh;
  std::unique_ptr<Facet> my_facet;
  mutable std::unique_ptr<ProlongationPointData> my_prolongation_point_data;
  mutable std::vector<const ProlongationNodeData *> my_prolong_nodes;
  mutable std::vector<unsigned> my_common_fields;

private:
  //: copy constructor not allowed
  ProlongationFacet(const ProlongationFacet & copy);
  static void communicate_shared_nodes( const CDMesh & mesh, EntityProlongationNodeMap & proc_prolong_nodes );
  static void communicate_facets( const CDMesh & mesh, ProlongFacetVec & proc_prolong_facets, const std::vector<BoundingBox> & proc_target_bboxes );
  static void communicate_facet_nodes( const CDMesh & mesh, const ProlongFacetVec & proc_prolong_facets, EntityProlongationNodeMap & proc_prolong_nodes, const std::vector<BoundingBox> & proc_target_bboxes );
};

} // namespace krino

#endif // Akri_ProlongationData_h
