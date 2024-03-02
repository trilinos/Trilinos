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
#include <Akri_MasterElement.hpp>
#include <stk_math/StkVector.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/Entity.hpp>

namespace krino {

class CDFEM_Support;
class CDMesh;
class FieldRef;
class ProlongationFacet;
class ProlongationNodeData;
typedef std::vector< const ProlongationFacet * > ProlongFacetVec;
typedef std::unordered_map<stk::mesh::EntityId, ProlongationNodeData *> EntityProlongationNodeMap;

class PartCollection {
public:
  PartCollection() = default;
  PartCollection(std::vector<unsigned> & parts) { myParts.swap(parts); }
  const std::vector<unsigned> & get_parts() const { return myParts; }
  bool operator<(const PartCollection & other) const{ return get_parts() < other.get_parts(); }
  bool operator==(const PartCollection & other) const { return get_parts() == other.get_parts(); }
private:
  std::vector<unsigned> myParts;
};

class FieldCollection {
public:
  FieldCollection() = default;
  FieldCollection(std::vector<unsigned> & fields) { myFields.swap(fields); }
  const std::vector<unsigned> & get_fields() const { return myFields; }
  bool operator<(const FieldCollection & other) const { return get_fields() < other.get_fields(); }
  bool operator==(const FieldCollection & other) const { return get_fields() == other.get_fields(); }
  const std::vector<int> & get_field_storage_indices() const {  return myFieldStorageIndices; }
  size_t get_field_storage_size() const {  return myFieldStorageSize; }
  void compute_field_storage_indices(const stk::mesh::BulkData & mesh) const;

private:
  std::vector<unsigned> myFields;
  mutable size_t myFieldStorageSize{0};
  mutable std::vector<int> myFieldStorageIndices;
};

class PartAndFieldCollections {
public:
  void build(const stk::mesh::BulkData & mesh);
  bool have_part_collection(const PartCollection & parts) const;
  bool have_field_collection(const FieldCollection & fields) const;
  int get_part_collection_id(const PartCollection & parts) const;
  int get_field_collection_id(const FieldCollection & fields) const;
  int get_part_collection_id(const stk::mesh::BulkData& mesh, const stk::mesh::Bucket & bucket) const;
  int get_field_collection_id(const stk::mesh::BulkData& mesh, const stk::mesh::Bucket & bucket) const;
  const std::vector<unsigned> & get_parts(const unsigned partCollectionId) const { return myPartCollections[partCollectionId].get_parts(); }
  const std::vector<unsigned> & get_fields(const unsigned fieldCollectionId) const { return myFieldCollections[fieldCollectionId].get_fields(); }
  const std::vector<int> & get_field_storage_indices(const unsigned fieldCollectionId) const { return myFieldCollections[fieldCollectionId].get_field_storage_indices(); }
  size_t get_field_storage_size(const unsigned fieldCollectionId) const { return myFieldCollections[fieldCollectionId].get_field_storage_size(); }

  static std::vector<unsigned> determine_io_parts(const stk::mesh::Bucket & bucket);
  static std::vector<unsigned> determine_fields(const stk::mesh::BulkData& mesh, const stk::mesh::Bucket & bucket);
private:
  void communicate(const stk::mesh::BulkData & mesh);
  std::vector<PartCollection> myPartCollections;
  std::vector<FieldCollection> myFieldCollections;
};

class ProlongationData {
public:
  ProlongationData() {}

  const double * get_field_data( const stk::mesh::FieldBase& state_field ) const;
  const std::vector<double> & get_field_data() const { return myFieldData; }

protected:
  void save_field_data(const stk::mesh::BulkData& stk_mesh, const PartAndFieldCollections & partAndFieldCollections, const int fieldCollectionId, const stk::mesh::Entity entity);
  void save_field_data(const std::vector<int> & fieldStorageIndices, std::vector<double> & fieldData) { myFieldStorageIndices = &fieldStorageIndices; myFieldData.swap(fieldData); }
  std::string debug_data_output(const stk::mesh::BulkData & mesh) const;

private:
  const std::vector<int> * myFieldStorageIndices{nullptr};
  std::vector<double> myFieldData;
};

class ProlongationPointData : public ProlongationData {
public:
  static void set_coords_fields(const int spatialDim, FieldRef coordsField, FieldRef snapDisplacementsField);
  ProlongationPointData() {}

  stk::math::Vector3d get_previous_coordinates() const;
  stk::math::Vector3d get_post_snap_coordinates() const;

protected:
  static int theSpatialDim;
  static FieldRef theSnapDisplacementsField;
  static FieldRef theCoordsField;
};

class ProlongationFacetPointData : public ProlongationPointData {
public:
  ProlongationFacetPointData(const CDMesh & mesh, const FacetDistanceQuery<Facet> & facetDistanceQuery, const std::vector<const ProlongationNodeData *> & facetNodes);
private:
  void interpolate_to_point(const stk::mesh::MetaData & meta, const FacetDistanceQuery<Facet> & facetDistanceQuery, const std::vector<const ProlongationNodeData *> & facetNodes);
  std::vector<int> myFieldStorageIndices;
};

class ProlongationElementData : public ProlongationData {
public:
  ProlongationElementData(const CDMesh & cdmesh, const stk::mesh::Entity element);
  virtual ~ProlongationElementData() {}
  ProlongationElementData (const ProlongationElementData&) = delete;
  ProlongationElementData& operator= (const ProlongationElementData&) = delete;
  void evaluate_prolongation_field(const CDFEM_Support & cdfemSupport, const FieldRef field, const unsigned field_length, const stk::math::Vector3d & p_coords, double * result) const;
  stk::math::Vector3d compute_parametric_coords_at_point(const stk::math::Vector3d & pointCoords) const;
  bool have_prolongation_data_stored_for_all_nodes() const;
  void fill_integration_weights(std::vector<double> & childIntgWeights) const;
  virtual void find_subelement_and_parametric_coordinates_at_point(const stk::math::Vector3d & pointCoordinates, const ProlongationElementData *& interpElem, stk::math::Vector3d & interpElemParamCoords) const = 0;
  virtual bool have_subelements() const = 0;
private:
  const MasterElement& myMasterElem;
  std::vector<const ProlongationNodeData *> myElemNodesData;
};

class ProlongationLeafElementData : public ProlongationElementData {
public:
  ProlongationLeafElementData(const CDMesh & cdmesh, const PartAndFieldCollections & partAndFieldCollections, const stk::mesh::Entity element);
  virtual ~ProlongationLeafElementData() {}
  virtual void find_subelement_and_parametric_coordinates_at_point(const stk::math::Vector3d & pointCoordinates, const ProlongationElementData *& interpElem, stk::math::Vector3d & interpElemParamCoords) const override;
  virtual bool have_subelements() const override { return false; }
private:
  std::string debug_output(const stk::mesh::BulkData & mesh, const stk::mesh::Entity element) const;
};

class ProlongationParentElementData : public ProlongationElementData {
public:
  ProlongationParentElementData(const CDMesh & cdmesh, const stk::mesh::Entity element, const std::vector<const ProlongationElementData *> & subelementsData, const bool doStoreElementFields);
  virtual ~ProlongationParentElementData() {}
  virtual void find_subelement_and_parametric_coordinates_at_point(const stk::math::Vector3d & pointCoordinates, const ProlongationElementData *& interpElem, stk::math::Vector3d & interpElemParamCoords) const override;
  virtual bool have_subelements() const override { return true; }
private:
  std::string debug_output(const stk::mesh::BulkData & mesh, const stk::mesh::Entity element) const;
  void homogenize_subelement_fields(const CDMesh & cdmesh,
    const std::vector<const ProlongationElementData *> & subelementsData);
  bool have_prolongation_data_stored_for_all_nodes_of_subelements() const;

  std::vector<const ProlongationElementData *> mySubelementsData;
  std::vector<int> myFieldStorageIndices;
};

class ProlongationNodeData : public ProlongationPointData {
public:
  ProlongationNodeData(const CDMesh & mesh, const PartAndFieldCollections & partAndFieldCollections, const stk::mesh::Entity node, bool communicate_me_to_all_sharers);
  ProlongationNodeData(const PartAndFieldCollections & partAndFieldCollections, const stk::mesh::EntityId in_entityId, const int partCollectionId, const int fieldCollectionId, std::vector<double> & fieldData);
  ProlongationNodeData (const ProlongationNodeData&) = delete;
  ProlongationNodeData& operator= (const ProlongationNodeData&) = delete;

  ~ProlongationNodeData() {}

  static ProlongationNodeData * unpack_from_buffer( stk::CommBuffer & b, const PartAndFieldCollections & partAndFieldCollections );
  void pack_into_buffer(stk::CommBuffer & b) const;
  bool communicate_me_to_all_sharers() const { return myCommunicateMeToAllSharersFlag; }
  int get_part_collection_id() const { return myPartCollectionId; }
  int get_field_collection_id() const { return myFieldCollectionId; }

  stk::mesh::EntityId entityId() const { return my_entityId; }

protected:
  std::string debug_output(const stk::mesh::BulkData & mesh, const stk::mesh::Entity node) const;

  stk::mesh::EntityId my_entityId;
  bool myCommunicateMeToAllSharersFlag;
  int myPartCollectionId;
  int myFieldCollectionId;
};

class ProlongationFacet {
public:
  ProlongationFacet(const CDMesh & mesh, const std::vector<const ProlongationNodeData *> & prolong_nodes);
  ProlongationFacet(const CDMesh & mesh, stk::mesh::Entity side);
  ProlongationFacet (const ProlongationFacet&) = delete;
  ProlongationFacet& operator= (const ProlongationFacet&) = delete;

  static void communicate( const CDMesh & mesh, ProlongFacetVec & proc_prolong_facets, EntityProlongationNodeMap & proc_prolong_nodes, const std::vector<BoundingBox> & proc_target_bboxes );
  static std::set<const ProlongationNodeData *> get_facet_nodes_to_communicate( const ProlongFacetVec & proc_prolong_facets, const BoundingBox & proc_target_bbox );

  void pack_into_buffer(stk::CommBuffer & b) const;
  static ProlongationFacet * unpack_from_buffer( const CDMesh & mesh, stk::CommBuffer & b );

  Facet * get_facet() const { return my_facet.get(); }
  std::vector<unsigned> compute_common_fields() const;
  std::unique_ptr<ProlongationFacetPointData> get_prolongation_point_data(const FacetDistanceQuery<Facet> & dist_query) const;
  const std::vector<const ProlongationNodeData *> & get_prolongation_nodes() const { return my_prolong_nodes; }
  bool communicate_me(const BoundingBox & proc_target_bbox) const;
  static void insert_into_bounding_box(const ProlongationFacet * prolongFacet, BoundingBox & bbox) { return prolongFacet->get_facet()->insert_into(bbox); }
  static stk::math::Vector3d get_centroid(const ProlongationFacet * prolongFacet) { return prolongFacet->get_facet()->centroid(); }

  void build_and_append_edge_facets(ProlongFacetVec & facetVec) const;

protected:
  const CDMesh & my_mesh;
  std::unique_ptr<Facet> my_facet;
  mutable std::vector<const ProlongationNodeData *> my_prolong_nodes;

private:
  static void communicate_shared_nodes( const CDMesh & mesh, EntityProlongationNodeMap & proc_prolong_nodes );
  static void communicate_facets( const CDMesh & mesh, ProlongFacetVec & proc_prolong_facets, const std::vector<BoundingBox> & proc_target_bboxes );
  static void communicate_facet_nodes( const CDMesh & mesh, const ProlongFacetVec & proc_prolong_facets, EntityProlongationNodeMap & proc_prolong_nodes, const std::vector<BoundingBox> & proc_target_bboxes );
};

class ProlongationQuery {
public:
  ProlongationQuery() = default;
  ProlongationQuery(const ProlongationFacet & facet, const FacetDistanceQuery<Facet> & distQuery) { myProlongFacetPointData = facet.get_prolongation_point_data(distQuery); }
  ProlongationQuery(const ProlongationNodeData * prolongNodeData) : myProlongNodeData(prolongNodeData) {}
  const ProlongationPointData * get_prolongation_point_data() const { return myProlongFacetPointData ? get_facet_point_data() : myProlongNodeData; }
private:
  const ProlongationPointData * get_facet_point_data() const { return myProlongFacetPointData.get(); }
  const ProlongationPointData * myProlongNodeData{nullptr};
  std::unique_ptr<ProlongationFacetPointData> myProlongFacetPointData;
};

} // namespace krino

#endif // Akri_ProlongationData_h
