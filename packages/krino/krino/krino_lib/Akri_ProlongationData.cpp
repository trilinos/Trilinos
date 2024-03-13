// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_CDFEM_Support.hpp>
#include <Akri_ProlongationData.hpp>
#include <Akri_DiagWriter.hpp>
#include <Akri_CDMesh.hpp>
#include <Akri_Element.hpp>
#include <Akri_MasterElementDeterminer.hpp>
#include <Akri_MathUtil.hpp>
#include <Akri_MeshHelpers.hpp>
#include <Akri_FieldRef.hpp>

#include <stk_mesh/base/CommunicateMeshTypes.hpp>
#include <stk_util/parallel/CommSparse.hpp>

#include <cmath>

namespace krino{

int ProlongationPointData::theSpatialDim = 0;
FieldRef ProlongationPointData::theCoordsField{};
FieldRef ProlongationPointData::theSnapDisplacementsField{};

void FieldCollection::compute_field_storage_indices(const stk::mesh::BulkData & mesh) const
{
  const stk::mesh::FieldVector & allFields = mesh.mesh_meta_data().get_fields();

  myFieldStorageIndices.clear();
  myFieldStorageIndices.resize(allFields.size(), -1);

  myFieldStorageSize = 0;
  for ( auto fieldIndex : myFields )
  {
    myFieldStorageIndices[fieldIndex] = myFieldStorageSize;

    const FieldRef field = allFields[fieldIndex];
    myFieldStorageSize += field.length();
  }
}

static const std::vector<unsigned> & get_data(const PartCollection & partCollection)
{
  return partCollection.get_parts();
}

static const std::vector<unsigned> & get_data(const FieldCollection & fieldCollection)
{
  return fieldCollection.get_fields();
}

template<typename CollectionType>
void pack_collections_for_all_other_procs(const std::vector<CollectionType> & partOrFieldCollections, stk::CommSparse &commSparse)
{
  stk::pack_and_communicate(commSparse,[&]()
  {
    for ( int procId=0; procId<commSparse.parallel_size(); ++procId )
    {
      if ( commSparse.parallel_rank() == procId ) continue;  // Don't talk to yourself, it's embarrassing
      stk::CommBuffer & buffer = commSparse.send_buffer(procId);

      for (auto && partOrFieldCollection : partOrFieldCollections)
      {
        const std::vector<unsigned> & partsOrFields = get_data(partOrFieldCollection);
        const size_t len = partsOrFields.size();
        buffer.pack(len);
        buffer.pack(partsOrFields.data(), len);
      }
    }
  });
}

template<typename CollectionType>
std::vector<CollectionType> receive_collections_that_this_proc_does_not_have_yet(const std::vector<CollectionType> & existingCollections, stk::CommSparse &commSparse)
{
  std::vector<CollectionType> missingCollections;
  std::vector<unsigned> partsOrFields;

  stk::unpack_communications(commSparse, [&](int procId)
  {
    stk::CommBuffer & buffer = commSparse.recv_buffer(procId);

    while ( buffer.remaining() )
    {
      size_t len = 0;
      buffer.unpack(len);
      partsOrFields.resize(len);
      buffer.unpack(partsOrFields.data(), len);
      CollectionType partOrFieldCollection(partsOrFields);
      auto iter = std::lower_bound(existingCollections.begin(), existingCollections.end(), partOrFieldCollection);
      const bool isMissing = (iter == existingCollections.end() || !(*iter == partOrFieldCollection));
      if (isMissing)
        missingCollections.push_back(partOrFieldCollection);
    }
  });

  return missingCollections;
}

void PartAndFieldCollections::communicate(const stk::mesh::BulkData & mesh)
{
  {
    stk::CommSparse commSparse(mesh.parallel());
    pack_collections_for_all_other_procs(myPartCollections, commSparse);
    const std::vector<PartCollection> missingPartCollections = receive_collections_that_this_proc_does_not_have_yet(myPartCollections, commSparse);
    for (auto && missingPartCollection : missingPartCollections)
      myPartCollections.push_back(missingPartCollection);
    stk::util::sort_and_unique(myPartCollections);
  }

  {
    stk::CommSparse commSparse(mesh.parallel());
    pack_collections_for_all_other_procs(myFieldCollections, commSparse);
    const std::vector<FieldCollection> missingFieldCollections = receive_collections_that_this_proc_does_not_have_yet(myFieldCollections, commSparse);
    for (auto && missingFieldCollection : missingFieldCollections)
      myFieldCollections.push_back(missingFieldCollection);
    stk::util::sort_and_unique(myFieldCollections);
  }
}

void PartAndFieldCollections::build(const stk::mesh::BulkData & mesh)
{
  for(const auto & bucketPtr : mesh.buckets(stk::topology::NODE_RANK))
  {
    std::vector<unsigned> parts = determine_io_parts(*bucketPtr);
    std::vector<unsigned> fields = determine_fields(mesh, *bucketPtr);
    myPartCollections.emplace_back(parts);
    myFieldCollections.emplace_back(fields);
  }

  for(const auto & bucketPtr : mesh.buckets(stk::topology::ELEMENT_RANK))
  {
    std::vector<unsigned> fields = determine_fields(mesh, *bucketPtr);
    myFieldCollections.emplace_back(fields);
  }

  stk::util::sort_and_unique(myPartCollections);
  stk::util::sort_and_unique(myFieldCollections);

  communicate(mesh);

  for (auto && fieldCollection : myFieldCollections)
    fieldCollection.compute_field_storage_indices(mesh);
}

int PartAndFieldCollections::get_part_collection_id(const PartCollection & partCollection) const
{
  const auto iter = std::lower_bound(myPartCollections.begin(), myPartCollections.end(), partCollection);
  STK_ThrowAssertMsg(iter != myPartCollections.end() && *iter == partCollection, "Failed to find collection of parts.");
  return std::distance(myPartCollections.begin(), iter);
}

int PartAndFieldCollections::get_field_collection_id(const FieldCollection & fieldCollection) const
{
  const auto iter = std::lower_bound(myFieldCollections.begin(), myFieldCollections.end(), fieldCollection);
  STK_ThrowAssertMsg(iter != myFieldCollections.end() && *iter == fieldCollection, "Failed to find collection of parts.");
  return std::distance(myFieldCollections.begin(), iter);
}

bool PartAndFieldCollections::have_part_collection(const PartCollection & partCollection) const
{
  const auto iter = std::lower_bound(myPartCollections.begin(), myPartCollections.end(), partCollection);
  return iter != myPartCollections.end() && *iter == partCollection;
}

bool PartAndFieldCollections::have_field_collection(const FieldCollection & fieldCollection) const
{
  const auto iter = std::lower_bound(myFieldCollections.begin(), myFieldCollections.end(), fieldCollection);
  return iter != myFieldCollections.end() && *iter == fieldCollection;
}

int PartAndFieldCollections::get_part_collection_id(const stk::mesh::BulkData& mesh, const stk::mesh::Bucket & bucket) const
{
  std::vector<unsigned> parts = determine_io_parts(bucket);
  PartCollection partCollection(parts);
  return get_part_collection_id(partCollection);
}

int PartAndFieldCollections::get_field_collection_id(const stk::mesh::BulkData& mesh, const stk::mesh::Bucket & bucket) const
{
  std::vector<unsigned> fields = determine_fields(mesh, bucket);
  FieldCollection fieldCollection(fields);
  return get_field_collection_id(fieldCollection);
}

std::vector<unsigned> PartAndFieldCollections::determine_io_parts(const stk::mesh::Bucket & bucket)
{
  // This list of parts is used to determine if a node needs to be ALE prolonged
  stk::mesh::PartVector node_parts = filter_non_io_parts(bucket.supersets());
  std::vector<unsigned> node_part_ids;
  node_part_ids.reserve(node_parts.size());
  for (auto * node_part : node_parts)
    node_part_ids.push_back(node_part->mesh_meta_data_ordinal());
  std::sort(node_part_ids.begin(), node_part_ids.end());
  return node_part_ids;
}

std::vector<unsigned> PartAndFieldCollections::determine_fields(const stk::mesh::BulkData& mesh, const stk::mesh::Bucket & bucket)
{
  const stk::mesh::EntityRank entity_rank = bucket.entity_rank();
  const stk::mesh::FieldVector & all_fields = mesh.mesh_meta_data().get_fields();
  std::vector<unsigned> entity_fields;
  entity_fields.reserve(all_fields.size());
  for ( auto && fieldPtr : all_fields )
  {
    const FieldRef field = *fieldPtr;
    if (field.field().entity_rank() == entity_rank && field.type_is<double>() && nullptr != field_data<double>(field, bucket))
      entity_fields.push_back(field.field().mesh_meta_data_ordinal());
  }

  std::sort(entity_fields.begin(), entity_fields.end());
  return entity_fields;
}

stk::math::Vector3d ProlongationPointData::get_previous_coordinates() const
{
  STK_ThrowAssertMsg(theCoordsField.valid(), "Static member coordinates field is not yet set.");

  stk::math::Vector3d coords(get_field_data(theCoordsField), theSpatialDim);
  if (theSnapDisplacementsField.valid())
  {
    FieldRef oldCdfemSnapDispField = theSnapDisplacementsField.field_state(stk::mesh::StateOld);
    const double * cdfemSnapDispPtr = get_field_data(theSnapDisplacementsField);
    const double * oldCdfemSnapDispPtr = get_field_data(oldCdfemSnapDispField);
    if (nullptr != cdfemSnapDispPtr)
      coords += stk::math::Vector3d(oldCdfemSnapDispPtr, theSpatialDim) - stk::math::Vector3d(cdfemSnapDispPtr, theSpatialDim);
  }
  return coords;
}

stk::math::Vector3d ProlongationPointData::get_post_snap_coordinates() const
{
  STK_ThrowAssertMsg(theCoordsField.valid(), "Static member coordinates field is not yet set.");
  stk::math::Vector3d coords(get_field_data(theCoordsField), theSpatialDim);
  return coords;
}

void ProlongationPointData::set_coords_fields(const int spatialDim, FieldRef coordsField, FieldRef snapDisplacementsField)
{
  STK_ThrowRequireMsg(coordsField.valid(), "Invalid coordinates field in ProlongationPointData::set_coords_fields()");
  theSpatialDim = spatialDim;
  theCoordsField = coordsField;
  theSnapDisplacementsField = snapDisplacementsField;
}

ProlongationElementData::ProlongationElementData(const CDMesh & cdmesh, const stk::mesh::Entity element)
: ProlongationData(),
  myMasterElem(MasterElementDeterminer::getMasterElement(cdmesh.stk_bulk().bucket(element).topology()))
{
  const stk::mesh::BulkData & mesh = cdmesh.stk_bulk();
  const StkMeshEntities elemNodes{mesh.begin_nodes(element), mesh.end_nodes(element)};
  myElemNodesData.reserve(elemNodes.size());
  for (auto node : elemNodes)
    myElemNodesData.push_back(cdmesh.fetch_prolong_node(mesh.identifier(node)));
}

void ProlongationElementData::fill_integration_weights(std::vector<double> & childIntgWeights) const
{
  STK_ThrowAssert(have_prolongation_data_stored_for_all_nodes());
  const unsigned dim = myMasterElem.topology_dimension();
  std::vector<double> flatCoords(myElemNodesData.size()*dim);
  for (size_t n=0; n<myElemNodesData.size(); ++n)
  {
    const stk::math::Vector3d nodeCoords = myElemNodesData[n]->get_post_snap_coordinates();
    for (unsigned d=0; d<dim; ++d)
      flatCoords[n*dim+d] = nodeCoords[d];
  }

  ElementObj::integration_weights(childIntgWeights, dim, flatCoords, myMasterElem, myMasterElem);
}

void
ProlongationElementData::evaluate_prolongation_field(const CDFEM_Support & cdfemSupport, const FieldRef field, const unsigned field_length, const stk::math::Vector3d & paramCoords, double * result) const
{
  // Figuring out the field master element here is actually quite hard since the entity may not exist any more.
  // We'll assume that the field master element is the master_elem or the one with topology master_elem->get_topology().base().
  // This will handle the Q2Q1 case.

  for (unsigned i=0; i<field_length; ++i) result[i] = 0.0;
  FieldRef initial_field;

  const MasterElement* calcMasterElem = &myMasterElem;
  const int elemNPE = myMasterElem.get_topology().num_nodes();
  std::vector<const double *> node_data(elemNPE, nullptr);
  const std::vector<double> zeros(field_length,0.0);

  for ( int n = 0; n < elemNPE; n++ )
  {
    const ProlongationNodeData * prolong_data = myElemNodesData[n];
    if (nullptr != prolong_data) node_data[n] = prolong_data->get_field_data(field);
    if (node_data[n] == nullptr)
    {
      if (!initial_field.valid())
      {
        initial_field = cdfemSupport.get_initial_prolongation_field( field );
      }
      if (initial_field.valid())
      {
        node_data[n] = prolong_data->get_field_data(initial_field);
      }
    }
    if (node_data[n] == nullptr)
    {
      calcMasterElem = &MasterElementDeterminer::getMasterElement(myMasterElem.get_topology().base());
      node_data[n] = zeros.data();
    }
  }

  const int fieldNPE = calcMasterElem->get_topology().num_nodes();
  std::vector<double> shapefcn (fieldNPE,0.);
  calcMasterElem->shape_fcn(1, paramCoords.data(), shapefcn.data());

  for ( int n = 0; n < fieldNPE; n++ )
  {
    STK_ThrowRequire(nullptr != node_data[n]);
    for (unsigned i=0; i<field_length; ++i) result[i] += shapefcn[n]*node_data[n][i];
  }
}

ProlongationLeafElementData::ProlongationLeafElementData(const CDMesh & cdmesh, const PartAndFieldCollections & partAndFieldCollections, const stk::mesh::Entity element)
: ProlongationElementData(cdmesh, element)
{
  const stk::mesh::BulkData& mesh = cdmesh.stk_bulk();
  const int fieldCollectionId = partAndFieldCollections.get_field_collection_id(mesh, mesh.bucket(element));
  save_field_data(mesh, partAndFieldCollections, fieldCollectionId, element);

  if (krinolog.shouldPrint(LOG_DEBUG))
    krinolog << debug_output(mesh, element) << stk::diag::dendl;
}

void ProlongationLeafElementData::find_subelement_and_parametric_coordinates_at_point(const stk::math::Vector3d & pointCoordinates, const ProlongationElementData *& interpElem, stk::math::Vector3d & interpElemParamCoords) const
{
  interpElem = this;
  interpElemParamCoords = compute_parametric_coords_at_point(pointCoordinates);
}

ProlongationParentElementData::ProlongationParentElementData(const CDMesh & cdmesh,
    const stk::mesh::Entity element,
    const std::vector<const ProlongationElementData *> & subelementsData,
    const bool doStoreElementFields)
: ProlongationElementData(cdmesh, element),
  mySubelementsData(subelementsData)
{
  STK_ThrowAssertMsg(have_prolongation_data_stored_for_all_nodes_of_subelements(), "Missing prolongation data for one or more nodes of subelement of parent element " << debug_entity_1line(cdmesh.stk_bulk(), element));

  if (doStoreElementFields)
    homogenize_subelement_fields(cdmesh, subelementsData);

  if (krinolog.shouldPrint(LOG_DEBUG))
    krinolog << debug_output(cdmesh.stk_bulk(), element) << stk::diag::dendl;
}

static std::vector<std::vector<double>> calculate_children_integration_weights(const std::vector<const ProlongationElementData *> & subelementsData)
{
  std::vector<std::vector<double>> childIntgWts(subelementsData.size());
  for (size_t iSub=0; iSub<subelementsData.size(); ++iSub)
    subelementsData[iSub]->fill_integration_weights(childIntgWts[iSub]);
  return childIntgWts;
}

bool ProlongationParentElementData::have_prolongation_data_stored_for_all_nodes_of_subelements() const
{
  for (auto && subelementData : mySubelementsData)
    if (!subelementData->have_prolongation_data_stored_for_all_nodes())
      return false;
  return true;
}

static bool any_child_has_field(const std::vector<const ProlongationElementData *> & subelementsData, const FieldRef field)
{
  for (auto && subelemData : subelementsData)
    if (nullptr != subelemData->get_field_data(field))
      return true;
  return false;
}

void ProlongationParentElementData::homogenize_subelement_fields(const CDMesh & cdmesh,
    const std::vector<const ProlongationElementData *> & subelementsData)
{
  const std::vector<std::vector<double>> childIntgWts = calculate_children_integration_weights(subelementsData);
  STK_ThrowAssert(!subelementsData.empty() && subelementsData.size() == childIntgWts.size());

  const unsigned num_intg_pts = childIntgWts[0].size();

  // homogenize fields
  const stk::mesh::FieldVector & allFields = cdmesh.stk_meta().get_fields();
  myFieldStorageIndices.resize(allFields.size(), -1);
  std::vector<double> fieldData;
  for ( auto && fieldPtr : allFields)
  {
    const FieldRef field = *fieldPtr;

    if (field.entity_rank()==stk::topology::ELEMENT_RANK && field.type_is<double>() && any_child_has_field(subelementsData, field))
    {
      const unsigned field_length = field.length();
      const unsigned field_data_index = fieldData.size();
      myFieldStorageIndices[field.field().mesh_meta_data_ordinal()] = field_data_index;
      fieldData.resize(field_data_index+field_length, 0.0);

      // TODO: Add a method to distinguish between vector fields and gauss point fields.
      const bool data_is_gauss_pt_field = (num_intg_pts == field_length);

      if (data_is_gauss_pt_field)
      {
        double tot_child_sum = 0.;
        double tot_child_vol = 0.;

        for (unsigned n = 0; n < subelementsData.size(); ++n)
        {
          const double * subelem_field_data = subelementsData[n]->get_field_data(field);
          if(NULL == subelem_field_data)
          {
            continue;
          }

          const std::vector<double> & child_intg_wts = childIntgWts[n];
          STK_ThrowAssertMsg(child_intg_wts.size() == num_intg_pts, "Children have different integration rules.");

          for (unsigned j=0; j<num_intg_pts; ++j)
          {
            tot_child_sum += subelem_field_data[j] * child_intg_wts[j];
            tot_child_vol += child_intg_wts[j];
          }
        }

        const double tot_child_avg = tot_child_sum / tot_child_vol;
        for (unsigned i=0; i<field_length; ++i)
        {
          fieldData[field_data_index+i] = tot_child_avg;
        }
      }
      else // vector field (includes scalar case)
      {
        double tot_child_vol = 0.;

        for (unsigned n = 0; n < subelementsData.size(); ++n )
        {
          const double * subelem_field_data = subelementsData[n]->get_field_data(field);
          if(NULL == subelem_field_data)
          {
            continue;
          }

          const std::vector<double> & child_intg_wts = childIntgWts[n];
          // We could relax this assertion if we had another way to distinguish gauss point fields from vector fields
          STK_ThrowAssertMsg(child_intg_wts.size() == num_intg_pts, "Children have different integration rules.");

          double child_vol = 0.;
          for (unsigned j=0; j<num_intg_pts; ++j)
          {
            child_vol += child_intg_wts[j];
          }
          tot_child_vol += child_vol;

          for (unsigned i=0; i<field_length; ++i)
          {
            fieldData[field_data_index+i] += child_vol * subelem_field_data[i];
          }
        }

        for (unsigned i=0; i<field_length; ++i)
        {
          fieldData[field_data_index+i] /= tot_child_vol;
        }
      }
    }
  }

  save_field_data(myFieldStorageIndices, fieldData);
}

bool ProlongationElementData::have_prolongation_data_stored_for_all_nodes() const
{
  for (auto && nodeData : myElemNodesData)
    if (!nodeData)
      return false;
  return true;
}

stk::math::Vector3d
ProlongationElementData::compute_parametric_coords_at_point(const stk::math::Vector3d & pointCoords) const
{
  stk::topology baseTopo = myMasterElem.get_topology().base();
  STK_ThrowAssertMsg(have_prolongation_data_stored_for_all_nodes(), "Missing prolongation data at node for prolongation at point " << pointCoords);

  std::vector<stk::math::Vector3d> baseElemNodeCoords;
  baseElemNodeCoords.reserve(baseTopo.num_nodes());
  for (unsigned n=0; n<baseTopo.num_nodes(); ++n)
    baseElemNodeCoords.push_back(myElemNodesData[n]->get_post_snap_coordinates());

  return get_parametric_coordinates_of_point(baseElemNodeCoords, pointCoords);
}

void ProlongationParentElementData::find_subelement_and_parametric_coordinates_at_point(const stk::math::Vector3d & pointCoordinates, const ProlongationElementData *& interpElem, stk::math::Vector3d & interpElemParamCoords) const
{
  STK_ThrowRequire(!mySubelementsData.empty());

  double minSqrDist = std::numeric_limits<double>::max();
  for (auto && subelemData : mySubelementsData)
  {
    const stk::math::Vector3d currentElemParamCoords = subelemData->compute_parametric_coords_at_point(pointCoordinates);
    const double currentChildSqrDist = compute_parametric_square_distance(currentElemParamCoords);
    if (currentChildSqrDist < minSqrDist)
    {
      minSqrDist = currentChildSqrDist;
      interpElem = subelemData;
      interpElemParamCoords = currentElemParamCoords;
    }
  }
}

void
ProlongationData::save_field_data(const stk::mesh::BulkData& stk_mesh, const PartAndFieldCollections & partAndFieldCollections, const int fieldCollectionId, const stk::mesh::Entity entity)
{
  const stk::mesh::FieldVector & allFields = stk_mesh.mesh_meta_data().get_fields();
  myFieldStorageIndices = &partAndFieldCollections.get_field_storage_indices(fieldCollectionId);
  myFieldData.clear();
  myFieldData.reserve(partAndFieldCollections.get_field_storage_size(fieldCollectionId));

  for ( auto fieldIndex : partAndFieldCollections.get_fields(fieldCollectionId) )
  {
    const FieldRef field = allFields[fieldIndex];

    STK_ThrowRequireMsg( field.entity_rank() == stk_mesh.entity_rank(entity) && field.type_is<double>(),
        "Error in prolongation field data storage.  Field " << field.name() << " has rank " << field.entity_rank() << " is double = " << field.type_is<double>() << " on " << stk_mesh.entity_key(entity));

    double * val = field_data<double>(field, entity);
    STK_ThrowRequireMsg( nullptr != val, "Error in prolongation field data storage.  Field " << field.name() << " is missing on " << stk_mesh.entity_key(entity));

    const unsigned field_length = field.length();
    STK_ThrowRequireMsg( static_cast<size_t>((*myFieldStorageIndices)[fieldIndex]) == myFieldData.size(), "Error in prolongation field data storage");
    for (unsigned i=0; i<field_length; ++i)
    {
      myFieldData.push_back(val[i]);
    }
  }
}

ProlongationNodeData::ProlongationNodeData(const CDMesh & mesh, const PartAndFieldCollections & partAndFieldCollections, const stk::mesh::Entity node, bool communicate_me_to_all_sharers)
  : ProlongationPointData(),
    my_entityId(mesh.stk_bulk().identifier(node)),
    myCommunicateMeToAllSharersFlag(communicate_me_to_all_sharers)
{
  const stk::mesh::BulkData& stk_mesh = mesh.stk_bulk();

  myPartCollectionId = partAndFieldCollections.get_part_collection_id(stk_mesh, stk_mesh.bucket(node));
  myFieldCollectionId = partAndFieldCollections.get_field_collection_id(stk_mesh, stk_mesh.bucket(node));
  save_field_data(stk_mesh, partAndFieldCollections, myFieldCollectionId, node);

  if (krinolog.shouldPrint(LOG_DEBUG))
    krinolog << debug_output(mesh.stk_bulk(), node) << stk::diag::dendl;
}

ProlongationNodeData::ProlongationNodeData( const PartAndFieldCollections & partAndFieldCollections, const stk::mesh::EntityId in_entityId, const int partCollectionId, const int fieldCollectionId, std::vector<double> & fieldData)
: ProlongationPointData(),
  my_entityId(in_entityId),
  myCommunicateMeToAllSharersFlag(false),
  myPartCollectionId(partCollectionId),
  myFieldCollectionId(fieldCollectionId)
{
  save_field_data(partAndFieldCollections.get_field_storage_indices(fieldCollectionId), fieldData);
}

static bool does_any_node_have_field(const FieldRef field, const std::vector<const ProlongationNodeData *> & nodes)
{
  for (auto && node : nodes)
    if (NULL != node->get_field_data(field))
      return true;
  return false;
}

ProlongationFacetPointData::ProlongationFacetPointData(const CDMesh & mesh,
    const FacetDistanceQuery<Facet> & facetDistanceQuery,
    const std::vector<const ProlongationNodeData *> & facetNodes)
{
  interpolate_to_point(mesh.stk_meta(), facetDistanceQuery, facetNodes);
}

void ProlongationFacetPointData::interpolate_to_point(const stk::mesh::MetaData & meta,
    const FacetDistanceQuery<Facet> & facetDistanceQuery,
    const std::vector<const ProlongationNodeData *> & facetNodes)
{
  const stk::math::Vector3d node_wts = facetDistanceQuery.closest_point_weights();

  // interpolate fields
  const stk::mesh::FieldVector & all_fields = meta.get_fields();
  myFieldStorageIndices.resize(all_fields.size(), -1);
  std::vector<double> fieldData;
  for ( stk::mesh::FieldVector::const_iterator it = all_fields.begin(); it != all_fields.end() ; ++it )
  {
    const FieldRef field = **it;

    if( field.entity_rank()!=stk::topology::NODE_RANK || !field.type_is<double>() ) continue;

    const unsigned field_length = field.length();

    if (does_any_node_have_field(field, facetNodes))
    {
      const unsigned faceDataIndex = fieldData.size();
      myFieldStorageIndices[field.field().mesh_meta_data_ordinal()] = faceDataIndex;
      fieldData.resize(faceDataIndex+field_length, 0.0);

      double node_wt_sum = 0.0;
      for (unsigned n = 0; n < facetNodes.size(); ++n)
      {
        const double * node_field_data = facetNodes[n]->get_field_data(field);

        if(node_field_data)
        {
          node_wt_sum += node_wts[n];
          for (unsigned i=0; i<field_length; ++i)
          {
            fieldData[faceDataIndex+i] += node_wts[n] * node_field_data[i];
          }
        }
      }

      for(unsigned i=0; i<field_length; ++i)
      {
        fieldData[faceDataIndex+i] /= node_wt_sum;
      }
    }
  }

  save_field_data(myFieldStorageIndices, fieldData);
}

static
void pack_nodes_that_need_to_be_communicated_to_sharers(const stk::mesh::BulkData & mesh, const EntityProlongationNodeMap & proc_prolong_nodes, stk::CommSparse &commSparse)
{
  std::vector<int> nodeSharedProcs;
  stk::pack_and_communicate(commSparse,[&]()
  {
    for (auto entry : proc_prolong_nodes)
    {
      const ProlongationNodeData & prolongNode = *entry.second;
      if (prolongNode.communicate_me_to_all_sharers())
      {
        stk::mesh::Entity node = mesh.get_entity(stk::topology::NODE_RANK, prolongNode.entityId());
        mesh.comm_shared_procs(node, nodeSharedProcs);
        for (int procId : nodeSharedProcs)
          prolongNode.pack_into_buffer(commSparse.send_buffer(procId));
      }
    }
  });
}

static
void pack_facet_prolong_nodes(const stk::mesh::BulkData & mesh, const ProlongFacetVec & proc_prolong_facets, const std::vector<BoundingBox> & proc_target_bboxes, stk::CommSparse &commSparse)
{
  stk::pack_and_communicate(commSparse,[&]()
  {
    for ( int procId=0; procId<commSparse.parallel_size(); ++procId )
    {
      if ( commSparse.parallel_rank() == procId ) continue;  // Don't talk to yourself, it's embarrassing
      stk::CommBuffer & buffer = commSparse.send_buffer(procId);

      for ( auto && node : ProlongationFacet::get_facet_nodes_to_communicate(proc_prolong_facets, proc_target_bboxes[procId]) )
        node->pack_into_buffer(buffer);
    }
  });
}

static
void receive_and_build_prolong_nodes(const PartAndFieldCollections & partAndFieldCollections, EntityProlongationNodeMap & proc_prolong_nodes, stk::CommSparse &commSparse)
{
  stk::unpack_communications(commSparse, [&](int procId)
  {
    stk::CommBuffer & buffer = commSparse.recv_buffer(procId);

    while ( buffer.remaining() )
    {
      ProlongationNodeData * node = ProlongationNodeData::unpack_from_buffer( buffer, partAndFieldCollections ); // This calls new to create a new ProlongationNodeData
      EntityProlongationNodeMap::iterator it = proc_prolong_nodes.find(node->entityId());
      if( it == proc_prolong_nodes.end() || it->second == nullptr )
      {
        proc_prolong_nodes[node->entityId()] = node;
      }
      else
      {
        delete node;
      }
    }
  });
}

void ProlongationFacet::communicate_shared_nodes( const CDMesh & mesh, EntityProlongationNodeMap & proc_prolong_nodes )
{
  stk::CommSparse commSparse(mesh.stk_bulk().parallel());
  pack_nodes_that_need_to_be_communicated_to_sharers(mesh.stk_bulk(), proc_prolong_nodes, commSparse);
  receive_and_build_prolong_nodes(mesh.get_prolong_part_and_field_collections(), proc_prolong_nodes, commSparse);
}


void
ProlongationFacet::communicate_facet_nodes( const CDMesh & mesh, const ProlongFacetVec & proc_prolong_facets, EntityProlongationNodeMap & proc_prolong_nodes, const std::vector<BoundingBox> & proc_target_bboxes )
{
  const int num_procs = mesh.stk_bulk().parallel_size();
  if ( num_procs == 1 ) return;  // Don't talk to yourself, it's embarrassing

  stk::CommSparse commSparse(mesh.stk_bulk().parallel());
  pack_facet_prolong_nodes(mesh.stk_bulk(), proc_prolong_facets, proc_target_bboxes, commSparse);
  receive_and_build_prolong_nodes(mesh.get_prolong_part_and_field_collections(), proc_prolong_nodes, commSparse);
}

void
ProlongationNodeData::pack_into_buffer(stk::CommBuffer & b) const
{
  b.pack(my_entityId);
  b.pack(myPartCollectionId);
  b.pack(myFieldCollectionId);

  const size_t numFieldData = get_field_data().size();
  b.pack(numFieldData);
  b.pack(get_field_data().data(), numFieldData);
}

ProlongationNodeData *
ProlongationNodeData::unpack_from_buffer( stk::CommBuffer & b, const PartAndFieldCollections & partAndFieldCollections )
{
  stk::mesh::EntityId global_id;
  b.unpack(global_id);

  int partCollectionId;
  b.unpack(partCollectionId);

  int fieldCollectionId;
  b.unpack(fieldCollectionId);

  size_t numFieldData = 0;
  b.unpack(numFieldData);
  std::vector<double> fieldData;
  fieldData.resize(numFieldData);
  b.unpack(fieldData.data(), numFieldData);

  ProlongationNodeData * node = new ProlongationNodeData(partAndFieldCollections, global_id, partCollectionId, fieldCollectionId, fieldData);

  return node;
}

const double * ProlongationData::get_field_data( const stk::mesh::FieldBase& state_field ) const
{
  const int field_index = (*myFieldStorageIndices)[state_field.mesh_meta_data_ordinal()]; return (field_index < 0) ? nullptr : &myFieldData[field_index];
}

static bool do_all_nodes_have_same_field_collection_id(const std::vector<const ProlongationNodeData *> & prolongNodes)
{
  for (unsigned iNode=1; iNode<prolongNodes.size(); ++iNode)
    if (prolongNodes[iNode]->get_field_collection_id() != prolongNodes[0]->get_field_collection_id())
      return false;
  return true;
}

static std::vector<unsigned> compute_fields_common_to_nodes(const PartAndFieldCollections & partAndFieldCollections, const std::vector<const ProlongationNodeData *> & prolongNodes)
{
  STK_ThrowAssert(!prolongNodes.empty());
  if (do_all_nodes_have_same_field_collection_id(prolongNodes))
    return partAndFieldCollections.get_fields(prolongNodes[0]->get_field_collection_id());

  std::vector<unsigned> commonFields;

  for (unsigned side_node_index=0; side_node_index<prolongNodes.size(); ++side_node_index)
  {
    const std::vector<unsigned> & nodeFields = partAndFieldCollections.get_fields(prolongNodes[side_node_index]->get_field_collection_id());
    if (0 == side_node_index)
    {
      commonFields = nodeFields;
    }
    else
    {
      std::vector<unsigned> working_set;
      working_set.swap(commonFields);
      std::set_intersection(working_set.begin(),working_set.end(),nodeFields.begin(),nodeFields.end(),std::back_inserter(commonFields));
    }
  }

  return commonFields;
}

std::vector<unsigned> ProlongationFacet::compute_common_fields() const
{
  return compute_fields_common_to_nodes(my_mesh.get_prolong_part_and_field_collections(), my_prolong_nodes);
}

static bool does_edge_have_field_not_on_facet(const PartAndFieldCollections & partAndFieldCollections, const std::vector<const ProlongationNodeData *> & edgeProlongNodes, const std::vector<unsigned> & facetFields)
{
  const auto edgeFields = compute_fields_common_to_nodes(partAndFieldCollections, edgeProlongNodes);
  for (auto && edgeField : edgeFields)
    if (!std::binary_search(facetFields.begin(), facetFields.end(), edgeField))
      return true;
  return false;
}

ProlongationFacet::ProlongationFacet(const CDMesh & mesh, stk::mesh::Entity side)
: my_mesh(mesh)
{
  const stk::mesh::BulkData & stk_mesh = my_mesh.stk_bulk();

  STK_ThrowAssert(stk_mesh.num_elements(side) > 0);
  stk::mesh::Entity elem0 = stk_mesh.begin_elements(side)[0];
  const PhaseTag elem0_phase = mesh.determine_entity_phase(elem0);

  const unsigned num_side_nodes = stk_mesh.bucket(side).topology().base().num_nodes();
  const stk::mesh::Entity* side_nodes = stk_mesh.begin_nodes(side);
  my_prolong_nodes.resize(num_side_nodes);

  for (unsigned side_node_index=0; side_node_index<num_side_nodes; ++side_node_index)
  {
    stk::mesh::Entity side_node = side_nodes[side_node_index];
    const ProlongationNodeData * prolong_node = mesh.fetch_prolong_node(stk_mesh.identifier(side_node));

    if (NULL == prolong_node)
    {
      krinolog << stk::diag::dendl;
      krinolog << "Failed to find node " << stk_mesh.identifier(side_node) << " while stashing facet " << debug_entity(mesh.stk_bulk(), side);
    }
    STK_ThrowRequire(NULL != prolong_node);
    my_prolong_nodes[side_node_index] = prolong_node;
  }

  STK_ThrowAssert((int)my_prolong_nodes.size() == my_mesh.spatial_dim());
  if (2 == my_prolong_nodes.size())
  {
    my_facet = std::make_unique<Facet2d>( my_prolong_nodes[0]->get_previous_coordinates(), my_prolong_nodes[1]->get_previous_coordinates());
  }
  else
  {
    STK_ThrowAssert(3 == my_prolong_nodes.size());
    my_facet = std::make_unique<Facet3d>( my_prolong_nodes[0]->get_previous_coordinates(), my_prolong_nodes[1]->get_previous_coordinates(), my_prolong_nodes[2]->get_previous_coordinates());
  }
}

ProlongationFacet::ProlongationFacet(const CDMesh & mesh, const std::vector<const ProlongationNodeData *> & prolong_nodes)
: my_mesh (mesh),
  my_prolong_nodes(prolong_nodes)
{
  STK_ThrowAssert((int)my_prolong_nodes.size() <= my_mesh.spatial_dim());
  if (2 == my_prolong_nodes.size())
  {
    my_facet = std::make_unique<Facet2d>( my_prolong_nodes[0]->get_previous_coordinates(), my_prolong_nodes[1]->get_previous_coordinates());
  }
  else
  {
    STK_ThrowAssert(3 == my_prolong_nodes.size());
    my_facet = std::make_unique<Facet3d>( my_prolong_nodes[0]->get_previous_coordinates(), my_prolong_nodes[1]->get_previous_coordinates(), my_prolong_nodes[2]->get_previous_coordinates());
  }
}

static void build_and_append_edge_facet_if_it_has_field_not_on_facet(const CDMesh & mesh, const std::vector<unsigned> & facetFields, const std::vector<const ProlongationNodeData *> & edgeProlongNodes, ProlongFacetVec & facetVec)
{
  if (does_edge_have_field_not_on_facet(mesh.get_prolong_part_and_field_collections(), edgeProlongNodes, facetFields))
  {
    ProlongationFacet * edgeFacet = new ProlongationFacet(mesh, edgeProlongNodes);
    facetVec.push_back(edgeFacet);
  }
}

void ProlongationFacet::build_and_append_edge_facets(ProlongFacetVec & facetVec) const
{
  if (my_prolong_nodes.size() == 3 && !do_all_nodes_have_same_field_collection_id(my_prolong_nodes))
  {
    const auto facetFields = compute_common_fields();
    build_and_append_edge_facet_if_it_has_field_not_on_facet(my_mesh, facetFields, {my_prolong_nodes[0], my_prolong_nodes[1]}, facetVec);
    build_and_append_edge_facet_if_it_has_field_not_on_facet(my_mesh, facetFields, {my_prolong_nodes[1], my_prolong_nodes[2]}, facetVec);
    build_and_append_edge_facet_if_it_has_field_not_on_facet(my_mesh, facetFields, {my_prolong_nodes[2], my_prolong_nodes[0]}, facetVec);
  }
}

std::unique_ptr<ProlongationFacetPointData> ProlongationFacet::get_prolongation_point_data(const FacetDistanceQuery<Facet> & dist_query) const
{
  STK_ThrowAssert(&dist_query.facet() == my_facet.get());
  return std::make_unique<ProlongationFacetPointData>(my_mesh, dist_query, my_prolong_nodes);
}

bool ProlongationFacet::communicate_me(const BoundingBox & proc_target_bbox) const
{
  for(auto && prolong_node : my_prolong_nodes)
  {
    if (!proc_target_bbox.contains(prolong_node->get_previous_coordinates()))
    {
      return false;
    }
  }
  return true;
}

std::set<const ProlongationNodeData *>
ProlongationFacet::get_facet_nodes_to_communicate( const ProlongFacetVec & proc_prolong_facets, const BoundingBox & proc_target_bbox )
{
  std::set<const ProlongationNodeData *> procFacetNodes;

  for ( auto && facet : proc_prolong_facets )
    if (facet->communicate_me(proc_target_bbox))
      for (auto node : facet->my_prolong_nodes)
        procFacetNodes.insert(node);

  return procFacetNodes;
}

void
ProlongationFacet::communicate( const CDMesh & mesh, ProlongFacetVec & proc_prolong_facets, EntityProlongationNodeMap & proc_prolong_nodes, const std::vector<BoundingBox> & proc_target_bboxes )
{
  communicate_shared_nodes(mesh, proc_prolong_nodes);
  communicate_facet_nodes(mesh, proc_prolong_facets, proc_prolong_nodes, proc_target_bboxes);
  communicate_facets(mesh, proc_prolong_facets, proc_target_bboxes);
}

static
void pack_facets_within_proc_bboxes(const stk::mesh::BulkData & mesh, const ProlongFacetVec & proc_prolong_facets, const std::vector<BoundingBox> & proc_target_bboxes, stk::CommSparse &commSparse)
{
  stk::pack_and_communicate(commSparse,[&]()
  {
    for ( int procId=0; procId<commSparse.parallel_size(); ++procId )
    {
      if ( commSparse.parallel_rank() == procId ) continue;  // Don't talk to yourself, it's embarrassing
      const BoundingBox & proc_target_bbox = proc_target_bboxes[procId];
      stk::CommBuffer & buffer = commSparse.send_buffer(procId);

      for ( auto && f_i : proc_prolong_facets )
      {
        const ProlongationFacet & facet = *f_i;
        if (facet.communicate_me(proc_target_bbox))
          facet.pack_into_buffer(buffer);
      }
    }
  });
}

static
void receive_and_build_prolong_facets(const CDMesh & mesh, ProlongFacetVec & proc_prolong_facets, stk::CommSparse &commSparse)
{
  stk::unpack_communications(commSparse, [&](int procId)
  {
    stk::CommBuffer & buffer = commSparse.recv_buffer(procId);

    while ( buffer.remaining() )
    {
      ProlongationFacet * facet = ProlongationFacet::unpack_from_buffer( mesh, buffer );
      proc_prolong_facets.push_back( facet );
    }
  });
}

void ProlongationFacet::communicate_facets( const CDMesh & mesh, ProlongFacetVec & proc_prolong_facets, const std::vector<BoundingBox> & proc_target_bboxes )
{
  stk::CommSparse commSparse(mesh.stk_bulk().parallel());
  pack_facets_within_proc_bboxes(mesh.stk_bulk(), proc_prolong_facets, proc_target_bboxes, commSparse);
  receive_and_build_prolong_facets(mesh, proc_prolong_facets, commSparse);
}

void
ProlongationFacet::pack_into_buffer(stk::CommBuffer & b) const
{
  const size_t num_prolong_nodes = my_prolong_nodes.size();
  b.pack(num_prolong_nodes);
  for(auto && prolong_node : my_prolong_nodes)
  {
    b.pack(prolong_node->entityId());
  }
}

ProlongationFacet *
ProlongationFacet::unpack_from_buffer(const CDMesh & mesh, stk::CommBuffer & b )
{
  size_t num_prolong_nodes = 0;
  b.unpack(num_prolong_nodes);
  std::vector<const ProlongationNodeData *> prolong_nodes(num_prolong_nodes);

  for(auto && prolong_node : prolong_nodes)
  {
    stk::mesh::EntityId node_id;
    b.unpack(node_id);
    prolong_node = mesh.fetch_prolong_node(node_id);
    STK_ThrowRequireMsg(prolong_node, "Communication error, missing prolongation node " << node_id << " on processor " << mesh.stk_bulk().parallel_rank());
  }

  ProlongationFacet * facet = new ProlongationFacet(mesh, prolong_nodes);

  return facet;
}

std::string ProlongationData::debug_data_output(const stk::mesh::BulkData & mesh) const
{
  std::ostringstream os;
  const stk::mesh::FieldVector & allFields = mesh.mesh_meta_data().get_fields();
  for ( auto && fieldPtr : allFields)
  {
    const int fieldIndex = (*myFieldStorageIndices)[fieldPtr->mesh_meta_data_ordinal()];
    if (fieldIndex >= 0)
    {
      const FieldRef field = *fieldPtr;
      os << "  Field: field_name=" << field.name() << ", field_state=" << static_cast<int>(field.state()) << ", ";
      const unsigned fieldLength = field.length();
      if (1 == fieldLength)
      {
        os << "value=" << myFieldData[fieldIndex] << "\n";
      }
      else
      {
        os << "values[] = ";
        for (unsigned i = 0; i < fieldLength; ++i) os << myFieldData[fieldIndex+i] << " ";
        os << "\n";
      }
    }
  }
  return os.str();
}

std::string ProlongationParentElementData::debug_output(const stk::mesh::BulkData & mesh, const stk::mesh::Entity element) const
{
  std::ostringstream os;
  os << "ProlongationParentElementData for " << mesh.identifier(element) << "\n";
  os << debug_data_output(mesh);
  return os.str();
}

std::string ProlongationLeafElementData::debug_output(const stk::mesh::BulkData & mesh, const stk::mesh::Entity element) const
{
  std::ostringstream os;
  os << "ProlongationLeafElementData for " << mesh.identifier(element) << "\n";
  os << debug_data_output(mesh);
  return os.str();
}

std::string ProlongationNodeData::debug_output(const stk::mesh::BulkData & mesh, const stk::mesh::Entity node) const
{
  std::ostringstream os;
  os << "ProlongationNodeData for " << mesh.identifier(node) << " at previous loc " << get_previous_coordinates() << " at post snap loc " << get_post_snap_coordinates() << "\n";
  os << debug_data_output(mesh);
  return os.str();
}

} // namespace krino
