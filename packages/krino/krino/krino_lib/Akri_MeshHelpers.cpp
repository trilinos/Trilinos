// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <stk_mesh/baseImpl/elementGraph/ElemElemGraph.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/MeshUtils.hpp>
#include <stk_mesh/base/FieldParallel.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>
#include <stk_io/IossBridge.hpp>
#include <stk_util/parallel/ParallelReduce.hpp> // Needed for all_reduce_max
#include <stk_util/parallel/CommSparse.hpp>

#include <Akri_DiagWriter.hpp>
#include <Akri_EntityIdPool.hpp>
#include <Akri_MeshHelpers.hpp>
#include <Akri_ParallelCommHelpers.hpp>
#include <Akri_ParallelErrorMessage.hpp>
#include <Akri_AuxMetaData.hpp>
#include <Akri_FieldRef.hpp>
#include <Akri_Vec.hpp>
#include <stk_util/environment/Env.hpp>
#include "Akri_ReportHandler.hpp"

namespace krino{

size_t get_global_num_entities(const stk::mesh::BulkData& mesh, stk::mesh::EntityRank entityRank)
{
  size_t numEntities = stk::mesh::count_selected_entities(mesh.mesh_meta_data().locally_owned_part(), mesh.buckets(entityRank));
  const size_t localNumEntities = numEntities;
  stk::all_reduce_sum(mesh.parallel(), &localNumEntities, &numEntities, 1);
  return numEntities;
}

template <class CONTAINER>
class ContainerResizer
{
public:
  ContainerResizer(CONTAINER & container) : myContainer(container) {}
  void resize(size_t size) { myContainer.resize(size); }
private:
  CONTAINER & myContainer;
};

template <class T, size_t SIZE>
class ContainerResizer<std::array<T, SIZE>>
{
public:
  ContainerResizer(std::array<T, SIZE> & container) : myContainer(container) {}
  void resize(size_t size) { ThrowRequire(SIZE == size); }
private:
  std::array<T, SIZE> & myContainer;
};

template <class CONTAINER>
void resize_container(CONTAINER & container, size_t size)
{
  ContainerResizer<CONTAINER> resizer(container);
  resizer.resize(size);
}

Vector3d get_side_normal(const stk::mesh::BulkData& mesh,
    const FieldRef coordsField,
    stk::mesh::Entity side)
{
  const auto * sideNodes = mesh.begin_nodes(side);
  const stk::topology sideTopology = mesh.bucket(side).topology();
  if (sideTopology == stk::topology::TRIANGLE_3 || sideTopology == stk::topology::TRIANGLE_6)
  {
    const Vector3d v0(field_data<double>(coordsField, sideNodes[0]));
    const Vector3d v1(field_data<double>(coordsField, sideNodes[1]));
    const Vector3d v2(field_data<double>(coordsField, sideNodes[2]));
    return Cross(v1-v0,v2-v0).unit_vector();
  }
  else if (sideTopology == stk::topology::LINE_2 || sideTopology == stk::topology::LINE_3)
  {
    const Vector3d v0(field_data<double>(coordsField, sideNodes[0]), 2);
    const Vector3d v1(field_data<double>(coordsField, sideNodes[1]), 2);
    return crossZ(v1-v0).unit_vector();
  }
  ThrowRequireMsg(false, "Unsupported topology " << sideTopology);

  return Vector3d::ZERO;
}

void fill_procs_owning_or_sharing_or_ghosting_node(const stk::mesh::BulkData& bulkData, stk::mesh::Entity node, std::vector<int> & procsOwningSharingOrGhostingNode)
{
    ThrowAssert(bulkData.parallel_owner_rank(node)==bulkData.parallel_rank());
    bulkData.comm_procs(node, procsOwningSharingOrGhostingNode);
    procsOwningSharingOrGhostingNode.push_back(bulkData.parallel_rank());
}

//--------------------------------------------------------------------------------

template<int DIM, class CONTAINER>
void fill_element_node_coordinates(const stk::mesh::BulkData & mesh, stk::mesh::Entity element, const FieldRef coordsField, CONTAINER & elementNodeCoords)
{
  StkMeshEntities elementNodes{mesh.begin_nodes(element), mesh.end_nodes(element)};
  resize_container(elementNodeCoords, elementNodes.size());
  for (size_t n=0; n<elementNodeCoords.size(); ++n)
  {
    double * coords = field_data<double>(coordsField, elementNodes[n]);
    for ( int d = 0; d < DIM; ++d )
      elementNodeCoords[n][d] = coords[d];
  }
}

void fill_element_node_coordinates(const stk::mesh::BulkData & mesh, stk::mesh::Entity element, const FieldRef coordsField, const int dim, std::vector<Vector3d> & elementNodeCoords)
{
  ThrowAssert(dim == 3 || dim == 2);
  if (dim == 2)
    fill_element_node_coordinates<2>(mesh, element, coordsField, elementNodeCoords);
  else
    fill_element_node_coordinates<3>(mesh, element, coordsField, elementNodeCoords);
}

void fill_element_node_coordinates(const stk::mesh::BulkData & mesh, stk::mesh::Entity element, const FieldRef coordsField, std::vector<Vector3d> & elementNodeCoords)
{
  const int dim = mesh.mesh_meta_data().spatial_dimension();
  elementNodeCoords.clear();
  for (auto node : StkMeshEntities{mesh.begin_nodes(element), mesh.end_nodes(element)})
    elementNodeCoords.emplace_back(field_data<double>(coordsField, node), dim);
}

static std::array<krino::Vector3d,4> gather_tet_coordinates(const stk::mesh::BulkData & mesh, stk::mesh::Entity element, const FieldRef coordsField)
{
  ThrowAssert(mesh.bucket(element).topology() == stk::topology::TETRAHEDRON_4);
  std::array<krino::Vector3d,4> elementNodeCoords;
  fill_element_node_coordinates<3>(mesh, element, coordsField, elementNodeCoords);
  return elementNodeCoords;
}

static std::array<krino::Vector3d,3> gather_tri_coordinates(const stk::mesh::BulkData & mesh, stk::mesh::Entity element, const FieldRef coordsField)
{
  ThrowAssert(mesh.bucket(element).topology() == stk::topology::TRIANGLE_3_2D);
  std::array<krino::Vector3d,3> elementNodeCoords;
  fill_element_node_coordinates<2>(mesh, element, coordsField, elementNodeCoords);
  return elementNodeCoords;
}

double compute_tri_volume(const std::array<krino::Vector3d,3> & elementNodeCoords)
{
  return 0.5*(Cross(elementNodeCoords[1]-elementNodeCoords[0], elementNodeCoords[2]-elementNodeCoords[0]).length());
}

double compute_tri_volume(const krino::Vector3d * elementNodeCoords)
{
  return 0.5*(Cross(elementNodeCoords[1]-elementNodeCoords[0], elementNodeCoords[2]-elementNodeCoords[0]).length());
}

double compute_tri_volume(const std::array<krino::Vector2d,3> & elementNodeCoords)
{
  return 0.5 * (elementNodeCoords[0][0]*(elementNodeCoords[1][1] - elementNodeCoords[2][1]) +
      elementNodeCoords[1][0]*(elementNodeCoords[2][1] - elementNodeCoords[0][1]) +
      elementNodeCoords[2][0]*(elementNodeCoords[0][1] - elementNodeCoords[1][1]));
}

double compute_tet_volume(const std::array<krino::Vector3d,4> & elementNodeCoords)
{
  return Dot(elementNodeCoords[3]-elementNodeCoords[0],Cross(elementNodeCoords[1]-elementNodeCoords[0], elementNodeCoords[2]-elementNodeCoords[0]))/6.0;
}

double compute_tet_volume(const krino::Vector3d * elementNodeCoords)
{
  return Dot(elementNodeCoords[3]-elementNodeCoords[0],Cross(elementNodeCoords[1]-elementNodeCoords[0], elementNodeCoords[2]-elementNodeCoords[0]))/6.0;
}

double compute_tri_or_tet_volume(const std::vector<krino::Vector3d> & elementNodeCoords)
{
  ThrowAssert(elementNodeCoords.size() == 4 || elementNodeCoords.size() == 3);
  if (elementNodeCoords.size() == 4)
    return compute_tet_volume({{elementNodeCoords[0],elementNodeCoords[1],elementNodeCoords[2],elementNodeCoords[3]}});
  return compute_tri_volume({{elementNodeCoords[0],elementNodeCoords[1],elementNodeCoords[2]}});
}

static double compute_tri_or_tet_volume(const stk::mesh::BulkData & mesh, stk::mesh::Entity element, const FieldRef coordsField)
{
  stk::topology elemTopology = mesh.bucket(element).topology();

  if (elemTopology == stk::topology::TETRAHEDRON_4)
  {
    const auto elementNodeCoords = gather_tet_coordinates(mesh, element, coordsField);
    return compute_tet_volume(elementNodeCoords);
  }

  ThrowRequireMsg(elemTopology == stk::topology::TRIANGLE_3_2D, "Topology " << elemTopology << " not supported in compute_tri_or_tet_volume.");
  const auto elementNodeCoords = gather_tri_coordinates(mesh, element, coordsField);
  return compute_tri_volume(elementNodeCoords);
}

static void update_min_max_values(const double currentValue, double & minValue, double & maxValue)
{
  minValue = std::min(minValue, currentValue);
  maxValue = std::max(maxValue, currentValue);
}

template<class CONTAINER>
void update_min_max_edge_lengths_squared(const CONTAINER & elementNodeCoords, double & minEdgeLengthSqr, double & maxEdgeLengthSqr)
{
  for ( size_t inode = 0; inode < elementNodeCoords.size(); ++inode )
    for ( size_t jnode = inode+1; jnode < elementNodeCoords.size(); ++jnode )
      update_min_max_values((elementNodeCoords[inode] - elementNodeCoords[jnode]).length_squared(), minEdgeLengthSqr, maxEdgeLengthSqr);
}

template<class CONTAINER>
void update_max_edge_lengths_squared(const CONTAINER & elementNodeCoords, double & maxEdgeLengthSqr)
{
  // This is a little strange for non-simplex elements since it goes from each node to every other node
  for ( size_t inode = 0; inode < elementNodeCoords.size(); ++inode )
    for ( size_t jnode = inode+1; jnode < elementNodeCoords.size(); ++jnode )
      maxEdgeLengthSqr = std::max(maxEdgeLengthSqr, (elementNodeCoords[inode] - elementNodeCoords[jnode]).length_squared());
}

//--------------------------------------------------------------------------------

double
compute_maximum_element_size(const stk::mesh::BulkData& mesh, const stk::mesh::Selector & selector)
{
  const unsigned ndim = mesh.mesh_meta_data().spatial_dimension();
  double max_sqr_edge_length = 0.0;

  const FieldRef coordsField(mesh.mesh_meta_data().coordinate_field());

  stk::mesh::Selector locally_owned_selector = selector & mesh.mesh_meta_data().locally_owned_part();

  const stk::mesh::BucketVector & buckets = mesh.get_buckets( stk::topology::ELEMENT_RANK, locally_owned_selector );
  std::vector<krino::Vector3d> elementNodeCoords;

  for ( auto && bucket : buckets )
  {
    for ( auto && elem : *bucket )
    {
      fill_element_node_coordinates(mesh, elem, coordsField, ndim, elementNodeCoords);
      update_max_edge_lengths_squared(elementNodeCoords, max_sqr_edge_length);
    }
  }

  const double local_max = max_sqr_edge_length;
  stk::all_reduce_max(mesh.parallel(), &local_max, &max_sqr_edge_length, 1);

  return std::sqrt(max_sqr_edge_length);
}

double
compute_maximum_size_of_selected_elements_using_node(const stk::mesh::BulkData& mesh, const stk::mesh::Selector & selector, const stk::mesh::Entity node)
{
  const unsigned ndim = mesh.mesh_meta_data().spatial_dimension();
  double maxSqrEdgeLength = 0.0;

  const FieldRef coordsField(mesh.mesh_meta_data().coordinate_field());

  std::vector<krino::Vector3d> elementNodeCoords;

  for (auto && element : StkMeshEntities{mesh.begin_elements(node), mesh.end_elements(node)})
  {
    if (selector(mesh.bucket(element)))
    {
      fill_element_node_coordinates(mesh, element, coordsField, ndim, elementNodeCoords);
      update_max_edge_lengths_squared(elementNodeCoords, maxSqrEdgeLength);
    }
  }

  return std::sqrt(maxSqrEdgeLength);
}

//--------------------------------------------------------------------------------

void compute_element_quality(const stk::mesh::BulkData & mesh, double & minEdgeLength, double & maxEdgeLength, double & minVolume, double & maxVolume)
{
  double minEdgeLengthSqr = std::numeric_limits<double>::max();
  double maxEdgeLengthSqr = -std::numeric_limits<double>::max();
  minVolume = std::numeric_limits<double>::max();
  maxVolume = std::numeric_limits<double>::lowest();

  const FieldRef coordsField(mesh.mesh_meta_data().coordinate_field());

  stk::mesh::Selector locally_owned_selector = mesh.mesh_meta_data().locally_owned_part();

  const stk::mesh::BucketVector & buckets = mesh.get_buckets( stk::topology::ELEMENT_RANK, locally_owned_selector );

  for ( auto && bucket : buckets )
  {
    stk::topology elem_topology = bucket->topology();
    const unsigned num_nodes = elem_topology.num_nodes();
    std::vector<krino::Vector3d> elem_node_coords(num_nodes);

    for ( auto && elem : *bucket )
    {
      if (elem_topology == stk::topology::TETRAHEDRON_4)
      {
        const auto elementNodeCoords = gather_tet_coordinates(mesh, elem, coordsField);
        update_min_max_values(compute_tet_volume(elementNodeCoords), minVolume, maxVolume);
        update_min_max_edge_lengths_squared(elementNodeCoords, minEdgeLengthSqr, maxEdgeLengthSqr);
      }
      else if (elem_topology == stk::topology::TRIANGLE_3_2D)
      {
        const auto elementNodeCoords = gather_tri_coordinates(mesh, elem, coordsField);
        update_min_max_values(compute_tri_volume(elementNodeCoords), minVolume, maxVolume);
        update_min_max_edge_lengths_squared(elementNodeCoords, minEdgeLengthSqr, maxEdgeLengthSqr);
      }
      else
      {
        ThrowRuntimeError("Topology " << elem_topology << " not supported in compute_element_quality.");
      }
    }
  }

  const double localMinEdgeLength = std::sqrt(minEdgeLengthSqr);
  stk::all_reduce_min(mesh.parallel(), &localMinEdgeLength, &minEdgeLength, 1);
  const double localMaxEdgeLength = std::sqrt(maxEdgeLengthSqr);
  stk::all_reduce_max(mesh.parallel(), &localMaxEdgeLength, &maxEdgeLength, 1);
  const double localMinVolume = minVolume;
  stk::all_reduce_min(mesh.parallel(), &localMinVolume, &minVolume, 1);
  const double localMaxVolume = maxVolume;
  stk::all_reduce_max(mesh.parallel(), &localMaxVolume, &maxVolume, 1);
}

static std::vector<stk::mesh::Entity> get_owned_nodes_with_nodal_volume_below_threshold(const stk::mesh::BulkData & mesh, const stk::mesh::Selector & blockSelector, const double threshold)
{
  ThrowRequireMsg(mesh.is_automatic_aura_on() || mesh.parallel_size() == 1, "Method requires automatic aura.");

  // This would be more efficient if a nodal field was used because it could compute the element volume only once.
  std::vector<stk::mesh::Entity> ownedNodesWithNodalVolBelowThreshold;

  const FieldRef coordsField(mesh.mesh_meta_data().coordinate_field());
  std::vector<krino::Vector3d> elementNodeCoords;

  stk::mesh::Selector ownedBlockSelector = mesh.mesh_meta_data().locally_owned_part() & blockSelector;

  const stk::mesh::BucketVector & buckets = mesh.get_buckets( stk::topology::NODE_RANK, ownedBlockSelector );

  for ( auto && bucket : buckets )
  {
    for ( auto && node : *bucket )
    {
      double nodalVolume = 0.;
      for (auto && element : StkMeshEntities{mesh.begin_elements(node), mesh.end_elements(node)})
      {
        if (blockSelector(mesh.bucket(element)))
          nodalVolume += compute_tri_or_tet_volume(mesh, element, coordsField);

        if (nodalVolume >= threshold)
          break;
      }
      if (nodalVolume < threshold)
        ownedNodesWithNodalVolBelowThreshold.push_back(node);
    }
  }
  return ownedNodesWithNodalVolBelowThreshold;
}

static std::vector<stk::mesh::Entity> get_nodes_with_no_attached_elements(const stk::mesh::BulkData & mesh)
{
  ThrowRequireMsg(mesh.is_automatic_aura_on() || mesh.parallel_size() == 1, "Method requires automatic aura.");

  std::vector<stk::mesh::Entity> nodesWithNoAttachedElements;

  for ( auto && bucket : mesh.buckets(stk::topology::NODE_RANK) )
    for ( auto && node : *bucket )
      if (mesh.num_elements(node) == 0)
        nodesWithNoAttachedElements.push_back(node);

  return nodesWithNoAttachedElements;
}

void pack_entities_for_sharing_procs(const stk::mesh::BulkData & mesh,
    const std::vector<stk::mesh::Entity> & entities,
    stk::CommSparse &commSparse)
{
  std::vector<int> sharingProcs;
  stk::pack_and_communicate(commSparse,[&]()
  {
    for (auto entity : entities)
    {
      if (mesh.bucket(entity).shared())
      {
        mesh.comm_shared_procs(entity, sharingProcs);
        for (int procId : sharingProcs)
          commSparse.send_buffer(procId).pack(mesh.entity_key(entity));
      }
    }
  });
}

void unpack_shared_entities(const stk::mesh::BulkData & mesh,
    std::vector<stk::mesh::Entity> & sharedEntities,
    stk::CommSparse &commSparse)
{
  stk::unpack_communications(commSparse, [&](int procId)
  {
    stk::CommBuffer & buffer = commSparse.recv_buffer(procId);

    while ( buffer.remaining() )
    {
      stk::mesh::EntityKey entityKey;
      commSparse.recv_buffer(procId).unpack(entityKey);
      stk::mesh::Entity entity = mesh.get_entity(entityKey);
      ThrowAssert(mesh.is_valid(entity));
      sharedEntities.push_back(entity);
    }
  });
}

static
void append_shared_entities_to_owned_ones(const stk::mesh::BulkData & mesh,
    std::vector<stk::mesh::Entity> & entities)
{
  stk::CommSparse commSparse(mesh.parallel());
  pack_entities_for_sharing_procs(mesh, entities, commSparse);

  std::vector<stk::mesh::Entity> sharedEntities;
  unpack_shared_entities(mesh, sharedEntities, commSparse);
  entities.insert(entities.end(), sharedEntities.begin(), sharedEntities.end());
}

//--------------------------------------------------------------------------------

void delete_node_and_all_entities_using_it(stk::mesh::BulkData & mesh, const stk::mesh::Entity node)
{
  std::vector<stk::mesh::Entity> relatives;

  const stk::mesh::EntityRank highestEntityRank = static_cast<stk::mesh::EntityRank>(mesh.mesh_meta_data().entity_rank_count()-1);
  for (stk::mesh::EntityRank irank = highestEntityRank; irank != stk::topology::NODE_RANK; --irank)
  {
    relatives.assign(mesh.begin(node, irank), mesh.end(node, irank));
    for (auto && relative : relatives)
      ThrowRequire(mesh.destroy_entity(relative));
  }
  ThrowRequire(mesh.destroy_entity(node));
}

static void delete_nodes_and_all_entities_using_them(stk::mesh::BulkData & mesh, const std::vector<stk::mesh::Entity> & nodesToDelete)
{
  mesh.modification_begin();
  for (auto && node : nodesToDelete)
    delete_node_and_all_entities_using_it(mesh, node);
  mesh.modification_end();
}

static size_t delete_nodes_with_nodal_volume_below_threshold_and_all_entities_using_them(stk::mesh::BulkData & mesh, const stk::mesh::Selector & blockSelector, const double threshold)
{
  std::vector<stk::mesh::Entity> nodesToDelete = get_owned_nodes_with_nodal_volume_below_threshold(mesh, blockSelector, threshold);
  const size_t globalNumNodesToDelete = stk::get_global_sum(mesh.parallel(), nodesToDelete.size());

  if (globalNumNodesToDelete > 0)
  {
    append_shared_entities_to_owned_ones(mesh, nodesToDelete);
    delete_nodes_and_all_entities_using_them(mesh, nodesToDelete);
  }
  return globalNumNodesToDelete;
}

static size_t delete_nodes_with_no_attached_elements(stk::mesh::BulkData & mesh)
{
  const std::vector<stk::mesh::Entity> nodesToDelete = get_nodes_with_no_attached_elements(mesh);
  const size_t globalNumNodesToDelete = stk::get_global_sum(mesh.parallel(), nodesToDelete.size());

  if (globalNumNodesToDelete > 0)
    delete_nodes_and_all_entities_using_them(mesh, nodesToDelete);

  return globalNumNodesToDelete;
}

void delete_all_entities_using_nodes_with_nodal_volume_below_threshold(stk::mesh::BulkData & mesh, const stk::mesh::Selector & blockSelector, const double threshold)
{
  const int maxIterations = 10;
  int iteration = 0;
  while (++iteration <= maxIterations)
  {
    const size_t numNodesDeletedWithSmallVolume = delete_nodes_with_nodal_volume_below_threshold_and_all_entities_using_them(mesh, blockSelector, threshold);
    if (numNodesDeletedWithSmallVolume > 0)
    {
      sierra::Env::outputP0() << "Iteration " << iteration << ":" << std::endl;
      sierra::Env::outputP0() << "  Deleted " << numNodesDeletedWithSmallVolume << " node(s) with a nodal volume less than " << threshold << " (and all attached entities)." << std::endl;
      const size_t numNodesDeletedWithNoElements = delete_nodes_with_no_attached_elements(mesh);
      sierra::Env::outputP0() << "  Deleted " << numNodesDeletedWithNoElements << " node(s) that have no attached elements." << std::endl;
    }
    else
      break;
  }
  if (iteration < maxIterations)
    sierra::Env::outputP0() << "Successfully deleted all nodes with a nodal volume less than " << threshold << "." << std::endl;
  else
    sierra::Env::outputP0() << "Terminating after performing max iterations.  There still may be nodes with a nodal volume less than " << threshold << "." << std::endl;
}

template <class FIRST_CONTAINER_TYPE>
inline
void append_1st_container_to_end_of_2nd_vector(
        const FIRST_CONTAINER_TYPE &firstContainer,
        std::vector<typename FIRST_CONTAINER_TYPE::value_type> &targetVec)
{
    targetVec.insert(targetVec.end(), firstContainer.begin(), firstContainer.end());
}

void get_nodes_of_subentity(const stk::mesh::BulkData & bulk,
    const stk::mesh::Entity entity,
    const stk::mesh::EntityRank subEntityRank,
    const stk::mesh::ConnectivityOrdinal subEntityOrdinal,
    std::vector<stk::mesh::Entity> & subEntityNodes)
{
  const stk::mesh::Entity *entityNodeEntities = bulk.begin_nodes(entity);
  const stk::topology entityTopology = bulk.bucket(entity).topology();
  subEntityNodes.resize(entityTopology.sub_topology(subEntityRank,subEntityOrdinal).num_nodes());
  entityTopology.sub_topology_nodes(entityNodeEntities, subEntityRank, subEntityOrdinal, subEntityNodes.data());
}

void get_nodes_of_element_side(const stk::mesh::BulkData & bulk,
    const stk::mesh::Entity element,
    const stk::mesh::ConnectivityOrdinal elementSideOrdinal,
    std::vector<stk::mesh::Entity> & elementSideNodes)
{
  get_nodes_of_subentity(bulk, element, bulk.mesh_meta_data().side_rank(), elementSideOrdinal, elementSideNodes);
}

bool does_element_side_exist(stk::mesh::BulkData& mesh, stk::mesh::Entity element, stk::mesh::ConnectivityOrdinal side_ordinal)
{
  stk::mesh::Entity side = stk::mesh::Entity();
  stk::mesh::EntityRank side_rank = mesh.mesh_meta_data().side_rank();

  unsigned elem_num_sides = mesh.num_connectivity(element, side_rank);
  const stk::mesh::Entity * elem_sides = mesh.begin(element, side_rank);
  const stk::mesh::ConnectivityOrdinal * elem_ord_it = mesh.begin_ordinals(element, side_rank);
  for (unsigned i=0 ; i<elem_num_sides ; ++i)
  {
    if (elem_ord_it[i] == side_ordinal)
    {
      side = elem_sides[i];
      break;
    }
  }

  return mesh.is_valid(side);
}

//--------------------------------------------------------------------------------

void
batch_create_sides(stk::mesh::BulkData & mesh, const std::vector<SideDescription> & sideDescriptions)
{
  if (!mesh.has_face_adjacent_element_graph())
  {
    mesh.initialize_face_adjacent_element_graph();
  }

  mesh.modification_begin();
  for (auto && sideDescription : sideDescriptions)
  {
    if (mesh.is_valid(sideDescription.element))
    {
      ThrowAssertMsg(mesh.bucket(sideDescription.element).owned(), "Expecting owned entity");
      if (!does_element_side_exist(mesh, sideDescription.element, sideDescription.elementSideOrdinal))
          mesh.declare_element_side(sideDescription.element, sideDescription.elementSideOrdinal, sideDescription.sideParts);
    }
  }
  mesh.modification_end();

  ParallelThrowAssert(mesh.parallel(), check_face_and_edge_ownership(mesh));
  ParallelThrowAssert(mesh.parallel(), check_face_and_edge_relations(mesh));
}

//--------------------------------------------------------------------------------

void
make_side_ids_consistent_with_stk_convention(stk::mesh::BulkData & mesh)
{
  const stk::mesh::MetaData & meta = mesh.mesh_meta_data();
  const stk::mesh::EntityRank side_rank = meta.side_rank();
  stk::mesh::Selector not_ghost_selector(meta.locally_owned_part() | meta.globally_shared_part());

  std::vector<stk::mesh::Entity> sides;
  stk::mesh::get_selected_entities( not_ghost_selector, mesh.buckets(side_rank), sides );

  std::vector<SideDescription> side_requests;

  mesh.modification_begin();
  for (auto&& side : sides)
  {
    const unsigned num_side_elems = mesh.num_elements(side);
    const stk::mesh::Entity* side_elems = mesh.begin_elements(side);
    const stk::mesh::ConnectivityOrdinal* side_elem_ordinals = mesh.begin_element_ordinals(side);
    ThrowRequire(num_side_elems > 0);
    stk::mesh::EntityId newId = 0;
    for (unsigned i=0; i<num_side_elems; ++i)
    {
      const stk::mesh::EntityId stk_elem_side_id = 10*mesh.identifier(side_elems[i]) + side_elem_ordinals[i] + 1;
      if (i==0 || stk_elem_side_id < newId) newId = stk_elem_side_id;
    }

    if (newId != mesh.identifier(side))
    {
      // Can't use because of locality assertions: mesh.change_entity_id(newId, side);
      const stk::mesh::PartVector side_parts = get_removable_parts(mesh, side);
      side_requests.emplace_back(side_elems[0], side_elem_ordinals[0], side_parts);
      disconnect_and_destroy_entity(mesh, side);
    }
  }
  mesh.modification_end();

  batch_create_sides(mesh, side_requests);
}

//--------------------------------------------------------------------------------

double
compute_element_volume_to_edge_ratio(stk::mesh::BulkData & mesh, stk::mesh::Entity element, const stk::mesh::Field<double> * const coordsField)
{
  stk::topology elem_topology = mesh.bucket(element).topology();
  if (elem_topology == stk::topology::TETRAHEDRON_4)
  {
    const std::array<krino::Vector3d,4> nodes = gather_tet_coordinates(mesh, element, coordsField);
    const double vol = compute_tet_volume(nodes);
    compute_tri_or_tet_volume(mesh, element, *coordsField);
    const double edge_rms = std::sqrt(
        ((nodes[1]-nodes[0]).length_squared() +
         (nodes[2]-nodes[0]).length_squared() +
         (nodes[3]-nodes[0]).length_squared() +
         (nodes[3]-nodes[1]).length_squared() +
         (nodes[3]-nodes[2]).length_squared() +
         (nodes[2]-nodes[1]).length_squared())/6.);
    return vol/(edge_rms*edge_rms*edge_rms);
  }
  ThrowRuntimeError("Topology " << elem_topology << " not supported in compute_element_volume_to_edge_ratio.");
}

//--------------------------------------------------------------------------------

template<typename T>
void debug_write_entity_field_if_present(std::ostream & output, const stk::mesh::BulkData & mesh, const FieldRef field, const stk::mesh::Entity entity)
{
  field.field().sync_to_host();
  if (field.entity_rank() == mesh.entity_rank(entity) && field.type_is<T>())
  {
    const auto * data = field_data<T>(field, entity);
    if (nullptr != data)
    {
      const unsigned field_length = field.length();
      if (1 == field_length)
      {
        output << "  Field: " << field.name() << ", state=" << static_cast<int>(field.state()) << ", value=" << *data << "\n";
      }
      else
      {
        output << "  Field: " << field.name() << ", state=" << static_cast<int>(field.state()) << ", values[] = ";
        for (unsigned i = 0; i < field_length; ++i) output << data[i] << " ";
        output << "\n";
      }
    }
  }
}

static void
debug_entity(std::ostream & output, const stk::mesh::BulkData & mesh, stk::mesh::Entity entity, const bool includeFields)
{
  if (!mesh.is_valid(entity))
  {
    output << "Invalid entity: " << mesh.entity_key(entity) << std::endl;
    return;
  }
  output << mesh.entity_key(entity) << ", parallel owner = " << mesh.parallel_owner_rank(entity) << " \n";
  output << "  Connectivity:" << "\n";
  const stk::mesh::EntityRank end_rank = static_cast<stk::mesh::EntityRank>(mesh.mesh_meta_data().entity_rank_count());
  for (stk::mesh::EntityRank r = stk::topology::BEGIN_RANK; r < end_rank; ++r) {
    unsigned num_rels = mesh.num_connectivity(entity, r);
    stk::mesh::Entity const *rel_entities = mesh.begin(entity, r);
    stk::mesh::ConnectivityOrdinal const *rel_ordinals = mesh.begin_ordinals(entity, r);
    stk::mesh::Permutation const *rel_permutations = mesh.begin_permutations(entity, r);
    for (unsigned i = 0; i < num_rels; ++i) {
      output << " " << mesh.entity_key(rel_entities[i])
         << " @" << rel_ordinals[i];
      if (rel_permutations) output << ":" << (int)rel_permutations[i];
      output << "\n";
    }
  }
  output << "  Parts: ";
  const stk::mesh::PartVector & parts = mesh.bucket(entity).supersets();
  for(stk::mesh::PartVector::const_iterator part_iter = parts.begin(); part_iter != parts.end(); ++part_iter)
  {
    const stk::mesh::Part * const part = *part_iter;
    output << part->name() << " ";
  }
  output << "\n";

  if (includeFields)
  {
    for ( auto && fieldPtr : mesh.mesh_meta_data().get_fields())
    {
      const FieldRef field = *fieldPtr;

      debug_write_entity_field_if_present<double>(output, mesh, field, entity);
      debug_write_entity_field_if_present<int>(output, mesh, field, entity);
      debug_write_entity_field_if_present<uint64_t>(output, mesh, field, entity);
    }
    output << "\n";
  }
}

std::string
debug_entity(const stk::mesh::BulkData & mesh, stk::mesh::Entity entity, const bool includeFields)
{
  std::ostringstream out;
  debug_entity(out, mesh, entity, includeFields);
  return out.str();
}

std::string
debug_entity(const stk::mesh::BulkData & mesh, stk::mesh::Entity entity)
{
  return debug_entity(mesh, entity, false);
}

static void
debug_entity_1line(std::ostream & output, const stk::mesh::BulkData & mesh, stk::mesh::Entity entity, const bool omitSideRank)
{
  if (!mesh.is_valid(entity))
  {
    output << "Invalid entity: " << mesh.entity_key(entity) << "\n";
    return;
  }
  output << mesh.entity_key(entity);
  output << ",  Owner: " << mesh.parallel_owner_rank(entity);
  output << ",  Connectivity: ";
  const stk::mesh::EntityRank end_rank = static_cast<stk::mesh::EntityRank>(mesh.mesh_meta_data().entity_rank_count());
  for (stk::mesh::EntityRank r = stk::topology::BEGIN_RANK; r < end_rank; ++r)
  {
    if (omitSideRank && r == mesh.mesh_meta_data().side_rank())
      continue;
    unsigned num_rels = mesh.num_connectivity(entity, r);
    stk::mesh::Entity const *rel_entities = mesh.begin(entity, r);
    stk::mesh::ConnectivityOrdinal const *rel_ordinals = mesh.begin_ordinals(entity, r);
    for (unsigned i = 0; i < num_rels; ++i) {
      output << " " << mesh.entity_key(rel_entities[i])
         << " @" << rel_ordinals[i] << " ";
    }
  }
  output << ",  Parts: ";
  for(auto&& part : mesh.bucket(entity).supersets())
  {
    output << part->name() << " ";
  }
}

std::string
debug_entity_1line(const stk::mesh::BulkData & mesh, stk::mesh::Entity entity, const bool omitSideRank)
{
  std::ostringstream out;
  debug_entity_1line(out, mesh, entity, omitSideRank);
  return out.str();
}

//--------------------------------------------------------------------------------

std::vector<unsigned>
get_side_permutation(stk::topology topology, stk::mesh::Permutation node_permutation)
{
  const unsigned perm = node_permutation;
  switch (topology)
  {
    case stk::topology::TRIANGLE_3_2D:
    case stk::topology::TRIANGLE_6_2D:
      switch (perm)
      {
        case 0:  return {0, 1, 2};
        case 1:  return {2, 0, 1};
        case 2:  return {1, 2, 0};
        default: ThrowRuntimeError("find_side_permutation error, invalid triangle permutation.");
      }
      break;
    case stk::topology::TETRAHEDRON_4:
    case stk::topology::TETRAHEDRON_10:
      switch (perm)
      {
        case  0:  return {0, 1, 2, 3};
        case  1:  return {1, 2, 0, 3};
        case  2:  return {2, 0, 1, 3};
        case  3:  return {2, 1, 3, 0};
        case  4:  return {1, 3, 2, 0};
        case  5:  return {3, 2, 1, 0};
        case  6:  return {3, 1, 0, 2};
        case  7:  return {1, 0, 3, 2};
        case  8:  return {0, 3, 1, 2};
        case  9:  return {0, 2, 3, 1};
        case 10:  return {2, 3, 0, 1};
        case 11:  return {3, 0, 2, 1};
        default: ThrowRuntimeError("find_side_permutation error, invalid tetrahedron permutation.");
      }
      break;
    default: ThrowRuntimeError("find_side_permutation error, unsupported topology.");
  }
}

//--------------------------------------------------------------------------------

const stk::mesh::Part &
find_element_part(const stk::mesh::BulkData& mesh, stk::mesh::Entity elem)
{
  ThrowAssert(mesh.entity_rank(elem) == stk::topology::ELEMENT_RANK);
  const stk::mesh::Part * elem_io_part = nullptr;

  const stk::mesh::PartVector & elem_parts = mesh.bucket(elem).supersets();
  for(stk::mesh::PartVector::const_iterator part_iter = elem_parts.begin(); part_iter != elem_parts.end(); ++part_iter)
  {
    const stk::mesh::Part * const part = *part_iter;
    if (part->primary_entity_rank() == stk::topology::ELEMENT_RANK && part->subsets().empty() && part->topology() != stk::topology::INVALID_TOPOLOGY)
    {
      // there should only be one element rank part without subsets with topology on the element
      ThrowRequireMsg(nullptr == elem_io_part, "For element " << mesh.identifier(elem) << ", more than one element rank part was found: " << elem_io_part->name() << " " << part->name());
      elem_io_part = part;
    }
  }
  ThrowRequire(NULL != elem_io_part);

  return *elem_io_part;
}

//--------------------------------------------------------------------------------

void
disconnect_entity(stk::mesh::BulkData & mesh, stk::mesh::Entity entity)
{
  stk::mesh::EntityRank entity_rank = mesh.entity_rank(entity);
  std::vector<stk::mesh::Entity> relatives;
  std::vector<stk::mesh::ConnectivityOrdinal> relative_ordinals;

  const stk::mesh::EntityRank highest_entity_rank = static_cast<stk::mesh::EntityRank>(mesh.mesh_meta_data().entity_rank_count()-1);
  for (stk::mesh::EntityRank irank = highest_entity_rank; irank != entity_rank; --irank)
  {
    // Previously this attempted to delete forward or backward and still the list got corrupted,
    // so just copy into vector and delete from there.
    relatives.assign(mesh.begin(entity, irank),mesh.end(entity, irank));
    relative_ordinals.assign(mesh.begin_ordinals(entity, irank), mesh.end_ordinals(entity, irank));

    for (size_t irel = 0; irel < relatives.size(); ++irel)
    {
      mesh.destroy_relation( relatives[irel], entity, relative_ordinals[irel]);
    }
  }
}

//--------------------------------------------------------------------------------

bool
disconnect_and_destroy_entity(stk::mesh::BulkData & mesh, stk::mesh::Entity entity)
{
  disconnect_entity(mesh, entity);
  return mesh.destroy_entity(entity);
}

//--------------------------------------------------------------------------------

bool
check_induced_parts(const stk::mesh::BulkData & mesh)
{

  // This method requires aura to work correctly.
  if (!mesh.is_automatic_aura_on() && mesh.parallel_size() > 1)
  {
    // Skip check if we don't have aura
    return true;
  }

  bool success = true;

  const stk::mesh::MetaData & meta = mesh.mesh_meta_data();
  stk::mesh::Selector not_ghost_selector(meta.locally_owned_part() | meta.globally_shared_part());

  std::vector< stk::mesh::Entity> entities;

  for (stk::mesh::EntityRank entity_rank = stk::topology::NODE_RANK; entity_rank <= stk::topology::ELEMENT_RANK; ++entity_rank)
  {
    stk::mesh::get_selected_entities( not_ghost_selector, mesh.buckets(entity_rank), entities );

    for (unsigned i=0; i<entities.size(); ++i)
    {
      stk::mesh::Entity entity = entities[i];
      const stk::mesh::PartVector & entity_parts = mesh.bucket(entity).supersets();

      // Make sure all entity_rank parts parts are found on the lower rank entities
      for (stk::mesh::EntityRank relative_rank = stk::topology::NODE_RANK; relative_rank < entity_rank; ++relative_rank)
      {
        for(stk::mesh::PartVector::const_iterator part_iter = entity_parts.begin(); part_iter != entity_parts.end(); ++part_iter)
        {
          const stk::mesh::Part * const entity_part = *part_iter;
          if (entity_part->primary_entity_rank() == entity_rank && !entity_part->force_no_induce())
          {
            bool have_relative_missing_part = false;
            const unsigned num_relatives = mesh.num_connectivity(entity, relative_rank);
            const stk::mesh::Entity* relatives = mesh.begin(entity, relative_rank);
            for (unsigned it_rel=0; it_rel<num_relatives && !have_relative_missing_part; ++it_rel)
            {
              stk::mesh::Entity relative = relatives[it_rel];
              if (!mesh.is_valid(relative)) continue;
              const stk::mesh::PartVector & relative_parts = mesh.bucket(relative).supersets();
              if (std::find(relative_parts.begin(), relative_parts.end(), entity_part) == relative_parts.end())
              {
                have_relative_missing_part = true;
              }
            }

            if (have_relative_missing_part)
            {
              success = false;
              krinolog << "Missing induction. Part " << entity_part->name() << " is found on " << mesh.entity_key(entity) << " but is missing from relatives: ";
              for (unsigned it_rel=0; it_rel<num_relatives; ++it_rel)
              {
                stk::mesh::Entity relative = relatives[it_rel];
                if (!mesh.is_valid(relative)) continue;
                const stk::mesh::PartVector & relative_parts = mesh.bucket(relative).supersets();
                if (std::find(relative_parts.begin(), relative_parts.end(), entity_part) == relative_parts.end())
                {
                  krinolog << mesh.entity_key(relative) << " ";
                }
              }
              krinolog << stk::diag::dendl;
            }
          }
        }
      }

      // Make sure higher ranked parts are found on an actual higher ranked entity
      for (stk::mesh::EntityRank relative_rank = stk::mesh::EntityRank(entity_rank+1); relative_rank <= stk::topology::ELEMENT_RANK; ++relative_rank)
      {
        for(stk::mesh::PartVector::const_iterator part_iter = entity_parts.begin(); part_iter != entity_parts.end(); ++part_iter)
        {
          const stk::mesh::Part * const entity_part = *part_iter;
          if (entity_part->primary_entity_rank() == relative_rank)
          {
            bool found_relative_with_part = false;
            const unsigned num_relatives = mesh.num_connectivity(entity, relative_rank);
            const stk::mesh::Entity* relatives = mesh.begin(entity, relative_rank);
            for (unsigned it_rel=0; it_rel<num_relatives && !found_relative_with_part; ++it_rel)
            {
              stk::mesh::Entity relative = relatives[it_rel];
              if (!mesh.is_valid(relative)) continue;
              const stk::mesh::PartVector & relative_parts = mesh.bucket(relative).supersets();
              if (std::find(relative_parts.begin(), relative_parts.end(), entity_part) != relative_parts.end())
              {
                found_relative_with_part = true;
              }
            }

            if (!found_relative_with_part)
            {
              success = false;
              krinolog << "Extraneous induction. Part " << entity_part->name() << " is found on " << mesh.entity_key(entity) << " but not on any of its relatives: ";
              for (unsigned it_rel=0; it_rel<num_relatives; ++it_rel)
              {
                stk::mesh::Entity relative = relatives[it_rel];
                if (!mesh.is_valid(relative)) continue;
                krinolog << mesh.entity_key(relative) << " ";
              }
              krinolog << stk::diag::dendl;
            }
          }
        }
      }
    }
  }

  return success;
}

//--------------------------------------------------------------------------------

bool
check_shared_entity_nodes(const stk::mesh::BulkData & mesh, stk::mesh::EntityKey remote_entity_key, std::vector<stk::mesh::EntityId> & remote_entity_node_ids)
{
  stk::mesh::Entity entity = mesh.get_entity(remote_entity_key);
  if (!mesh.is_valid(entity))
  {
    krinolog << "Shared entity error, local entity does not exist, remote entity: " << remote_entity_key << stk::diag::dendl;
    return false;
  }

  std::vector<stk::mesh::Entity> entity_nodes(mesh.begin_nodes(entity), mesh.end_nodes(entity));
  if (entity_nodes.size() != remote_entity_node_ids.size())
  {
    krinolog << "Shared entity error, number_of nodes don't match, number of remote nodes = " << remote_entity_node_ids.size() << stk::diag::dendl;
    krinolog << "Local entity: " << debug_entity(mesh, entity) << stk::diag::dendl;
    return false;
  }

  bool nodes_match = true;
  for (size_t node_index=0;node_index<entity_nodes.size();++node_index)
  {
    stk::mesh::Entity remote_entity_node = mesh.get_entity(stk::topology::NODE_RANK, remote_entity_node_ids[node_index]);
    if (remote_entity_node != entity_nodes[node_index]) nodes_match = false;
  }
  if (!nodes_match)
  {
    krinolog << "Shared entity error, node mismatch: " << stk::diag::dendl;
    krinolog << "Remote entity nodes: ";
    for (size_t node_index=0;node_index<entity_nodes.size();++node_index)
    {
      krinolog << remote_entity_node_ids[node_index] << " ";
    }
    krinolog << stk::diag::dendl;
    krinolog << "Local entity: " << debug_entity(mesh, entity) << stk::diag::dendl;
    return false;
  }

  return true;
}

bool
check_shared_entity_nodes(const stk::mesh::BulkData & mesh, std::vector<stk::mesh::Entity> & entities)
{
  bool success = true;

  std::vector<int> sharing_procs;
  stk::CommSparse comm_spec(mesh.parallel());

  for (int phase=0;phase<2;++phase)
  {
    for (std::vector<stk::mesh::Entity>::iterator it_entity = entities.begin(); it_entity != entities.end(); ++it_entity)
    {
      stk::mesh::Entity entity = *it_entity;
      std::vector<stk::mesh::Entity> entity_nodes(mesh.begin_nodes(entity), mesh.end_nodes(entity));
      stk::mesh::EntityKey entity_key = mesh.entity_key(entity);
      ThrowRequire(mesh.bucket(entity).shared());

      mesh.shared_procs_intersection(entity_nodes, sharing_procs);
      for (size_t proc_index=0;proc_index<sharing_procs.size();++proc_index)
      {
        const int other_proc = sharing_procs[proc_index];
        if (other_proc != mesh.parallel_rank())
        {
          comm_spec.send_buffer(other_proc).pack(entity_key);
          const size_t num_nodes = entity_nodes.size();
          comm_spec.send_buffer(other_proc).pack(num_nodes);
          for (size_t node_index=0; node_index<entity_nodes.size(); ++node_index )
          {
            comm_spec.send_buffer(other_proc).pack(mesh.identifier(entity_nodes[node_index]));
          }
        }
      }
    }

    if ( phase == 0 )
    {
      comm_spec.allocate_buffers();
    }
    else
    {
      comm_spec.communicate();
    }
  }

  for(int i = 0; i < mesh.parallel_size(); ++i)
  {
    if(i != mesh.parallel_rank())
    {
      while(comm_spec.recv_buffer(i).remaining())
      {
        stk::mesh::EntityKey entity_key;
        comm_spec.recv_buffer(i).unpack(entity_key);
        size_t num_nodes;
        comm_spec.recv_buffer(i).unpack(num_nodes);
        std::vector<stk::mesh::EntityId> entity_node_ids(num_nodes);
        for (size_t node_index=0; node_index<num_nodes; ++node_index )
        {
          comm_spec.recv_buffer(i).unpack(entity_node_ids[node_index]);
        }

        if (!check_shared_entity_nodes(mesh, entity_key, entity_node_ids)) success = false;
      }
    }
  }
  return success;
}


bool
check_shared_entity_nodes(const stk::mesh::BulkData & mesh)
{
  const stk::mesh::MetaData & meta = mesh.mesh_meta_data();
  std::vector<int> sharing_procs;
  stk::CommSparse comm_spec(mesh.parallel());

  bool success = true;
  stk::mesh::Selector shared_selector = meta.globally_shared_part();
  std::vector<stk::mesh::Entity> entities;

  for (stk::mesh::EntityRank entity_rank = stk::topology::EDGE_RANK; entity_rank < stk::topology::ELEMENT_RANK; ++entity_rank)
  {
    stk::mesh::get_selected_entities( shared_selector, mesh.buckets( entity_rank ), entities );

    if (!check_shared_entity_nodes(mesh, entities)) success = false;
  }
  return success;
}

//--------------------------------------------------------------------------------

bool
check_face_and_edge_relations(const stk::mesh::BulkData & mesh)
{
  const stk::mesh::MetaData & meta = mesh.mesh_meta_data();

  bool success = true;
  stk::mesh::Selector not_ghost_selector = meta.locally_owned_part() | meta.globally_shared_part();

  for (stk::mesh::EntityRank entity_rank = stk::topology::EDGE_RANK; entity_rank <= stk::topology::FACE_RANK; ++entity_rank)
  {
    const stk::mesh::BucketVector & buckets = mesh.get_buckets( entity_rank, not_ghost_selector );

    stk::mesh::BucketVector::const_iterator ib = buckets.begin();
    stk::mesh::BucketVector::const_iterator ib_end = buckets.end();

    for ( ; ib != ib_end; ++ib )
    {
      const stk::mesh::Bucket & b = **ib;
      const size_t length = b.size();
      for (size_t it_entity = 0; it_entity < length; ++it_entity)
      {
        stk::mesh::Entity entity = b[it_entity];
        stk::topology entity_topology = mesh.bucket(entity).topology();
        std::vector<stk::mesh::Entity> entity_nodes(mesh.begin_nodes(entity), mesh.end_nodes(entity));
        std::vector<stk::mesh::Entity> entity_elems;
        stk::mesh::get_entities_through_relations(mesh, entity_nodes, stk::topology::ELEMENT_RANK, entity_elems);

        if (entity_elems.empty())
        {
          krinolog << "Relation error, entity not attached to any elements: " << stk::diag::dendl;
          krinolog << "Entity: " << debug_entity(mesh, entity) << stk::diag::dendl;
          success = false;
        }

        bool have_coincident_shell = false;
        std::pair<stk::mesh::ConnectivityOrdinal, stk::mesh::Permutation> shell_relationship(stk::mesh::INVALID_CONNECTIVITY_ORDINAL, stk::mesh::INVALID_PERMUTATION);
        for (auto&& elem : entity_elems)
        {
          stk::topology elem_topology = mesh.bucket(elem).topology();
          if (elem_topology.is_shell() && elem_topology.num_nodes() == entity_topology.num_nodes())
          {
            have_coincident_shell = true;
            shell_relationship = determine_shell_side_ordinal_and_permutation(mesh, elem, entity);
            break;
          }
        }

        for (auto&& elem : entity_elems)
        {
          stk::topology elem_topology = mesh.bucket(elem).topology();
          const bool is_coincident_shell = (elem_topology.is_shell() && elem_topology.num_nodes() == entity_topology.num_nodes());
          bool should_be_attached = true;
          if (have_coincident_shell && !is_coincident_shell)
          {
            // Volume elements should only be attached to inward pointing faces when the surface has a shell.
            std::pair<stk::mesh::ConnectivityOrdinal, stk::mesh::Permutation> relationship = determine_ordinal_and_permutation(mesh, elem, entity);
            const bool elem_polarity = entity_topology.is_positive_polarity(relationship.second);
            const bool shell_polarity = entity_topology.is_positive_polarity(shell_relationship.second);
            should_be_attached = elem_polarity != shell_polarity;
          }

          const unsigned num_elem_entities = mesh.num_connectivity(elem, entity_rank);
          const stk::mesh::Entity* elem_entities = mesh.begin(elem, entity_rank);
          const stk::mesh::ConnectivityOrdinal * elem_ordinals = mesh.begin_ordinals(elem, entity_rank);
          const stk::mesh::Permutation * elem_permutations = mesh.begin_permutations(elem, entity_rank);
          bool already_attached = false;
          for (unsigned it_s=0; it_s<num_elem_entities; ++it_s)
          {
            if (elem_entities[it_s] == entity)
            {
              if (!should_be_attached)
              {
                krinolog << "Relation error, entity is attached to element, but should not be (due to coincident shell on side): " << stk::diag::dendl;
                krinolog << "Entity: " << debug_entity(mesh, entity) << stk::diag::dendl;
                krinolog << "Element: " << debug_entity(mesh, elem) << stk::diag::dendl;
                success = false;
              }
              if (already_attached)
              {
                krinolog << "Relation error, entity attached to element more than once: " << stk::diag::dendl;
                krinolog << "Entity: " << debug_entity(mesh, entity) << stk::diag::dendl;
                krinolog << "Element: " << debug_entity(mesh, elem) << stk::diag::dendl;
                success = false;
              }
              else
              {
                already_attached = true;
                std::pair<stk::mesh::ConnectivityOrdinal, stk::mesh::Permutation> relationship =
                    is_coincident_shell ?
                    shell_relationship :
                    determine_ordinal_and_permutation(mesh, elem, entity);
                if (relationship.first != elem_ordinals[it_s])
                {
                  krinolog << "Relation error, ordinal is incorrect: " << relationship.first << "!=" << elem_ordinals[it_s] << stk::diag::dendl;
                  krinolog << "Entity: " << debug_entity(mesh, entity) << stk::diag::dendl;
                  krinolog << "Element: " << debug_entity(mesh, elem) << stk::diag::dendl;
                  success = false;
                }
                if (relationship.second != elem_permutations[it_s])
                {
                  krinolog << "Relation error, permutation is incorrect: " << relationship.second << "!=" << elem_permutations[it_s] << stk::diag::dendl;
                  krinolog << "Entity: " << debug_entity(mesh, entity) << stk::diag::dendl;
                  krinolog << "Element: " << debug_entity(mesh, elem) << stk::diag::dendl;
                  success = false;
                }
              }
            }
          }
          if (!already_attached && should_be_attached)
          {

            krinolog << "Relation error, entity is not attached to element: " << stk::diag::dendl;
            krinolog << "Entity: " << debug_entity(mesh, entity) << stk::diag::dendl;
            krinolog << "Element: " << debug_entity(mesh, elem) << stk::diag::dendl;
            std::pair<stk::mesh::ConnectivityOrdinal, stk::mesh::Permutation> relationship = determine_ordinal_and_permutation(mesh, elem, entity);
            for (unsigned it_s=0; it_s<num_elem_entities; ++it_s)
            {
              if (elem_ordinals[it_s] == relationship.first)
              {
                krinolog << "Another side is already attached to this element with ordinal " << relationship.first << ": " << debug_entity(mesh, elem_entities[it_s]) << stk::diag::dendl;
              }
            }
            success = false;
          }
        }
      }
    }
  }
  return success;
}

void
attach_sides_to_elements(stk::mesh::BulkData & mesh)
{
  stk::mesh::MetaData & meta = mesh.mesh_meta_data();

  std::vector<stk::mesh::Entity> entities;
  stk::mesh::get_entities( mesh, meta.side_rank(), entities );

  mesh.modification_begin();
  for (auto&& entity : entities)
  {
    attach_entity_to_elements(mesh, entity);
  }
  mesh.modification_end();
}

static bool is_entity_attached_to_element(const stk::mesh::BulkData & mesh, const stk::mesh::EntityRank entityRank, const stk::mesh::Entity entity, const stk::mesh::Entity element)
{
  for (auto && elemEntity : StkMeshEntities{mesh.begin(element, entityRank), mesh.end(element, entityRank)})
    if (elemEntity == entity)
      return true;
  return false;
}

void
attach_entity_to_element(stk::mesh::BulkData & mesh, const stk::mesh::EntityRank entityRank, const stk::mesh::Entity entity, const stk::mesh::Entity element)
{
  //Sorry! Passing these scratch vectors into stk's declare_relation function is
  //a performance improvement (fewer allocations). But stk will try to clean up
  //this ugliness soon. (i.e., find a better way to get the performance.)
  stk::mesh::OrdinalVector scratch1, scratch2, scratch3;

  //const auto & [ordinal, permutation] = determine_ordinal_and_permutation(mesh, element, entity);
  const auto & ordinalAndPermutations = determine_ordinal_and_permutation(mesh, element, entity);
  const auto & ordinal = ordinalAndPermutations.first;
  const auto & permutation = ordinalAndPermutations.second;

  mesh.declare_relation( element, entity, ordinal, permutation, scratch1, scratch2, scratch3 );
  ThrowRequireMsg(is_entity_attached_to_element(mesh, entityRank, entity, element),  "Could not attach " << debug_entity_1line(mesh,entity) << "  to element " << debug_entity_1line(mesh,element));
}

void
attach_entity_to_elements(stk::mesh::BulkData & mesh, stk::mesh::Entity entity)
{
  //Sorry! Passing these scratch vectors into stk's declare_relation function is
  //a performance improvement (fewer allocations). But stk will try to clean up
  //this ugliness soon. (i.e., find a better way to get the performance.)
  stk::mesh::OrdinalVector scratch1, scratch2, scratch3;

  stk::mesh::MetaData & meta = mesh.mesh_meta_data();
  stk::topology entity_topology = mesh.bucket(entity).topology();
  stk::mesh::EntityRank entity_rank = entity_topology.rank();
  std::vector<stk::mesh::Entity> entity_nodes(mesh.begin_nodes(entity), mesh.end_nodes(entity));
  std::vector<stk::mesh::Entity> entity_elems;
  stk::mesh::get_entities_through_relations(mesh, entity_nodes, stk::topology::ELEMENT_RANK, entity_elems);

  bool have_coincident_shell = false;
  std::pair<stk::mesh::ConnectivityOrdinal, stk::mesh::Permutation> shell_relationship(stk::mesh::INVALID_CONNECTIVITY_ORDINAL, stk::mesh::INVALID_PERMUTATION);
  for (auto&& elem : entity_elems)
  {
    stk::topology elem_topology = mesh.bucket(elem).topology();
    if (elem_topology.is_shell() && elem_topology.num_nodes() == entity_topology.num_nodes())
    {
      have_coincident_shell = true;
      shell_relationship = determine_shell_side_ordinal_and_permutation(mesh, elem, entity);
      break;
    }
  }

  for (auto&& elem : entity_elems)
  {
    if (!mesh.bucket(elem).member(meta.locally_owned_part()))
    {
      continue;
    }
    const bool alreadyAttached = is_entity_attached_to_element(mesh, entity_rank, entity, elem);
    if (!alreadyAttached)
    {
      std::pair<stk::mesh::ConnectivityOrdinal, stk::mesh::Permutation> relationship(stk::mesh::INVALID_CONNECTIVITY_ORDINAL, stk::mesh::INVALID_PERMUTATION);
      if (!have_coincident_shell)
      {
        relationship = determine_ordinal_and_permutation(mesh, elem, entity);
      }
      else
      {
        stk::topology elem_topology = mesh.bucket(elem).topology();
        if (elem_topology.is_shell() && elem_topology.num_nodes() == entity_topology.num_nodes())
        {
          ThrowAssertMsg(shell_relationship.second == determine_permutation(mesh, elem, entity, shell_relationship.first), "All shells should have same permutation for side.");
          relationship = shell_relationship;
        }
        else
        {
          relationship = determine_ordinal_and_permutation(mesh, elem, entity);
          const bool elem_polarity = entity_topology.is_positive_polarity(relationship.second);
          const bool shell_polarity = entity_topology.is_positive_polarity(shell_relationship.second);
          if (elem_polarity == shell_polarity)
          {
            // Side does not touch volume element;
            continue;
          }
        }
      }

      mesh.declare_relation( elem, entity, relationship.first, relationship.second, scratch1, scratch2, scratch3 );
      const bool successfullyAttached = (find_entity_by_ordinal(mesh, elem, entity_rank, relationship.first) == entity);
      if (!successfullyAttached)
      {
        krinolog << "Could not attach " << debug_entity_1line(mesh,entity) << "\n";
        krinolog << "  to element " << debug_entity_1line(mesh,elem) << stk::diag::dendl;
      }
      else if (krinolog.shouldPrint(LOG_DEBUG))
      {
        krinolog << "Attached " << debug_entity_1line(mesh,entity) << "\n";
        krinolog << "  to element " << debug_entity_1line(mesh,elem) << stk::diag::dendl;
      }
      ThrowRequire(successfullyAttached);
    }
  }
}

void unpack_entities_from_other_procs(const stk::mesh::BulkData & mesh,
    std::set<stk::mesh::Entity> & entities,
    stk::CommSparse &commSparse)
{
  entities.clear();
  stk::unpack_communications(commSparse, [&](int procId)
  {
    stk::CommBuffer & buffer = commSparse.recv_buffer(procId);

    while ( buffer.remaining() )
    {
      stk::mesh::EntityKey entityKey;
      commSparse.recv_buffer(procId).unpack(entityKey);
      entities.insert(mesh.get_entity(entityKey));
    }
  });
}

void
update_node_activation(stk::mesh::BulkData & mesh, stk::mesh::Part & active_part)
{
  stk::mesh::MetaData & meta = mesh.mesh_meta_data();

  stk::mesh::PartVector active_part_vec(1, &active_part);
  stk::mesh::PartVector inactive_part_vec;

  std::vector<stk::mesh::Entity> entities;
  stk::mesh::Selector locally_owned(meta.locally_owned_part());
  stk::mesh::get_selected_entities( locally_owned, mesh.buckets( stk::topology::NODE_RANK ), entities );

  for (std::vector<stk::mesh::Entity>::iterator i_node = entities.begin(); i_node != entities.end(); ++i_node)
  {
    stk::mesh::Entity node = *i_node;

    const unsigned num_node_elems = mesh.num_elements(node);
    const stk::mesh::Entity* node_elems = mesh.begin_elements(node);
    bool have_active_elems = false;
    for (unsigned node_elem_index=0; node_elem_index<num_node_elems && !have_active_elems; ++node_elem_index)
    {
      stk::mesh::Entity elem = node_elems[node_elem_index];

      if (mesh.bucket(elem).member(active_part))
      {
        have_active_elems = true;
      }
    }

    if (have_active_elems)
    {
      if (!mesh.bucket(node).member(active_part))
      {
        mesh.change_entity_parts(node, active_part_vec, inactive_part_vec);
      }
    }
    else
    {
      if (mesh.bucket(node).member(active_part))
      {
        mesh.change_entity_parts(node, inactive_part_vec, active_part_vec);
      }
    }
  }
}

void
activate_all_entities(stk::mesh::BulkData & mesh, stk::mesh::Part & active_part)
{
  std::vector<stk::mesh::PartVector> add_parts;
  std::vector<stk::mesh::PartVector> remove_parts;
  std::vector<stk::mesh::Entity> entities;

  stk::mesh::Selector inactive_locally_owned = mesh.mesh_meta_data().locally_owned_part() & !active_part;

  for (stk::mesh::EntityRank entity_rank = stk::topology::NODE_RANK; entity_rank <= stk::topology::ELEMENT_RANK; ++entity_rank)
  {
    const stk::mesh::BucketVector & buckets = mesh.get_buckets(entity_rank, inactive_locally_owned);
    for (auto&& bucket_ptr : buckets)
    {
      entities.insert(entities.end(), bucket_ptr->begin(), bucket_ptr->end());
    }
  }
  add_parts.assign(entities.size(), {&active_part});
  remove_parts.resize(entities.size());

  mesh.batch_change_entity_parts(entities, add_parts, remove_parts);
}

//--------------------------------------------------------------------------------

void destroy_custom_ghostings(stk::mesh::BulkData & mesh)
{
  const std::vector<stk::mesh::Ghosting *> & ghostings = mesh.ghostings();
  for(unsigned i = stk::mesh::BulkData::AURA+1; i < ghostings.size(); ++i)
  {
    mesh.destroy_ghosting(*ghostings[i]);
  }
}

//--------------------------------------------------------------------------------

void
delete_mesh_entities(stk::mesh::BulkData & mesh, const std::vector<stk::mesh::Entity> & child_elems)
{
  if (child_elems.empty())
    return;

  stk::mesh::MetaData & meta = mesh.mesh_meta_data();

  stk::mesh::Selector universal_selector = meta.universal_part();

  std::vector<stk::mesh::Entity> child_sides;
  std::vector<stk::mesh::Entity> child_edges;
  std::vector<stk::mesh::Entity> child_nodes;

  for (unsigned i=0; i<child_elems.size(); ++i)
  {
    stk::mesh::Entity child = child_elems[i];
    if (!mesh.is_valid(child)) continue;

    const unsigned num_child_sides = mesh.num_connectivity(child,meta.side_rank());
    const stk::mesh::Entity* this_child_sides = mesh.begin(child,meta.side_rank());
    for (unsigned child_side_index=0; child_side_index<num_child_sides; ++child_side_index)
    {
      stk::mesh::Entity side = this_child_sides[child_side_index];
      child_sides.push_back(side);
    }

    if (meta.spatial_dimension() == 3)
    {
      const unsigned num_child_edges = mesh.num_edges(child);
      const stk::mesh::Entity* this_child_edges = mesh.begin_edges(child);
      for (unsigned child_edge_index=0; child_edge_index<num_child_edges; ++child_edge_index)
      {
        stk::mesh::Entity edge = this_child_edges[child_edge_index];
        child_edges.push_back(edge);
      }
    }

    const unsigned num_child_nodes = mesh.num_nodes(child);
    const stk::mesh::Entity* this_child_nodes = mesh.begin_nodes(child);
    for (unsigned child_node_index=0; child_node_index<num_child_nodes; ++child_node_index)
    {
      stk::mesh::Entity node = this_child_nodes[child_node_index];
      child_nodes.push_back(node);
    }

    if(krinolog.shouldPrint(LOG_DEBUG))
      krinolog << "Destroying child, elem#" << mesh.identifier(child) << ":" << stk::diag::dendl;

    ThrowRequireMsg(mesh.destroy_entity(child), "Failed to destroy entity" << mesh.entity_key(child));
  }

  for (unsigned i=0; i<child_sides.size(); ++i)
  {
    stk::mesh::Entity side = child_sides[i];
    if (mesh.is_valid(side))
    {
      const unsigned num_side_elements = mesh.num_elements(side);
      if (0 == num_side_elements)
      {
        if(krinolog.shouldPrint(LOG_DEBUG))
          krinolog << "Destroying child side#" << mesh.identifier(side) << ":" << stk::diag::dendl;

        ThrowRequireMsg(mesh.destroy_entity(side), "Failed to destroy entity" << mesh.entity_key(side));
      }
      else if(krinolog.shouldPrint(LOG_DEBUG))
      {
        krinolog << "Can't destroy child side#" << mesh.identifier(side) << ":" << stk::diag::dendl;
      }
    }
  }

  for (unsigned i=0; i<child_edges.size(); ++i)
  {
    stk::mesh::Entity edge = child_edges[i];
    if (mesh.is_valid(edge))
    {
      const unsigned num_edge_elements = mesh.num_elements(edge);
      if (0 == num_edge_elements)
      {
        if(krinolog.shouldPrint(LOG_DEBUG))
          krinolog << "Destroying child edge#" << mesh.identifier(edge) << ":" << stk::diag::dendl;

        ThrowRequireMsg(mesh.destroy_entity(edge), "Failed to destroy entity" << mesh.entity_key(edge));
      }
    }
  }

  for (unsigned i=0; i<child_nodes.size(); ++i)
  {
    stk::mesh::Entity node = child_nodes[i];
    if (mesh.is_valid(node))
    {
      const unsigned num_node_elements = mesh.num_elements(node);
      const unsigned num_node_faces = mesh.num_faces(node);
      const unsigned num_node_edges = mesh.num_edges(node);
      if (0 == num_node_elements && 0 == num_node_faces && 0 == num_node_edges)
      {
        if(krinolog.shouldPrint(LOG_DEBUG))
          krinolog << "Destroying child node#" << mesh.identifier(node) << ":" << stk::diag::dendl;

        ThrowRequireMsg(mesh.destroy_entity(node), "Failed to destroy entity" << mesh.entity_key(node));
      }
    }
  }
}

double compute_child_position(const stk::mesh::BulkData & mesh, stk::mesh::Entity child, stk::mesh::Entity parent0, stk::mesh::Entity parent1)
{
  // if this function is not sufficiently precise or too expensive, then this position should be stored elsewhere
  const stk::mesh::Field<double> * const coords_field = reinterpret_cast<const stk::mesh::Field<double>*>(mesh.mesh_meta_data().coordinate_field());
  ThrowRequireMsg(nullptr != coords_field, "Coordinates must be defined.");

  double * child_coords = stk::mesh::field_data(*coords_field, child);
  double * parent0_coords = stk::mesh::field_data(*coords_field, parent0);
  double * parent1_coords = stk::mesh::field_data(*coords_field, parent1);

  const unsigned ndim = mesh.mesh_meta_data().spatial_dimension();
  unsigned best_dim = 0;
  double best_extent = std::abs(parent1_coords[0] - parent0_coords[0]);
  for (unsigned dim = 1; dim < ndim; ++dim)
  {
    const double extent = std::abs(parent1_coords[dim] - parent0_coords[dim]);
    if (extent > best_extent)
    {
      best_dim = dim;
      best_extent = extent;
    }
  }
  return std::abs(child_coords[best_dim] - parent0_coords[best_dim])/best_extent;
}

//--------------------------------------------------------------------------------

void
store_edge_node_parent_ids(const stk::mesh::BulkData & mesh,
    const FieldRef & parent_id_field,
    stk::mesh::Entity edge_node_entity,
    stk::mesh::EntityId parent0_id,
    stk::mesh::EntityId parent1_id)
{
  if (parent_id_field.type_is<unsigned>())
  {
    auto * stored_parent_ids = field_data<unsigned>(parent_id_field, edge_node_entity);
    ThrowAssertMsg(stored_parent_ids, "Node " << mesh.identifier(edge_node_entity)
        << " does not have the parent_ids field suggesting it is a mesh node.");
    stored_parent_ids[0] = parent0_id;
    stored_parent_ids[1] = parent1_id;
  }
  else if (parent_id_field.type_is<uint64_t>())
  {
    auto * stored_parent_ids = field_data<uint64_t>(parent_id_field, edge_node_entity);
    ThrowAssertMsg(stored_parent_ids, "Node " << mesh.identifier(edge_node_entity)
        << " does not have the parent_ids field suggesting it is a mesh node.");
    stored_parent_ids[0] = parent0_id;
    stored_parent_ids[1] = parent1_id;
  }
  else
  {
    ThrowRequireMsg(false, "Unsupported field type for parent_node_ids_field");
  }
}

//--------------------------------------------------------------------------------

std::array<stk::mesh::EntityId, 2>
get_edge_node_parent_ids(const stk::mesh::BulkData & mesh,
    const FieldRef & parent_id_field,
    const stk::mesh::Entity edge_node_entity)
{
  std::array<stk::mesh::EntityId, 2> parent_ids;

  if (parent_id_field.type_is<unsigned>())
  {
    const auto * stored_parent_ids = field_data<unsigned>(parent_id_field, edge_node_entity);
    ThrowRequireMsg(stored_parent_ids, "No SubElementNode found for node " << mesh.identifier(edge_node_entity)
        << ", but it does not have the parent_ids field suggesting it is a mesh node.");
    parent_ids[0] = stored_parent_ids[0];
    parent_ids[1] = stored_parent_ids[1];
  }
  else if (parent_id_field.type_is<uint64_t>())
  {
    const auto * stored_parent_ids = field_data<uint64_t>(parent_id_field, edge_node_entity);
    ThrowRequireMsg(stored_parent_ids, "No SubElementNode found for node " << mesh.identifier(edge_node_entity)
        << ", but it does not have the parent_ids field suggesting it is a mesh node.");
    parent_ids[0] = stored_parent_ids[0];
    parent_ids[1] = stored_parent_ids[1];
  }
  else
  {
    ThrowRequireMsg(false, "Unsupported field type for parent_node_ids_field");
  }
  return parent_ids;
}

//--------------------------------------------------------------------------------

void get_parent_nodes_from_child(const stk::mesh::BulkData & mesh,
    stk::mesh::Entity child, const FieldRef & parent_id_field,
    std::set<stk::mesh::Entity> & parent_nodes)
{
  if (has_field_data(parent_id_field, child))
  {
    auto parent_ids = get_edge_node_parent_ids(mesh, parent_id_field, child);
    const stk::mesh::Entity parent0 = mesh.get_entity(stk::topology::NODE_RANK, parent_ids[0]);
    const stk::mesh::Entity parent1 = mesh.get_entity(stk::topology::NODE_RANK, parent_ids[1]);
    ThrowRequire(mesh.is_valid(parent0) && mesh.is_valid(parent1));
    get_parent_nodes_from_child(mesh, parent0, parent_id_field, parent_nodes);
    get_parent_nodes_from_child(mesh, parent1, parent_id_field, parent_nodes);
  }
  else
  {
    parent_nodes.insert(child);
  }
}

//--------------------------------------------------------------------------------

void debug_print_selector_parts(const stk::mesh::Selector & selector)
{
  stk::mesh::PartVector parts;
  selector.get_parts(parts);
  krinolog << "Selector contains parts: " << stk::diag::push << stk::diag::dendl;
  for(stk::mesh::PartVector::const_iterator it = parts.begin(); it != parts.end(); ++it)
  {
    krinolog << (*it)->name() << stk::diag::dendl;
  }
  krinolog << stk::diag::pop << stk::diag::dendl;
}

//--------------------------------------------------------------------------------

stk::mesh::PartVector filter_non_io_parts(const stk::mesh::PartVector & all_parts)
{
  stk::mesh::PartVector io_parts;

  for(stk::mesh::PartVector::const_iterator it = all_parts.begin(); it != all_parts.end(); ++it)
  {
    if( stk::io::is_part_io_part(**it) )
    {
      io_parts.push_back(*it);
    }
  }

  return io_parts;
}

void
activate_selected_sides_touching_active_elements(stk::mesh::BulkData & mesh, const stk::mesh::Selector & side_selector, stk::mesh::Part & active_part)
{
  // This method requires AURA
  ThrowRequire(mesh.is_automatic_aura_on());

  mesh.modification_begin();
  stk::mesh::PartVector active_part_vec(1, &active_part);
  stk::mesh::PartVector inactive_part_vec;
  stk::mesh::Selector select_locally_owned = side_selector & mesh.mesh_meta_data().locally_owned_part();

  std::vector<stk::mesh::Entity> sides;
  stk::mesh::get_selected_entities( select_locally_owned, mesh.buckets( mesh.mesh_meta_data().side_rank() ), sides );
  for (auto && side : sides)
  {
    bool have_active_elem = false;
    const stk::mesh::Entity* side_elems = mesh.begin_elements(side);
    const unsigned num_side_elems = mesh.num_elements(side);
    for (unsigned ielem=0; ielem<num_side_elems; ++ielem)
    {
      if (mesh.bucket(side_elems[ielem]).member(active_part))
      {
        have_active_elem = true;
        break;
      }
    }

    if (have_active_elem) mesh.change_entity_parts(side, active_part_vec, inactive_part_vec);
    else  mesh.change_entity_parts(side, inactive_part_vec, active_part_vec);
  }
  mesh.modification_end();
}

void
get_partially_and_fully_coincident_elements(const stk::mesh::BulkData & mesh, stk::mesh::Entity elem, std::vector<stk::mesh::Entity> & coincident_elems)
{
  stk::topology elem_topology = mesh.bucket(elem).topology();
  const stk::mesh::Entity* elem_nodes = mesh.begin_nodes(elem);
  std::vector<stk::mesh::Entity> elem_side_nodes;
  std::vector<stk::mesh::Entity> elem_nbrs;
  std::vector<stk::mesh::Entity> nbr_side_nodes;

  coincident_elems.clear();
  const unsigned num_sides = elem_topology.num_sides();
  for (unsigned iside=0; iside<num_sides; ++iside)
  {
    stk::topology side_topology = elem_topology.side_topology(iside);
    elem_side_nodes.resize(side_topology.num_nodes());
    elem_topology.side_nodes(elem_nodes, iside, elem_side_nodes.data());

    std::vector<stk::mesh::Entity> sorted_elem_side_nodes = elem_side_nodes;
    std::sort(sorted_elem_side_nodes.begin(), sorted_elem_side_nodes.end(), stk::mesh::EntityLess(mesh));
    const unsigned unique_len = std::distance(sorted_elem_side_nodes.begin(), std::unique( sorted_elem_side_nodes.begin(), sorted_elem_side_nodes.end() ));
    const bool is_degenerate_side = unique_len < mesh.mesh_meta_data().spatial_dimension();
    if (is_degenerate_side) continue;

    stk::mesh::get_entities_through_relations(mesh, elem_side_nodes, stk::topology::ELEMENT_RANK, elem_nbrs);
    ThrowRequire(!elem_nbrs.empty());

    for (auto && nbr : elem_nbrs)
    {
      if (nbr == elem) continue;

      stk::topology nbr_topology = mesh.bucket(nbr).topology();
      const stk::mesh::Entity* nbr_nodes = mesh.begin_nodes(nbr);
      const unsigned num_nbr_sides = nbr_topology.num_sides();
      for (unsigned inbr_side=0; inbr_side<num_nbr_sides; ++inbr_side)
      {
        if (nbr_topology.side_topology(inbr_side) != side_topology)
        {
          // This doesn't handle a degenerate quad being coincident with tri,
          // maybe that will be needed eventually too.
          continue;
        }
        nbr_side_nodes.resize(side_topology.num_nodes());
        nbr_topology.side_nodes(nbr_nodes, inbr_side, nbr_side_nodes.data());

        stk::EquivalentPermutation result = side_topology.is_equivalent(nbr_side_nodes.data(), elem_side_nodes.data());
        if (result.is_equivalent && side_topology.is_positive_polarity(result.permutation_number))
        {
          coincident_elems.push_back(nbr);
          break;
        }
      }
    }
  }
}

bool connectivity_of_side_with_nodes_is_good(const stk::mesh::BulkData & mesh,
    const stk::mesh::Part & active_part,
    const std::vector<stk::mesh::Entity> & sideNodes,
    std::vector<stk::mesh::Entity> & sideEntities,
    std::vector<stk::mesh::Entity> & activeSideElements)
{
  stk::mesh::get_entities_through_relations(mesh, sideNodes, stk::topology::ELEMENT_RANK, sideEntities);
  ThrowRequire(!sideEntities.empty());

  activeSideElements.clear();
  for (auto && sideElement : sideEntities)
    if (mesh.bucket(sideElement).member(active_part) && !mesh.bucket(sideElement).topology().is_shell())
      activeSideElements.push_back(sideElement);

  if (activeSideElements.size() != 2)
  {
    const auto sideRank = mesh.mesh_meta_data().side_rank();
    krinolog << stk::diag::dendl;
    krinolog << "Found mismatched internal side with " << activeSideElements.size()
        << " elements, from analysis of side with nodes ";
    for (auto && side_node : sideNodes) krinolog << mesh.identifier(side_node) << " ";
    krinolog << stk::diag::dendl;

    stk::mesh::get_entities_through_relations(mesh, sideNodes, mesh.mesh_meta_data().side_rank(), sideEntities);
    if (!sideEntities.empty())
    {
      krinolog << "Sides using side nodes: " << stk::diag::dendl;
      for (auto && side : sideEntities)
        krinolog << debug_entity(mesh, side) << stk::diag::dendl;
    }

    krinolog << "Elements using side nodes: " << stk::diag::dendl;
    for (auto && sideElement : activeSideElements)
    {
      krinolog << debug_entity(mesh, sideElement) << stk::diag::dendl;
      for (auto && elemSide : StkMeshEntities{mesh.begin(sideElement, sideRank), mesh.end(sideElement, sideRank)})
        krinolog << debug_entity(mesh, elemSide) << stk::diag::dendl;
      for (auto && elemNode : StkMeshEntities{mesh.begin_nodes(sideElement), mesh.end_nodes(sideElement)})
        krinolog << debug_entity(mesh, elemNode) << stk::diag::dendl;
    }

    if (activeSideElements.size() != 2)
    return false;
  }
  return true;
}

bool
check_element_side_connectivity(const stk::mesh::BulkData & mesh, const stk::mesh::Part & exterior_boundary_part, const stk::mesh::Part & active_part)
{
  // This method exploits aura to make sure that there are exactly two active volume elements that support every side
  // that is not on the exterior boundary.
  if (!mesh.is_automatic_aura_on() && mesh.parallel_size() > 1)
  {
    // Skip check if we don't have aura
    return true;
  }

  bool found_mismatched_side = false;

  std::vector<stk::mesh::Entity> elementNodes;
  std::vector<stk::mesh::Entity> sideNodes;
  std::vector<stk::mesh::Entity> sideEntities;
  std::vector<stk::mesh::Entity> activeSideElements;

  const stk::mesh::BucketVector & elemBuckets = mesh.get_buckets( stk::topology::ELEMENT_RANK, active_part & mesh.mesh_meta_data().locally_owned_part() );

  for ( auto && bucket : elemBuckets )
  {
    stk::topology element_topology = bucket->topology();
    if (element_topology.is_shell()) continue;
    const unsigned num_sides = element_topology.num_sides();

    for ( auto && element : *bucket )
    {
      elementNodes.assign(mesh.begin_nodes(element), mesh.end_nodes(element));
      for (unsigned iside=0; iside<num_sides; ++iside)
      {
        stk::mesh::Entity existing_element_side = find_entity_by_ordinal(mesh, element, mesh.mesh_meta_data().side_rank(), iside);
        if (mesh.is_valid(existing_element_side) && mesh.bucket(existing_element_side).member(exterior_boundary_part))
        {
          continue;
        }

        stk::topology side_topology = element_topology.side_topology(iside);
        sideNodes.resize(side_topology.num_nodes());
        element_topology.side_nodes(elementNodes, iside, sideNodes.data());

        if (!connectivity_of_side_with_nodes_is_good(mesh, active_part, sideNodes, sideEntities, activeSideElements))
          found_mismatched_side = true;
      }
    }
  }

  const stk::mesh::BucketVector & sideBuckets = mesh.get_buckets( mesh.mesh_meta_data().side_rank(), active_part & mesh.mesh_meta_data().locally_owned_part() & !exterior_boundary_part );

  for ( auto && bucket : sideBuckets )
  {
    for ( auto && side : *bucket )
    {
      sideNodes.assign(mesh.begin_nodes(side), mesh.end_nodes(side));
      if (!connectivity_of_side_with_nodes_is_good(mesh, active_part, sideNodes, sideEntities, activeSideElements))
        found_mismatched_side = true;
    }
  }

  const int local_found_mismatched_side = found_mismatched_side;
  int global_found_mismatched_side = false;
  stk::all_reduce_max(mesh.parallel(), &local_found_mismatched_side, &global_found_mismatched_side, 1);
  return !global_found_mismatched_side;
}

bool
check_coincident_elements(const stk::mesh::BulkData & mesh, const stk::mesh::Part & active_part)
{
  // For this method to successfully detect issues, coincident elements must already exist on the same processor.

  std::vector< stk::mesh::Entity> elements;
  std::vector< stk::mesh::Entity> coincident_elements;
  std::vector< stk::mesh::Entity> element_nodes;

  ParallelErrorMessage err(mesh.parallel());

  stk::mesh::get_entities( mesh, stk::topology::ELEMENT_RANK, elements );
  for (auto && element : elements)
  {
    element_nodes.assign(mesh.begin_nodes(element), mesh.end_nodes(element));
    stk::mesh::get_entities_through_relations(mesh, element_nodes, stk::topology::ELEMENT_RANK, coincident_elements);
    ThrowRequire(!coincident_elements.empty());
    stk::topology element_topology = mesh.bucket(element).topology();

    bool coincident_element_error = false;
    for (auto && coincident : coincident_elements)
    {
      if (coincident == element || element_topology != mesh.bucket(coincident).topology()) continue;
      if (!element_topology.is_shell())
      {
        err << "Non-shell elements " << mesh.entity_key(element) << " and "
            << mesh.entity_key(coincident) << " are fully coincident (overlapping).\n";
        coincident_element_error = true;
      }
      else
      {
        stk::EquivalentPermutation result = element_topology.is_equivalent(mesh.begin_nodes(coincident), element_nodes.data());
        ThrowRequire(result.is_equivalent);
        if (result.permutation_number != stk::mesh::DEFAULT_PERMUTATION)
        {
          err << "Elements " << mesh.entity_key(element) << " and " << mesh.entity_key(coincident)
              << " are fully coincident shell elements but have different node order.\n";
          coincident_element_error = true;
        }
      }
    }

    // Detect partially coincident, non-shell elements that are both active
    if (!coincident_element_error && mesh.bucket(element).member(active_part) && !element_topology.is_shell() && element_topology.base() != stk::topology::BEAM_2)
    {
      get_partially_and_fully_coincident_elements(mesh, element, coincident_elements);

      for (auto && coincident : coincident_elements)
      {
        if (mesh.bucket(coincident).member(active_part) && !mesh.bucket(coincident).topology().is_shell() &&
            mesh.bucket(coincident).topology().base() != stk::topology::BEAM_2)
        {
          err << "Non-shell elements " << mesh.entity_key(element) << " with topology "
              << element_topology.name() << " and " << mesh.entity_key(coincident)
              << " with topology " << mesh.bucket(coincident).topology().name()
              << " are partially coincident (overlapping) and active.\n";
          coincident_element_error = true;
        }
      }
    }
  }
  auto err_msg = err.gather_message();
  krinolog << err_msg.second;
  return !err_msg.first;
}

bool
fix_coincident_element_ownership(stk::mesh::BulkData & mesh)
{
  // This method exploits aura to choose which processor the faces and edges should be owned by
  if (!mesh.is_automatic_aura_on())
  {
    // Make no changes, hope for the best.
    return false;
  }

  stk::mesh::MetaData & meta = mesh.mesh_meta_data();
  stk::mesh::Selector locally_owned_selector(meta.locally_owned_part());

  std::vector<stk::mesh::Entity> elems;
  std::vector<stk::mesh::Entity> coincident_elems;

  std::vector<stk::mesh::EntityProc> entities_to_move;

  stk::mesh::get_selected_entities( locally_owned_selector, mesh.buckets(stk::topology::ELEMENT_RANK), elems );

  for (auto && elem : elems)
  {
    get_partially_and_fully_coincident_elements(mesh, elem, coincident_elems);

    int new_owner = mesh.parallel_owner_rank(elem);
    for (auto && nbr : coincident_elems)
    {
      const int elem_owner = mesh.parallel_owner_rank(nbr);
      if (elem_owner < new_owner) new_owner = elem_owner;
    }

    if (new_owner != mesh.parallel_owner_rank(elem))
    {
      entities_to_move.push_back(stk::mesh::EntityProc(elem, new_owner));
    }
  }

  const int local_made_moves = !entities_to_move.empty();
  int global_made_moves = false;
  stk::all_reduce_max(mesh.parallel(), &local_made_moves, &global_made_moves, 1);

  if (global_made_moves)
  {
    mesh.change_entity_owner(entities_to_move);
  }
  return global_made_moves;
}

static int determine_new_owner_for_owned_face_or_edge(const stk::mesh::BulkData & mesh, const stk::mesh::Entity faceOrEdge)
{
  const int currentOwner = mesh.parallel_rank();
  int newOwner = mesh.parallel_size();
  for (auto elem : StkMeshEntities{mesh.begin_elements(faceOrEdge), mesh.end_elements(faceOrEdge)})
  {
    const int elemOwner = mesh.parallel_owner_rank(elem);
    if (elemOwner == currentOwner)
      return currentOwner;
    if (elemOwner < newOwner)
      newOwner = elemOwner;
  }
  return newOwner;
}

bool
fix_face_and_edge_ownership(stk::mesh::BulkData & mesh)
{
  // This method exploits aura to choose which processor the faces and edges should be owned by
  if (!mesh.is_automatic_aura_on())
  {
    // Make no changes, hope for the best.
    return false;
  }

  const int parallelRank = mesh.parallel_rank();
  stk::mesh::Selector locallyOwnedSelector(mesh.mesh_meta_data().locally_owned_part());

  std::vector<stk::mesh::EntityProc> entitiesToMove;

  for (auto && edgeOrFaceRank : {stk::topology::EDGE_RANK, stk::topology::FACE_RANK})
  {
    for (auto && bucket : mesh.get_buckets( edgeOrFaceRank, locallyOwnedSelector ))
    {
      for ( auto && faceOrEdge : *bucket )
      {
        const int newOwner = determine_new_owner_for_owned_face_or_edge(mesh, faceOrEdge);
        if (newOwner != parallelRank)
          entitiesToMove.push_back(stk::mesh::EntityProc(faceOrEdge, newOwner));
      }
    }
  }

  if(stk::is_true_on_any_proc(mesh.parallel(), !entitiesToMove.empty()))
  {
    mesh.change_entity_owner(entitiesToMove);
    return true;
  }

  return false;
}

static int determine_new_owner_for_owned_node_to_assure_active_owned_element(const stk::mesh::BulkData & mesh, const stk::mesh::Entity node, const stk::mesh::Part & activePart)
{
  const int currentOwner = mesh.parallel_rank();
  int newOwner = mesh.parallel_size();
  for (auto elem : StkMeshEntities{mesh.begin_elements(node), mesh.end_elements(node)})
  {
    if (mesh.bucket(elem).member(activePart))
    {
      const int elemOwner = mesh.parallel_owner_rank(elem);
      if (elemOwner == currentOwner)
        return currentOwner;
      if (elemOwner < newOwner)
        newOwner = elemOwner;
    }
  }
  if (newOwner == mesh.parallel_size())
    return currentOwner;
  return newOwner;
}

bool
fix_node_owners_to_assure_active_owned_element_for_node(stk::mesh::BulkData & mesh, const stk::mesh::Part & activePart)
{
  // This method exploits aura to choose which processor the faces and edges should be owned by
  if (!mesh.is_automatic_aura_on())
  {
    // Make no changes, hope for the best.
    return false;
  }

  const int parallelRank = mesh.parallel_rank();
  stk::mesh::Selector locallyOwnedSelector(mesh.mesh_meta_data().locally_owned_part());

  std::vector<stk::mesh::EntityProc> entitiesToMove;

  for (auto && bucket : mesh.get_buckets( stk::topology::NODE_RANK, locallyOwnedSelector ))
  {
    for ( auto && node : *bucket )
    {
      const int newOwner = determine_new_owner_for_owned_node_to_assure_active_owned_element(mesh, node, activePart);
      if (newOwner != parallelRank)
        entitiesToMove.push_back(stk::mesh::EntityProc(node, newOwner));
    }
  }

  if(stk::is_true_on_any_proc(mesh.parallel(), !entitiesToMove.empty()))
  {
    mesh.change_entity_owner(entitiesToMove);
    return true;
  }

  return false;
}

static int does_face_or_edge_have_an_element_with_the_same_owning_proc(const stk::mesh::BulkData & mesh, const stk::mesh::Entity faceOrEdge)
{
  int entityOwner = mesh.parallel_owner_rank(faceOrEdge);
  for (auto elem : StkMeshEntities{mesh.begin_elements(faceOrEdge), mesh.end_elements(faceOrEdge)})
    if (mesh.parallel_owner_rank(elem) == entityOwner)
      return true;
  return false;
}

bool
check_face_and_edge_ownership(const stk::mesh::BulkData & mesh)
{
  // This method exploits aura to choose which processor the faces and edges should be owned by
  if (!mesh.is_automatic_aura_on())
  {
    // Skip check if we don't have aura
    return true;
  }

  stk::mesh::Selector locallyOwnedSelector(mesh.mesh_meta_data().locally_owned_part());

  bool error = false;
  for (auto && edgeOrFaceRank : {stk::topology::EDGE_RANK, stk::topology::FACE_RANK})
  {
    for (auto && bucket : mesh.get_buckets( edgeOrFaceRank, locallyOwnedSelector ))
    {
      for ( auto && faceOrEdge : *bucket )
      {
        if (!does_face_or_edge_have_an_element_with_the_same_owning_proc(mesh, faceOrEdge))
        {
          error = true;
          krinolog << "Error: " << mesh.entity_key(faceOrEdge)
            << " is owned on processor " << mesh.parallel_owner_rank(faceOrEdge)
            << " but does not have any elements that are owned on this processor.\n";
          krinolog << debug_entity_1line(mesh, faceOrEdge) << "\n";
          for (auto elem : StkMeshEntities{mesh.begin_elements(faceOrEdge), mesh.end_elements(faceOrEdge)})
            krinolog << "  " << debug_entity_1line(mesh, elem) << "\n";
        }
      }
    }
  }
  return !(stk::is_true_on_any_proc(mesh.parallel(), error));
}

bool
set_region_id_on_unset_entity_and_neighbors(stk::mesh::BulkData & mesh, stk::mesh::Entity & entity, stk::mesh::EntityRank entity_rank, stk::mesh::FieldBase & region_id_field, const unsigned region_id)
{
  unsigned * region_id_data = stk::mesh::field_data(reinterpret_cast<stk::mesh::Field<unsigned>&>(region_id_field), entity);
  ThrowAssert(NULL != region_id_data);
  if (*region_id_data != 0)
  {
    ThrowAssert(*region_id_data == region_id);
    return false;
  }

  *region_id_data = region_id;

  // now visit neighbors
  const stk::mesh::Entity* begin_nodes = mesh.begin_nodes(entity);
  const stk::mesh::Entity* end_nodes = mesh.end_nodes(entity);

  for (const stk::mesh::Entity* it_node = begin_nodes; it_node != end_nodes; ++it_node)
  {
    stk::mesh::Entity node = *it_node;

    const stk::mesh::Entity* node_entities_begin = mesh.begin(node, entity_rank);
    const stk::mesh::Entity* node_entities_end = mesh.end(node, entity_rank);

    for (const stk::mesh::Entity* it_nbr = node_entities_begin; it_nbr != node_entities_end; ++it_nbr)
    {
      stk::mesh::Entity neighbor = *it_nbr;
      set_region_id_on_unset_entity_and_neighbors(mesh, neighbor, entity_rank, region_id_field, region_id);
    }
  }

  return true;
}

void
identify_isolated_regions(stk::mesh::Part & part, stk::mesh::FieldBase & region_id_field)
{
  stk::mesh::BulkData & mesh = part.mesh_bulk_data();
  stk::mesh::MetaData & meta = mesh.mesh_meta_data();

  stk::mesh::EntityRank entity_rank = part.primary_entity_rank();
  stk::mesh::Selector selector = part & meta.locally_owned_part();

  std::vector<stk::mesh::Entity> entities;
  stk::mesh::get_selected_entities( selector, mesh.buckets(entity_rank), entities );

  // initialize region id
  const stk::mesh::BucketVector & buckets = mesh.get_buckets(entity_rank, selector);
  for (stk::mesh::BucketVector::const_iterator ib = buckets.begin(); ib != buckets.end(); ++ib )
  {
    const stk::mesh::Bucket & b = **ib;
    const size_t length = b.size();
    unsigned * region_id = stk::mesh::field_data(reinterpret_cast<stk::mesh::Field<unsigned>&>(region_id_field), b);
    ThrowAssert(NULL != region_id);
    std::fill(region_id, region_id+length, 0);
  }

  unsigned local_region_counter = 0;
  unsigned region_id = 1 + local_region_counter*mesh.parallel_size() + mesh.parallel_rank();

  for (std::vector<stk::mesh::Entity>::iterator it_entity=entities.begin(); it_entity!=entities.end(); ++it_entity)
  {
    stk::mesh::Entity entity = *it_entity;

    const bool made_changes = set_region_id_on_unset_entity_and_neighbors(mesh, entity, entity_rank, region_id_field, region_id);
    if (made_changes)
    {
      ++local_region_counter;
      region_id = 1 + local_region_counter*mesh.parallel_size() + mesh.parallel_rank();
    }
  }
}

const unsigned * get_side_node_ordinals(stk::topology topology, unsigned side_ordinal)
{
  static std::vector< std::vector< std::vector<unsigned> > > all_side_node_ordinals(stk::topology::BEGIN_TOPOLOGY + stk::topology::NUM_TOPOLOGIES);
  std::vector< std::vector<unsigned> > & topology_side_node_ordinals =  all_side_node_ordinals[topology.value()];
  if (topology_side_node_ordinals.empty())
  {
    const size_t num_sides = topology.num_sides();
    topology_side_node_ordinals.resize(num_sides);
    for (size_t side_index = 0; side_index < num_sides; ++side_index)
    {
      std::vector<unsigned> & side_node_ordinals =  topology_side_node_ordinals[side_index];
      side_node_ordinals.resize(topology.side_topology(side_index).num_nodes());
      topology.side_node_ordinals(side_index, side_node_ordinals.begin());
    }
  }
  return &topology_side_node_ordinals[side_ordinal][0];
}

const unsigned * get_edge_node_ordinals(stk::topology topology, unsigned edge_ordinal)
{
  static std::vector< std::vector< std::vector<unsigned> > > all_edge_node_ordinals(stk::topology::BEGIN_TOPOLOGY + stk::topology::NUM_TOPOLOGIES);
  std::vector< std::vector<unsigned> > & topology_edge_node_ordinals =  all_edge_node_ordinals[topology.value()];
  if (topology_edge_node_ordinals.empty())
  {
    const size_t num_edges = topology.num_edges();
    topology_edge_node_ordinals.resize(num_edges);
    for (size_t edge_index = 0; edge_index < num_edges; ++edge_index)
    {
      std::vector<unsigned> & edge_node_ordinals =  topology_edge_node_ordinals[edge_index];
      edge_node_ordinals.resize(topology.edge_topology(edge_index).num_nodes());
      topology.edge_node_ordinals(edge_index, edge_node_ordinals.begin());
    }
  }
  return &topology_edge_node_ordinals[edge_ordinal][0];
}

stk::mesh::PartVector get_common_io_parts(const stk::mesh::BulkData & mesh, const std::vector<stk::mesh::Entity> entities)
{
  stk::mesh::PartVector common_io_parts;

  bool first = true;
  for (auto&& entity : entities)
  {
    stk::mesh::PartVector entity_io_parts = filter_non_io_parts(mesh.bucket(entity).supersets());
    std::sort(entity_io_parts.begin(), entity_io_parts.end());
    if (first)
    {
      first = false;
      common_io_parts.swap(entity_io_parts);
    }
    else
    {
      stk::mesh::PartVector working_set;
      working_set.swap(common_io_parts);
      std::set_intersection(working_set.begin(),working_set.end(),entity_io_parts.begin(),entity_io_parts.end(),std::back_inserter(common_io_parts));
    }
  }
  return common_io_parts;
}

stk::mesh::PartVector get_removable_parts(const stk::mesh::BulkData & mesh, const stk::mesh::Bucket & bucket)
{
  stk::mesh::PartVector removable_parts;
  for ( auto&& part : bucket.supersets() )
  {
    stk::mesh::EntityRank part_rank = part->primary_entity_rank();
    if ((part_rank == stk::topology::INVALID_RANK || part_rank == bucket.entity_rank()) &&
        (!stk::mesh::is_auto_declared_part(*part) || stk::mesh::is_topology_root_part(*part)))
    {
      removable_parts.push_back(part);
    }
  }
  return removable_parts;
}

stk::mesh::PartVector get_removable_parts(const stk::mesh::BulkData & mesh, const stk::mesh::Entity entity)
{
  return get_removable_parts(mesh, mesh.bucket(entity));
}

void set_relation_permutation(stk::mesh::BulkData & mesh, stk::mesh::Entity from, stk::mesh::Entity to, stk::mesh::ConnectivityOrdinal to_ord, stk::mesh::Permutation to_permutation)
{
  const stk::mesh::EntityRank from_rank = mesh.entity_rank(from);
  const stk::mesh::EntityRank to_rank   = mesh.entity_rank(to);

  stk::mesh::Entity const*              fwd_rels  = mesh.begin(from, to_rank);
  stk::mesh::ConnectivityOrdinal const* fwd_ords  = mesh.begin_ordinals(from, to_rank);
  stk::mesh::Permutation *              fwd_perms = const_cast<stk::mesh::Permutation*>(mesh.begin_permutations(from, to_rank));
  const int                  num_fwd   = mesh.num_connectivity(from, to_rank);

  stk::mesh::Entity const*              back_rels  = mesh.begin(to, from_rank);
  stk::mesh::ConnectivityOrdinal const* back_ords  = mesh.begin_ordinals(to, from_rank);
  stk::mesh::Permutation *              back_perms = const_cast<stk::mesh::Permutation*>(mesh.begin_permutations(to, from_rank));
  const int                  num_back   = mesh.num_connectivity(to,from_rank);

  // Find and change fwd connectivity
  for (int i = 0; i < num_fwd; ++i, ++fwd_rels, ++fwd_ords, ++fwd_perms) {
    // Allow clients to make changes to permutation
    // Permutations do not affect Relation ordering, so this is safe.
    if (*fwd_rels == to && *fwd_ords == to_ord) {
      *fwd_perms = to_permutation;
    }
  }

  // Find and change back connectivity
  for (int i = 0; i < num_back; ++i, ++back_rels, ++back_ords, ++back_perms) {
    // Allow clients to make changes to permutation
    // Permutations do not affect Relation ordering, so this is safe.
    if (*back_rels == from && *back_ords == to_ord) {
      *back_perms = to_permutation;
    }
  }
}

std::pair<stk::mesh::ConnectivityOrdinal, stk::mesh::Permutation>
determine_shell_side_ordinal_and_permutation(const stk::mesh::BulkData & mesh, stk::mesh::Entity shell, stk::mesh::Entity side)
{
  // The input shell may or may not be attached to side, but should be coincident with the side.
  // We will figure out what the ordinal for this shell-side relation based on the existing
  // connectivity of the side.
  const unsigned num_side_elems = mesh.num_elements(side);
  ThrowRequireMsg(num_side_elems > 0, "Cannot determine shell_side_ordinal for completely disconnected side.");

  stk::mesh::EntityRank side_rank = mesh.mesh_meta_data().side_rank();
  stk::topology side_topology = mesh.bucket(side).topology();
  const stk::mesh::Entity * side_elems = mesh.begin_elements(side);
  const stk::mesh::ConnectivityOrdinal * side_elem_ordinals = mesh.begin_element_ordinals(side);
  const stk::mesh::Permutation * side_elem_permutations = mesh.begin_element_permutations(side);

  std::pair<stk::mesh::ConnectivityOrdinal, stk::mesh::Permutation> result(stk::mesh::INVALID_CONNECTIVITY_ORDINAL, stk::mesh::INVALID_PERMUTATION);
  for (unsigned it = 0; it<num_side_elems; ++it)
  {
    stk::mesh::Entity elem = side_elems[it];
    if (mesh.bucket(elem).topology().is_shell())
    {
      const auto shell_relationship = std::make_pair(side_elem_ordinals[it],side_elem_permutations[it]);
      ThrowRequireMsg(0==it || shell_relationship == result, "Side is inconsistently attached to a shell.");
      result = shell_relationship;
    }
    else
    {
      const stk::mesh::ConnectivityOrdinal ordinal = side_elem_ordinals[it];
      stk::EquivalentPermutation equiv = stk::mesh::sub_rank_equivalent(mesh, elem, ordinal, side_rank, mesh.begin_nodes(shell));
      ThrowRequireMsg(equiv.is_equivalent, "Logic error in figuring out element-side connectivity.");
      const stk::mesh::ConnectivityOrdinal shell_side_ordinal = (side_topology.is_positive_polarity(equiv.permutation_number)) ? stk::mesh::ConnectivityOrdinal(1) : stk::mesh::ConnectivityOrdinal(0);
      const stk::mesh::Permutation shell_side_perm = determine_permutation(mesh, shell, side, shell_side_ordinal);
      const auto shell_relationship = std::make_pair(shell_side_ordinal, shell_side_perm);
      ThrowRequireMsg(0==it || shell_relationship == result, "Side is inconsistently attached to a volume element.");
      result = shell_relationship;
    }
  }
  return result;
}

stk::mesh::Permutation
determine_permutation(const stk::mesh::BulkData & mesh, const stk::mesh::Entity entity, const stk::mesh::Entity relative, const stk::mesh::ConnectivityOrdinal ordinal)
{
  const stk::mesh::EntityRank relative_rank = mesh.entity_rank(relative);
  ThrowAssert(mesh.entity_rank(entity) > relative_rank);

  const stk::mesh::Entity * relative_nodes = mesh.begin_nodes(relative);

  const stk::EquivalentPermutation equiv = stk::mesh::sub_rank_equivalent(mesh, entity, ordinal, relative_rank, relative_nodes);
  if(!equiv.is_equivalent)
  {
    ThrowErrorMsg("Could not find connection between " << mesh.entity_key(entity) <<" and "
        << mesh.entity_key(relative) << " with ordinal " << ordinal
        << debug_entity(mesh, entity) << debug_entity(mesh, relative));
  }

  return static_cast<stk::mesh::Permutation>(equiv.permutation_number);
}

std::pair<stk::mesh::ConnectivityOrdinal, stk::mesh::Permutation>
determine_ordinal_and_permutation(const stk::mesh::BulkData & mesh, const stk::mesh::Entity entity, const stk::mesh::Entity relative)
{
  const stk::mesh::EntityRank relative_rank = mesh.entity_rank(relative);
  ThrowAssert(mesh.entity_rank(entity) > relative_rank);
  stk::topology relative_topology = mesh.bucket(relative).topology();
  ThrowAssert(relative_topology.num_nodes() == mesh.num_nodes(relative));

  const stk::mesh::Entity * relative_nodes = mesh.begin_nodes(relative);

  stk::topology entity_topology = mesh.bucket(entity).topology();
  const bool looking_for_shell_side = entity_topology.is_shell() && relative_rank == mesh.mesh_meta_data().side_rank();

  for(size_t i = 0; i < entity_topology.num_sub_topology(relative_rank); ++i)
  {
    if (entity_topology.sub_topology(relative_rank, i) == relative_topology)
    {
      const stk::EquivalentPermutation equiv = stk::mesh::sub_rank_equivalent(mesh, entity, stk::mesh::ConnectivityOrdinal(i), relative_rank, relative_nodes);
      const bool match = equiv.is_equivalent && (!looking_for_shell_side || equiv.permutation_number < relative_topology.num_positive_permutations());
      if(match)
      {
        return std::pair<stk::mesh::ConnectivityOrdinal, stk::mesh::Permutation>(static_cast<stk::mesh::ConnectivityOrdinal>(i), static_cast<stk::mesh::Permutation>(equiv.permutation_number));
      }
    }
  }
  ThrowRuntimeError("Could not find connection between " << mesh.entity_key(entity) << " and " << mesh.entity_key(relative) << debug_entity(mesh, entity) << debug_entity(mesh, relative));
}

static
void pack_owned_entities_for_ghosting_procs(const stk::mesh::BulkData & mesh,
    const std::vector<stk::mesh::Entity> & entities,
    stk::CommSparse &commSparse)
{
  std::vector<int> elemCommProcs;
  stk::pack_and_communicate(commSparse,[&]()
  {
    for (auto entity : entities)
    {
      if (mesh.bucket(entity).owned())
      {
        mesh.comm_procs(entity, elemCommProcs);
        for (int procId : elemCommProcs)
          if (procId != commSparse.parallel_rank())
            commSparse.send_buffer(procId).pack(mesh.entity_key(entity));
      }
    }
  });
}

static
void unpack_ghosted_entities_from_owners(const stk::mesh::BulkData & mesh,
    std::vector<stk::mesh::Entity> & entities,
    stk::CommSparse &commSparse)
{
  stk::unpack_communications(commSparse, [&](int procId)
  {
    stk::CommBuffer & buffer = commSparse.recv_buffer(procId);

    while ( buffer.remaining() )
    {
      stk::mesh::EntityKey entityKey;
      commSparse.recv_buffer(procId).unpack(entityKey);
      stk::mesh::Entity entity = mesh.get_entity(entityKey);
      ThrowRequire(mesh.is_valid(entity));
      entities.push_back(entity);
    }
  });
}

void communicate_owned_entities_to_ghosting_procs(const stk::mesh::BulkData & mesh, std::vector<stk::mesh::Entity> & entities)
{
  stk::CommSparse commSparse(mesh.parallel());
  pack_owned_entities_for_ghosting_procs(mesh, entities, commSparse);
  unpack_ghosted_entities_from_owners(mesh, entities, commSparse);
}

} // namespace krino
