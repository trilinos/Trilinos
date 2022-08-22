#include "Akri_StkMeshBuilder.hpp"
#include <Akri_AuxMetaData.hpp>
#include <Akri_Phase_Support.hpp>
#include <stk_util/util/ReportHandler.hpp>
#include <stk_util/util/SortAndUnique.hpp>
#include <stk_math/StkVector.hpp>
#include <stk_mesh/base/FieldBase.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_io/IossBridge.hpp>
#include "../../../stk/stk_mesh/stk_mesh/base/SkinBoundary.hpp"

namespace krino
{

template <int DIM>
StkMeshBuilder<DIM>::StkMeshBuilder(stk::mesh::BulkData & mesh, const stk::ParallelMachine comm)
: mMesh(mesh), mAuxMeta(AuxMetaData::create(mesh.mesh_meta_data())), mPhaseSupport(Phase_Support::get(mesh.mesh_meta_data())), mComm(comm)
{
  declare_coordinates();
  mMesh.mesh_meta_data().use_simple_fields();
}

template <int DIM>
void StkMeshBuilder<DIM>::declare_coordinates()
{
  stk::mesh::Field<double> & coordsField = mMesh.mesh_meta_data().template declare_field<double>(
      stk::topology::NODE_RANK, "coordinates", 1u);
  stk::mesh::put_field_on_entire_mesh(coordsField, DIM);
  stk::io::set_field_role(coordsField, Ioss::Field::MESH);
}

template <int DIM>
void StkMeshBuilder<DIM>::create_block_parts(const std::vector<unsigned> &elementBlockIDs)
{
  stk::topology simplexTopology = ((DIM == 2) ? stk::topology::TRIANGLE_3_2D : stk::topology::TETRAHEDRON_4);

  for (unsigned blockId : elementBlockIDs)
  {
    const std::string blockName = "block_"+std::to_string(blockId);
    stk::mesh::Part &part = mMesh.mesh_meta_data().declare_part_with_topology(blockName, simplexTopology);
    mMesh.mesh_meta_data().set_part_id(part, blockId);
    stk::io::put_io_part_attribute(part);
  }
}

std::string get_surface_name(const unsigned sidesetId)
{
  const std::string surfaceName = "surface_"+std::to_string(sidesetId);
  return surfaceName;
}

template <int DIM>
void StkMeshBuilder<DIM>::create_sideset_parts(const std::vector<unsigned> &sidesetIds)
{
    for (unsigned sidesetId : sidesetIds)
    {
      stk::mesh::Part &sidesetPart = mMesh.mesh_meta_data().declare_part(get_surface_name(sidesetId), mMesh.mesh_meta_data().side_rank());
      mMesh.mesh_meta_data().set_part_id(sidesetPart, sidesetId);
      stk::io::put_io_part_attribute(sidesetPart);
    }
}

std::vector<stk::mesh::PartVector> convert_vector_of_vector_of_sideset_ids_to_parts(const stk::mesh::MetaData & meta, const std::vector<std::vector<unsigned>>& vectorOfVectorsOfSidesetIds)
{
    std::vector<stk::mesh::PartVector> addParts {vectorOfVectorsOfSidesetIds.size()};
    for(size_t i {0}; i < vectorOfVectorsOfSidesetIds.size(); ++i)
    {
        addParts[i].reserve(vectorOfVectorsOfSidesetIds[i].size());
        for(const size_t sidesetId : vectorOfVectorsOfSidesetIds[i])
            addParts[i].push_back(meta.get_part(get_surface_name(sidesetId)));
    }
    return addParts;
}

template <int DIM>
void StkMeshBuilder<DIM>::add_sides_to_sidesets(const std::vector<stk::mesh::Entity> &sides, const std::vector<std::vector<unsigned>> &sidesetIdsPerSide)
{
    ThrowRequireWithSierraHelpMsg(sides.size() == sidesetIdsPerSide.size());
    const std::vector<stk::mesh::PartVector> addParts = convert_vector_of_vector_of_sideset_ids_to_parts(mMesh.mesh_meta_data(), sidesetIdsPerSide);
    const std::vector<stk::mesh::PartVector> remParts(sidesetIdsPerSide.size(), stk::mesh::PartVector{});
    mMesh.batch_change_entity_parts(sides, addParts, remParts);
}

template <int DIM>
stk::mesh::Entity StkMeshBuilder<DIM>::get_side_with_nodes(const std::vector<stk::mesh::Entity> &nodesOfSide) const
{
  std::vector<stk::mesh::Entity> sidesWithNodes;

  stk::mesh::get_entities_through_relations(mMesh, nodesOfSide, mMesh.mesh_meta_data().side_rank(), sidesWithNodes);
  ThrowRequireMsg(sidesWithNodes.size() == 1, "Expected to find one side with nodes, but found " << sidesWithNodes.size());
  return sidesWithNodes[0];
}

template <int DIM>
void StkMeshBuilder<DIM>::set_node_coordinates(const stk::mesh::Entity node, const stk::math::Vector3d &newLoc)
{
    double* node_coords = (double*)stk::mesh::field_data(*mMesh.mesh_meta_data().coordinate_field(), node);
    node_coords[0] = newLoc[0];
    node_coords[1] = newLoc[1];
    if (mMesh.mesh_meta_data().spatial_dimension() == 3) node_coords[2] = newLoc[2];
}

template <int DIM>
stk::mesh::Entity StkMeshBuilder<DIM>::create_node(const stk::math::Vector3d &loc, const std::vector<int> &sharingProcs, stk::mesh::EntityId nodeId)
{
    stk::mesh::Entity node = mMesh.declare_node(nodeId);

    int proc = mMesh.parallel_rank();
    for(int sharingProc : sharingProcs)
    {
        if ( sharingProc != proc)
            mMesh.add_node_sharing(node, sharingProc);
    }

    set_node_coordinates(node, loc);
    return node;
}

stk::mesh::Part * get_block_part(const stk::mesh::MetaData &meta, const unsigned blockId)
{
    stk::mesh::Part *blockPart{nullptr};
    for (stk::mesh::Part * part : meta.get_parts())
    {
        if (part->primary_entity_rank() == stk::topology::ELEM_RANK && (unsigned)part->id() == blockId)
        {
            blockPart = part;
            break;
        }
    }
    ThrowRequireMsg(blockPart!=nullptr, "Can't find a block with id " << blockId);
    return blockPart;
}

template <int DIM>
void StkMeshBuilder<DIM>::create_boundary_sides()
{
    stk::mesh::create_exposed_block_boundary_sides(mMesh, mMesh.mesh_meta_data().universal_part(), {&mAuxMeta.exposed_boundary_part()});
}

template <int DIM>
bool StkMeshBuilder<DIM>::check_boundary_sides() const
{
  return stk::mesh::check_exposed_block_boundary_sides(mMesh, mMesh.mesh_meta_data().universal_part(), mAuxMeta.exposed_boundary_part());
}

template <int DIM>
void StkMeshBuilder<DIM>::create_block_boundary_sides()
{
  stk::mesh::create_exposed_block_boundary_sides(mMesh, mMesh.mesh_meta_data().universal_part(), {&mAuxMeta.block_boundary_part()});
}

template <int DIM>
bool StkMeshBuilder<DIM>::check_block_boundary_sides() const
{
  return stk::mesh::check_interior_block_boundary_sides(mMesh, mMesh.mesh_meta_data().universal_part(), mAuxMeta.block_boundary_part());
}

template <int DIM>
stk::mesh::Entity StkMeshBuilder<DIM>::create_element(const std::vector<stk::mesh::Entity> &nodes, stk::mesh::EntityId elementId, unsigned blockId)
{
    const stk::mesh::Part *blockPart = get_block_part(mMesh.mesh_meta_data(), blockId);
    stk::mesh::Entity element = mMesh.declare_element(elementId, stk::mesh::ConstPartVector{blockPart});
    unsigned idx = 0;
    for (auto nd : nodes)
        mMesh.declare_relation(element, stk::mesh::Entity(nd), idx++);
    return element;
}

template <int DIM>
std::vector<stk::mesh::Entity>
StkMeshBuilder<DIM>::create_parallel_nodes(const std::vector<stk::math::Vec<double,DIM>>& nodeLocs,
    const std::map<unsigned,std::vector<int>> &nodeIndicesWithSharingProcs,
    const std::vector<stk::mesh::EntityId> & assignedGlobalNodeIdsforAllNodes)
{
    std::vector<stk::mesh::Entity> nodesWhichAreValidIfTheyExistOnProc(nodeLocs.size(), stk::mesh::Entity());
    int curProc = stk::parallel_machine_rank(mComm);
    for(auto &nodeIndexWithSharingProcs : nodeIndicesWithSharingProcs)
    {
        unsigned nodeIndex = nodeIndexWithSharingProcs.first;
        stk::mesh::EntityId nodeGlobalId = assignedGlobalNodeIdsforAllNodes[nodeIndex];
        const std::vector<int> &sharingProcs = nodeIndexWithSharingProcs.second;

        if (std::find(sharingProcs.begin(), sharingProcs.end(), curProc) != sharingProcs.end() )
            nodesWhichAreValidIfTheyExistOnProc[nodeIndex] = create_node(stk::math::Vector3d{nodeLocs[nodeIndex].data(),DIM},
                sharingProcs,
                nodeGlobalId);
    }
    return nodesWhichAreValidIfTheyExistOnProc;
}

std::vector<stk::mesh::EntityId> get_ids_available_for_rank(stk::mesh::BulkData & mesh, stk::mesh::EntityRank rank, size_t numRequested)
{
    stk::mesh::EntityIdVector requestedIds;
    mesh.generate_new_ids(rank, numRequested, requestedIds);
    std::vector<stk::mesh::EntityId> idsToReturn(requestedIds.begin(), requestedIds.end());
    std::reverse(idsToReturn.begin(), idsToReturn.end());
    return idsToReturn;
}

template <int DIM>
std::vector<stk::mesh::Entity>
StkMeshBuilder<DIM>::create_parallel_elements(const std::vector<std::array<unsigned, NPE>> &elementConn,
    const std::vector<unsigned> &elementBlockIDs,
    const std::vector<int> &elementProcOwners,
    const std::vector<stk::mesh::Entity>& nodesWhichAreValidIfTheyExistOnProc)
{
    const int proc = stk::parallel_machine_rank(mComm);

    size_t numOwnedElements = 0;
    for (int elemProc : elementProcOwners)
      if (elemProc == proc) ++numOwnedElements;

    std::vector<stk::mesh::EntityId> elementIds = get_ids_available_for_rank(mMesh, stk::topology::ELEM_RANK, numOwnedElements);

    std::vector<stk::mesh::Entity> ownedElems;
    for (size_t iElem=0; iElem<elementConn.size(); ++iElem)
    {
        if (elementProcOwners[iElem] == proc)
        {
            ThrowRequireWithSierraHelpMsg((size_t)NPE == elementConn[iElem].size());

            stk::mesh::EntityId elementId = elementIds.back();
            elementIds.pop_back();

            std::vector<stk::mesh::Entity> oneElementConnWithLocalIds(NPE);
            for(unsigned i = 0; i < NPE; i++)
              oneElementConnWithLocalIds[i] = nodesWhichAreValidIfTheyExistOnProc[elementConn[iElem][i]];

            stk::mesh::Entity elem = create_element(oneElementConnWithLocalIds, elementId, elementBlockIDs[iElem]);
            ownedElems.push_back(elem);
        }
    }

    return ownedElems;
}

template <int DIM>
std::map<unsigned,std::vector<int>>
StkMeshBuilder<DIM>::build_node_sharing_procs(const std::vector<std::array<unsigned, NPE>> &elementConn,
    const std::vector<int> &elementProcOwners) const
{
  std::map<unsigned,std::vector<int>> nodeIndicesWithSharingProcs;
  for (size_t iElem=0; iElem<elementConn.size(); ++iElem)
      for (auto && node : elementConn[iElem])
          nodeIndicesWithSharingProcs[node].push_back(elementProcOwners[iElem]);

  for (auto && entry : nodeIndicesWithSharingProcs)
  {
      auto & procs = entry.second;
      stk::util::sort_and_unique(procs);
  }
  return nodeIndicesWithSharingProcs;
}

template <int DIM>
std::map<unsigned,std::vector<int>>
StkMeshBuilder<DIM>::build_node_sharing_procs_for_all_nodes_on_all_procs(const unsigned numNodes, const unsigned numProcs) const
{
  std::map<unsigned,std::vector<int>> nodeIndicesWithSharingProcs;
  for (unsigned iNode{0}; iNode<numNodes; ++iNode)
    for (unsigned proc=0; proc<numProcs; ++proc)
      nodeIndicesWithSharingProcs[iNode].push_back(proc);
  return nodeIndicesWithSharingProcs;
}

template <int DIM>
void StkMeshBuilder<DIM>::build_mesh(const std::vector<stk::math::Vec<double,DIM>> &nodeLocs,
    const std::vector<std::vector<std::array<unsigned, NPE>>> &elementConnPerProc,
    const unsigned blockId)
{
    ThrowRequireWithSierraHelpMsg(elementConnPerProc.size() == (size_t)stk::parallel_machine_size(mComm));
    std::vector<std::array<unsigned, NPE>> elementConn;
    std::vector<unsigned> elementBlockIDs;
    std::vector<int> elementProcOwners;
    for (unsigned proc=0; proc<elementConnPerProc.size(); ++proc)
    {
      for (auto && elemConn : elementConnPerProc[proc])
      {
        elementConn.push_back(elemConn);
        elementBlockIDs.push_back(blockId);
        elementProcOwners.push_back(proc);
      }
    }

    build_mesh_with_all_needed_block_ids(nodeLocs, elementConn, elementBlockIDs, {blockId}, elementProcOwners);
}

template <int DIM>
void StkMeshBuilder<DIM>::build_mesh(const std::vector<stk::math::Vec<double,DIM>> &nodeLocs,
    const std::vector<std::array<unsigned, NPE>> &elementConn,
    const std::vector<unsigned> &elementBlockIDs,
    const std::vector<int> &specifiedElementProcOwners)
{
    std::vector<unsigned> allBlockIDs = elementBlockIDs;
    stk::util::sort_and_unique(allBlockIDs);

    build_mesh_with_all_needed_block_ids(nodeLocs, elementConn, elementBlockIDs, allBlockIDs, specifiedElementProcOwners);
}

template <int DIM>
void StkMeshBuilder<DIM>::build_mesh_nodes_and_elements(
    const std::vector<stk::math::Vec<double,DIM>> &nodeLocs,
    const std::vector<std::array<unsigned, NPE>> &elementConn,
    const std::vector<unsigned> &elementBlockIDs,
    const std::vector<int> &specifiedElementProcOwners
)
{
    create_block_parts(elementBlockIDs);

    const size_t numGlobalElems = elementConn.size();
    std::vector<int> elementProcOwners = specifiedElementProcOwners;
    if (elementProcOwners.empty()) // Put all elements on proc 0 if called with empty specifiedElementProcOwners
      elementProcOwners.assign(numGlobalElems, 0);

    ThrowRequireWithSierraHelpMsg(elementBlockIDs.size() == numGlobalElems);
    ThrowRequireWithSierraHelpMsg(elementProcOwners.size() == numGlobalElems);

    mAssignedGlobalNodeIdsforAllNodes.resize(nodeLocs.size());
    for (unsigned iNode=0; iNode<mAssignedGlobalNodeIdsforAllNodes.size(); ++iNode)
      mAssignedGlobalNodeIdsforAllNodes[iNode] = iNode+101;

    const std::map<unsigned,std::vector<int>> nodeIndicesWithSharingProcs =
        (0 == numGlobalElems) ?
        build_node_sharing_procs_for_all_nodes_on_all_procs(nodeLocs.size(), stk::parallel_machine_size(mComm)) :
        build_node_sharing_procs(elementConn, elementProcOwners);

    mMesh.modification_begin();
    const auto nodeHandlesWhichAreValidForNodesThatExistOnProc = create_parallel_nodes(nodeLocs, nodeIndicesWithSharingProcs, mAssignedGlobalNodeIdsforAllNodes);
    mOwnedElems = create_parallel_elements(elementConn, elementBlockIDs, elementProcOwners, nodeHandlesWhichAreValidForNodesThatExistOnProc);
    mMesh.modification_end();

    create_boundary_sides();
    create_block_boundary_sides();
}

template <int DIM>
void StkMeshBuilder<DIM>::build_mesh_with_all_needed_block_ids
(
    const std::vector<stk::math::Vec<double,DIM>> &nodeLocs,
    const std::vector<std::array<unsigned, NPE>> &elementConn,
    const std::vector<unsigned> &elementBlockIDs,
    const std::vector<unsigned> &allBlocksIncludingThoseThatDontHaveElements,
    const std::vector<int> &specifiedElementProcOwners
)
{
    build_mesh_nodes_and_elements(nodeLocs, elementConn, elementBlockIDs, specifiedElementProcOwners);
}

// Explicit template instantiation
template class StkMeshBuilder<2>;
template class StkMeshBuilder<3>;

}
