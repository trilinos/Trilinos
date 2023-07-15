#include <Akri_VolumePreservingSnappingLimiter.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <Akri_FieldRef.hpp>
#include <Akri_MeshHelpers.hpp>
#include <stk_math/StkVector.hpp>
#include <Akri_AuxMetaData.hpp>
#include <Akri_DiagWriter.hpp>

namespace krino {

static void replace_coordinates_of_node_with_new_location(const StkMeshEntities & elemNodes, const stk::mesh::Entity node, const stk::math::Vector3d & newNodeLoc, std::vector<stk::math::Vector3d> & elemNodeCoords)
{
  for (size_t n=0; n<elemNodes.size(); ++n)
  {
    if (elemNodes[n] == node)
    {
      elemNodeCoords[n] = newNodeLoc;
      return;
    }
  }
  STK_ThrowRequireMsg(false, "Did not find the expected node in replace_coordinates_of_node_with_new_location");
}

static double compute_relative_volume_change(const stk::mesh::BulkData & mesh,
    const FieldRef coordsField,
    const VolumePreservingSnappingLimiter::ElementToBlockConverter & elementToBlockConverter,
    const stk::mesh::Part & blockPart,
    const stk::mesh::Entity node,
    const stk::math::Vector3d & newNodeLoc,
    std::vector<stk::math::Vector3d> & elemNodeCoordsWorkspace)
{
  const int dim = mesh.mesh_meta_data().spatial_dimension();

  StkMeshEntities nodeElements{mesh.begin_elements(node), mesh.end_elements(node)};

  double volumeBefore = 0.;
  double volumeAfter = 0.;
  unsigned numElements = 0.;
  for (auto elem : nodeElements)
  {
    if (elementToBlockConverter(mesh, elem) == &blockPart)
    {
      ++numElements;
      StkMeshEntities elemNodes{mesh.begin_nodes(elem), mesh.end_nodes(elem)};
      fill_node_locations(dim, coordsField, elemNodes, elemNodeCoordsWorkspace);

      volumeBefore += compute_tri_or_tet_volume(elemNodeCoordsWorkspace);
      replace_coordinates_of_node_with_new_location(elemNodes, node, newNodeLoc, elemNodeCoordsWorkspace);
      volumeAfter += compute_tri_or_tet_volume(elemNodeCoordsWorkspace);
    }
  }

  if (0 == numElements)
    return 0.;

  const double elemAverageVol = std::max(volumeBefore,volumeAfter)/numElements;
  return std::abs(volumeAfter-volumeBefore)/elemAverageVol;
}

VolumePreservingSnappingLimiter::VolumePreservingSnappingLimiter(
    const stk::mesh::BulkData & mesh,
    const FieldRef coordsField,
    const ElementToBlockConverter & elementToBlockConverter,
    const double volumeConservationTol)
  : myMesh(mesh),
    myAuxMeta(AuxMetaData::get(mesh.mesh_meta_data())),
    myElementToBlockConverter(elementToBlockConverter),
    myCoordsField(coordsField),
    myVolumeConservationTol(volumeConservationTol)
{
}

std::set<stk::mesh::Part*> VolumePreservingSnappingLimiter::get_blocks_to_consider(const stk::mesh::Entity node) const
{
  std::set<stk::mesh::Part*> blocksToConsider;
  for (auto && elem : StkMeshEntities{myMesh.begin_elements(node), myMesh.end_elements(node)})
  {
    stk::mesh::Part * blockPart = myElementToBlockConverter(myMesh, elem);
    if (nullptr != blockPart)
      blocksToConsider.insert(blockPart);
  }
  return blocksToConsider;
}

bool VolumePreservingSnappingLimiter::is_snap_allowed(const stk::mesh::Entity node, const stk::math::Vector3d & snapLocation) const
{
  const std::set<stk::mesh::Part*> blocksToConsider = get_blocks_to_consider(node);
  if (blocksToConsider.size() == 1 && !myMesh.bucket(node).member(myAuxMeta.exposed_boundary_part()))
    return true;

  std::vector<stk::math::Vector3d> elemNodeCoords;
  for (auto && blockPart : blocksToConsider)
  {
    const double volChange = compute_relative_volume_change(myMesh, myCoordsField,  myElementToBlockConverter, *blockPart, node, snapLocation, elemNodeCoords);
    if (volChange > myVolumeConservationTol)
      return false;
  }
  return true;
}

}
