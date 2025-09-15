#ifndef KRINO_KRINO_KRINO_LIB_AKRI_LEVELSETSHAPESENSITIVITIES_HPP_
#define KRINO_KRINO_KRINO_LIB_AKRI_LEVELSETSHAPESENSITIVITIES_HPP_

#include <stk_mesh/base/Types.hpp>
#include <stk_math/StkVector.hpp>

#include <Akri_Phase_Support.hpp>

namespace krino {

struct LevelSetShapeSensitivity
{
  LevelSetShapeSensitivity(const stk::mesh::EntityId inInterfaceNodeId, const std::vector<stk::mesh::EntityId> & inParentNodeIds, const std::vector<stk::math::Vector3d> & inDCoordsdParentLevelSets)
  : interfaceNodeId(inInterfaceNodeId), parentNodeIds(inParentNodeIds), dCoordsdParentLevelSets(inDCoordsdParentLevelSets) {}
  stk::mesh::EntityId interfaceNodeId;
  std::vector<stk::mesh::EntityId> parentNodeIds;
  std::vector<stk::math::Vector3d> dCoordsdParentLevelSets;
};

std::vector<LS_Field> create_levelset_copies_and_set_to_use_as_snap_fields(const stk::mesh::MetaData & meta, const std::vector<LS_Field> & origLSFields);

std::vector<LS_Field> get_levelset_copy_fields(const stk::mesh::MetaData & meta, const std::vector<LS_Field> & origLSFields);

std::vector<LS_Field> update_levelset_copies_to_prepare_for_snapping(const stk::mesh::MetaData & meta, const std::vector<LS_Field> & origLSFields);

std::string output_sensitivity(const LevelSetShapeSensitivity & sens);

std::vector<LevelSetShapeSensitivity> get_levelset_shape_sensitivities(const stk::mesh::BulkData & mesh, const std::vector<LS_Field> & lsFields, const bool doComputeClosestPointSensitivities);

using Facet2dSensitivityIndices = std::array<std::size_t,2>;
using Facet3dSensitivityIndices = std::array<std::size_t,3>;

template <size_t NNODES>
void fill_levelset_facets_and_shape_sensitivities(const stk::mesh::BulkData & mesh,
   const LS_Field & lsField,
   const bool doComputeClosestPointSensitivities,
   std::vector<std::array<size_t,NNODES>> & facetSensitivityIndices,
   std::vector<LevelSetShapeSensitivity> & sortedShapeSensitivities);

template <size_t NNODES>
std::array<const LevelSetShapeSensitivity*,NNODES> get_facet_node_sensitivities(const std::vector<LevelSetShapeSensitivity> & sensitivities, const std::array<size_t,NNODES> & facetSensIndices);

template <size_t NNODES>
std::array<stk::mesh::Entity,NNODES> get_facet_nodes(const stk::mesh::BulkData & mesh, const std::array<const LevelSetShapeSensitivity*,NNODES> & facetNodeSens);

}



#endif /* KRINO_KRINO_KRINO_LIB_AKRI_LEVELSETSHAPESENSITIVITIES_HPP_ */
