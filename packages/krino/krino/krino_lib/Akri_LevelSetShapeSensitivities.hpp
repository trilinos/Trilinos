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

std::string output_sensitivity(const LevelSetShapeSensitivity & sens);

std::vector<LevelSetShapeSensitivity> get_levelset_shape_sensitivities(const stk::mesh::BulkData & mesh, const std::vector<LS_Field> & lsFields);

}



#endif /* KRINO_KRINO_KRINO_LIB_AKRI_LEVELSETSHAPESENSITIVITIES_HPP_ */
