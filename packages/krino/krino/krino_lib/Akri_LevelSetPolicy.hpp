#ifndef KRINO_KRINO_KRINO_LIB_AKRI_LEVELSETPOLICY_HPP_
#define KRINO_KRINO_KRINO_LIB_AKRI_LEVELSETPOLICY_HPP_

#include <stk_mesh/base/MetaData.hpp>
#include <Akri_Phase_Support.hpp>

namespace krino {


struct LSPerPhasePolicy
{
  static std::vector<LS_Field> setup_levelsets_on_blocks(stk::mesh::MetaData & meta,
      const unsigned numLevelSets,
      const stk::mesh::PartVector & blocks,
      Block_Surface_Connectivity & blockSurfaceInfo,
      const bool registerFields,
      CDFEM_Inequality_Spec* deathSpec = nullptr);
  static std::vector<LS_Field> setup_levelsets_on_all_blocks(stk::mesh::MetaData & meta, const unsigned numLevelSets);
};

struct LSPerInterfacePolicy
{
  static std::vector<LS_Field> setup_levelsets_on_blocks(stk::mesh::MetaData & meta,
      const unsigned numLevelSets,
      const stk::mesh::PartVector & blocks,
      Block_Surface_Connectivity & blockSurfaceInfo,
      const bool registerFields,
      CDFEM_Inequality_Spec* deathSpec = nullptr);
  static std::vector<LS_Field> setup_levelsets_on_all_blocks(stk::mesh::MetaData & meta, const unsigned numLevelSets);
  static std::vector<LS_Field> setup_levelsets_on_all_blocks_with_void_phase_for_any_negative_levelset(stk::mesh::MetaData & meta, const unsigned numLevelSets);
};

}

#endif /* KRINO_KRINO_KRINO_LIB_AKRI_LEVELSETPOLICY_HPP_ */
