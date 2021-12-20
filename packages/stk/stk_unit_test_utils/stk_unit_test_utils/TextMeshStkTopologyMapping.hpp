#ifndef STK_STK_UNIT_TEST_UTILS_STK_UNIT_TEST_UTILS_TEXTMESHSTKTOPOLOGYMAPPING_HPP_
#define STK_STK_UNIT_TEST_UTILS_STK_UNIT_TEST_UTILS_TEXTMESHSTKTOPOLOGYMAPPING_HPP_

#include "TextMeshUtils.hpp"

#include "stk_mesh/base/Types.hpp"                   // for EntityId, etc
#include "stk_topology/topology.hpp"                 // for topology, etc

#include <string>                                    // for basic_string, etc
#include <utility>                                   // for pair
#include <vector>                                    // for vector

namespace stk
{
namespace unit_test_util
{
class StkTopologyMapping : public text_mesh::TopologyMapping<stk::topology>
{
 public:
  stk::topology invalid_topology() const override { return stk::topology::INVALID_TOPOLOGY; }

  void initialize_topology_map() override
  {
    m_nameToTopology = {
        {"NODE", stk::topology::NODE},
        {"LINE_2", stk::topology::LINE_2},
        {"LINE_3", stk::topology::LINE_3},
        {"TRI_3", stk::topology::TRI_3},
        {"TRI_4", stk::topology::TRI_4},
        {"TRI_6", stk::topology::TRI_6},
        {"QUAD_4", stk::topology::QUAD_4},
        {"QUAD_6", stk::topology::QUAD_6},
        {"QUAD_8", stk::topology::QUAD_8},
        {"QUAD_9", stk::topology::QUAD_9},
        {"PARTICLE", stk::topology::PARTICLE},
        {"LINE_2_1D", stk::topology::LINE_2_1D},
        {"LINE_3_1D", stk::topology::LINE_3_1D},
        {"BEAM_2", stk::topology::BEAM_2},
        {"BEAM_3", stk::topology::BEAM_3},
        {"SHELL_LINE_2", stk::topology::SHELL_LINE_2},
        {"SHELL_LINE_3", stk::topology::SHELL_LINE_3},
        {"SPRING_2", stk::topology::SPRING_2},
        {"SPRING_3", stk::topology::SPRING_3},
        {"TRI_3_2D", stk::topology::TRI_3_2D},
        {"TRI_4_2D", stk::topology::TRI_4_2D},
        {"TRI_6_2D", stk::topology::TRI_6_2D},
        {"QUAD_4_2D", stk::topology::QUAD_4_2D},
        {"QUAD_8_2D", stk::topology::QUAD_8_2D},
        {"QUAD_9_2D", stk::topology::QUAD_9_2D},
        {"SHELL_TRI_3", stk::topology::SHELL_TRI_3},
        {"SHELL_TRI_4", stk::topology::SHELL_TRI_4},
        {"SHELL_TRI_6", stk::topology::SHELL_TRI_6},
        {"SHELL_QUAD_4", stk::topology::SHELL_QUAD_4},
        {"SHELL_QUAD_8", stk::topology::SHELL_QUAD_8},
        {"SHELL_QUAD_9", stk::topology::SHELL_QUAD_9},
        {"TET_4", stk::topology::TET_4},
        {"TET_8", stk::topology::TET_8},
        {"TET_10", stk::topology::TET_10},
        {"TET_11", stk::topology::TET_11},
        {"PYRAMID_5", stk::topology::PYRAMID_5},
        {"PYRAMID_13", stk::topology::PYRAMID_13},
        {"PYRAMID_14", stk::topology::PYRAMID_14},
        {"WEDGE_6", stk::topology::WEDGE_6},
        {"WEDGE_12", stk::topology::WEDGE_12},
        {"WEDGE_15", stk::topology::WEDGE_15},
        {"WEDGE_18", stk::topology::WEDGE_18},
        {"HEX_8", stk::topology::HEX_8},
        {"HEX_20", stk::topology::HEX_20},
        {"HEX_27", stk::topology::HEX_27},
    };
  }
};

} // namespace unit_test_util
} // namespace stk

#endif /* STK_STK_UNIT_TEST_UTILS_STK_UNIT_TEST_UTILS_TEXTMESHSTKTOPOLOGYMAPPING_HPP_ */
