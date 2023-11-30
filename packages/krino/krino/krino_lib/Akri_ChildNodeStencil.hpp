#ifndef KRINO_KRINO_KRINO_LIB_AKRI_CHILDNODESTENCIL_HPP_
#define KRINO_KRINO_KRINO_LIB_AKRI_CHILDNODESTENCIL_HPP_

#include <vector>

#include <stk_mesh/base/Entity.hpp>

namespace stk { namespace mesh { class BulkData; } }
namespace stk { namespace mesh { class Part; } }
namespace krino { class FieldRef; }

namespace krino {

struct ChildNodeStencil
{
  ChildNodeStencil(const stk::mesh::Entity child, const std::vector<stk::mesh::Entity> & parents, const std::vector<double> & parentWts)
  : childNode(child), parentNodes(parents), parentWeights(parentWts) {}
  stk::mesh::Entity childNode;
  std::vector<stk::mesh::Entity> parentNodes;
  std::vector<double> parentWeights;
};

void fill_child_node_stencils(const stk::mesh::BulkData & mesh,
    const stk::mesh::Part & childNodePart,
    const FieldRef & parentIdsField,
    const FieldRef & parentWtsField,
    std::vector<ChildNodeStencil> & childNodeStencils);

}



#endif /* KRINO_KRINO_KRINO_LIB_AKRI_CHILDNODESTENCIL_HPP_ */
