#ifndef PACKAGES_STK_STK_MESH_STK_MESH_BASEIMPL_ELEMENTGRAPH_GRAPHTYPES_HPP_
#define PACKAGES_STK_STK_MESH_STK_MESH_BASEIMPL_ELEMENTGRAPH_GRAPHTYPES_HPP_

namespace stk
{
namespace mesh
{
namespace impl
{

typedef int64_t LocalId;

typedef std::vector<stk::mesh::PartOrdinal> PartOrdinals;
typedef std::map<LocalId, PartOrdinals> ParallelPartInfo;
typedef std::map<LocalId, bool> ParallelSelectedInfo;

}
}
}

#endif /* PACKAGES_STK_STK_MESH_STK_MESH_BASEIMPL_ELEMENTGRAPH_GRAPHTYPES_HPP_ */
