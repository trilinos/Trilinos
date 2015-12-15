#ifndef BULK_DATA_ID_MAPPER_HPP
#define BULK_DATA_ID_MAPPER_HPP

#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/ElemElemGraphImpl.hpp>
#include <stk_mesh/base/ElemGraphCoincidentElems.hpp>

namespace stk
{
namespace mesh
{
namespace impl
{

class BulkDataIdMapper : public IdMapper
{
public:
    BulkDataIdMapper(const stk::mesh::BulkData &b,
                     const stk::mesh::EntityVector &e,
                     const std::vector<stk::mesh::impl::LocalId> &localIds) :
            bulk(b),
            elements(e),
            elemToLocalIds(localIds)
    {
    }
    virtual stk::mesh::EntityId localToGlobal(stk::mesh::impl::LocalId local) const
    {
        return bulk.identifier(elements[local]);
    }
    virtual stk::mesh::impl::LocalId globalToLocal(stk::mesh::EntityId global) const
    {
        stk::mesh::Entity elem = bulk.get_entity(stk::mesh::EntityKey(stk::topology::ELEM_RANK, global));
        return elemToLocalIds[elem.local_offset()];
    }
private:
    const stk::mesh::BulkData &bulk;
    const stk::mesh::EntityVector &elements;
    const std::vector<stk::mesh::impl::LocalId> &elemToLocalIds;
};

}
}
}

#endif
