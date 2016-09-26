#ifndef EXODUSTRANSLATOR_HPP_
#define EXODUSTRANSLATOR_HPP_

#include <stk_topology/topology.hpp>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/CoordinateSystems.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/TopologyDimensions.hpp>
#include <stk_mesh/base/FindRestriction.hpp>

namespace stk
{
namespace mesh
{

class ExodusTranslator
{
public:
    ExodusTranslator(const stk::mesh::BulkData& b) : mBulkData(b)
    {

    }

    size_t get_number_global_entities(stk::mesh::EntityRank rank) const
    {
        std::vector<size_t> counts;
        stk::mesh::comm_mesh_counts(mBulkData, counts);
        return counts[rank];
    }

    size_t get_number_global_elements() const
    {
        return get_number_global_entities(stk::topology::ELEMENT_RANK);
    }

private:
    const stk::mesh::BulkData& mBulkData;
};

}
}

#endif /* EXODUSTRANSLATOR_HPP_ */
