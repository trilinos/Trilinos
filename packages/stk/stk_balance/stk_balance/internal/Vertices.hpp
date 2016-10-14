#ifndef BALANCE_VERTEX_INFO_HPP
#define BALANCE_VERTEX_INFO_HPP

#include <stk_mesh/base/Types.hpp>
#include <string>
#include <vector>
#include "balanceTypes.hpp"
#include <stk_util/environment/ReportHandler.hpp>  // for ThrowRequireMsg

namespace stk {namespace mesh {class BulkData;}}
namespace stk {namespace mesh {class Selector;}}
namespace stk { namespace balance {class BalanceSettings;}}

namespace stk {
namespace balance {
namespace internal {

class Vertices
{
public:
    virtual ~Vertices() {}

    size_t num_vertices() const { return mVertexIds.size(); }
    const std::vector<BalanceGlobalNumber> & get_vertex_ids() const { return mVertexIds; }
    const std::vector<double> & get_vertex_coords() const { return mVertexCoordinates; }
    const std::vector<double> & get_vertex_weights() const { return mVertexWeights; }

    void set_vertex_weight( size_t idx, double weight)
    {
        ThrowRequireMsg(mNumFieldCriteria==1 && idx < mVertexWeights.size(), "invalid index for " << __PRETTY_FUNCTION__);
        mVertexWeights[idx] = weight;
    }

    void set_vertex_weights( const std::vector<double> &weights )
    {
        mVertexWeights = weights;
    }

    unsigned get_spatial_dim() const { return mSpatialDim; }
    void set_spatial_dim(unsigned dim) { mSpatialDim = dim; }

    unsigned get_num_field_criteria() const { return mNumFieldCriteria; }
    void set_num_field_criteria(unsigned num) {mNumFieldCriteria = num; }

protected:
    void fillVertexIds(const stk::mesh::BulkData& bulkData, const stk::mesh::EntityVector &entities);

    void fillCoordinates(const stk::mesh::BulkData& bulkData, const std::string& coords_field_name, const stk::mesh::EntityVector &entities);

    void fillVertexWeights(const stk::mesh::BulkData& bulkData, const stk::balance::BalanceSettings& balanceSettings, const stk::mesh::EntityVector &entities, const std::vector<stk::mesh::Selector> &selectors);

    void fillFieldVertexWeights(const stk::balance::BalanceSettings& balanceSettings,
                                const stk::mesh::BulkData& stkMeshBulkData,
                                const std::vector<stk::mesh::Selector>& selectors,
                                const stk::mesh::EntityVector &entitiesToBalance);

    std::vector<double> mVertexCoordinates;
    std::vector<BalanceGlobalNumber> mVertexIds;
    std::vector<double> mVertexWeights;
    unsigned mSpatialDim = 0;
    unsigned mNumFieldCriteria = 1;
};

}
}
}
#endif
