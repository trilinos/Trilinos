#include "GeometricVertices.hpp"
#include "stk_balance/balanceUtils.hpp"

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>

namespace stk {
namespace balance {
namespace internal {

GeometricVertices::GeometricVertices(const stk::balance::BalanceSettings& balanceSettings,
                                     const stk::mesh::BulkData& bulkData,
                                     const stk::mesh::EntityVector& entities_to_balance,
                                     const std::vector<stk::mesh::Selector> criteria)
  : mEntitiesToBalance(entities_to_balance),
    mSelectors(criteria)
{
  set_num_field_criteria( balanceSettings.getNumCriteria() );
  set_spatial_dim( bulkData.mesh_meta_data().spatial_dimension() );
  fillVertexIds(bulkData, mEntitiesToBalance);
  fillCoordinates(bulkData, balanceSettings.getCoordinateFieldName(), mEntitiesToBalance);
  fillVertexWeights(bulkData, balanceSettings, mEntitiesToBalance, mSelectors);
}

size_t GeometricVertices::getNumWeightsPerVertex() const
{
  return mSelectors.size()*get_num_field_criteria();
}

GeometricVertices::~GeometricVertices() {}

}
}
}
