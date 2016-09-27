#ifndef GEOMETRICVERTICES_HPP_
#define GEOMETRICVERTICES_HPP_

#include "Vertices.hpp"
#include <vector>

namespace stk {namespace mesh {class BulkData;}}
namespace stk {namespace mesh {class Selector;}}
namespace stk { namespace balance {class BalanceSettings;}}

namespace stk {
namespace balance {
namespace internal {

class GeometricVertices : public stk::balance::internal::Vertices
{
public:
    GeometricVertices(const stk::balance::BalanceSettings& balanceSettings,
                      const stk::mesh::BulkData& bulkData,
                      const stk::mesh::EntityVector& entities_to_balance,
                      const std::vector<stk::mesh::Selector> criteria);

    virtual ~GeometricVertices();

    size_t getNumWeightsPerVertex() const;

private:
    const stk::mesh::EntityVector& mEntitiesToBalance;
    std::vector<stk::mesh::Selector> mSelectors;
};

}
}
}
#endif /* GEOMETRICVERTICES_HPP_ */
