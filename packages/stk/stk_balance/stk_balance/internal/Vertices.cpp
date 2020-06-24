#include "Vertices.hpp"
#include "privateDeclarations.hpp"
#include "balanceUtils.hpp"

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>

namespace stk {
namespace balance {
namespace internal {

void Vertices::fillVertexIds(const stk::mesh::BulkData& bulkData,
                             const stk::mesh::EntityVector &entities)
{
    mVertexIds.resize(entities.size(), 0);
    for(size_t i=0;i<entities.size();++i) {
        mVertexIds[i] = bulkData.identifier(entities[i]);
    }
}

void Vertices::fillCoordinates(const stk::mesh::BulkData& bulkData,
                               const std::string& coords_field_name,
                               const stk::mesh::EntityVector &entities)
{
    mVertexCoordinates.resize(entities.size()*mSpatialDim, 0);
    const stk::mesh::FieldBase * coord = bulkData.mesh_meta_data().get_field(stk::topology::NODE_RANK, coords_field_name);

    for(size_t i=0;i<entities.size();++i)
        stk::balance::internal::fillEntityCentroid(bulkData, coord, entities[i], &mVertexCoordinates[mSpatialDim*i]);
}

void Vertices::fillVertexWeights(const stk::mesh::BulkData& bulkData,
                                 const stk::balance::BalanceSettings& balanceSettings,
                                 const stk::mesh::EntityVector &entities,
                                 const std::vector<stk::mesh::Selector> &selectors)
{
    if(balanceSettings.areVertexWeightsProvidedViaFields())
    {
        fillFieldVertexWeights(balanceSettings, bulkData, selectors, entities);
    }
    else
    {
        mVertexWeights.resize(entities.size(), 0.0);
        for(size_t i=0;i<entities.size();++i)
            mVertexWeights[i] = balanceSettings.getGraphVertexWeight(bulkData.bucket(entities[i]).topology());
    }
}

void Vertices::fillFieldVertexWeights(const stk::balance::BalanceSettings& balanceSettings,
                            const stk::mesh::BulkData& stkMeshBulkData,
                            const std::vector<stk::mesh::Selector>& selectors,
                            const stk::mesh::EntityVector &entitiesToBalance)
{
    unsigned numSelectors = selectors.size();
    unsigned numEntities = entitiesToBalance.size();
    unsigned numCriteria = balanceSettings.getNumCriteria();

    mVertexWeights.resize(numSelectors*numEntities*numCriteria);
    for (double &weight : mVertexWeights) { weight = 0.0; }

    for(size_t i=0;i<entitiesToBalance.size();++i)
    {
        for (unsigned sel=0 ; sel<numSelectors ; ++sel)
        {
            stk::mesh::Selector selector = selectors[sel];
            if (selector(stkMeshBulkData.bucket(entitiesToBalance[i])))
            {
                for(size_t weight_index=0;weight_index<numCriteria;weight_index++)
                {
                    unsigned index = stk::balance::internal::get_index(numSelectors, numCriteria, i, sel, weight_index);
                    mVertexWeights[index] = balanceSettings.getGraphVertexWeight(entitiesToBalance[i], weight_index);
                }
            }
        }
    }
}

}
}
}
