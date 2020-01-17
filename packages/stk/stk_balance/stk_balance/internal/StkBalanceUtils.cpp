#include "StkBalanceUtils.hpp"
#include <stk_balance/balanceUtils.hpp>
#include <stk_io/StkIoUtils.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include "stk_mesh/base/SkinMeshUtil.hpp"
#include "stk_mesh/base/FEMHelpers.hpp"

#include <stk_search/SearchMethod.hpp>
#include <stk_search/CoarseSearch.hpp>

#include <iostream>
#include <iomanip>


namespace stk { namespace balance { namespace internal {

int getNumSharedNodesBetweenElements(const ::stk::mesh::BulkData& stkMeshBulkData,
                                     const ::stk::mesh::Entity element1,
                                     const ::stk::mesh::Entity element2)
{
    const stk::mesh::Entity* nodes1 = stkMeshBulkData.begin_nodes(element1);
    int numNodes1 = stkMeshBulkData.num_nodes(element1);
    const stk::mesh::Entity* nodes2 = stkMeshBulkData.begin_nodes(element2);
    int numNodes2 = stkMeshBulkData.num_nodes(element2);

    int numShared = 0;
    for(int i1=0; i1<numNodes1; ++i1) {
        for(int i2=0; i2<numNodes2; ++i2) {
            if (nodes1[i1] == nodes2[i2]) {
                ++numShared;
                break;
            }
        }
    }

    return numShared;
}

std::string get_parallel_filename(int subdomainIndex, int numSubdomains, const std::string& baseFilename)
{
    return stk::io::construct_parallel_filename(baseFilename, numSubdomains, subdomainIndex);
}

void addBoxForNodes(stk::mesh::BulkData& stkMeshBulkData,
                    unsigned numNodes,
                    const stk::mesh::Entity* nodes,
                    const stk::mesh::FieldBase* coordField,
                    const double eps,
                    stk::mesh::EntityId elementId,
                    stk::balance::internal::SearchBoxIdentProcs& faceBoxes)
{
    unsigned dim = stkMeshBulkData.mesh_meta_data().spatial_dimension();

    std::vector<double> x(numNodes, 0);
    std::vector<double> y(numNodes, 0);
    std::vector<double> z(numNodes, 0);

    for (unsigned nodeIndex = 0; nodeIndex < numNodes; nodeIndex++) {
        double* xyz = static_cast<double*>(stk::mesh::field_data(*coordField, nodes[nodeIndex]));
        x[nodeIndex] = xyz[0];
        y[nodeIndex] = xyz[1];
        if (dim == 3) z[nodeIndex] = xyz[2];
    }

    double maxX = *std::max_element(x.begin(), x.end());
    double maxY = *std::max_element(y.begin(), y.end());
    double maxZ = *std::max_element(z.begin(), z.end());

    double minX = *std::min_element(x.begin(), x.end());
    double minY = *std::min_element(y.begin(), y.end());
    double minZ = *std::min_element(z.begin(), z.end());

    stk::balance::internal::StkBox faceBox(minX - eps, minY - eps, minZ - eps, maxX + eps, maxY + eps, maxZ + eps);
    stk::balance::internal::SearchIdentProc id(elementId, stkMeshBulkData.parallel_rank());
    faceBoxes.emplace_back(faceBox, id);
}

void fillFaceBoxesWithIds(stk::mesh::BulkData &stkMeshBulkData, const BalanceSettings & balanceSettings, const stk::mesh::FieldBase* coord, stk::balance::internal::SearchBoxIdentProcs &faceBoxes, const stk::mesh::Selector& searchSelector)
{
    std::vector<stk::mesh::SideSetEntry> skinnedSideSet = stk::mesh::SkinMeshUtil::get_skinned_sideset_excluding_region(stkMeshBulkData, searchSelector, !searchSelector);
    stk::mesh::EntityVector sideNodes;
    for (stk::mesh::SideSetEntry sidesetEntry : skinnedSideSet)
    {
        stk::mesh::Entity sidesetElement = sidesetEntry.element;
        stk::mesh::ConnectivityOrdinal sidesetSide = sidesetEntry.side;
        stk::mesh::get_subcell_nodes(stkMeshBulkData, sidesetElement,
                                     stkMeshBulkData.mesh_meta_data().side_rank(),
                                     sidesetSide, sideNodes);
        const double eps = balanceSettings.getToleranceForFaceSearch(stkMeshBulkData, *coord,
                                                                     sideNodes.data(), sideNodes.size());
        addBoxForNodes(stkMeshBulkData, sideNodes.size(), &sideNodes[0], coord, eps,
                       stkMeshBulkData.identifier(sidesetElement), faceBoxes);
    }
}

void fillParticleBoxesWithIds(stk::mesh::BulkData &stkMeshBulkData,
                              const BalanceSettings & balanceSettings,
                              const stk::mesh::FieldBase* coord,
                              stk::balance::internal::SearchBoxIdentProcs &boxes)
{
    const stk::mesh::BucketVector &elementBuckets = stkMeshBulkData.buckets(stk::topology::ELEMENT_RANK);

    for (size_t i=0;i<elementBuckets.size();i++)
    {
        stk::mesh::Bucket &bucket = *elementBuckets[i];
        if( bucket.owned() && (bucket.topology() == stk::topology::PARTICLE ) )
        {
            for(size_t j = 0; j < bucket.size(); j++)
            {
                const stk::mesh::Entity *node = stkMeshBulkData.begin_nodes(bucket[j]);
                double *xyz = static_cast<double *>(stk::mesh::field_data(*coord, *node));
                double eps = balanceSettings.getAbsoluteToleranceForParticleSearch(bucket[j]);

                stk::balance::internal::StkBox box(xyz[0] - eps, xyz[1] - eps, xyz[2] - eps, xyz[0] + eps, xyz[1] + eps, xyz[2] + eps);

                unsigned int val1 = stkMeshBulkData.identifier(bucket[j]);
                int val2 = stkMeshBulkData.parallel_rank();
                stk::balance::internal::SearchIdentProc id(val1, val2);

                boxes.emplace_back(box, id);
            }
        }
    }
}

SearchElemPairs getBBIntersectionsForFacesParticles(stk::mesh::BulkData& stkMeshBulkData,
                                                    const BalanceSettings &balanceSettings,
                                                    const stk::mesh::Selector& searchSelector)
{
    bool useLocalIds = balanceSettings.getGraphOption() == BalanceSettings::COLOR_MESH ||
                       balanceSettings.getGraphOption() == BalanceSettings::COLOR_MESH_BY_TOPOLOGY ||
                       balanceSettings.getGraphOption() == BalanceSettings::COLOR_MESH_AND_OUTPUT_COLOR_FIELDS;
    ThrowRequireWithSierraHelpMsg(useLocalIds != true);

    const stk::mesh::FieldBase * coord = stkMeshBulkData.mesh_meta_data().coordinate_field();

    stk::balance::internal::SearchBoxIdentProcs faceBoxes;
    fillFaceBoxesWithIds(stkMeshBulkData, balanceSettings, coord, faceBoxes, searchSelector);

    if ( balanceSettings.getEdgesForParticlesUsingSearch() )
    {
        fillParticleBoxesWithIds(stkMeshBulkData, balanceSettings, coord, faceBoxes);
    }

    stk::balance::internal::SearchElemPairs searchResults;

    stk::search::coarse_search(faceBoxes, faceBoxes, stk::search::KDTREE, stkMeshBulkData.parallel(), searchResults);

    stk::balance::internal::SearchElemPairs::iterator iter = std::unique(searchResults.begin(), searchResults.end());
    searchResults.resize(iter - searchResults.begin());
    return searchResults;
}

}}}
