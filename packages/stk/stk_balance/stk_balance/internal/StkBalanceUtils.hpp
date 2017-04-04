#ifndef STK_BALANCE_UTILITIES_HPP
#define STK_BALANCE_UTILITIES_HPP

#include <string>
#include <vector>

// stk search (order matters!)
#include <stk_search/IdentProc.hpp>
#include <stk_search/BoundingBox.hpp>
#include <stk_mesh/base/Types.hpp>

namespace stk { namespace mesh { class BulkData; } }
namespace stk { namespace mesh { class MetaData; } }
namespace stk { namespace mesh { class Selector; } }
namespace stk { namespace mesh { class Entity; } }
namespace stk { namespace mesh { class FieldBase; } }

namespace stk { namespace balance { class BalanceSettings; } }

namespace stk { namespace balance { namespace internal {

typedef stk::search::Box<float> StkBox;
typedef stk::search::IdentProc<stk::mesh::EntityId, int> StkMeshIdent;
typedef std::pair<StkBox, StkMeshIdent> BoxWithStkId;
typedef std::vector< BoxWithStkId > BoxVectorWithStkId;
typedef std::vector<std::pair<StkMeshIdent, StkMeshIdent> > StkSearchResults;

std::string get_parallel_filename(int subdomainIndex, int numSubdomains, const std::string& baseFilename);

int getNumSharedNodesBetweenElements(const ::stk::mesh::BulkData& stkMeshBulkData,
                                     const ::stk::mesh::Entity element1,
                                     const ::stk::mesh::Entity element2);

StkSearchResults getSearchResultsForFacesParticles(stk::mesh::BulkData& stkMeshBulkData, const BalanceSettings &balanceSettings, const stk::mesh::Selector& searchSelector);

void addBoxForNodes(stk::mesh::BulkData& stkMeshBulkData,
                    unsigned numNodes,
                    const stk::mesh::Entity* nodes,
                    const stk::mesh::FieldBase* coord,
                    const double eps,
                    stk::mesh::EntityId elementId,
                    BoxVectorWithStkId& faceBoxes);


void fillFaceBoxesWithIds(stk::mesh::BulkData &stkMeshBulkData, const double eps, const stk::mesh::FieldBase* coord, BoxVectorWithStkId &faceBoxes, const stk::mesh::Selector& searchSelector);

const stk::mesh::FieldBase * get_coordinate_field(const stk::mesh::MetaData& meta_data, const std::string& coordinateFieldName);

}}}

#endif
