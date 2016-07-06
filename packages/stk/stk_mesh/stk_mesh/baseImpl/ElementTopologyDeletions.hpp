#include "stk_mesh/base/BulkData.hpp"
#include "stk_mesh/base/Entity.hpp"
#include "stk_mesh/base/MetaData.hpp"
#include "stk_mesh/base/Selector.hpp"
#include "stk_mesh/base/Types.hpp"

namespace stk
{
namespace mesh
{
namespace impl
{

struct RelationEntityToNode
{
    Entity entity;
    Entity node;
    ConnectivityOrdinal ordinal;
};

class ElementTopologyDeletions
{
public:
    ElementTopologyDeletions(stk::mesh::BulkData &bulkData, stk::topology topologyToDelete);
    void find_buckets_to_destroy_and_relations_to_sever();
    const std::vector<RelationEntityToNode> & get_relations_to_sever() const;
    const stk::mesh::BucketVector & get_buckets_to_delete() const;

private:
    stk::mesh::Selector get_all_other_element_topologies_selector(const stk::mesh::MetaData &meta, stk::topology topology);
    void append_upward_relations_to_sever(stk::mesh::Bucket &bucket);
private:
    stk::mesh::BulkData & bulk;
    stk::mesh::Selector topoSelector;
    stk::mesh::Selector selectsOtherElemTopologies;
    stk::mesh::BucketVector bucketsToDelete;
    std::vector<RelationEntityToNode> relationsToDestroy;
};

}
}
}
