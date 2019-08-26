#include "stk_tools/mesh_tools/BucketStatistics.hpp"
#include "stk_mesh/base/BulkData.hpp"
#include "stk_mesh/base/MetaData.hpp"
#include "stk_topology/topology.hpp"
#include <cmath>
#include <map>

namespace {
struct BucketStatistics
{
  unsigned numBuckets;
  unsigned totalNumElems;
};
}

namespace stk
{
namespace tools
{

void dump_bucket_statistics(const stk::mesh::BulkData& bulk, std::ostream& os)
{
  std::map< stk::topology, BucketStatistics > topologyToBucketMap;
  std::set< stk::topology > topologies;

  const stk::mesh::BucketVector& buckets = bulk.get_buckets(stk::topology::ELEMENT_RANK, bulk.mesh_meta_data().locally_owned_part());
  for (const stk::mesh::Bucket* bucket : buckets)
  {
    topologies.insert(bucket->topology());
  }

  for (stk::topology topology : topologies)
  {
    BucketStatistics stats{0, 0};
    topologyToBucketMap[topology] = stats;
    for (const stk::mesh::Bucket* bucket : buckets)
    {
      if (bucket->topology() == topology)
      {
        topologyToBucketMap[topology].numBuckets += 1;
        topologyToBucketMap[topology].totalNumElems += bucket->size();
      }
    }
  }

  for (auto& mapIt : topologyToBucketMap)
  {
    int avgElems = std::nearbyint((double)mapIt.second.totalNumElems / mapIt.second.numBuckets);
    os << mapIt.first << ":\t NumBuckets = " << mapIt.second.numBuckets << ", avgNumElems = " << avgElems << std::endl;
  }
}

}
}
