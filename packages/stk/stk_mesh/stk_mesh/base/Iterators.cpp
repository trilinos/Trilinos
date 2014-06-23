#include <stk_mesh/base/Iterators.hpp>
#include "stk_mesh/base/Selector.hpp"   // for Selector
namespace stk { namespace mesh { class Bucket; } }

namespace stk {
namespace mesh {

SelectedBucketVectorIteratorRange get_selected_bucket_range(const BucketVector& buckets, const Selector& selector)
{
  return std::make_pair(SelectedBucketVectorIterator(selector, buckets.begin(), buckets.end()),
                        SelectedBucketVectorIterator(buckets.end()));
}

} //namespace mesh
} //namespace stk
