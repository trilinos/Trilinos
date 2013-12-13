#include <stk_mesh/base/Iterators.hpp>

namespace stk {
namespace mesh {

SelectedBucketVectorIteratorRange get_selected_bucket_range(const std::vector<Bucket*>& buckets, const Selector& selector)
{
  return std::make_pair(SelectedBucketVectorIterator(selector, buckets.begin(), buckets.end()),
                        SelectedBucketVectorIterator(buckets.end()));
}

} //namespace mesh
} //namespace stk
