#include <stk_mesh/base/Iterators.hpp>

namespace stk {
namespace mesh {

BucketVectorEntityIteratorRange get_entity_range(const std::vector<Bucket*>& buckets)
{
  if ( buckets.empty() ) {
    BucketVectorEntityIterator end_itr(buckets.end());
    return std::make_pair(end_itr, end_itr);
  }
  else {
    return std::make_pair(BucketVectorEntityIterator(buckets.begin(), BucketIterator((*buckets.begin())->begin()),   buckets.end()),
                          BucketVectorEntityIterator(buckets.end()));
  }
}

SelectedBucketVectorEntityIteratorRange get_entity_range(const std::vector<Bucket*>& buckets, const Selector& selector)
{
  if ( buckets.empty() ) {
    SelectedBucketVectorIterator end_select_itr(buckets.end());
    //SelectedBucketVectorIterator end_select_itr(selector, buckets.end(), buckets.end());
    SelectedBucketVectorEntityIterator end_itr(end_select_itr);
    return std::make_pair(end_itr, end_itr);
  }
  else {
    //SelectedBucketVectorIterator end_select_itr(selector, buckets.end(), buckets.end());
    SelectedBucketVectorIterator end_select_itr(buckets.end());
    return std::make_pair(SelectedBucketVectorEntityIterator(SelectedBucketVectorIterator(selector, buckets.begin(), buckets.end()),
                                                             BucketIterator((*buckets.begin())->begin()),
                                                             end_select_itr),
                          SelectedBucketVectorEntityIterator(end_select_itr));
  }
}

SelectedBucketVectorIteratorRange get_selected_bucket_range(const std::vector<Bucket*>& buckets, const Selector& selector)
{
  return std::make_pair(SelectedBucketVectorIterator(selector, buckets.begin(), buckets.end()),
                        SelectedBucketVectorIterator(buckets.end()));
}

} //namespace mesh
} //namespace stk
