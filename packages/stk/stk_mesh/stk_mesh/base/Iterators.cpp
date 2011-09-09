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
    return std::make_pair(BucketVectorEntityIterator(buckets.begin(), BucketPtrIterator((*buckets.begin())->begin()),   buckets.end()),
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
                                                             BucketPtrIterator((*buckets.begin())->begin()),
                                                             end_select_itr),
                          SelectedBucketVectorEntityIterator(end_select_itr));
  }
}

AllBucketsRange get_bucket_range(const std::vector<std::vector<Bucket*> >& buckets)
{
  if (buckets.empty()) {
    AllBucketsIterator end_itr(buckets.end());
    return std::make_pair(end_itr, end_itr);
  }
  else {
    return std::make_pair(AllBucketsIterator(buckets.begin(), buckets.begin()->begin(), buckets.end()),
                          AllBucketsIterator(buckets.end()));
  }
}

AllBucketsRange get_bucket_range(const std::vector<std::vector<Bucket*> >& buckets,
                                 std::vector<std::vector<Bucket*> >::const_iterator itr)
{
  AllBucketsIterator end_itr(itr+1);
  return std::make_pair(AllBucketsIterator(itr, itr->begin(), itr+1),
                        end_itr);
}

AllSelectedBucketsRange get_selected_bucket_range(const AllBucketsRange& bucket_range, const Selector& selector)
{
  return std::make_pair(AllSelectedBucketsIterator(selector, bucket_range.first, bucket_range.second),
                        AllSelectedBucketsIterator(selector, bucket_range.second, bucket_range.second));
}

SelectedBucketRangeEntityIteratorRange get_selected_bucket_entity_range(const AllBucketsRange& bucket_range, const Selector& selector)
{
  if (bucket_range.first == bucket_range.second) {
    AllSelectedBucketsIterator end_select_itr(bucket_range.second);
    //AllSelectedBucketsIterator end_select_itr(selector, bucket_range.second, bucket_range.second);
    SelectedBucketRangeEntityIterator end_itr(end_select_itr);
    return std::make_pair(end_itr, end_itr);
  }
  else {
    AllSelectedBucketsRange selected_bucket_range = get_selected_bucket_range(bucket_range, selector);
    return std::make_pair(SelectedBucketRangeEntityIterator(selected_bucket_range.first,
                                                            BucketPtrIterator( (*selected_bucket_range.first)->begin()),
                                                            selected_bucket_range.second),
                          SelectedBucketRangeEntityIterator(selected_bucket_range.second));
  }
}

} //namespace mesh
} //namespace stk
