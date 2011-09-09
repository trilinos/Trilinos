#ifndef Toolkit_Iterators_hpp
#define Toolkit_Iterators_hpp

#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/base/Selector.hpp>

#include <boost/iterator/filter_iterator.hpp>

#include <vector>
#include <algorithm>

namespace stk {
namespace mesh {

// Requirements:
//   Dereference of HighLevelItrType gives a container (OR pointer to container) with iterators compatible with LowLevelItrType
//   Dereference of LowLevelItrType gives ValueType
//   LowLevelItrType must default-construct to a deterministic invalid value
//
// Incrementing this iterator will take us to the next valid low-level iterator, skipping past
// empty high-level containers, until the end iterator.
template<typename HighLevelItrType, typename LowLevelItrType, typename ValueType>
class TwoLevelIterator : public std::iterator<std::forward_iterator_tag, ValueType>
{
 public:
  typedef TwoLevelIterator<HighLevelItrType, LowLevelItrType, ValueType>  self;

  // Construct an iterator from a starting point specified by high_itr and low_itr
  TwoLevelIterator(HighLevelItrType high_itr, LowLevelItrType low_itr, HighLevelItrType high_end_itr) :
    m_high_itr(high_itr),
    m_low_itr(low_itr),
    m_high_end_itr(high_end_itr)
  {
    if (high_itr != high_end_itr) {
      find_next_valid_item();
    }
  }

  // Construct the "end" iterator
  TwoLevelIterator(HighLevelItrType high_end_itr) :
    m_high_itr(high_end_itr),
    m_low_itr(),
    m_high_end_itr(high_end_itr)
  {}

  TwoLevelIterator(const self& rhs) :
    m_high_itr(rhs.m_high_itr),
    m_low_itr(rhs.m_low_itr),
    m_high_end_itr(rhs.m_high_end_itr)
  {}

  self& operator=(const self& rhs)
  {
    m_high_itr     = rhs.m_high_itr;
    m_low_itr      = rhs.m_low_itr;
    m_high_end_itr = rhs.m_high_end_itr;

    return *this;
  }

  bool operator==(const self& rhs) const
  {
    return (m_high_itr == rhs.m_high_itr && m_low_itr == rhs.m_low_itr);
  }

  bool operator!=(const self& rhs) const
  {
    return !(*this == rhs);
  }

  // ++x
  self& operator++(int)
  {
    increment();
    return *this;
  }

  // x++
  self operator++()
  {
    self copy = *this;
    increment();
    return copy;
  }

  ValueType operator*() const
  {
    return *m_low_itr;
  }

  ValueType* operator->() const
  {
    return &*m_low_itr;
  }

 private:
  void find_next_valid_item()
  {
    // if low_itr is at the end of current container, go to next container
    while (m_low_itr == get_end(*m_high_itr)) {
      ++m_high_itr;
      if (m_high_itr == m_high_end_itr) {
        // We reached the end! Set low_itr to invalid and return
        m_low_itr = LowLevelItrType();
        return;
      }
      m_low_itr = get_begin(*m_high_itr);
    }
  }

  // The 4 methods below are needed for supporting the notion that high_itr
  // can dereference to a container or a pointer to a container.

  template <typename Cont>
  LowLevelItrType
  get_begin(Cont& container)
  {
    return LowLevelItrType(container.begin());
  }

  template <typename Cont>
  LowLevelItrType
  get_end(Cont& container)
  {
    return LowLevelItrType(container.end());
  }

  template <class Cont>
  LowLevelItrType
  get_begin(Cont* container_ptr)
  {
    return LowLevelItrType(container_ptr->begin());
  }

  template <class Cont>
  LowLevelItrType
  get_end(Cont* container_ptr)
  {
    return LowLevelItrType(container_ptr->end());
  }

  void increment()
  {
    ++m_low_itr;
    find_next_valid_item();
  }

  HighLevelItrType   m_high_itr;
  LowLevelItrType    m_low_itr;
  HighLevelItrType   m_high_end_itr;
};

// Requirements:
//   BucketIteratorType must dereference to a Bucket*
//
// Incrementing this iterator will take us to the next *selected* bucket, skipping past
// unselected buckets, until the end.
template <typename BucketIteratorType>
class SelectedBucketIterator : public std::iterator<std::forward_iterator_tag, Bucket*>
{
 public:
  typedef SelectedBucketIterator<BucketIteratorType> self;

  SelectedBucketIterator(const Selector& selector, BucketIteratorType bucket_itr, BucketIteratorType bucket_end_itr) :
    m_bucket_itr(bucket_itr),
    m_bucket_end_itr(bucket_end_itr),
    m_selector(selector)
  {
    if (bucket_itr != bucket_end_itr) {
      find_next_valid_item();
    }
  }

  // Construct the "end" iterator
  SelectedBucketIterator(BucketIteratorType bucket_end_itr) :
    m_bucket_itr(bucket_end_itr),
    m_bucket_end_itr(bucket_end_itr),
    m_selector()
  {}

  SelectedBucketIterator(const self& rhs) :
    m_bucket_itr(rhs.m_bucket_itr),
    m_bucket_end_itr(rhs.m_bucket_end_itr),
    m_selector(rhs.m_selector)
  {}

  self& operator=(const self& rhs)
  {
    m_bucket_itr     = rhs.m_bucket_itr;
    m_bucket_end_itr = rhs.m_bucket_end_itr;
    m_selector       = rhs.m_selector;

    return *this;
  }

  bool operator==(const self& rhs) const
  {
    return (m_bucket_itr == rhs.m_bucket_itr);
  }

  bool operator==(const BucketIteratorType& rhs) const
  {
    return (m_bucket_itr == rhs);
  }

  bool operator!=(const self& rhs) const
  {
    return !(*this == rhs);
  }

  // ++x
  self& operator++(int)
  {
    increment();
    return *this;
  }

  // x++
  self operator++()
  {
    self copy = *this;
    increment();
    return copy;
  }

  // The method below is why boost::filter_iterator won't work for us. filter_iterator
  // deferences to a reference, tranform iterator dereferences to a copy, making them
  // incompatible.
  Bucket* operator*() const
  {
    return *m_bucket_itr;
  }

  Bucket** operator->() const
  {
    return &*m_bucket_itr;
  }

 private:
  void find_next_valid_item()
  {
    while (m_bucket_itr != m_bucket_end_itr && !m_selector(**m_bucket_itr)) {
      ++m_bucket_itr;
    }
  }

  void increment()
  {
    ++m_bucket_itr;
    find_next_valid_item();
  }

  BucketIteratorType m_bucket_itr;
  BucketIteratorType m_bucket_end_itr;
  Selector           m_selector;
};

// Iterator for iterating over all entities within each bucket of a vector of buckets
typedef TwoLevelIterator<std::vector<Bucket*>::const_iterator, BucketPtrIterator, Entity* const> BucketVectorEntityIterator;
typedef std::pair<BucketVectorEntityIterator, BucketVectorEntityIterator>                        BucketVectorEntityIteratorRange;

// Iterator for iterating over selected buckets within a vector of buckets
typedef SelectedBucketIterator<std::vector<Bucket*>::const_iterator>                      SelectedBucketVectorIterator;
//typedef boost::filter_iterator<Selector, std::vector<Bucket*>::const_iterator>            SelectedBucketVectorIterator;

// Iterator for iterating over all entities within each *selected* bucket of a vector of buckets
typedef TwoLevelIterator<SelectedBucketVectorIterator, BucketPtrIterator, Entity* const>  SelectedBucketVectorEntityIterator;
typedef std::pair<SelectedBucketVectorEntityIterator, SelectedBucketVectorEntityIterator> SelectedBucketVectorEntityIteratorRange;

// Iterator for iterating over all buckets in a vector of vectors of buckets
typedef TwoLevelIterator<std::vector<std::vector<Bucket*> >::const_iterator, std::vector<Bucket*>::const_iterator, Bucket* const> AllBucketsIterator;
typedef std::pair<AllBucketsIterator, AllBucketsIterator>                                                                         AllBucketsRange;

// Iterator for iterating over all *selected* buckets in a bucket range
typedef SelectedBucketIterator<AllBucketsIterator>                         AllSelectedBucketsIterator;
//typedef boost::filter_iterator<Selector, AllBucketsIterator>               AllSelectedBucketsIterator;
typedef std::pair<AllSelectedBucketsIterator, AllSelectedBucketsIterator>  AllSelectedBucketsRange;

// Iterator for iterating over all entities within each bucket of a bucket range
typedef TwoLevelIterator<AllBucketsIterator, BucketPtrIterator, Entity* const>  BucketRangeEntityIterator;
typedef std::pair<BucketRangeEntityIterator, BucketRangeEntityIterator>         BucketRangeEntityIteratorRange;

// Iterator for iterating over all *selected* entities withing a bucket range
typedef TwoLevelIterator<AllSelectedBucketsIterator, BucketPtrIterator, Entity* const>    SelectedBucketRangeEntityIterator;
typedef std::pair<SelectedBucketRangeEntityIterator, SelectedBucketRangeEntityIterator>   SelectedBucketRangeEntityIteratorRange;

//
// API - Convert collections into ranges. For internal use only. Clients should use
//       GetBuckets.hpp, GetEntities.hpp or their BulkData object.
//

// Get a range allowing you to iterate over all entities withing a collection of buckets
BucketVectorEntityIteratorRange get_entity_range(const std::vector<Bucket*>& buckets);

// Get a range allowing you to iterate over all *selected* entities withing a collection of buckets
SelectedBucketVectorEntityIteratorRange get_entity_range(const std::vector<Bucket*>& buckets, const Selector& selector);

// Get a range allowing you to iterate over all buckets within a collection of collections of buckets
AllBucketsRange get_bucket_range(const std::vector<std::vector<Bucket*> >& buckets);

// Get a range allowing you to iterate over a single collection of buckets within a collection of collections of buckets;
// the single collection is specified by the itr argument.
AllBucketsRange get_bucket_range(const std::vector<std::vector<Bucket*> >& buckets,
                                 std::vector<std::vector<Bucket*> >::const_iterator itr);

// Get a range allowing you to iterate over all *selected* buckets within a collection of collections of buckets
AllSelectedBucketsRange get_selected_bucket_range(const AllBucketsRange& bucket_range, const Selector& selector);

// Get a range allowing you to iterate over all *selected* buckets within a collection of collections of buckets
SelectedBucketRangeEntityIteratorRange get_selected_bucket_entity_range(const AllBucketsRange& bucket_range, const Selector& selector);

} //namespace mesh
} //namespace stk

#endif
