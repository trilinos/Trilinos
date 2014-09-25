// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
// 
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 

#ifndef Toolkit_Iterators_hpp
#define Toolkit_Iterators_hpp

#include <iterator>                     // for iterator_traits, etc
#include <stk_mesh/base/Bucket.hpp>     // for BucketIterator
#include <stk_mesh/base/Selector.hpp>   // for Selector
#include <utility>                      // for pair
#include <vector>                       // for vector<>::const_iterator, etc
#include "stk_mesh/base/Types.hpp"      // for BucketVector



namespace stk {
namespace mesh {

// Requirements:
//   Dereference of HighLevelItrType gives a container (OR pointer to container) with iterators compatible with LowLevelItrType
//   Dereference of LowLevelItrType gives ValueType
//   LowLevelItrType must default-construct to a deterministic invalid value
//
// Incrementing this iterator will take us to the next valid low-level iterator, skipping past
// empty high-level containers, until the end iterator.
template<typename HighLevelItrType, typename LowLevelItrType>
class TwoLevelIterator : public std::iterator<std::forward_iterator_tag, typename std::iterator_traits<LowLevelItrType>::value_type>
{
 public:
  typedef TwoLevelIterator<HighLevelItrType, LowLevelItrType>  self;
  typedef typename std::iterator_traits<LowLevelItrType>::reference reference;
  typedef typename std::iterator_traits<LowLevelItrType>::pointer pointer;

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

  TwoLevelIterator() :
    m_high_itr(),
    m_low_itr(),
    m_high_end_itr()
  {}

  bool operator==(const self& rhs) const
  {
    return (m_high_itr == rhs.m_high_itr && m_low_itr == rhs.m_low_itr);
  }

  bool operator!=(const self& rhs) const
  {
    return !(*this == rhs);
  }

  // x++
  self operator++(int)
  {
    self copy = *this;
    increment();
    return copy;
  }

  // ++x
  self& operator++()
  {
    increment();
    return *this;
  }

  reference operator*() const
  {
    return *m_low_itr;
  }

  pointer operator->() const
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
//
// As long as we're using a pointer as the value type, we need to
// specify the reference type to be the value_type in order for this
// class to work with boost
template <typename BucketIteratorType>
class SelectedBucketIterator : public std::iterator<std::forward_iterator_tag,
                                                    typename BucketIteratorType::value_type,
                                                    typename BucketIteratorType::difference_type,
                                                    typename BucketIteratorType::pointer,
                                                    typename BucketIteratorType::value_type>
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

  SelectedBucketIterator() :
    m_bucket_itr(),
    m_bucket_end_itr(),
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

  // x++
  self operator++(int)
  {
    self copy = *this;
    increment();
    return copy;
  }

  // ++x
  self& operator++()
  {
    increment();
    return *this;
  }

  // The method below is why boost::filter_iterator won't work for us. filter_iterator
  // deferences to a reference, tranform iterator dereferences to a copy, making them
  // incompatible.
  typename BucketIteratorType::value_type operator*() const
  {
    return *m_bucket_itr;
  }

  typename BucketIteratorType::pointer operator->() const
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
typedef TwoLevelIterator<BucketVector::const_iterator, BucketIterator> BucketVectorEntityIterator;
typedef std::pair<BucketVectorEntityIterator, BucketVectorEntityIterator>      BucketVectorEntityIteratorRange;

// Iterator for iterating over selected buckets within a vector of buckets
typedef SelectedBucketIterator<BucketVector::const_iterator>                      SelectedBucketVectorIterator;
//typedef boost::filter_iterator<Selector, BucketVector::const_iterator>            SelectedBucketVectorIterator;
typedef std::pair<SelectedBucketVectorIterator, SelectedBucketVectorIterator>             SelectedBucketVectorIteratorRange;

// Iterator for iterating over all entities within each *selected* bucket of a vector of buckets
typedef TwoLevelIterator<SelectedBucketVectorIterator, BucketIterator>                 SelectedBucketVectorEntityIterator;
typedef std::pair<SelectedBucketVectorEntityIterator, SelectedBucketVectorEntityIterator> SelectedBucketVectorEntityIteratorRange;

// Iterator for iterating over all buckets in a vector of vectors of buckets
typedef TwoLevelIterator<std::vector<BucketVector >::const_iterator, BucketVector::const_iterator> AllBucketsIterator;
typedef std::pair<AllBucketsIterator, AllBucketsIterator>                                                          AllBucketsRange;

// Iterator for iterating over all *selected* buckets in a bucket range
typedef SelectedBucketIterator<AllBucketsIterator>                         AllSelectedBucketsIterator;
//typedef boost::filter_iterator<Selector, AllBucketsIterator>               AllSelectedBucketsIterator;
typedef std::pair<AllSelectedBucketsIterator, AllSelectedBucketsIterator>  AllSelectedBucketsRange;

// Iterator for iterating over all entities within each bucket of a bucket range
typedef TwoLevelIterator<AllBucketsIterator, BucketIterator>            BucketRangeEntityIterator;
typedef std::pair<BucketRangeEntityIterator, BucketRangeEntityIterator> BucketRangeEntityIteratorRange;

// Iterator for iterating over all *selected* entities withing a bucket range
typedef TwoLevelIterator<AllSelectedBucketsIterator, BucketIterator>                     SelectedBucketRangeEntityIterator;
typedef std::pair<SelectedBucketRangeEntityIterator, SelectedBucketRangeEntityIterator>  SelectedBucketRangeEntityIteratorRange;

//
// API - Convert collections into ranges. For internal use only. Clients should use
//       GetBuckets.hpp, GetEntities.hpp or their BulkData object.
//

// Get a range allowing you iterate over selected buckets in a vector
SelectedBucketVectorIteratorRange get_selected_bucket_range(const BucketVector& buckets, const Selector& selector);

} //namespace mesh
} //namespace stk

#endif
