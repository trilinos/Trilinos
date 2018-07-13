/*
// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
// @HEADER
*/

#ifndef TPETRA_DETAILS_FIXEDHASHTABLE_DEF_HPP
#define TPETRA_DETAILS_FIXEDHASHTABLE_DEF_HPP

#include "Tpetra_Details_FixedHashTable_decl.hpp"
#include "Tpetra_Details_computeOffsets.hpp"
#ifdef TPETRA_USE_MURMUR_HASH
#  include "Kokkos_Functional.hpp" // hash function used by Kokkos::UnorderedMap
#endif // TPETRA_USE_MURMUR_HASH
#include "Kokkos_ArithTraits.hpp"
#include "Teuchos_TypeNameTraits.hpp"
#include <type_traits>

namespace Tpetra {
namespace Details {
//
// This namespace stores utility functions and Kokkos
// functors for use in FixedHashTable construction.
//
namespace FHT {


// Is it worth actually using building the FixedHashTable using
// parallel threads, instead of just counting in a sequential loop?
//
// The parallel version of FixedHashTable construction isn't just a
// parallelization of the sequential loops.  It incurs additional
// overheads.  For example, the CountBuckets kernel uses atomic update
// instructions to count the number of "buckets" per offsets array
// (ptr) entry.  Atomic updates have overhead, even if only one thread
// issues them.  The Kokkos kernels are still correct in that case,
// but I would rather not incur overhead then.  It might make sense to
// set the minimum number of threads to something greater than 1, but
// we would need experiments to find out.
template<class ExecSpace>
bool worthBuildingFixedHashTableInParallel () {
  return ExecSpace::concurrency() > 1;
}

// If the input kokkos::View<const KeyType*, ArrayLayout,
// InputExecSpace, Kokkos::MemoryUnamanged> is NOT accessible from the
// OutputExecSpace execution space, make and return a deep copy.
// Otherwise, just return the original input.
//
// The point of this is to avoid unnecessary copies, when the input
// array of keys comes in as a Teuchos::ArrayView (which we wrap in an
// unmanaged Kokkos::View).
template<class KeyType,
         class ArrayLayout,
         class InputExecSpace,
         class OutputExecSpace,
         const bool mustDeepCopy =
           ! std::is_same<typename InputExecSpace::memory_space,
                          typename OutputExecSpace::memory_space>::value>
struct DeepCopyIfNeeded {
  // The default implementation is trivial; all the work happens in
  // partial specializations.
};

// Specialization for when a deep copy is actually needed.
template<class KeyType,
         class ArrayLayout,
         class InputExecSpace,
         class OutputExecSpace>
struct DeepCopyIfNeeded<KeyType, ArrayLayout, InputExecSpace, OutputExecSpace, true>
{
  typedef Kokkos::View<const KeyType*, ArrayLayout,
                       InputExecSpace, Kokkos::MemoryUnmanaged> input_view_type;
  // In this case, a deep copy IS needed.  As a result, the output
  // type is a managed Kokkos::View, which differs from the input
  // type.  Clients must get the correct return type from this struct,
  // either from the typedef below or from 'auto'.  Assigning an
  // unmanaged View to a managed View is a syntax error.
  typedef Kokkos::View<const KeyType*, ArrayLayout, OutputExecSpace> output_view_type;

  static output_view_type copy (const input_view_type& src) {
    typedef typename output_view_type::non_const_type NC;

    NC dst (Kokkos::ViewAllocateWithoutInitializing (src.tracker ().label ()),
            src.extent (0));
    Kokkos::deep_copy (dst, src);
    return output_view_type (dst);
  }
};

// Specialization if no need to make a deep copy.
template<class KeyType,
         class ArrayLayout,
         class InputExecSpace,
         class OutputExecSpace>
struct DeepCopyIfNeeded<KeyType, ArrayLayout, InputExecSpace, OutputExecSpace, false> {
  typedef Kokkos::View<const KeyType*, ArrayLayout,
                       InputExecSpace, Kokkos::MemoryUnmanaged> input_view_type;
  typedef Kokkos::View<const KeyType*, ArrayLayout, OutputExecSpace,
                       Kokkos::MemoryUnmanaged> output_view_type;

  static output_view_type copy (const input_view_type& src) {
    return output_view_type (src);
  }
};

//
// Functors for FixedHashTable initialization
//
// NOTE (mfh 23 May 2015): Once we can use lambdas with CUDA, we
// should consider replacing all of these functors with in-line
// lambdas.  The only issue is that we would need to bake the
// execution space into the policy, since the default execution space
// might differ from the one Tpetra wants to use.

/// \brief Parallel for functor for counting "buckets" in the FixedHashTable.
///
/// \tparam CountsViewType Type of the Kokkos::View specialization
///   used to store the bucket counts; the output of this functor.
/// \tparam KeysViewType Type of the Kokkos::View specialization
///   used to store the keys; the input of this functor.
/// \tparam SizeType The parallel loop index type; a built-in integer
///   type.  Defaults to the type of the input View's dimension.  You
///   may use a shorter type to improve performance.
template<class CountsViewType,
         class KeysViewType,
         class SizeType = typename KeysViewType::size_type>
class CountBuckets {
public:
  typedef CountsViewType counts_view_type;
  typedef KeysViewType keys_view_type;
  typedef typename CountsViewType::execution_space execution_space;
  typedef typename CountsViewType::memory_space memory_space;
  typedef SizeType size_type;
  typedef typename keys_view_type::non_const_value_type key_type;
  // mfh 21 May 2015: Having a device_type typedef in the functor
  // along with an execution_space typedef causes compilation issues.
  // This is because one of Kokkos' partial specializations picks up
  // on the device_type typedef, and another picks up on the
  // execution_space typedef.  The former is a legacy of a previous
  // design iteration of Kokkos, which did not separate memory and
  // execution spaces.
  typedef Tpetra::Details::Hash<key_type, Kokkos::Device<execution_space, memory_space> > hash_type;

  /// \brief Constructor
  ///
  /// \param counts [out] (Preallocated) View of the bucket counts
  /// \param keys [in] View of the keys
  /// \param size [in] Number of buckets; length of \c counts
  CountBuckets (const counts_view_type& counts,
                const keys_view_type& keys,
                const size_type size) :
    counts_ (counts),
    keys_ (keys),
    size_ (size)
  {}

  /// \brief Do this for every entry of \c keys_.
  ///
  /// Count the number of keys in \c keys_ that hash to the same
  /// value.
  KOKKOS_INLINE_FUNCTION void
  operator () (const size_type& i) const
  {
    typedef typename hash_type::result_type hash_value_type;

    const hash_value_type hashVal = hash_type::hashFunc (keys_[i], size_);
    Kokkos::atomic_increment (&counts_[hashVal]);
  }

private:
  //! Bucket counts (output argument).
  counts_view_type counts_;
  //! Keys for the FixedHashTable to construct (input argument).
  keys_view_type keys_;
  //! Number of buckets plus 1 (or 0, if no buckets).
  size_type size_;
};

/// \brief Reduction result for FillPairs functor below.
///
/// The reduction result finds the min and max keys, and reports
/// whether FillPairs succeeded (it should always succeed, unless
/// there a bug in earlier code manifested).  The default values
/// (minKey_ set to max KeyType value, and maxKey_ set to min KeyType
/// value) ensure correct behavior even if there is only one key.
///
/// \tparam KeyType Type of each input of the hash function.
///   It must be an integer type.
template<class KeyType>
struct FillPairsResult {
  KOKKOS_INLINE_FUNCTION
  FillPairsResult () :
    minKey_ (::Kokkos::Details::ArithTraits<KeyType>::max ()),
    // min() for a floating-point type returns the minimum _positive_
    // normalized value.  This is different than for integer types.
    // lowest() is new in C++11 and returns the least value, always
    // negative for signed finite types.
    //
    // mfh 23 May 2015: I have heard reports that
    // std::numeric_limits<int>::lowest() does not exist with the
    // Intel compiler.  I'm not sure if the users in question actually
    // enabled C++11.  However, it's easy enough to work around this
    // issue.  The standard floating-point types are signed and have a
    // sign bit, so lowest() is just -max().  For integer types, we
    // can use min() instead.
    maxKey_ (::Kokkos::Details::ArithTraits<KeyType>::is_integer ?
             ::Kokkos::Details::ArithTraits<KeyType>::min () :
             -::Kokkos::Details::ArithTraits<KeyType>::max ()),
    success_ (true)
  {
    static_assert (std::is_arithmetic<KeyType>::value, "FillPairsResult: "
                   "KeyType must be some kind of number type.");
  }

  KOKKOS_INLINE_FUNCTION
  FillPairsResult (const KeyType& initMinKey,
                   const KeyType& initMaxKey) :
    minKey_ (initMinKey),
    maxKey_ (initMaxKey),
    success_ (true)
  {
    static_assert (std::is_arithmetic<KeyType>::value, "FillPairsResult: "
                   "KeyType must be some kind of number type.");
  }

  KeyType minKey_; //!< The current minimum key
  KeyType maxKey_; //!< The current maximum key
  bool success_;   //!< Whether fill succeeded (it can only fail on a bug)
};

/// \brief Parallel reduce functor for filling the FixedHashTable, and
///   computing the min and max keys.
///
/// This is also a parallel reduce functor in order to check for
/// failure.  Failure should only happen on a bug (in CountBuckets or
/// computeOffsetsFromCounts), but checking for it is cheap and easy.
///
/// \tparam PairsViewType Type of the Kokkos::View specialization used
///   to store the (key,value) pairs in the FixedHashTable; output of
///   this functor.
/// \tparam KeysViewType Type of the Kokkos::View specialization
///   used to store the keys; input of this functor.
/// \tparam CountsViewType Type of the Kokkos::View specialization
///   used to store the bucket counts; input of this functor, used as
///   scratch space (so it must be nonconst; offsets_view_type is the
///   const version of this).
/// \tparam SizeType The parallel loop index type; a built-in integer
///   type.  Defaults to the type of the input View's dimension.  You
///   may use a shorter type to improve performance.
///
/// \note This functor CANNOT also find duplicate keys.  This is
///   because different threads might have different keys that hash to
///   the same bucket, concurrently.  Threads resolve this by using
///   atomic updates to "reserve" a position in the 'vals' output
///   array.  Thus, two threads might concurrently check for existing
///   duplicates and find none, but then concurrently insert
///   (key,value) pairs with duplicate keys.  The only thread-scalable
///   way to check for duplicate keys is to wait until after the table
///   has been filled.
template<class PairsViewType,
         class KeysViewType,
         class CountsViewType,
         class SizeType = typename KeysViewType::size_type>
class FillPairs {
public:
  typedef typename CountsViewType::non_const_type counts_view_type;
  typedef typename counts_view_type::const_type offsets_view_type;

  typedef PairsViewType pairs_view_type;
  typedef typename KeysViewType::const_type keys_view_type;
  typedef typename offsets_view_type::execution_space execution_space;
  typedef typename offsets_view_type::memory_space memory_space;
  typedef SizeType size_type;

  typedef typename keys_view_type::non_const_value_type key_type;
  typedef typename pairs_view_type::non_const_value_type pair_type;

  typedef FillPairsResult<key_type> value_type;

  // mfh 23 May 2015: Having a device_type typedef in the functor
  // along with an execution_space typedef causes compilation issues.
  // This is because one of Kokkos' partial specializations picks up
  // on the device_type typedef, and another picks up on the
  // execution_space typedef.  The former is a legacy of a previous
  // design iteration of Kokkos, which did not separate memory and
  // execution spaces.
  typedef Tpetra::Details::Hash<key_type, Kokkos::Device<execution_space, memory_space> > hash_type;

  /// \brief Constructor
  ///
  /// \param pairs [out] (Preallocated) View of (key,value) pairs
  /// \param counts [in/out] View of bucket counts; overwritten as
  ///   scratch space
  /// \param ptr [in] View of offsets
  /// \param keys [in] View of the keys
  /// \param startingValue [in] Starting value.  For each key keys[i],
  ///   the corresponding value (in the (key,value) pair) is
  ///   startingValue + i.
  FillPairs (const pairs_view_type& pairs,
             const counts_view_type& counts,
             const offsets_view_type& ptr,
             const keys_view_type& keys,
             const typename pair_type::second_type startingValue) :
    pairs_ (pairs),
    counts_ (counts),
    ptr_ (ptr),
    keys_ (keys),
    size_ (counts.extent (0)),
    startingValue_ (startingValue),
    initMinKey_ (::Kokkos::Details::ArithTraits<key_type>::max ()),
    initMaxKey_ (::Kokkos::Details::ArithTraits<key_type>::is_integer ?
                 ::Kokkos::Details::ArithTraits<key_type>::min () :
                 -::Kokkos::Details::ArithTraits<key_type>::max ())
  {}

  /// \brief Constructor that takes initial min and max key values.
  ///
  /// This constructor is useful for Tpetra::Map's noncontiguous
  /// constructor.  That constructor first harvests an initial
  /// sequence of contiguous global indices, then puts any remaining
  /// global indices that follow into the hash table.  That initial
  /// sequence defines initial min and max keys.
  ///
  /// \param pairs [out] (Preallocated) View of (key,value) pairs
  /// \param counts [in/out] View of bucket counts; overwritten as
  ///   scratch space
  /// \param ptr [in] View of offsets
  /// \param keys [in] View of the keys
  /// \param startingValue [in] Starting value.  For each key keys[i],
  ///   the corresponding value (in the (key,value) pair) is
  ///   startingValue + i.
  /// \param initMinKey [in] Initial min key value
  /// \param initMaxKey [in] Initial max key value
  FillPairs (const pairs_view_type& pairs,
             const counts_view_type& counts,
             const offsets_view_type& ptr,
             const keys_view_type& keys,
             const typename pair_type::second_type startingValue,
             const key_type initMinKey,
             const key_type initMaxKey) :
    pairs_ (pairs),
    counts_ (counts),
    ptr_ (ptr),
    keys_ (keys),
    size_ (counts.extent (0)),
    startingValue_ (startingValue),
    initMinKey_ (initMinKey),
    initMaxKey_ (initMaxKey)
  {}

  //! Set the initial value of the reduction result.
  KOKKOS_INLINE_FUNCTION void init (value_type& dst) const
  {
    dst.minKey_ = initMinKey_;
    dst.maxKey_ = initMaxKey_;
    dst.success_ = true;
  }

  KOKKOS_INLINE_FUNCTION void
  join (volatile value_type& dst,
        const volatile value_type& src) const
  {
    if (src.maxKey_ > dst.maxKey_) {
      dst.maxKey_ = src.maxKey_;
    }
    if (src.minKey_ < dst.minKey_) {
      dst.minKey_ = src.minKey_;
    }
    dst.success_ = dst.success_ && src.success_;
  }

  /// \brief Parallel loop body; do this for every entry of \c keys_.
  ///
  /// Add (key = keys_[i], value = startingValue_ + i) pair to the
  /// hash table.  Compute min and max key value.
  KOKKOS_INLINE_FUNCTION void
  operator () (const size_type& i, value_type& dst) const
  {
    typedef typename hash_type::result_type hash_value_type;
    typedef typename offsets_view_type::non_const_value_type offset_type;
    typedef typename pair_type::second_type val_type;
    typedef typename std::remove_reference< decltype( counts_[0] ) >::type atomic_incr_type;

    const key_type key = keys_[i];
    if (key > dst.maxKey_) {
      dst.maxKey_ = key;
    }
    if (key < dst.minKey_) {
      dst.minKey_ = key;
    }
    const val_type theVal = startingValue_ + static_cast<val_type> (i);
    const hash_value_type hashVal = hash_type::hashFunc (key, size_);

    // Return the old count; decrement afterwards.
    const offset_type count = Kokkos::atomic_fetch_add (&counts_[hashVal], atomic_incr_type(-1));
    if (count == 0) {
      dst.success_ = false; // FAILURE!
    }
    else {
      const offset_type curPos = ptr_[hashVal+1] - count;

      pairs_[curPos].first = key;
      pairs_[curPos].second = theVal;
    }
  }

private:
  pairs_view_type pairs_;
  counts_view_type counts_;
  offsets_view_type ptr_;
  keys_view_type keys_;
  size_type size_;
  typename pair_type::second_type startingValue_;
  //! Initial minimum key.
  key_type initMinKey_;
  //! Initial maximum key.
  key_type initMaxKey_;
};

/// \brief Functor for checking whether a FixedHashTable has one or
///   more duplicate entries.
///
/// \tparam OffsetsViewType Type of the Kokkos::View specialization
///   used to store the "row" offsets; input of this functor.
/// \tparam PairsViewType Type of the Kokkos::View specialization used
///   to store the (key,value) pairs in the FixedHashTable; input of
///   this functor.
/// \tparam SizeType The parallel loop index type; a built-in integer
///   type.  Defaults to the type of the input View's dimension.  You
///   may use a shorter type to improve performance.
///
/// This functor works by iterating over all the hash "buckets" (the
/// entries of \c ptr).  For each i, all the keys in the (key,value)
/// pairs in the half-exclusive range ptr[i], ptr[i+1] must have the
/// same hash value.  Thus, any duplicate keys must be in that range;
/// if the hash function is actually a function, it's impossible to
/// find duplicate keys elsewhere.
///
/// If the hash function isn't actually a function, the above
/// algorithm won't work.  However, users are not (currently) allowed
/// to supply arbitrary hash functions to Tpetra, so we don't have to
/// worry about checking the hash function here.
template<class OffsetsViewType,
         class PairsViewType,
         class SizeType = typename OffsetsViewType::size_type>
class CheckForDuplicateKeys {
public:
  typedef typename OffsetsViewType::const_type offsets_view_type;
  typedef typename PairsViewType::const_type pairs_view_type;
  typedef typename offsets_view_type::execution_space execution_space;
  typedef typename offsets_view_type::memory_space memory_space;
  typedef SizeType size_type;

  // The result of the check is whether the table has one or more duplicates.
  typedef int value_type;

  /// \brief Constructor
  ///
  /// \param pairs [in] View of the FixedHashTable's (key,value) pairs
  /// \param ptr [in] View of the FixedHashTable's "bucket" offsets
  CheckForDuplicateKeys (const pairs_view_type& pairs,
                         const offsets_view_type& ptr) :
    pairs_ (pairs),
    ptr_ (ptr),
    size_ (ptr_.extent (0) == 0 ?
           size_type (0) :
           ptr_.extent (0) - 1)
  {}

  //! Set the initial value of the reduction result.
  KOKKOS_INLINE_FUNCTION void init (value_type& dst) const
  {
    dst = 0;
  }

  //! Combine two intermediate reduction results.
  KOKKOS_INLINE_FUNCTION void
  join (volatile value_type& dst,
        const volatile value_type& src) const
  {
    dst = dst + src > 0?1:0;
  }

  //! Parallel loop body.
  KOKKOS_INLINE_FUNCTION void
  operator () (const size_type& i, value_type& dst) const
  {
    typedef typename offsets_view_type::non_const_value_type offset_type;
    typedef typename pairs_view_type::non_const_value_type pair_type;
    typedef typename pair_type::first_type key_type;

    if (dst>0) {
      return; // we've already found duplicate keys elsewhere
    }
    else {
      const offset_type beg = ptr_[i];
      const offset_type end = ptr_[i+1];
      bool foundDuplicateKey = false;
      // This is an ~ n^2 algorithm in the worst case, where n is the
      // max number of keys that hash to the same bucket.  However, if
      // the hash function is reasonable, n should be much less than
      // the total number of keys.
      for (offset_type j = beg + 1; j < end; ++j) {
        const key_type curKey = pairs_[j].first;
        for (offset_type k = beg; k < j; ++k) {
          if (pairs_[k].first == curKey) {
            foundDuplicateKey = true;
            break;
          }
        }
      }
      dst = (dst>0) || foundDuplicateKey?1:0;
    }
  }

private:
  pairs_view_type pairs_;
  offsets_view_type ptr_;
  size_type size_;
};

} // namespace FHT

//
// Here begins the actual implementation of FixedHashTable.
//

template<class KeyType, class ValueType, class DeviceType>
bool
FixedHashTable<KeyType, ValueType, DeviceType>::
hasKeys () const {
  // This also works if FixedHashTable has no entries.  getKey()
  // works in that case, but always returns the flag (invalid).
  //
  // FIXME (31 May 2015) This only works because vals_ contains no
  // padding.  If we ever pad within a "row" of vals_, we'll have to
  // change this.
  return keys_.extent (0) == val_.extent (0);
}

template<class KeyType, class ValueType, class DeviceType>
void
FixedHashTable<KeyType, ValueType, DeviceType>::
check () const
{
  // const char prefix[] = "Tpetra::Details::FixedHashTable: ";
  // const char suffix[] = "  Please report this bug to the Tpetra developers.";
}

template<class KeyType, class ValueType, class DeviceType>
FixedHashTable<KeyType, ValueType, DeviceType>::
FixedHashTable () :
  minKey_ (::Kokkos::Details::ArithTraits<KeyType>::max ()),
  maxKey_ (::Kokkos::Details::ArithTraits<KeyType>::is_integer ?
           ::Kokkos::Details::ArithTraits<KeyType>::min () :
           -::Kokkos::Details::ArithTraits<KeyType>::max ()),
  minVal_ (::Kokkos::Details::ArithTraits<ValueType>::max ()),
  maxVal_ (::Kokkos::Details::ArithTraits<ValueType>::is_integer ?
           ::Kokkos::Details::ArithTraits<ValueType>::min () :
           -::Kokkos::Details::ArithTraits<ValueType>::max ()),
  firstContigKey_ (::Kokkos::Details::ArithTraits<KeyType>::max ()),
  lastContigKey_ (::Kokkos::Details::ArithTraits<KeyType>::is_integer ?
                  ::Kokkos::Details::ArithTraits<KeyType>::min () :
                  -::Kokkos::Details::ArithTraits<KeyType>::max ()),
  contiguousValues_ (true), // trivially
  checkedForDuplicateKeys_ (true), // it's an empty table; no need to check
  hasDuplicateKeys_ (false)
{
#ifdef HAVE_TPETRA_DEBUG
  check ();
#endif // HAVE_TPETRA_DEBUG
}

template<class KeyType, class ValueType, class DeviceType>
FixedHashTable<KeyType, ValueType, DeviceType>::
FixedHashTable (const keys_type& keys) :
  keys_ (keys),
  minKey_ (::Kokkos::Details::ArithTraits<KeyType>::max ()), // to be set in init()
  maxKey_ (::Kokkos::Details::ArithTraits<KeyType>::is_integer ?
           ::Kokkos::Details::ArithTraits<KeyType>::min () :
           -::Kokkos::Details::ArithTraits<KeyType>::max ()), // to be set in init()
  minVal_ (0),
  maxVal_ (keys.size () == 0 ?
           static_cast<ValueType> (0) :
           static_cast<ValueType> (keys.size () - 1)),
  firstContigKey_ (::Kokkos::Details::ArithTraits<KeyType>::max ()),
  lastContigKey_ (::Kokkos::Details::ArithTraits<KeyType>::is_integer ?
                  ::Kokkos::Details::ArithTraits<KeyType>::min () :
                  -::Kokkos::Details::ArithTraits<KeyType>::max ()),
  contiguousValues_ (true),
  checkedForDuplicateKeys_ (false),
  hasDuplicateKeys_ (false) // to revise in hasDuplicateKeys()
{
  const ValueType startingValue = static_cast<ValueType> (0);
  const KeyType initMinKey = this->minKey_;
  const KeyType initMaxKey = this->maxKey_;
  this->init (keys, startingValue, initMinKey, initMaxKey,
              initMinKey, initMinKey, false);

#ifdef HAVE_TPETRA_DEBUG
  check ();
#endif // HAVE_TPETRA_DEBUG
}

template<class KeyType, class ValueType, class DeviceType>
FixedHashTable<KeyType, ValueType, DeviceType>::
FixedHashTable (const Teuchos::ArrayView<const KeyType>& keys,
                const bool keepKeys) :
  minKey_ (::Kokkos::Details::ArithTraits<KeyType>::max ()), // to be set in init()
  maxKey_ (::Kokkos::Details::ArithTraits<KeyType>::is_integer ?
           ::Kokkos::Details::ArithTraits<KeyType>::min () :
           -::Kokkos::Details::ArithTraits<KeyType>::max ()), // to be set in init()
  minVal_ (0),
  maxVal_ (keys.size () == 0 ?
           static_cast<ValueType> (0) :
           static_cast<ValueType> (keys.size () - 1)),
  firstContigKey_ (::Kokkos::Details::ArithTraits<KeyType>::max ()),
  lastContigKey_ (::Kokkos::Details::ArithTraits<KeyType>::is_integer ?
                  ::Kokkos::Details::ArithTraits<KeyType>::min () :
                  -::Kokkos::Details::ArithTraits<KeyType>::max ()),
  contiguousValues_ (true),
  checkedForDuplicateKeys_ (false),
  hasDuplicateKeys_ (false) // to revise in hasDuplicateKeys()
{
  typedef typename keys_type::non_const_type nonconst_keys_type;

  // mfh 01 May 2015: I don't trust that
  // Teuchos::ArrayView::getRawPtr() returns NULL when the size is 0,
  // so I ensure this manually.
  const ValueType startingValue = static_cast<ValueType> (0);
  host_input_keys_type keys_k (keys.size () == 0 ? NULL : keys.getRawPtr (),
                               keys.size ());
  using Kokkos::ViewAllocateWithoutInitializing;
  nonconst_keys_type keys_d (ViewAllocateWithoutInitializing ("FixedHashTable::keys"),
                             keys_k.extent (0));
  Kokkos::deep_copy (keys_d, keys_k);
  const KeyType initMinKey = this->minKey_;
  const KeyType initMaxKey = this->maxKey_;
  this->init (keys_d, startingValue, initMinKey, initMaxKey,
              initMinKey, initMinKey, false);
  if (keepKeys) {
    keys_ = keys_d;
#ifdef HAVE_TPETRA_DEBUG
    typedef typename keys_type::size_type size_type;
    TEUCHOS_TEST_FOR_EXCEPTION
      (keys_.extent (0) != static_cast<size_type> (keys.size ()),
       std::logic_error, "Tpetra::Details::FixedHashTable constructor: "
       "keepKeys is true, but on return, keys_.extent(0) = " <<
       keys_.extent (0) << " != keys.size() = " << keys.size () <<
       ".  Please report this bug to the Tpetra developers.");
#endif // HAVE_TPETRA_DEBUG
  }

#ifdef HAVE_TPETRA_DEBUG
  check ();
#endif // HAVE_TPETRA_DEBUG
}

template<class KeyType, class ValueType, class DeviceType>
FixedHashTable<KeyType, ValueType, DeviceType>::
FixedHashTable (const Teuchos::ArrayView<const KeyType>& keys,
                const ValueType startingValue,
                const bool keepKeys) :
  minKey_ (::Kokkos::Details::ArithTraits<KeyType>::max ()),
  maxKey_ (::Kokkos::Details::ArithTraits<KeyType>::is_integer ?
           ::Kokkos::Details::ArithTraits<KeyType>::min () :
           -::Kokkos::Details::ArithTraits<KeyType>::max ()),
  minVal_ (startingValue),
  maxVal_ (keys.size () == 0 ?
           startingValue :
           static_cast<ValueType> (startingValue + keys.size () - 1)),
  firstContigKey_ (::Kokkos::Details::ArithTraits<KeyType>::max ()),
  lastContigKey_ (::Kokkos::Details::ArithTraits<KeyType>::is_integer ?
                  ::Kokkos::Details::ArithTraits<KeyType>::min () :
                  -::Kokkos::Details::ArithTraits<KeyType>::max ()),
  contiguousValues_ (true),
  checkedForDuplicateKeys_ (false),
  hasDuplicateKeys_ (false) // to revise in hasDuplicateKeys()
{
  typedef typename keys_type::non_const_type nonconst_keys_type;

  // mfh 01 May 2015: I don't trust that
  // Teuchos::ArrayView::getRawPtr() returns NULL when the size is 0,
  // so I ensure this manually.
  host_input_keys_type keys_k (keys.size () == 0 ? NULL : keys.getRawPtr (),
                               keys.size ());
  using Kokkos::ViewAllocateWithoutInitializing;
  nonconst_keys_type keys_d (ViewAllocateWithoutInitializing ("FixedHashTable::keys"),
                             keys_k.extent (0));
  Kokkos::deep_copy (keys_d, keys_k);

  const KeyType initMinKey = ::Kokkos::Details::ArithTraits<KeyType>::max ();
  // min() for a floating-point type returns the minimum _positive_
  // normalized value.  This is different than for integer types.
  // lowest() is new in C++11 and returns the least value, always
  // negative for signed finite types.
  //
  // mfh 23 May 2015: I have heard reports that
  // std::numeric_limits<int>::lowest() does not exist with the Intel
  // compiler.  I'm not sure if the users in question actually enabled
  // C++11.  However, it's easy enough to work around this issue.  The
  // standard floating-point types are signed and have a sign bit, so
  // lowest() is just -max().  For integer types, we can use min()
  // instead.
  const KeyType initMaxKey = ::Kokkos::Details::ArithTraits<KeyType>::is_integer ?
    ::Kokkos::Details::ArithTraits<KeyType>::min () :
    -::Kokkos::Details::ArithTraits<KeyType>::max ();
  this->init (keys_d, startingValue, initMinKey, initMaxKey,
              initMinKey, initMinKey, false);
  if (keepKeys) {
    keys_ = keys_d;
#ifdef HAVE_TPETRA_DEBUG
    typedef typename keys_type::size_type size_type;
    TEUCHOS_TEST_FOR_EXCEPTION
      (keys_.extent (0) != static_cast<size_type> (keys.size ()),
       std::logic_error, "Tpetra::Details::FixedHashTable constructor: "
       "keepKeys is true, but on return, keys_.extent(0) = " <<
       keys_.extent (0) << " != keys.size() = " << keys.size () <<
       ".  Please report this bug to the Tpetra developers.");
#endif // HAVE_TPETRA_DEBUG
  }

#ifdef HAVE_TPETRA_DEBUG
  check ();
#endif // HAVE_TPETRA_DEBUG
}

template<class KeyType, class ValueType, class DeviceType>
FixedHashTable<KeyType, ValueType, DeviceType>::
FixedHashTable (const keys_type& keys,
                const KeyType firstContigKey,
                const KeyType lastContigKey,
                const ValueType startingValue,
                const bool keepKeys) :
  minKey_ (::Kokkos::Details::ArithTraits<KeyType>::max ()),
  maxKey_ (::Kokkos::Details::ArithTraits<KeyType>::is_integer ?
           ::Kokkos::Details::ArithTraits<KeyType>::min () :
           -::Kokkos::Details::ArithTraits<KeyType>::max ()),
  minVal_ (startingValue),
  maxVal_ (keys.size () == 0 ?
           startingValue :
           static_cast<ValueType> (startingValue + keys.size () - 1)),
  firstContigKey_ (firstContigKey),
  lastContigKey_ (lastContigKey),
  contiguousValues_ (true),
  checkedForDuplicateKeys_ (false),
  hasDuplicateKeys_ (false) // to revise in hasDuplicateKeys()
{
  const KeyType initMinKey = ::Kokkos::Details::ArithTraits<KeyType>::max ();
  // min() for a floating-point type returns the minimum _positive_
  // normalized value.  This is different than for integer types.
  // lowest() is new in C++11 and returns the least value, always
  // negative for signed finite types.
  //
  // mfh 23 May 2015: I have heard reports that
  // std::numeric_limits<int>::lowest() does not exist with the Intel
  // compiler.  I'm not sure if the users in question actually enabled
  // C++11.  However, it's easy enough to work around this issue.  The
  // standard floating-point types are signed and have a sign bit, so
  // lowest() is just -max().  For integer types, we can use min()
  // instead.
  const KeyType initMaxKey = ::Kokkos::Details::ArithTraits<KeyType>::is_integer ?
    ::Kokkos::Details::ArithTraits<KeyType>::min () :
    -::Kokkos::Details::ArithTraits<KeyType>::max ();
  this->init (keys, startingValue, initMinKey, initMaxKey,
              firstContigKey, lastContigKey, true);
  if (keepKeys) {
    keys_ = keys;
#ifdef HAVE_TPETRA_DEBUG
    typedef typename keys_type::size_type size_type;
    TEUCHOS_TEST_FOR_EXCEPTION
      (keys_.extent (0) != static_cast<size_type> (keys.size ()),
       std::logic_error, "Tpetra::Details::FixedHashTable constructor: "
       "keepKeys is true, but on return, keys_.extent(0) = " <<
       keys_.extent (0) << " != keys.size() = " << keys.size () <<
       ".  Please report this bug to the Tpetra developers.");
#endif // HAVE_TPETRA_DEBUG
  }

#ifdef HAVE_TPETRA_DEBUG
  check ();
#endif // HAVE_TPETRA_DEBUG
}

template<class KeyType, class ValueType, class DeviceType>
FixedHashTable<KeyType, ValueType, DeviceType>::
FixedHashTable (const Teuchos::ArrayView<const KeyType>& keys,
                const KeyType firstContigKey,
                const KeyType lastContigKey,
                const ValueType startingValue,
                const bool keepKeys) :
  minKey_ (::Kokkos::Details::ArithTraits<KeyType>::max ()),
  maxKey_ (::Kokkos::Details::ArithTraits<KeyType>::is_integer ?
           ::Kokkos::Details::ArithTraits<KeyType>::min () :
           -::Kokkos::Details::ArithTraits<KeyType>::max ()),
  minVal_ (startingValue),
  maxVal_ (keys.size () == 0 ?
           startingValue :
           static_cast<ValueType> (startingValue + keys.size () - 1)),
  firstContigKey_ (firstContigKey),
  lastContigKey_ (lastContigKey),
  contiguousValues_ (true),
  checkedForDuplicateKeys_ (false),
  hasDuplicateKeys_ (false) // to revise in hasDuplicateKeys()
{
  typedef typename keys_type::non_const_type nonconst_keys_type;

  // mfh 01 May 2015: I don't trust that
  // Teuchos::ArrayView::getRawPtr() returns NULL when the size is 0,
  // so I ensure this manually.
  host_input_keys_type keys_k (keys.size () == 0 ? NULL : keys.getRawPtr (),
                               keys.size ());
  using Kokkos::ViewAllocateWithoutInitializing;
  nonconst_keys_type keys_d (ViewAllocateWithoutInitializing ("FixedHashTable::keys"),
                             keys_k.extent (0));
  Kokkos::deep_copy (keys_d, keys_k);

  const KeyType initMinKey = ::Kokkos::Details::ArithTraits<KeyType>::max ();
  // min() for a floating-point type returns the minimum _positive_
  // normalized value.  This is different than for integer types.
  // lowest() is new in C++11 and returns the least value, always
  // negative for signed finite types.
  //
  // mfh 23 May 2015: I have heard reports that
  // std::numeric_limits<int>::lowest() does not exist with the Intel
  // compiler.  I'm not sure if the users in question actually enabled
  // C++11.  However, it's easy enough to work around this issue.  The
  // standard floating-point types are signed and have a sign bit, so
  // lowest() is just -max().  For integer types, we can use min()
  // instead.
  const KeyType initMaxKey = ::Kokkos::Details::ArithTraits<KeyType>::is_integer ?
    ::Kokkos::Details::ArithTraits<KeyType>::min () :
    -::Kokkos::Details::ArithTraits<KeyType>::max ();
  this->init (keys_d, startingValue, initMinKey, initMaxKey,
              firstContigKey, lastContigKey, true);
  if (keepKeys) {
    keys_ = keys_d;
#ifdef HAVE_TPETRA_DEBUG
    typedef typename keys_type::size_type size_type;
    TEUCHOS_TEST_FOR_EXCEPTION
      (keys_.extent (0) != static_cast<size_type> (keys.size ()),
       std::logic_error, "Tpetra::Details::FixedHashTable constructor: "
       "keepKeys is true, but on return, keys_.extent(0) = " <<
       keys_.extent (0) << " != keys.size() = " << keys.size () <<
       ".  Please report this bug to the Tpetra developers.");
#endif // HAVE_TPETRA_DEBUG
  }

#ifdef HAVE_TPETRA_DEBUG
  check ();
#endif // HAVE_TPETRA_DEBUG
}

template<class KeyType, class ValueType, class DeviceType>
FixedHashTable<KeyType, ValueType, DeviceType>::
FixedHashTable (const keys_type& keys,
                const ValueType startingValue) :
  keys_ (keys),
  minKey_ (::Kokkos::Details::ArithTraits<KeyType>::max ()),
  maxKey_ (::Kokkos::Details::ArithTraits<KeyType>::is_integer ?
           ::Kokkos::Details::ArithTraits<KeyType>::min () :
           -::Kokkos::Details::ArithTraits<KeyType>::max ()),
  minVal_ (startingValue),
  maxVal_ (keys.size () == 0 ?
           startingValue :
           static_cast<ValueType> (startingValue + keys.size () - 1)),
  firstContigKey_ (::Kokkos::Details::ArithTraits<KeyType>::max ()),
  lastContigKey_ (::Kokkos::Details::ArithTraits<KeyType>::is_integer ?
                  ::Kokkos::Details::ArithTraits<KeyType>::min () :
                  -::Kokkos::Details::ArithTraits<KeyType>::max ()),
  contiguousValues_ (true),
  checkedForDuplicateKeys_ (false),
  hasDuplicateKeys_ (false) // to revise in hasDuplicateKeys()
{
  const KeyType initMinKey = ::Kokkos::Details::ArithTraits<KeyType>::max ();
  // min() for a floating-point type returns the minimum _positive_
  // normalized value.  This is different than for integer types.
  // lowest() is new in C++11 and returns the least value, always
  // negative for signed finite types.
  //
  // mfh 23 May 2015: I have heard reports that
  // std::numeric_limits<int>::lowest() does not exist with the Intel
  // compiler.  I'm not sure if the users in question actually enabled
  // C++11.  However, it's easy enough to work around this issue.  The
  // standard floating-point types are signed and have a sign bit, so
  // lowest() is just -max().  For integer types, we can use min()
  // instead.
  const KeyType initMaxKey = ::Kokkos::Details::ArithTraits<KeyType>::is_integer ?
    ::Kokkos::Details::ArithTraits<KeyType>::min () :
    -::Kokkos::Details::ArithTraits<KeyType>::max ();
  this->init (keys, startingValue, initMinKey, initMaxKey,
              initMinKey, initMinKey, false);

#ifdef HAVE_TPETRA_DEBUG
  check ();
#endif // HAVE_TPETRA_DEBUG
}

template<class KeyType, class ValueType, class DeviceType>
FixedHashTable<KeyType, ValueType, DeviceType>::
FixedHashTable (const Teuchos::ArrayView<const KeyType>& keys,
                const Teuchos::ArrayView<const ValueType>& vals) :
  minKey_ (::Kokkos::Details::ArithTraits<KeyType>::max ()),
  maxKey_ (::Kokkos::Details::ArithTraits<KeyType>::is_integer ?
           ::Kokkos::Details::ArithTraits<KeyType>::min () :
           -::Kokkos::Details::ArithTraits<KeyType>::max ()),
  minVal_ (::Kokkos::Details::ArithTraits<ValueType>::max ()),
  maxVal_ (::Kokkos::Details::ArithTraits<ValueType>::is_integer ?
           ::Kokkos::Details::ArithTraits<ValueType>::min () :
           -::Kokkos::Details::ArithTraits<ValueType>::max ()),
  firstContigKey_ (::Kokkos::Details::ArithTraits<KeyType>::max ()),
  lastContigKey_ (::Kokkos::Details::ArithTraits<KeyType>::is_integer ?
                  ::Kokkos::Details::ArithTraits<KeyType>::min () :
                  -::Kokkos::Details::ArithTraits<KeyType>::max ()),
  contiguousValues_ (false),
  checkedForDuplicateKeys_ (false),
  hasDuplicateKeys_ (false) // to revise in hasDuplicateKeys()
{
  // mfh 01 May 2015: I don't trust that
  // Teuchos::ArrayView::getRawPtr() returns NULL when the size is 0,
  // so I ensure this manually.
  host_input_keys_type keys_k (keys.size () == 0 ? NULL : keys.getRawPtr (),
                               keys.size ());
  host_input_vals_type vals_k (vals.size () == 0 ? NULL : vals.getRawPtr (),
                               vals.size ());
  const KeyType initMinKey = ::Kokkos::Details::ArithTraits<KeyType>::max ();
  // min() for a floating-point type returns the minimum _positive_
  // normalized value.  This is different than for integer types.
  // lowest() is new in C++11 and returns the least value, always
  // negative for signed finite types.
  //
  // mfh 23 May 2015: I have heard reports that
  // std::numeric_limits<int>::lowest() does not exist with the Intel
  // compiler.  I'm not sure if the users in question actually enabled
  // C++11.  However, it's easy enough to work around this issue.  The
  // standard floating-point types are signed and have a sign bit, so
  // lowest() is just -max().  For integer types, we can use min()
  // instead.
  const KeyType initMaxKey = ::Kokkos::Details::ArithTraits<KeyType>::is_integer ?
    ::Kokkos::Details::ArithTraits<KeyType>::min () :
    -::Kokkos::Details::ArithTraits<KeyType>::max ();
  this->init (keys_k, vals_k, initMinKey, initMaxKey);

#ifdef HAVE_TPETRA_DEBUG
  check ();
#endif // HAVE_TPETRA_DEBUG
}

template<class KeyType, class ValueType, class DeviceType>
void
FixedHashTable<KeyType, ValueType, DeviceType>::
init (const keys_type& keys,
      ValueType startingValue,
      KeyType initMinKey,
      KeyType initMaxKey,
      KeyType firstContigKey,
      KeyType lastContigKey,
      const bool computeInitContigKeys)
{
  using Kokkos::subview;
  using Kokkos::ViewAllocateWithoutInitializing;
  using Teuchos::TypeNameTraits;
  typedef typename std::decay<decltype (keys.extent (0)) >::type size_type;
  const char prefix[] = "Tpetra::Details::FixedHashTable: ";

  const offset_type numKeys = static_cast<offset_type> (keys.extent (0));
  {
    const offset_type theMaxVal = ::Kokkos::Details::ArithTraits<offset_type>::max ();
    const size_type maxValST = static_cast<size_type> (theMaxVal);
    TEUCHOS_TEST_FOR_EXCEPTION
      (keys.extent (0) > maxValST, std::invalid_argument, prefix << "The "
       "number of keys " << keys.extent (0) << " does not fit in "
       "offset_type = " << TypeNameTraits<offset_type>::name () << ", whose "
       "max value is " << theMaxVal << ".  This means that it is not possible to "
       "use this constructor.");
  }
  TEUCHOS_TEST_FOR_EXCEPTION
    (static_cast<unsigned long long> (numKeys) >
     static_cast<unsigned long long> (::Kokkos::Details::ArithTraits<ValueType>::max ()),
     std::invalid_argument, "Tpetra::Details::FixedHashTable: The number of "
     "keys " << numKeys << " is greater than the maximum representable "
     "ValueType value " << ::Kokkos::Details::ArithTraits<ValueType>::max () << ".  "
     "This means that it is not possible to use this constructor.");
  TEUCHOS_TEST_FOR_EXCEPTION
    (numKeys > static_cast<offset_type> (INT_MAX), std::logic_error, prefix <<
     "This class currently only works when the number of keys is <= INT_MAX = "
     << INT_MAX << ".  If this is a problem for you, please talk to the Tpetra "
     "developers.");

  const bool buildInParallel =
    FHT::worthBuildingFixedHashTableInParallel<execution_space> ();

  // NOTE (mfh 14 May 2015) This method currently assumes UVM.  We
  // could change that by setting up ptr and val as Kokkos::DualView
  // instances.  If we do that, since we are filling on host for now,
  // we want to make sure that we only zero-fill ptr on host
  // initially, and that we don't fill val at all.  Once we finish
  // Kokkos-izing all the set-up kernels, we won't need DualView for
  // either ptr or val.

  if (computeInitContigKeys) {
    // Find the first and last initial contiguous keys.  If we find a
    // long sequence of initial contiguous keys, we can save space by
    // not storing them explicitly as pairs in the hash table.
    //
    // NOTE (mfh 01 Jun 2015) Doing this in parallel requires a scan
    // ("min index such that the difference between the current key and
    // the next != 1"), which takes multiple passes over the data.  We
    // could fuse it with CountBuckets (only update counts on 'final'
    // pass).  However, we're really just moving this sequential search
    // out of Map's constructor here, so there is no loss in doing it
    // sequentially for now.  Later, we can work on parallelization.
    //
    // NOTE (mfh 01 Jun 2015, 28 Mar 2016) The code below assumes UVM.
    // It is rational to assume UVM here, because this is "sparse"
    // access -- we might not need to look at all the entries of keys,
    // so it doesn't necessarily make sense to copy the whole thing
    // back to host.
    if (numKeys > 0) {
      firstContigKey_ = keys[0];
      // Start with one plus, then decrement at the end.  That lets us do
      // only one addition per loop iteration, rather than two (if we test
      // against lastContigKey + 1 and then increment lastContigKey).
      lastContigKey_ = firstContigKey_ + 1;

      // We will only store keys in the table that are not part of the
      // initial contiguous sequence.  It's possible for the initial
      // contiguous sequence to be trivial, which for a nonzero number of
      // keys means that the "sequence" has length 1.
      for (offset_type k = 1; k < numKeys; ++k) {
        if (lastContigKey_ != keys[k]) {
          break;
        }
        ++lastContigKey_;
      }
      --lastContigKey_;
    }
  }
  else {
    firstContigKey_ = firstContigKey;
    lastContigKey_ = lastContigKey;
  }

  offset_type startIndex;
  if (numKeys > 0) {
    initMinKey = std::min (initMinKey, firstContigKey_);
    initMaxKey = std::max (initMaxKey, lastContigKey_);
    startIndex = static_cast<offset_type> (lastContigKey_ - firstContigKey_);
  } else {
    startIndex = 0;
  }

  const offset_type theNumKeys = numKeys - startIndex;
  const offset_type size = hash_type::getRecommendedSize (theNumKeys);
#ifdef HAVE_TPETRA_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(
    size == 0 && numKeys != 0, std::logic_error,
    "Tpetra::Details::FixedHashTable constructor: "
    "getRecommendedSize(" << numKeys << ") returned zero, "
    "even though the number of keys " << numKeys << " is nonzero.  "
    "Please report this bug to the Tpetra developers.");
#endif // HAVE_TPETRA_DEBUG
  keys_type theKeys =
    subview (keys, std::pair<offset_type, offset_type> (startIndex, numKeys));

  // FIXME (mfh 28 Mar 2016) For some reason, we don't seem to need a
  // fence here, but we do before other allocations.

  // The array of counts must be separate from the array of offsets,
  // in order for parallel_scan to work correctly.
  typedef typename ptr_type::non_const_type counts_type;
  counts_type counts ("FixedHashTable::counts", size);

  //
  // Count the number of "buckets" per offsets array (ptr) entry.
  //

  // The Kokkos kernel uses atomic update instructions to count the
  // number of "buckets" per offsets array (ptr) entry.  Atomic
  // updates incur overhead, even in the sequential case.  The Kokkos
  // kernel is still correct in that case, but I would rather not
  // incur overhead then.
  if (buildInParallel) {
    FHT::CountBuckets<counts_type, keys_type> functor (counts, theKeys, size);

    typedef typename counts_type::execution_space execution_space;
    typedef Kokkos::RangePolicy<execution_space, offset_type> range_type;
    Kokkos::parallel_for (range_type (0, theNumKeys), functor);
  }
  else {
    // Access to counts is not necessarily contiguous, but is
    // irregular and likely touches all pages of the array.  Thus, it
    // probably makes sense to use a host copy explicitly, rather than
    // assume UVM.
    auto countsHost = Kokkos::create_mirror_view (counts);
    // FIXME (mfh 28 Mar 2016) Does create_mirror_view zero-fill?
    Kokkos::deep_copy (countsHost, static_cast<offset_type> (0));

    for (offset_type k = 0; k < theNumKeys; ++k) {
      typedef typename hash_type::result_type hash_value_type;
      const hash_value_type hashVal = hash_type::hashFunc (theKeys[k], size);
      ++countsHost[hashVal];
    }
    Kokkos::deep_copy (counts, countsHost);
  }

  // FIXME (mfh 28 Mar 2016) Need a fence here, otherwise SIGSEGV w/
  // CUDA when ptr is filled.
  execution_space::fence ();

  // Kokkos::View fills with zeros by default.
  typename ptr_type::non_const_type ptr ("FixedHashTable::ptr", size+1);

  // Compute row offsets via prefix sum:
  //
  // ptr[i+1] = \sum_{j=0}^{i} counts[j].
  //
  // Thus, ptr[i+1] - ptr[i] = counts[i], so that ptr[i+1] = ptr[i] +
  // counts[i].  If we stored counts[i] in ptr[i+1] on input, then the
  // formula is ptr[i+1] += ptr[i].
  //
  // parallel_scan does not incur overhead with Kokkos::Serial, but
  // with actual parallel execution spaces, it does require multiple
  // passes over the data.  Thus, it still makes sense to have a
  // sequential fall-back.
  if (buildInParallel) {
    ::Tpetra::Details::computeOffsetsFromCounts (ptr, counts);
  }
  else {
    // mfh 28 Mar 2016: We could use UVM here, but it's pretty easy to
    // use a host mirror too, so I'll just do that.
    typename decltype (ptr)::HostMirror ptrHost = Kokkos::create_mirror_view (ptr);

    ptrHost[0] = 0;
    for (offset_type i = 0; i < size; ++i) {
      //ptrHost[i+1] += ptrHost[i];
      ptrHost[i+1] = ptrHost[i] + counts[i];
    }
    Kokkos::deep_copy (ptr, ptrHost);
  }

  // FIXME (mfh 28 Mar 2016) Need a fence here, otherwise SIGSEGV w/
  // CUDA when val is filled.
  execution_space::fence ();

  // Allocate the array of (key,value) pairs.  Don't fill it with
  // zeros, because we will fill it with actual data below.
  using Kokkos::ViewAllocateWithoutInitializing;
  typedef typename val_type::non_const_type nonconst_val_type;
  nonconst_val_type val (ViewAllocateWithoutInitializing ("FixedHashTable::pairs"),
                         theNumKeys);

  // Fill in the hash table's "values" (the (key,value) pairs).
  typedef FHT::FillPairs<typename val_type::non_const_type, keys_type,
    typename ptr_type::non_const_type> functor_type;
  typename functor_type::value_type result (initMinKey, initMaxKey);

  const ValueType newStartingValue = startingValue + static_cast<ValueType> (startIndex);
  if (buildInParallel) {
    functor_type functor (val, counts, ptr, theKeys, newStartingValue,
                          initMinKey, initMaxKey);
    typedef Kokkos::RangePolicy<execution_space, offset_type> range_type;
    Kokkos::parallel_reduce (range_type (0, theNumKeys), functor, result);
  }
  else {
    for (offset_type k = 0; k < theNumKeys; ++k) {
      typedef typename hash_type::result_type hash_value_type;
      const KeyType key = theKeys[k];
      if (key > result.maxKey_) {
        result.maxKey_ = key;
      }
      if (key < result.minKey_) {
        result.minKey_ = key;
      }
      const ValueType theVal = newStartingValue + static_cast<ValueType> (k);
      const hash_value_type hashVal = hash_type::hashFunc (key, size);

      // Return the old count; decrement afterwards.
      //
      // NOTE (mfh 28 Mar 2016) This assumes UVM.  It might make more
      // sense to use a host mirror here, due to the access pattern.
      const offset_type count = counts[hashVal];
      --counts[hashVal];
      if (count == 0) {
        result.success_ = false; // FAILURE!
        break;
      }
      else {
        // NOTE (mfh 28 Mar 2016) This assumes UVM.
        const offset_type curPos = ptr[hashVal+1] - count;

        // NOTE (mfh 28 Mar 2016) This assumes UVM.
        val[curPos].first = key;
        val[curPos].second = theVal;
      }
    }
  }

  // FIXME (mfh 01 Jun 2015) Temporarily commented out because of
  // reports of exceptions being thrown in Albany.
  //
  // TEUCHOS_TEST_FOR_EXCEPTION
  //   (! result.success_, std::logic_error, "Tpetra::Details::FixedHashTable::"
  //    "init: Filling the hash table failed!  Please report this bug to the "
  //    "Tpetra developers.");
  (void) result; // prevent build warning (set but unused variable)

  // "Commit" the computed arrays and other computed quantities.
  ptr_ = ptr;
  val_ = val;
  minKey_ = result.minKey_;
  maxKey_ = result.maxKey_;
  // We've already set firstContigKey_ and lastContigKey_ above.
}


template<class KeyType, class ValueType, class DeviceType>
void
FixedHashTable<KeyType, ValueType, DeviceType>::
init (const host_input_keys_type& keys,
      const host_input_vals_type& vals,
      KeyType initMinKey,
      KeyType initMaxKey)
{
  const offset_type numKeys = static_cast<offset_type> (keys.extent (0));
  TEUCHOS_TEST_FOR_EXCEPTION
    (static_cast<unsigned long long> (numKeys) > static_cast<unsigned long long> (::Kokkos::Details::ArithTraits<ValueType>::max ()),
     std::invalid_argument, "Tpetra::Details::FixedHashTable: The number of "
     "keys " << numKeys << " is greater than the maximum representable "
     "ValueType value " << ::Kokkos::Details::ArithTraits<ValueType>::max () << ".");
  TEUCHOS_TEST_FOR_EXCEPTION
    (numKeys > static_cast<offset_type> (INT_MAX), std::logic_error, "Tpetra::"
     "Details::FixedHashTable: This class currently only works when the number "
     "of keys is <= INT_MAX = " << INT_MAX << ".  If this is a problem for you"
     ", please talk to the Tpetra developers.");

  // There's no need to find the first and last initial contiguous
  // keys in this case, because if we reach this init() function, then
  // hasContiguousValues() is false and so get() doesn't use the
  // initial contiguous sequence of keys.

  const offset_type size = hash_type::getRecommendedSize (numKeys);
#ifdef HAVE_TPETRA_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(
    size == 0 && numKeys != 0, std::logic_error,
    "Tpetra::Details::FixedHashTable constructor: "
    "getRecommendedSize(" << numKeys << ") returned zero, "
    "even though the number of keys " << numKeys << " is nonzero.  "
    "Please report this bug to the Tpetra developers.");
#endif // HAVE_TPETRA_DEBUG

  // NOTE (mfh 14 May 2015) This method currently assumes UVM.  We
  // could change that by setting up ptr and val as Kokkos::DualView
  // instances.  If we do that, since we are filling on host for now,
  // we want to make sure that we only zero-fill ptr on host
  // initially, and that we don't fill val at all.  Once we finish
  // Kokkos-izing all the set-up kernels, we won't need DualView for
  // either ptr or val.

  typename ptr_type::non_const_type ptr ("FixedHashTable::ptr", size + 1);

  // Allocate the array of key,value pairs.  Don't waste time filling
  // it with zeros, because we will fill it with actual data below.
  using Kokkos::ViewAllocateWithoutInitializing;
  typedef typename val_type::non_const_type nonconst_val_type;
  nonconst_val_type val (ViewAllocateWithoutInitializing ("FixedHashTable::pairs"),
                         numKeys);

  // Compute number of entries in each hash table position.
  for (offset_type k = 0; k < numKeys; ++k) {
    const typename hash_type::result_type hashVal =
      hash_type::hashFunc (keys[k], size);
    // Shift over one, so that counts[j] = ptr[j+1].  See below.
    ++ptr[hashVal+1];
  }

  // Compute row offsets via prefix sum:
  //
  // ptr[i+1] = \sum_{j=0}^{i} counts[j].
  //
  // Thus, ptr[i+1] - ptr[i] = counts[i], so that ptr[i+1] = ptr[i] +
  // counts[i].  If we stored counts[i] in ptr[i+1] on input, then the
  // formula is ptr[i+1] += ptr[i].
  for (offset_type i = 0; i < size; ++i) {
    ptr[i+1] += ptr[i];
  }
  //ptr[0] = 0; // We've already done this when initializing ptr above.

  // curRowStart[i] is the offset of the next element in row i.
  typename ptr_type::non_const_type curRowStart ("FixedHashTable::curRowStart", size);

  // Fill in the hash table.
  FHT::FillPairsResult<KeyType> result (initMinKey, initMaxKey);
  for (offset_type k = 0; k < numKeys; ++k) {
    typedef typename hash_type::result_type hash_value_type;
    const KeyType key = keys[k];
    if (key > result.maxKey_) {
      result.maxKey_ = key;
    }
    if (key < result.minKey_) {
      result.minKey_ = key;
    }
    const ValueType theVal = vals[k];
    if (theVal > maxVal_) {
      maxVal_ = theVal;
    }
    if (theVal < minVal_) {
      minVal_ = theVal;
    }
    const hash_value_type hashVal = hash_type::hashFunc (key, size);

    const offset_type offset = curRowStart[hashVal];
    const offset_type curPos = ptr[hashVal] + offset;
    if (curPos >= ptr[hashVal+1]) {
      result.success_ = false; // FAILURE!
    }
    else {
      val[curPos].first = key;
      val[curPos].second = theVal;
      ++curRowStart[hashVal];
    }
  }

  TEUCHOS_TEST_FOR_EXCEPTION
    (! result.success_, std::logic_error, "Tpetra::Details::FixedHashTable::"
     "init: Filling the hash table failed!  Please report this bug to the "
     "Tpetra developers.");

  // "Commit" the computed arrays and other computed quantities.
  ptr_ = ptr;
  val_ = val;
  minKey_ = result.minKey_;
  maxKey_ = result.maxKey_;
  // We've already assigned to minVal_ and maxVal_ above.
}

template <class KeyType, class ValueType, class DeviceType>
bool
FixedHashTable<KeyType, ValueType, DeviceType>::
hasDuplicateKeys ()
{
  if (! checkedForDuplicateKeys_) {
    hasDuplicateKeys_ = checkForDuplicateKeys ();
    checkedForDuplicateKeys_ = true;
  }
  return hasDuplicateKeys_;
}

template <class KeyType, class ValueType, class DeviceType>
bool
FixedHashTable<KeyType, ValueType, DeviceType>::
checkForDuplicateKeys () const
{
  const offset_type size = this->getSize ();
  // It's allowed for the hash table to have a positive number of
  // buckets (getSize()), but still store no entries (numPairs()).
  // Both cases trivially entail no duplicates.
  if (size == 0 || this->numPairs () == 0) {
    return false;
  }
  else {
    typedef FHT::CheckForDuplicateKeys<ptr_type, val_type> functor_type;
    functor_type functor (val_, ptr_);
    int hasDupKeys = 0;
    typedef Kokkos::RangePolicy<execution_space, offset_type> range_type;
    Kokkos::parallel_reduce (range_type (0, size), functor, hasDupKeys);
    return hasDupKeys > 0;
  }
}

template <class KeyType, class ValueType, class DeviceType>
std::string
FixedHashTable<KeyType, ValueType, DeviceType>::
description () const
{
  std::ostringstream oss;
  oss << "FixedHashTable<"
      << Teuchos::TypeNameTraits<KeyType>::name () << ","
      << Teuchos::TypeNameTraits<ValueType>::name () << ">: "
      << "{ numKeys: " << val_.extent (0)
      << ", tableSize: " << this->getSize () << " }";
  return oss.str();
}

template <class KeyType, class ValueType, class DeviceType>
void
FixedHashTable<KeyType, ValueType, DeviceType>::
describe (Teuchos::FancyOStream& out,
          const Teuchos::EVerbosityLevel verbLevel) const
{
  using std::endl;
  using std::setw;
  using Teuchos::OSTab;
  using Teuchos::rcpFromRef;
  using Teuchos::TypeNameTraits;
  using Teuchos::VERB_DEFAULT;
  using Teuchos::VERB_NONE;
  using Teuchos::VERB_LOW;
  using Teuchos::VERB_EXTREME;

  // NOTE (mfh 14 May 2015) This method currently assumes UVM for
  // access to ptr_ and val_ from the host.

  Teuchos::EVerbosityLevel vl = verbLevel;
  if (vl == VERB_DEFAULT) vl = VERB_LOW;

  if (vl == VERB_NONE) {
    // do nothing
  }
  else if (vl == VERB_LOW) {
    out << this->description() << endl;
  }
  else {  // MEDIUM, HIGH or EXTREME
    out << "FixedHashTable:" << endl;
    {
      OSTab tab1 (rcpFromRef (out));

      // const std::string label = this->getObjectLabel ();
      // if (label != "") {
      //   out << "label: " << label << endl;
      // }
      out << "Template parameters:" << endl;
      {
        OSTab tab2 (rcpFromRef (out));
        out << "KeyType: " << TypeNameTraits<KeyType>::name () << endl
            << "ValueType: " << TypeNameTraits<ValueType>::name () << endl;
      }

      const offset_type tableSize = this->getSize ();
      const offset_type numKeys = val_.extent (0);

      out << "Table parameters:" << endl;
      {
        OSTab tab2 (rcpFromRef (out));
        out << "numKeys: " << numKeys << endl
            << "tableSize: " << tableSize << endl;
      }

      if (vl >= VERB_EXTREME) {
        out << "Contents: ";
        if (tableSize == 0 || numKeys == 0) {
          out << "[]" << endl;
        } else {
          out << "[ " << endl;
          {
            OSTab tab2 (rcpFromRef (out));
            for (offset_type i = 0; i < tableSize; ++i) {
              OSTab tab3 (rcpFromRef (out));
              out << "[";
              for (offset_type k = ptr_[i]; k < ptr_[i+1]; ++k) {
                out << "(" << val_[k].first << "," << val_[k].second << ")";
                if (k + 1 < ptr_[i+1]) {
                  out << ", ";
                }
              }
              out << "]" << endl;
            } // for each table position i
          }
          out << "]" << endl;
        } // The table contains entries
      } // vl >= VERB_EXTREME
    }
    out << endl;
  } // if vl > VERB_LOW
}

} // namespace Details
} // namespace Tpetra

// Macro that explicitly instantiates FixedHashTable for the given local
// ordinal (LO) and global ordinal (GO) types.  Note that FixedHashTable's
// template parameters occur in the opposite order of most Tpetra
// classes.  This is because FixedHashTable performs global-to-local
// lookup, and the convention in templated C++ lookup tables (such as
// std::map) is <KeyType, ValueType>.
//
// This macro must be explanded within the Tpetra::Details namespace.
#define TPETRA_DETAILS_FIXEDHASHTABLE_INSTANT_DEFAULTNODE(LO,GO) \
  template class FixedHashTable< GO , LO >;

// Macro that explicitly instantiates FixedHashTable for the given
// local ordinal (LO), global ordinal (GO), and Kokkos device (DEVICE)
// types.  Note that FixedHashTable's first two template parameters
// occur in the opposite order of most Tpetra classes.  This is
// because FixedHashTable performs global-to-local lookup, and the
// convention in templated C++ lookup tables (such as std::map) is
// <KeyType, ValueType>.
//
// This macro must be explanded within the Tpetra::Details namespace.
#define TPETRA_DETAILS_FIXEDHASHTABLE_INSTANT(LO, GO, DEVICE) \
  template class FixedHashTable< GO , LO , DEVICE >;

#endif // TPETRA_DETAILS_FIXEDHASHTABLE_DEF_HPP
