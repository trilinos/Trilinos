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
#ifdef TPETRA_USE_MURMUR_HASH
#  include <Kokkos_Functional.hpp> // hash function used by Kokkos::UnorderedMap
#endif // TPETRA_USE_MURMUR_HASH
#include <limits> // std::numeric_limits
#include <type_traits>


//
// This anonymous namespace stores utility functions and Kokkos
// functors for use in FixedHashTable construction.
//
namespace { // (anonymous)

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
//
// FixedHashTable code should call the nonmember function below, that
// has the same name but starts with a lower-case w.
template<class ExecSpace>
struct WorthBuildingFixedHashTableInParallel {
  typedef typename ExecSpace::execution_space execution_space;

  static bool isWorth () {
    // NOTE: Kokkos::Cuda does NOT have this method.  That's why we
    // need the partial specialization below.
    return execution_space::max_hardware_threads () > 1;
  }
};

#ifdef KOKKOS_HAVE_CUDA
template<>
struct WorthBuildingFixedHashTableInParallel<Kokkos::Cuda> {
  // There could be more complicated expressions for whether this is
  // actually worthwhile, but for now I'll just say that with Cuda, we
  // will ALWAYS count buckets in parallel (that is, run a Kokkos
  // parallel kernel).
  static bool isWorth () {
    return true;
  }
};
#endif // KOKKOS_HAVE_CUDA

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
  return WorthBuildingFixedHashTableInParallel<ExecSpace>::isWorth ();
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
            src.dimension_0 ());
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
    Kokkos::atomic_fetch_add (&counts_[hashVal], 1);
  }

private:
  //! Bucket counts (output argument).
  counts_view_type counts_;
  //! Keys for the FixedHashTable to construct (input argument).
  keys_view_type keys_;
  //! Number of buckets plus 1 (or 0, if no buckets).
  size_type size_;
};

/// \brief Parallel scan functor for computing "row" offsets.
///
/// Kokkos::parallel_scan functor for computing the row offsets array
/// from the array of counts (which the above functor CountBuckets
/// computes).
///
/// \tparam OffsetsViewType Type of the Kokkos::View specialization
///   used to store the "row" offsets; the output of this functor.
/// \tparam SizeType The parallel loop index type; a built-in integer
///   type.  Defaults to the type of the input View's dimension.  You
///   may use a shorter type to improve performance.
template<class OffsetsViewType,
         class SizeType = typename OffsetsViewType::size_type>
class ComputeRowOffsets {
public:
  typedef OffsetsViewType offsets_view_type;
  typedef typename OffsetsViewType::const_type counts_view_type;
  typedef typename OffsetsViewType::execution_space execution_space;
  typedef typename OffsetsViewType::memory_space memory_space;
  typedef SizeType size_type;
  typedef typename OffsetsViewType::non_const_value_type value_type;

  /// \brief Constructor
  ///
  /// \param offsets [out] (Preallocated) offsets; one entry longer
  ///   than \c counts
  /// \param counts [in] View of bucket counts
  ComputeRowOffsets (const offsets_view_type& offsets,
                     const counts_view_type& counts) :
    offsets_ (offsets),
    counts_ (counts),
    size_ (counts.dimension_0 ())
  {}

  //! Set the initial value of the reduction result.
  KOKKOS_INLINE_FUNCTION void init (value_type& dst) const
  {
    dst = 0;
  }

  KOKKOS_INLINE_FUNCTION void
  join (volatile value_type& dst,
        const volatile value_type& src) const
  {
    dst += src;
  }

  KOKKOS_INLINE_FUNCTION void
  operator () (const size_type& i, value_type& update, const bool final) const
  {
    if (final) {
      offsets_[i] = update;
    }
    if (i < size_) {
      update += counts_[i];
    }
  }

private:
  //! Offsets (output argument)
  offsets_view_type offsets_;
  //! Bucket counts (input argument).
  counts_view_type counts_;
  //! Number of entries in counts_.
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
  FillPairsResult () :
    minKey_ (std::numeric_limits<KeyType>::max ()),
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
    maxKey_ (std::numeric_limits<KeyType>::is_integer ?
             std::numeric_limits<KeyType>::min () :
             -std::numeric_limits<KeyType>::max ()),
    success_ (true)
  {
    static_assert (std::is_arithmetic<KeyType>::value, "FillPairsResult: "
                   "KeyType must be some kind of number type.");
  }

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
/// ComputeRowOffsets), but checking for it is cheap and easy.
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

  typedef FillPairsResult<typename pair_type::second_type> value_type;

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
    size_ (counts.dimension_0 ()),
    startingValue_ (startingValue),
    initMinKey_ (std::numeric_limits<key_type>::max ()),
    initMaxKey_ (std::numeric_limits<key_type>::is_integer ?
                 std::numeric_limits<key_type>::min () :
                 -std::numeric_limits<key_type>::max ())
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
    size_ (counts.dimension_0 ()),
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
    const offset_type count = Kokkos::atomic_fetch_add (&counts_[hashVal], -1);
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
  typedef bool value_type;

  /// \brief Constructor
  ///
  /// \param pairs [in] View of the FixedHashTable's (key,value) pairs
  /// \param ptr [in] View of the FixedHashTable's "bucket" offsets
  CheckForDuplicateKeys (const pairs_view_type& pairs,
                         const offsets_view_type& ptr) :
    pairs_ (pairs),
    ptr_ (ptr),
    size_ (ptr_.dimension_0 () == 0 ?
           size_type (0) :
           ptr_.dimension_0 () - 1)
  {}

  //! Set the initial value of the reduction result.
  KOKKOS_INLINE_FUNCTION void init (value_type& dst) const
  {
    dst = false;
  }

  //! Combine two intermediate reduction results.
  KOKKOS_INLINE_FUNCTION void
  join (volatile value_type& dst,
        const volatile value_type& src) const
  {
    dst = dst || src;
  }

  //! Parallel loop body.
  KOKKOS_INLINE_FUNCTION void
  operator () (const size_type& i, value_type& dst) const
  {
    typedef typename offsets_view_type::non_const_value_type offset_type;
    typedef typename pairs_view_type::non_const_value_type pair_type;
    typedef typename pair_type::first_type key_type;

    if (dst) {
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
      dst = dst || foundDuplicateKey;
    }
  }

private:
  pairs_view_type pairs_;
  offsets_view_type ptr_;
  size_type size_;
};

} // namespace (anonymous)

//
// Here begins the actual implementation of FixedHashTable.
//

namespace Tpetra {
namespace Details {

template<class KeyType, class ValueType, class DeviceType>
void
FixedHashTable<KeyType, ValueType, DeviceType>::
check () const
{
  const char prefix[] = "Tpetra::Details::FixedHashTable: ";
  const char suffix[] = "  Please report this bug to the Tpetra developers.";

#if ! defined(TPETRA_HAVE_KOKKOS_REFACTOR)
  TEUCHOS_TEST_FOR_EXCEPTION
    (ptr_.ptr_on_device () != rawPtr_, std::logic_error,
     prefix << "ptr_.ptr_on_device () != rawPtr_." << suffix);
  TEUCHOS_TEST_FOR_EXCEPTION
    (val_.ptr_on_device () != rawVal_, std::logic_error,
     prefix << "val_.ptr_on_device () != rawVal_." << suffix);
#endif // ! defined(TPETRA_HAVE_KOKKOS_REFACTOR)

  TEUCHOS_TEST_FOR_EXCEPTION
    (invalidValue_ != Teuchos::OrdinalTraits<ValueType>::invalid (),
     std::logic_error, prefix << "invalidValue_ == " << invalidValue_
     << " != Teuchos::OrdinalTraits<ValueType>::invalid() == "
     << Teuchos::OrdinalTraits<ValueType>::invalid () << "." << suffix);
}

template<class KeyType, class ValueType, class DeviceType>
FixedHashTable<KeyType, ValueType, DeviceType>::
FixedHashTable () :
#if ! defined(TPETRA_HAVE_KOKKOS_REFACTOR)
  rawPtr_ (NULL),
  rawVal_ (NULL),
#endif // ! defined(TPETRA_HAVE_KOKKOS_REFACTOR)
  invalidValue_ (Teuchos::OrdinalTraits<ValueType>::invalid ()),
  checkedForDuplicateKeys_ (true), // it's an empty table; no need to check
  hasDuplicateKeys_ (false)
{
#ifdef HAVE_TPETRA_DEBUG
  check ();
#endif // HAVE_TPETRA_DEBUG
}

template<class KeyType, class ValueType, class DeviceType>
FixedHashTable<KeyType, ValueType, DeviceType>::
FixedHashTable (const Teuchos::ArrayView<const KeyType>& keys) :
#if ! defined(TPETRA_HAVE_KOKKOS_REFACTOR)
  rawPtr_ (NULL),
  rawVal_ (NULL),
#endif // ! defined(TPETRA_HAVE_KOKKOS_REFACTOR)
  invalidValue_ (Teuchos::OrdinalTraits<ValueType>::invalid ()),
  checkedForDuplicateKeys_ (false),
  hasDuplicateKeys_ (false) // to revise in hasDuplicateKeys()
{
  // mfh 01 May 2015: I don't trust that
  // Teuchos::ArrayView::getRawPtr() returns NULL when the size is 0,
  // so I ensure this manually.
  const ValueType startingValue = static_cast<ValueType> (0);
  host_input_keys_type keys_k (keys.size () == 0 ? NULL : keys.getRawPtr (),
                               keys.size ());
  const KeyType initMinKey = std::numeric_limits<KeyType>::max ();
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
  const KeyType initMaxKey = std::numeric_limits<KeyType>::is_integer ?
    std::numeric_limits<KeyType>::min () :
    -std::numeric_limits<KeyType>::max ();
  this->init (keys_k, startingValue, initMinKey, initMaxKey);

#ifdef HAVE_TPETRA_DEBUG
  check ();
#endif // HAVE_TPETRA_DEBUG
}

template<class KeyType, class ValueType, class DeviceType>
FixedHashTable<KeyType, ValueType, DeviceType>::
FixedHashTable (const Teuchos::ArrayView<const KeyType>& keys,
                const ValueType startingValue) :
#if ! defined(TPETRA_HAVE_KOKKOS_REFACTOR)
  rawPtr_ (NULL),
  rawVal_ (NULL),
#endif // ! defined(TPETRA_HAVE_KOKKOS_REFACTOR)
  invalidValue_ (Teuchos::OrdinalTraits<ValueType>::invalid ()),
  checkedForDuplicateKeys_ (false),
  hasDuplicateKeys_ (false) // to revise in hasDuplicateKeys()
{
  // mfh 01 May 2015: I don't trust that
  // Teuchos::ArrayView::getRawPtr() returns NULL when the size is 0,
  // so I ensure this manually.
  host_input_keys_type keys_k (keys.size () == 0 ? NULL : keys.getRawPtr (),
                               keys.size ());
  const KeyType initMinKey = std::numeric_limits<KeyType>::max ();
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
  const KeyType initMaxKey = std::numeric_limits<KeyType>::is_integer ?
    std::numeric_limits<KeyType>::min () :
    -std::numeric_limits<KeyType>::max ();
  this->init (keys_k, startingValue, initMinKey, initMaxKey);

#ifdef HAVE_TPETRA_DEBUG
  check ();
#endif // HAVE_TPETRA_DEBUG
}

template<class KeyType, class ValueType, class DeviceType>
FixedHashTable<KeyType, ValueType, DeviceType>::
FixedHashTable (const Teuchos::ArrayView<const KeyType>& keys,
                const Teuchos::ArrayView<const ValueType>& vals) :
#if ! defined(TPETRA_HAVE_KOKKOS_REFACTOR)
  rawPtr_ (NULL),
  rawVal_ (NULL),
#endif // ! defined(TPETRA_HAVE_KOKKOS_REFACTOR)
  invalidValue_ (Teuchos::OrdinalTraits<ValueType>::invalid ()),
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
  const KeyType initMinKey = std::numeric_limits<KeyType>::max ();
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
  const KeyType initMaxKey = std::numeric_limits<KeyType>::is_integer ?
    std::numeric_limits<KeyType>::min () :
    -std::numeric_limits<KeyType>::max ();
  this->init (keys_k, vals_k, initMinKey, initMaxKey);

#ifdef HAVE_TPETRA_DEBUG
  check ();
#endif // HAVE_TPETRA_DEBUG
}

template<class KeyType, class ValueType, class DeviceType>
void
FixedHashTable<KeyType, ValueType, DeviceType>::
init (const host_input_keys_type& keys,
      const ValueType startingValue,
      const KeyType initMinKey,
      const KeyType initMaxKey)
{
  const offset_type numKeys = static_cast<offset_type> (keys.dimension_0 ());
  TEUCHOS_TEST_FOR_EXCEPTION
    (numKeys > static_cast<offset_type> (INT_MAX), std::logic_error, "Tpetra::"
     "Details::FixedHashTable: This class currently only works when the number "
     "of keys is <= INT_MAX = " << INT_MAX << ".  If this is a problem for you"
     ", please talk to the Tpetra developers.");

  const offset_type size = hash_type::getRecommendedSize (numKeys);
#ifdef HAVE_TPETRA_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(
    size == 0 && numKeys != 0, std::logic_error,
    "Tpetra::Details::FixedHashTable constructor: "
    "getRecommendedSize(" << numKeys << ") returned zero, "
    "even though the number of keys " << numKeys << " is nonzero.  "
    "Please report this bug to the Tpetra developers.");
#endif // HAVE_TPETRA_DEBUG

  const bool buildInParallel =
    worthBuildingFixedHashTableInParallel<execution_space> ();

  // NOTE (mfh 14 May 2015) This method currently assumes UVM.  We
  // could change that by setting up ptr and val as Kokkos::DualView
  // instances.  If we do that, since we are filling on host for now,
  // we want to make sure that we only zero-fill ptr on host
  // initially, and that we don't fill val at all.  Once we finish
  // Kokkos-izing all the set-up kernels, we won't need DualView for
  // either ptr or val.

  // The array of counts must be separate from the array of offsets,
  // in order for parallel_scan to work correctly.
  typedef typename ptr_type::non_const_type counts_type;
  counts_type counts ("counts", size);

  // Only make a device copy of the input array 'keys' if the input
  // array lives in a different memory space.  Remember that with UVM,
  // host code can access CUDA device memory, but not the other way
  // around.
  typedef DeepCopyIfNeeded<KeyType, typename host_input_keys_type::array_layout,
                           typename host_input_keys_type::execution_space,
                           execution_space> copier_type;
  typedef typename copier_type::output_view_type keys_type;
  keys_type keys_d = copier_type::copy (keys);

  //
  // Count the number of "buckets" per offsets array (ptr) entry.
  //

  // The Kokkos kernel uses atomic update instructions to count the
  // number of "buckets" per offsets array (ptr) entry.  Atomic
  // updates incur overhead, even in the sequential case.  The Kokkos
  // kernel is still correct in that case, but I would rather not
  // incur overhead then.
  if (buildInParallel) {
    for (offset_type k = 0; k < numKeys; ++k) {
      typedef typename hash_type::result_type hash_value_type;

      const hash_value_type hashVal = hash_type::hashFunc (keys[k], size);
      ++counts[hashVal];
    }
  }
  else {
    CountBuckets<counts_type, keys_type> functor (counts, keys_d, size);
    Kokkos::parallel_for (numKeys, functor);
  }

  // Kokkos::View fills with zeros by default.
  typename ptr_type::non_const_type ptr ("ptr", size + 1);

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
    typedef ComputeRowOffsets<typename ptr_type::non_const_type> functor_type;
    functor_type functor (ptr, counts);
    Kokkos::parallel_scan (size+1, functor);
  }
  else {
    for (offset_type i = 0; i < size; ++i) {
      //ptr[i+1] += ptr[i];
      ptr[i+1] = ptr[i] + counts[i];
    }
    //ptr[0] = 0; // We've already done this when initializing ptr above.
  }

  // Allocate the array of (key,value) pairs.  Don't fill it with
  // zeros, because we will fill it with actual data below.
  typename val_type::non_const_type val (Kokkos::ViewAllocateWithoutInitializing ("val"), numKeys);

  // Fill in the hash table's "values" (the (key,value) pairs).
  typedef FillPairs<typename val_type::non_const_type, keys_type,
                    typename ptr_type::non_const_type> functor_type;
  typename functor_type::value_type result (initMinKey, initMaxKey);
  if (buildInParallel) {
    functor_type functor (val, counts, ptr, keys_d, startingValue,
                          initMinKey, initMaxKey);
    Kokkos::parallel_reduce (numKeys, functor, result);
  }
  else {
    for (offset_type k = 0; k < numKeys; ++k) {
      typedef typename hash_type::result_type hash_value_type;
      const KeyType key = keys[k];
      if (key > result.maxKey_) {
        result.maxKey_ = key;
      }
      if (key < result.minKey_) {
        result.minKey_ = key;
      }
      const ValueType theVal = startingValue + static_cast<ValueType> (k);
      const hash_value_type hashVal = hash_type::hashFunc (key, size);

      // Return the old count; decrement afterwards.
      const offset_type count = (counts[hashVal])--;
      if (count == 0) {
        result.success_ = false; // FAILURE!
        break;
      }
      else {
        const offset_type curPos = ptr[hashVal+1] - count;

        val[curPos].first = key;
        val[curPos].second = theVal;
      }
    }
  }

  TEUCHOS_TEST_FOR_EXCEPTION
    (! result.success_, std::logic_error, "Tpetra::Details::FixedHashTable::"
     "init: Filling the hash table failed!  Please report this bug to the "
     "Tpetra developers.");

  // "Commit" the computed arrays.
  ptr_ = ptr;
  val_ = val;
#if ! defined(TPETRA_HAVE_KOKKOS_REFACTOR)
  rawPtr_ = ptr.ptr_on_device ();
  rawVal_ = val.ptr_on_device ();
#endif // ! defined(TPETRA_HAVE_KOKKOS_REFACTOR)
}


template<class KeyType, class ValueType, class DeviceType>
void
FixedHashTable<KeyType, ValueType, DeviceType>::
init (const host_input_keys_type& keys,
      const host_input_vals_type& vals,
      const KeyType initMinKey,
      const KeyType initMaxKey)
{
  const offset_type numKeys = static_cast<offset_type> (keys.dimension_0 ());
  TEUCHOS_TEST_FOR_EXCEPTION
    (numKeys > static_cast<offset_type> (INT_MAX), std::logic_error, "Tpetra::"
     "Details::FixedHashTable: This class currently only works when the number "
     "of keys is <= INT_MAX = " << INT_MAX << ".  If this is a problem for you"
     ", please talk to the Tpetra developers.");

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

  typename ptr_type::non_const_type ptr ("ptr", size + 1);

  // Allocate the array of key,value pairs.  Don't waste time filling
  // it with zeros, because we will fill it with actual data below.
  typename val_type::non_const_type val (Kokkos::ViewAllocateWithoutInitializing ("val"), numKeys);

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
  typename ptr_type::non_const_type curRowStart ("curRowStart", size);

  // Fill in the hash table.
  FillPairsResult<KeyType> result (initMinKey, initMaxKey);
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

  // "Commit" the computed arrays.
  ptr_ = ptr;
  val_ = val;
#if ! defined(TPETRA_HAVE_KOKKOS_REFACTOR)
  rawPtr_ = ptr.ptr_on_device ();
  rawVal_ = val.ptr_on_device ();
#endif // ! defined(TPETRA_HAVE_KOKKOS_REFACTOR)
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
  if (size == 0) {
    return false;
  }
  else {
    typedef CheckForDuplicateKeys<ptr_type, val_type> functor_type;
    functor_type functor (val_, ptr_);
    bool hasDupKeys = false;
    Kokkos::parallel_reduce (size, functor, hasDupKeys);
    return hasDupKeys;
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
      << "{ numKeys: " << val_.dimension_0 ()
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

      const std::string label = this->getObjectLabel ();
      if (label != "") {
        out << "label: " << label << endl;
      }
      out << "Template parameters:" << endl;
      {
        OSTab tab2 (rcpFromRef (out));
        out << "KeyType: " << TypeNameTraits<KeyType>::name () << endl
            << "ValueType: " << TypeNameTraits<ValueType>::name () << endl;
      }

      const offset_type tableSize = this->getSize ();
      const offset_type numKeys = val_.dimension_0 ();

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
