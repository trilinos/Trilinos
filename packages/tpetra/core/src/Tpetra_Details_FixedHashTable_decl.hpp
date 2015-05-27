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

#ifndef TPETRA_DETAILS_FIXEDHASHTABLE_DECL_HPP
#define TPETRA_DETAILS_FIXEDHASHTABLE_DECL_HPP

#include <Tpetra_Details_Hash.hpp>
#include <Teuchos_Describable.hpp>
#include <Kokkos_Core.hpp>

namespace Tpetra {
namespace Details {

/// \class FixedHashTable
/// \tparam KeyType The type of the hash table's keys.  This must be a
///   built-in signed or unsigned integer type.
/// \tparam ValueType The type of the hash table's values.  This must
///   be a built-in signed or unsigned integer type.
/// \tparam DeviceType Specialization of Kokkos::Device.
///
/// This class implements a look-up table from integer keys to integer
/// values.  All the (key,value) pairs must be added at once, and
/// pairs may not be changed or removed.  Keys and values may have
/// different types.  Tpetra::Map may use this to implement
/// global-to-local index lookup.
///
/// The hash table uses a "compressed sparse row" storage strategy.
/// The hash function maps a key to its "row" in the table, and then
/// we search within that row to find the corresponding value.  In
/// each row, we store a key and its value adjacent to each other.
/// This strategy puts (key,value) pairs in a single contiguous array,
/// rather than in separately allocated buckets (as in a conventional
/// dynamically allocated hash table).  This saves initialization
/// time, as long as the hash function takes less than half the time
/// of a system call to allocate memory.  This is because there are
/// only \f$O(1)\f$ memory allocation calls, rather than one for each
/// (key,value) pair or hash bucket.  The compressed sparse row
/// strategy may also improve locality for hash table lookups.
template<class KeyType,
         class ValueType,
         class DeviceType>
class FixedHashTable : public Teuchos::Describable {
private:
  typedef typename DeviceType::execution_space execution_space;
  typedef typename DeviceType::memory_space memory_space;
  typedef Kokkos::Device<execution_space, memory_space> device_type;

  typedef Hash<KeyType, device_type> hash_type;
  typedef typename hash_type::offset_type offset_type;

  /// \brief Type of the array of hash table "buckets" (a.k.a. "row"
  ///   offsets).
  ///
  /// We specify LayoutLeft explicitly so that the layout is the same
  /// on all Kokkos devices.  It's a 1-D View so LayoutLeft and
  /// LayoutRight mean the same thing, but specifying the layout
  /// explicitly makes Kokkos::deep_copy work.
  typedef typename Kokkos::View<const offset_type*, Kokkos::LayoutLeft,
                                device_type> ptr_type;
  /// \brief Type of the array of (key, value) pairs in the hash table.
  ///
  /// We specify LayoutLeft explicitly so that the layout is the same
  /// on all Kokkos devices.  It's a 1-D View so LayoutLeft and
  /// LayoutRight mean the same thing, but specifying the layout
  /// explicitly makes Kokkos::deep_copy work.
  typedef typename Kokkos::View<const Kokkos::pair<KeyType, ValueType>*,
                                Kokkos::LayoutLeft, device_type> val_type;

public:
  //! Default constructor; makes an empty table.
  FixedHashTable ();

  /// \brief Constructor for arbitrary keys and contiguous values
  ///   starting with zero.
  ///
  /// Add <tt>(keys[i], i)</tt> to the table,
  /// for i = 0, 1, ..., <tt>keys.size()</tt>.
  FixedHashTable (const Teuchos::ArrayView<const KeyType>& keys);

  /// \brief Constructor for arbitrary keys and contiguous values
  ///   starting with \c startingValue.
  ///
  /// Add <tt>(keys[i], startingValue + i)</tt> to the table, for i =
  /// 0, 1, ..., <tt>keys.size()</tt>.  This version is useful if Map
  /// wants to exclude an initial sequence of contiguous GIDs from the
  /// table, and start with a given LID.
  FixedHashTable (const Teuchos::ArrayView<const KeyType>& keys,
                  const ValueType startingValue);

  /// \brief Constructor for arbitrary keys and arbitrary values.
  ///
  /// Add <tt>(keys[i], vals[i])</tt> to the table, for i = 0, 1, ...,
  /// <tt>keys.size()</tt>.  This version is useful for applications
  /// other than Map's GID-to-LID lookup table.
  FixedHashTable (const Teuchos::ArrayView<const KeyType>& keys,
                  const Teuchos::ArrayView<const ValueType>& vals);

  template<class K, class V, class D>
  friend class FixedHashTable;

  /// \brief "Copy" constructor that takes a FixedHashTable with the
  ///   same KeyType and ValueType, but a different DeviceType.
  ///
  /// This constructor makes a deep copy of the input's data if
  /// necessary.  Anything that it doesn't need to deep copy, it
  /// shallow copies.
  template<class InDeviceType>
  FixedHashTable (const FixedHashTable<KeyType, ValueType, InDeviceType>& src,
                  typename std::enable_if<! std::is_same<DeviceType, InDeviceType>::value, int>::type* = NULL)
  {
    using Kokkos::ViewAllocateWithoutInitializing;
    typedef typename ptr_type::non_const_type nonconst_ptr_type;
    typedef typename val_type::non_const_type nonconst_val_type;

    nonconst_ptr_type ptr (ViewAllocateWithoutInitializing ("ptr"),
                           src.ptr_.dimension_0 ());

    // NOTE (mfh 01 May 2015) deep_copy works here, because regardless
    // of the DeviceType, all FixedHashTable types use the same array
    // layout for their internal 1-D Views.
    Kokkos::deep_copy (ptr, src.ptr_);
    nonconst_val_type val (ViewAllocateWithoutInitializing ("val"),
                           src.val_.dimension_0 ());
    Kokkos::deep_copy (val, src.val_);

    this->ptr_ = ptr;
    this->val_ = val;
#if ! defined(TPETRA_HAVE_KOKKOS_REFACTOR)
    this->rawPtr_ = ptr.ptr_on_device ();
    this->rawVal_ = val.ptr_on_device ();
#endif // ! defined(TPETRA_HAVE_KOKKOS_REFACTOR)
    this->invalidValue_ = src.invalidValue_;
    this->minKey_ = src.minKey_;
    this->maxKey_ = src.maxKey_;
    this->checkedForDuplicateKeys_ = src.checkedForDuplicateKeys_;
    this->hasDuplicateKeys_ = src.hasDuplicateKeys_;

#if defined(HAVE_TPETRA_DEBUG)
    this->check ();
#endif // defined(HAVE_TPETRA_DEBUG)
  }

  //! Get the value corresponding to the given key.
  KOKKOS_INLINE_FUNCTION ValueType get (const KeyType& key) const {
    const offset_type size = this->getSize ();
    if (size == 0) {
      return invalidValue_;
    }
    else {
      const typename hash_type::result_type hashVal =
        hash_type::hashFunc (key, size);
#if defined(TPETRA_HAVE_KOKKOS_REFACTOR) || defined(HAVE_TPETRA_DEBUG)
      const offset_type start = ptr_[hashVal];
      const offset_type end = ptr_[hashVal+1];
      for (offset_type k = start; k < end; ++k) {
        if (val_[k].first == key) {
          return val_[k].second;
        }
      }
#else
      const offset_type start = rawPtr_[hashVal];
      const offset_type end = rawPtr_[hashVal+1];
      for (offset_type k = start; k < end; ++k) {
        if (rawVal_[k].first == key) {
          return rawVal_[k].second;
        }
      }
#endif // HAVE_TPETRA_DEBUG
      return invalidValue_;
    }
  }

  /// \brief Number of (key, value) pairs in the table.
  ///
  /// This counts duplicate keys separately.
  KOKKOS_INLINE_FUNCTION offset_type numPairs () const {
    // NOTE (mfh 26 May 2015) This only works because the table
    // _stores_ pairs with duplicate keys separately.  If the table
    // didn't do that, we would have to keep a separate numPairs_
    // field (remembering the size of the input array of keys).
    return val_.dimension_0 ();
  }

  /// \brief The minimum key.
  ///
  /// This function does not throw, and always returns a value.  If
  /// the table is empty, the value is undefined.  Furthermore, if the
  /// table is empty, we do not promise that minKey() <= maxKey().
  ///
  /// This class assumes that both keys and values are numbers.
  /// Therefore, keys are less-than comparable.
  KOKKOS_INLINE_FUNCTION KeyType minKey () const {
    return minKey_;
  }

  /// \brief The maximum key.
  ///
  /// This function does not throw, and always returns a value.  If
  /// the table is empty, the value is undefined.  Furthermore, if the
  /// table is empty, we do not promise that minKey() <= maxKey().
  ///
  /// This class assumes that both keys and values are numbers.
  /// Therefore, keys are less-than comparable.
  KOKKOS_INLINE_FUNCTION KeyType maxKey () const {
    return maxKey_;
  }

  /// \brief Whether the table has any duplicate keys.
  ///
  /// This is a nonconst function because it requires running a Kokkos
  /// kernel to search the keys.  The result of the first call is
  /// cached and reused on subsequent calls.
  bool hasDuplicateKeys ();

  //! Implementation of Teuchos::Describable
  //@{
  //! Return a simple one-line description of this object.
  std::string description () const;

  //! Print this object with the given verbosity to the output stream.
  void
  describe (Teuchos::FancyOStream &out,
            const Teuchos::EVerbosityLevel verbLevel=
            Teuchos::Describable::verbLevel_default) const;
  //@}

private:
  //! Array of "row" offsets.
  ptr_type ptr_;
  //! Array of hash table entries.
  val_type val_;

#if ! defined(TPETRA_HAVE_KOKKOS_REFACTOR)
  /// \brief <tt>rawPtr_ == ptr_.ptr_on_device()</tt>.
  ///
  /// This is redundant, but we keep it around as a fair performance
  /// comparison against the "classic" version of Tpetra.
  const offset_type* rawPtr_;
  /// \brief <tt>rawVal_ == val_.ptr_on_device()</tt>.
  ///
  /// This is redundant, but we keep it around as a fair performance
  /// comparison against the "classic" version of Tpetra.
  const Kokkos::pair<KeyType, ValueType>* rawVal_;
#endif // ! defined(TPETRA_HAVE_KOKKOS_REFACTOR)

  /// \brief "Invalid" value of ValueType, used as a flag.
  ///
  /// We store this here because
  /// Teuchos::OrdinalTraits<ValueType>::invalid() is not a CUDA
  /// device function.  This is nonconst because otherwise the
  /// compiler deletes the implicit copy constructor.
  ValueType invalidValue_;

  /// \brief Minimum key (computed in init()).
  ///
  /// This class assumes that keys are less-than comparable.
  KeyType minKey_;

  /// \brief Maximum key (computed in init()).
  ///
  /// This class assumes that keys are less-than comparable.
  KeyType maxKey_;

  /// \brief Whether the table has checked for duplicate keys.
  ///
  /// This is set at the end of the first call to hasDuplicateKeys().
  /// The results of that method are cached in hasDuplicateKeys_ (see
  /// below).
  bool checkedForDuplicateKeys_;

  /// \brief Whether the table noticed any duplicate keys.
  ///
  /// This is only valid if checkedForDuplicateKeys_ (above) is true.
  bool hasDuplicateKeys_;

  /// \brief Whether the table has duplicate keys.
  ///
  /// This method doesn't cache anything (and is therefore marked
  /// const); hasDuplicateKeys() (which see) caches this result.
  bool checkForDuplicateKeys () const;

  //! The number of "buckets" in the bucket array.
  KOKKOS_INLINE_FUNCTION offset_type getSize () const {
    return ptr_.dimension_0 () == 0 ?
      static_cast<offset_type> (0) :
      static_cast<offset_type> (ptr_.dimension_0 () - 1);
  }

  //! Sanity checks; throw std::logic_error if any of them fail.
  void check () const;

  typedef Kokkos::View<const KeyType*,
                       typename ptr_type::HostMirror::array_layout,
                       typename ptr_type::HostMirror::execution_space,
                       Kokkos::MemoryUnmanaged> host_input_keys_type;

  typedef Kokkos::View<const ValueType*,
                       typename ptr_type::HostMirror::array_layout,
                       typename ptr_type::HostMirror::execution_space,
                       Kokkos::MemoryUnmanaged> host_input_vals_type;

  /// \brief Allocate storage and initialize the table; use given
  ///   initial min and max keys.
  ///
  /// Add <tt>(keys[i], startingValue + i)</tt> to the table,
  /// for i = 0, 1, ..., <tt>keys.size()</tt>.
  void
  init (const host_input_keys_type& keys,
        const ValueType startingValue,
        const KeyType initMinKey,
        const KeyType initMaxKey);

  /// \brief Allocate storage and initialize the table; use given
  ///   initial min and max keys.
  ///
  /// Add <tt>(keys[i], vals[i])</tt> to the table, for i = 0, 1, ...,
  /// <tt>keys.size()</tt>.  This is called by the version of the
  /// constructor that takes the same arguments.
  void
  init (const host_input_keys_type& keys,
        const host_input_vals_type& vals,
        const KeyType initMinKey,
        const KeyType initMaxKey);
};

} // Details namespace
} // Tpetra namespace

#endif // TPETRA_DETAILS_FIXEDHASHTABLE_DECL_HPP
