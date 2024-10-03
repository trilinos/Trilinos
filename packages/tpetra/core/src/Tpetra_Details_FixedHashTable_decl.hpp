// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_DETAILS_FIXEDHASHTABLE_DECL_HPP
#define TPETRA_DETAILS_FIXEDHASHTABLE_DECL_HPP

#include "Tpetra_Details_Hash.hpp"
#include "Tpetra_Details_OrdinalTraits.hpp"
#include "Tpetra_Details_copyOffsets.hpp"
#include "Tpetra_Details_Profiling.hpp"

#include "Teuchos_Describable.hpp"
#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_VerbosityLevel.hpp"

#include "Kokkos_Core.hpp"
#include "Kokkos_ArithTraits.hpp"

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
class FixedHashTable {
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

  /// \brief Whether the table was created using one of the
  ///   constructors that assume contiguous values.
  ///
  /// \return false if this object was created using the two-argument
  ///   (keys, vals) constructor (that takes lists of both keys and
  ///   values), else true.
  KOKKOS_INLINE_FUNCTION bool hasContiguousValues () const {
    return contiguousValues_;
  }

public:
  /// \brief Type of a 1-D Kokkos::View (array) used to store keys.
  ///
  /// This is the type preferred by FixedHashTable's constructors.
  typedef Kokkos::View<const KeyType*, Kokkos::LayoutLeft, device_type> keys_type;

  //! Default constructor; makes an empty table.
  KOKKOS_DEFAULTED_FUNCTION FixedHashTable() = default;

  /// \brief Constructor for arbitrary keys and contiguous values
  ///   starting with zero.
  ///
  /// Add <tt>(keys[i], i)</tt> to the table,
  /// for i = 0, 1, ..., <tt>keys.extent(0)</tt>.
  ///
  /// \param keys [in] The keys in the hash table.  The table
  ///   <i>always</i> keeps a (shallow) copy, and thus hasKeys() is
  ///   true on return.
  FixedHashTable (const keys_type& keys);

  /// \brief Constructor for arbitrary keys and contiguous values
  ///   starting with zero.
  ///
  /// Add <tt>(keys[i], i)</tt> to the table,
  /// for i = 0, 1, ..., <tt>keys.size()</tt>.
  ///
  /// \param keys [in] The keys in the hash table.
  FixedHashTable (const Teuchos::ArrayView<const KeyType>& keys);

  /// \brief Constructor for arbitrary keys and contiguous values
  ///   starting with \c startingValue.
  ///
  /// Add <tt>(keys[i], startingValue + i)</tt> to the table, for i =
  /// 0, 1, ..., <tt>keys.extent(0)</tt>.  This version is useful
  /// if Map wants to exclude an initial sequence of contiguous GIDs
  /// from the table, and start with a given LID.
  ///
  /// \param keys [in] The keys in the hash table.  The table
  ///   <i>always</i> keeps a (shallow) copy, and thus hasKeys() is
  ///   true on return.
  /// \param startingValue [in] First value in the contiguous sequence
  ///   of values.
  FixedHashTable (const keys_type& keys,
                  const ValueType startingValue);

  /// \brief Constructor for arbitrary keys and contiguous values
  ///   starting with \c startingValue.
  ///
  /// Add <tt>(keys[i], startingValue + i)</tt> to the table, for i =
  /// 0, 1, ..., <tt>keys.size()</tt>.  This version is useful if Map
  /// wants to exclude an initial sequence of contiguous GIDs from the
  /// table, and start with a given LID.
  ///
  /// \param keys [in] The keys in the hash table.
  /// \param startingValue [in] First value in the contiguous sequence
  ///   of values.
  FixedHashTable (const Teuchos::ArrayView<const KeyType>& keys,
                  const ValueType startingValue);

  /// \brief Constructor for arbitrary keys and contiguous values
  ///   starting with \c startingValue, that takes an initial
  ///   contiguous sequence of keys, stored as a Kokkos::View.
  ///
  /// Add <tt>(keys[i], startingValue + i)</tt> to the table, for i =
  /// 0, 1, ..., <tt>keys.size()</tt>.
  ///
  /// \param keys [in] The keys in the hash table.  We reserve the
  ///   right to keep the input View without a deep copy.  If you
  ///   intend to modify this View after calling this constructor,
  ///   please make a deep copy yourself and give the copy to this
  ///   constructor.
  /// \param firstContigKey [in] First key in the initial contiguous
  ///   sequence of keys.
  /// \param lastContigKey [in] Last key (inclusive!) in the initial
  ///   contiguous sequence of keys.
  /// \param startingValue [in] First value in the contiguous sequence
  ///   of values.
  FixedHashTable (const keys_type& keys,
                  const KeyType firstContigKey,
                  const KeyType lastContigKey,
                  const ValueType startingValue);

  /// \brief Constructor for arbitrary keys and contiguous values
  ///   starting with \c startingValue, that takes an initial
  ///   contiguous sequence of keys, stored as a Teuchos::ArrayView
  ///   (host pointer).
  ///
  /// Add <tt>(keys[i], startingValue + i)</tt> to the table, for i =
  /// 0, 1, ..., <tt>keys.size()</tt>.
  ///
  /// \param keys [in] The keys in the hash table.
  /// \param firstContigKey [in] First key in the initial contiguous
  ///   sequence of keys.
  /// \param lastContigKey [in] Last key (inclusive!) in the initial
  ///   contiguous sequence of keys.
  /// \param startingValue [in] First value in the contiguous sequence
  ///   of values.
  FixedHashTable (const Teuchos::ArrayView<const KeyType>& keys,
                  const KeyType firstContigKey,
                  const KeyType lastContigKey,
                  const ValueType startingValue);

  /// \brief Constructor for arbitrary keys and arbitrary values.
  ///
  /// Add <tt>(keys[i], vals[i])</tt> to the table, for i = 0, 1, ...,
  /// <tt>keys.size()</tt>.  This version is useful for applications
  /// other than Map's GID-to-LID lookup table.
  ///
  /// \param keys [in] The keys in the hash table.
  /// \param vals [in] The values in the hash table.
  FixedHashTable (const Teuchos::ArrayView<const KeyType>& keys,
                  const Teuchos::ArrayView<const ValueType>& vals);

  template<class K, class V, class D>
  friend class FixedHashTable;

  /// \brief "Copy" constructor that takes a FixedHashTable with the
  ///   same KeyType and ValueType, but a different DeviceType.
  ///
  /// This constructor makes a deep copy of the input's data if
  /// necessary.
  template<class InDeviceType>
  FixedHashTable (const FixedHashTable<KeyType, ValueType, InDeviceType>& src,
                  typename std::enable_if<! std::is_same<DeviceType, InDeviceType>::value, int>::type* = NULL)
  {
    using Kokkos::ViewAllocateWithoutInitializing;
    typedef typename ptr_type::non_const_type nonconst_ptr_type;
    typedef typename val_type::non_const_type nonconst_val_type;

    Tpetra::Details::ProfilingRegion pr("Tpetra::Details::FixedHashTable::ctor(InDeviceType)");

    // FIXME (mfh 28 May 2015) The code below _always_ copies.  This
    // shouldn't be necessary if the input and output memory spaces
    // are the same.  However, it is always correct.

    // Different Devices may have different offset_type, because
    // offset_type comes from the memory space's size_type typedef.
    // That's why we use a specialized deep copy function here instead
    // of Kokkos::deep_copy.
    nonconst_ptr_type ptr (ViewAllocateWithoutInitializing ("Tpetra::FixedHashTable::ptr"),
                           src.ptr_.extent (0));
    ::Tpetra::Details::copyOffsets (ptr, src.ptr_);
    nonconst_val_type val (ViewAllocateWithoutInitializing ("Tpetra::FixedHashTable::val"),
                           src.val_.extent (0));
    // val and src.val_ have the same entry types, unlike (possibly)
    // ptr and src.ptr_.  Thus, we can use Kokkos::deep_copy here.
    // DEEP_COPY REVIEW - DEVICE-TO-DEVICE
    Kokkos::deep_copy (execution_space(), val, src.val_);

    this->ptr_ = ptr;
    this->val_ = val;
    this->minKey_ = src.minKey_;
    this->maxKey_ = src.maxKey_;
    this->minVal_ = src.minVal_;
    this->maxVal_ = src.maxVal_;
    this->firstContigKey_ = src.firstContigKey_;
    this->lastContigKey_ = src.lastContigKey_;
    this->contiguousValues_ = src.contiguousValues_;
    this->checkedForDuplicateKeys_ = src.checkedForDuplicateKeys_;
    this->hasDuplicateKeys_ = src.hasDuplicateKeys_;
  }

  //! Get the value corresponding to the given key.
  KOKKOS_INLINE_FUNCTION ValueType get (const KeyType& key) const {
    const offset_type size = this->getSize ();
    if (size == 0) {
      // Don't use Teuchos::OrdinalTraits or std::numeric_limits here,
      // because neither have the right device function markings.
      return Tpetra::Details::OrdinalTraits<ValueType>::invalid ();
    }

    // If this object assumes contiguous values, then it doesn't store
    // the initial sequence of >= 1 contiguous keys in the table.
    if (this->hasContiguousValues () &&
        key >= firstContigKey_ && key <= lastContigKey_) {
      return static_cast<ValueType> (key - firstContigKey_) + this->minVal ();
    }

    const typename hash_type::result_type hashVal =
      hash_type::hashFunc (key, size);

    const offset_type start = ptr_[hashVal];
    const offset_type end = ptr_[hashVal+1];
    for (offset_type k = start; k < end; ++k) {
      if (val_[k].first == key) {
        return val_[k].second;
      }
    }

    // Don't use Teuchos::OrdinalTraits or std::numeric_limits here,
    // because neither have the right device function markings.
    return Tpetra::Details::OrdinalTraits<ValueType>::invalid ();
  }

  /// \brief Number of (key, value) pairs in the table.
  ///
  /// This counts pairs with the same key value as separate pairs.
  KOKKOS_INLINE_FUNCTION offset_type numPairs () const {
    // NOTE (mfh 26 May 2015) Using val_.extent(0) only works
    // because the table stores pairs with duplicate keys separately.
    // If the table didn't do that, we would have to keep a separate
    // numPairs_ field (remembering the size of the input array of
    // keys).
    if (this->hasContiguousValues ()) {
      return val_.extent (0) + static_cast<offset_type> (lastContigKey_ - firstContigKey_);
    }
    else {
      return val_.extent (0);
    }
  }

  /// \brief The minimum key in the table.
  ///
  /// This function does not throw.  If the table is empty, the return
  /// value is undefined.  Furthermore, if the table is empty, we do
  /// not promise that minKey() <= maxKey().
  ///
  /// This class assumes that both keys and values are numbers.
  /// Therefore, keys are less-than comparable.
  KOKKOS_INLINE_FUNCTION KeyType minKey () const {
    return minKey_;
  }

  /// \brief The maximum key in the table.
  ///
  /// This function does not throw.  If the table is empty, the return
  /// value is undefined.  Furthermore, if the table is empty, we do
  /// not promise that minKey() <= maxKey().
  ///
  /// This class assumes that both keys and values are numbers.
  /// Therefore, keys are less-than comparable.
  KOKKOS_INLINE_FUNCTION KeyType maxKey () const {
    return maxKey_;
  }

  /// \brief The minimum value in the table.
  ///
  /// A "value" is the result of calling get() on a key.
  ///
  /// This function does not throw.  If the table is empty, the return
  /// value is undefined.  Furthermore, if the table is empty, we do
  /// not promise that minVal() <= maxVal().
  KOKKOS_INLINE_FUNCTION ValueType minVal () const {
    return minVal_;
  }

  /// \brief The maximum value in the table.
  ///
  /// A "value" is the result of calling get() on a key.
  ///
  /// This function does not throw.  If the table is empty, the return
  /// value is undefined.  Furthermore, if the table is empty, we do
  /// not promise that minVal() <= maxVal().
  KOKKOS_INLINE_FUNCTION ValueType maxVal () const {
    return maxVal_;
  }

  /// \brief Whether the table has any duplicate keys.
  ///
  /// This is a nonconst function because it requires running a Kokkos
  /// kernel to search the keys.  The result of the first call is
  /// cached and reused on subsequent calls.
  ///
  /// This function is the "local" (to an MPI process) version of
  /// Tpetra::Map::isOneToOne.  If a Tpetra::Map has duplicate keys
  /// (global indices) on any one MPI process, then it is most
  /// certainly not one to one.  The opposite may not necessarily be
  /// true, because a Tpetra::Map might have duplicate global indices
  /// that occur on different MPI processes.
  bool hasDuplicateKeys ();

  /// \brief Implementation of Teuchos::Describable interface.
  ///
  /// FixedHashTable can't actually inherit from Teuchos::Describable,
  /// because that would prevent instances of this class from living
  /// in GPU __device__ code.  See GitHub issue #1623.
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

  /// \brief Minimum key (computed in init()).
  ///
  /// In Tpetra::Map, this corresponds to the minimum global index
  /// (local to the MPI process).
  /// @remark It will be set in @ref init.
  KeyType minKey_ = ::Kokkos::ArithTraits<KeyType>::max();

  /// \brief Maximum key (computed in init()).
  ///
  /// In Tpetra::Map, this corresponds to the maximum global index
  /// (local to the MPI process).
  /// @remark It will be set in @ref init.
  KeyType maxKey_ = ::Kokkos::ArithTraits<KeyType>::max();

  /// \brief Minimum value.
  ///
  /// In Tpetra::Map, this corresponds to the minimum local index
  /// (local to the MPI process).
  ValueType minVal_ = ::Kokkos::ArithTraits<ValueType>::max();

  /// \brief Maximum value.
  ///
  /// In Tpetra::Map, this corresponds to the maximum local index
  /// (local to the MPI process).
  ValueType maxVal_ = ::Kokkos::ArithTraits<ValueType>::is_integer ?
                      ::Kokkos::ArithTraits<ValueType>::min() :
                     -::Kokkos::ArithTraits<ValueType>::max();

  /// \brief First key in any initial contiguous sequence.
  ///
  /// This only has a defined value if the number of keys is nonzero.
  /// In that case, the initial contiguous sequence of keys may have
  /// length 1 or more.  Length 1 means that the sequence is trivial
  /// (there are no initial contiguous keys).
  KeyType firstContigKey_ = ::Kokkos::ArithTraits<KeyType>::max();

  /// \brief Last key in any initial contiguous sequence.
  ///
  /// This only has a defined value if the number of keys is nonzero.
  /// In that case, the initial contiguous sequence of keys may have
  /// length 1 or more.  Length 1 means that the sequence is trivial
  /// (there are no initial contiguous keys).
  KeyType lastContigKey_ = ::Kokkos::ArithTraits<KeyType>::max();

  /// \brief Whether the table was created using one of the
  ///   constructors that assume contiguous values.
  ///
  /// This is false if this object was created using the two-argument
  /// (keys, vals) constructor (that takes lists of both keys and
  /// values), else true.
  bool contiguousValues_ = true;

  /// \brief Whether the table has checked for duplicate keys.
  ///
  /// This is set at the end of the first call to hasDuplicateKeys().
  /// The results of that method are cached in hasDuplicateKeys_ (see
  /// below).
  /// @remark It will be revised in @ref hasDuplicateKeys.
  bool checkedForDuplicateKeys_ = true;

  /// \brief Whether the table noticed any duplicate keys.
  ///
  /// This is only valid if checkedForDuplicateKeys_ (above) is true.
  bool hasDuplicateKeys_ = false;

  /// \brief Whether the table has duplicate keys.
  ///
  /// This method doesn't cache anything (and is therefore marked
  /// const); hasDuplicateKeys() (which see) caches this result.
  bool checkForDuplicateKeys () const;

  //! The number of "buckets" in the bucket array.
  KOKKOS_INLINE_FUNCTION offset_type getSize () const {
    return ptr_.extent (0) == 0 ?
      static_cast<offset_type> (0) :
      static_cast<offset_type> (ptr_.extent (0) - 1);
  }

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
  /// Add <tt>(keys[i], startingValue + i)</tt> to the table, for i =
  /// 0, 1, ..., <tt>keys.size()</tt>.  Optionally, compute any
  /// initial contiguous sequence of keys.
  void
  init (const keys_type& keys,
        const ValueType startingValue,
        KeyType initMinKey,
        KeyType initMaxKey,
        KeyType firstContigKey,
        KeyType lastContigKey,
        const bool computeInitContigKeys);

  /// \brief Allocate storage and initialize the table; use given
  ///   initial min and max keys.
  ///
  /// Add <tt>(keys[i], vals[i])</tt> to the table, for i = 0, 1, ...,
  /// <tt>keys.size()</tt>.  This is called by the version of the
  /// constructor that takes the same arguments.
  void
  init (const host_input_keys_type& keys,
        const host_input_vals_type& vals,
        KeyType initMinKey,
        KeyType initMaxKey);
};

} // Details namespace
} // Tpetra namespace

#endif // TPETRA_DETAILS_FIXEDHASHTABLE_DECL_HPP
