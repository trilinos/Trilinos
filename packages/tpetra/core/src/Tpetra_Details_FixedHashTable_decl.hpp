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

#include <Tpetra_ConfigDefs.hpp>
#include <Teuchos_Describable.hpp>
#include <Kokkos_Core.hpp>

namespace Tpetra {
namespace Details {

/// \class FixedHashTable
/// \tparam KeyType The type of the hash table's keys.  This must be a
///   built-in signed or unsinged integer type.
/// \tparam ValueType The type of the hash table's values.  This must
///   be a built-in signed or unsigned integer type.
///
/// This class implements a hash table from signed integer keys to
/// signed integer values, where all the (key,value) pairs must be
/// added at once, and may not be changed or removed after being
/// added.  Keys and values may have different types.  Tpetra::Map may
/// use this to implement global-to-local index lookup.
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
         class DeviceType = Kokkos::Device<Kokkos::Serial, Kokkos::HostSpace> >
class FixedHashTable : public Teuchos::Describable {
public:
  //! Default constructor; makes an empty table.
  FixedHashTable ();

  /// \brief Constructor for arbitrary keys and contiguous values
  ///   starting with zero.
  ///
  /// Add <tt>(keys[i], i)</tt> to the table,
  /// for i = 0, 1, ..., <tt>keys.size()</tt>.
  FixedHashTable (const ArrayView<const KeyType>& keys);

  /// \brief Constructor for arbitrary keys and contiguous values
  ///   starting with \c startingValue.
  ///
  /// Add <tt>(keys[i], startingValue + i)</tt> to the table, for i =
  /// 0, 1, ..., <tt>keys.size()</tt>.  This version is useful if Map
  /// wants to exclude an initial sequence of contiguous GIDs from the
  /// table, and start with a given LID.
  FixedHashTable (const ArrayView<const KeyType>& keys,
                  const ValueType startingValue);

  /// \brief Constructor for arbitrary keys and arbitrary values.
  ///
  /// Add <tt>(keys[i], vals[i])</tt> to the table, for i = 0, 1, ...,
  /// <tt>keys.size()</tt>.  This version is useful for applications
  /// other than Map's GID-to-LID lookup table.
  FixedHashTable (const ArrayView<const KeyType>& keys,
                  const ArrayView<const ValueType>& vals);

  template<class K, class V, class D>
  friend struct FixedHashTable;

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
    typedef typename Kokkos::View<size_type*, DeviceType>::HostMirror ptr_type;
    typedef typename Kokkos::View<Kokkos::pair<KeyType, ValueType>*, DeviceType>::HostMirror val_type;

    ptr_type ptr (Kokkos::ViewAllocateWithoutInitializing ("ptr"), src.ptr_.dimension_0 ());
    // FIXME (mfh 01 May 2015) deep_copy won't work once we switch
    // from using host mirrors to storing the data in device memory.
    // In that case, we'll either need to fix the layout, or use a
    // copy kernel.  Fixing the layout would probably be easier, and
    // it doesn't matter since all Views here are 1-D.
    Kokkos::deep_copy (ptr, src.ptr_);
    val_type val (Kokkos::ViewAllocateWithoutInitializing ("val"), src.val_.dimension_0 ());
    Kokkos::deep_copy (val, src.val_);

    this->ptr_ = ptr;
    this->val_ = val;
#if ! defined(TPETRA_HAVE_KOKKOS_REFACTOR)
    this->rawPtr_ = ptr.ptr_on_device ();
    this->rawVal_ = val.ptr_on_device ();
#endif // ! defined(TPETRA_HAVE_KOKKOS_REFACTOR)
    this->hasDuplicateKeys_ = src.hasDuplicateKeys_;

#if defined(HAVE_TPETRA_DEBUG)
    this->check ();
#endif // defined(HAVE_TPETRA_DEBUG)
  }

  //! Get the value corresponding to the given key.
  ValueType get (const KeyType key) const;

  //! Whether the table noticed any duplicate keys on construction.
  bool hasDuplicateKeys () const {
    return hasDuplicateKeys_;
  }

  //! Implementation of Teuchos::Describable
  //@{
  //! Return a simple one-line description of this object.
  std::string description() const;

  //! Print this object with the given verbosity to the output stream.
  void
  describe (Teuchos::FancyOStream &out,
            const Teuchos::EVerbosityLevel verbLevel=
            Teuchos::Describable::verbLevel_default) const;
  //@}

private:
  typedef typename Kokkos::View<const KeyType*, DeviceType>::size_type size_type;

  //! Array of "row" offsets.
  typename Kokkos::View<const size_type*, DeviceType>::HostMirror ptr_;
  //! Array of hash table entries.
  typename Kokkos::View<const Kokkos::pair<KeyType, ValueType>*, DeviceType>::HostMirror val_;

#if ! defined(TPETRA_HAVE_KOKKOS_REFACTOR)
  /// \brief <tt>rawPtr_ == ptr_.ptr_on_device()</tt>.
  ///
  /// This is redundant, but we keep it around as a fair performance
  /// comparison against the "classic" version of Tpetra.
  const size_type* rawPtr_;
  /// \brief <tt>rawVal_ == val_.ptr_on_device()</tt>.
  ///
  /// This is redundant, but we keep it around as a fair performance
  /// comparison against the "classic" version of Tpetra.
  const Kokkos::pair<KeyType, ValueType>* rawVal_;
#endif // ! defined(TPETRA_HAVE_KOKKOS_REFACTOR)

  //! Whether the table noticed any duplicate keys on construction.
  bool hasDuplicateKeys_;

  //! The number of "buckets" in the bucket array.
  size_type getSize () const {
    return ptr_.dimension_0 () == 0 ? size_type (0) : ptr_.dimension_0 () - 1;
  }

  //! Sanity checks; throw std::logic_error if any of them fail.
  void check () const;

  /// \brief Allocate storage and initialize the table.
  ///
  /// Add <tt>(keys[i], startingValue + i)</tt> to the table,
  /// for i = 0, 1, ..., <tt>keys.size()</tt>.
  void
  init (const ArrayView<const KeyType>& keys,
        const ValueType startingValue);

  /// \brief Allocate storage and initialize the table.
  ///
  /// Add <tt>(keys[i], vals[i])</tt> to the table, for i = 0, 1, ...,
  /// <tt>keys.size()</tt>.  This is called by the version of the
  /// constructor that takes the same arguments.
  void
  init (const ArrayView<const KeyType>& keys,
        const ArrayView<const ValueType>& vals);

  //! The hash function; it returns \c int no matter the value type.
  int hashFunc (const KeyType key, const size_type size) const;

  //! Number of "buckets" that the constructor should allocate.
  int getRecommendedSize (const int size);
};

} // Details namespace

} // Tpetra namespace

#endif // TPETRA_DETAILS_FIXEDHASHTABLE_DECL_HPP
