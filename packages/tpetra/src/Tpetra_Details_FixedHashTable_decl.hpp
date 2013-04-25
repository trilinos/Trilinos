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

#ifndef TPETRA_HASHTABLE_DECL_HPP
#define TPETRA_HASHTABLE_DECL_HPP

#include <Teuchos_Describable.hpp>
#include <Tpetra_ConfigDefs.hpp>
#include <Teuchos_ArrayRCP.hpp>

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
template<typename KeyType, typename ValueType>
class FixedHashTable : public Teuchos::Describable {
public:
  /// \brief Constructor: Add <tt>(keys[i], i)</tt> to the table, 
  ///   for i = 0, 1, ..., <tt>keys.size()</tt>.
  FixedHashTable (const ArrayView<const KeyType>& keys);

  /// \brief Constructor: Add <tt>(keys[i], startingValue + i)</tt> to
  ///   the table, for i = 0, 1, ..., <tt>keys.size()</tt>.
  ///
  /// This version is useful if Map wants to exclude an initial
  /// sequence of contiguous GIDs from the table, and start with a
  /// given LID.
  FixedHashTable (const ArrayView<const KeyType>& keys, 
		  const ValueType startingValue);

  //! Copy constructor: Make a shallow copy of the data.
  FixedHashTable (const FixedHashTable& obj);

  //! Get the value corresponding to the given key.
  ValueType get (const KeyType key) const;

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
  typedef Array<int>::size_type size_type;

  /// \brief <tt>ptr_.size() == size_ + 1</tt>.
  ///
  /// This is redundant, but we keep it around to avoid the function
  /// call (for <tt>ptr_.size()</tt>) in the hash table.
  KeyType size_;
  //! Array of "row" offsets.
  ArrayRCP<const size_type> ptr_;
  //! Array of hash table entries.
  ArrayRCP<const std::pair<KeyType, ValueType> > val_;
  /// \brief <tt>rawPtr_ == ptr_.getRawPtr()</tt>.
  ///
  /// This is redundant, but we keep it around to speed up get().
  const size_type* rawPtr_;
  /// \brief <tt>rawVal_ == val_.getRawPtr()</tt>.
  ///
  /// This is redundant, but we keep it around to speed up get().
  const std::pair<KeyType, ValueType>* rawVal_;

  /// \brief Initialize: Add <tt>(keys[i], startingValue + i)</tt> to
  ///   the table, for i = 0, 1, ..., <tt>keys.size()</tt>.
  void
  init (const ArrayView<const KeyType>& keys, 
	const ValueType startingValue);

  //! The hash function; it returns \c int no matter the value type.
  int hashFunc (const KeyType key) const;

  int getRecommendedSize (const int size);
};

} // Details namespace

} // Tpetra namespace

#endif
