// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// ***********************************************************************
// @HEADER


#ifndef TEUCHOS_STRING_INDEXED_ORDERED_VALUE_OBJECT_CONTAINER_HPP
#define TEUCHOS_STRING_INDEXED_ORDERED_VALUE_OBJECT_CONTAINER_HPP


#include "Teuchos_Array.hpp"
#include "Teuchos_FilteredIterator.hpp"


namespace Teuchos {



/** \brief Base types for StringIndexedOrderedValueObjectContainer.
 */
class StringIndexedOrderedValueObjectContainerBase {
public:

  /** \brief Ordinal used for the index. */
  typedef Teuchos_Ordinal Ordinal;

  /** \brief Return the value for invalid ordinal. */
  static Ordinal getInvalidOrdinal() { return -1; }

  /** \brief Destructor */
  virtual ~StringIndexedOrderedValueObjectContainerBase() {};

public: // Public members (but only for unit testing purposes!)

  /** \brief A safe ordinal index type that default initializes to a special
   * value.
   */
  class OrdinalIndex {
  public:
    /** \brief . */
    typedef StringIndexedOrderedValueObjectContainerBase::Ordinal Ordinal;
    /** \brief . */
    Ordinal idx;
    /** \brief . */
    OrdinalIndex()
      : idx(StringIndexedOrderedValueObjectContainerBase::getInvalidOrdinal())
      {}
    /** \brief . */
    OrdinalIndex(const Ordinal idx_in)
      : idx(idx_in)
      {}
  };

  /** \brief A simple aggregate type to bind a key string and and objects value.
   *
   * This is meant to be a drop-in replacement for std::pair<std::string,
   * ObjType>.  That is why the key and the object are called 'first' and
   * 'second'.
   *
   * Note that there is no invariant on the key and the object value so there is
   * no need to encapsulate them.  This is good because we can't since we need
   * to provide 'first' and 'second' as raw data members to match std::pair.
   * However, we don't want to allow users to be able to change the 'first' key
   * string (like you can't with the std::map::value_type returned in the
   * iterator) so we use a const reference for 'first'.
   */
  template<class ObjType>
  class KeyObjectPair {
  public:
    /**  \brief . */
    const std::string &first;
    /**  \brief . */
    ObjType second;
    /** \brief . */
    std::string key;
    /** \brief . */
    KeyObjectPair() : first(key), second(ObjType()), key(""), isActive_(true) {}
    /** \brief . */
    KeyObjectPair(const std::string &key_in, const ObjType &obj_in, bool isActive_in = true)
      : first(key), second(obj_in), key(key_in), isActive_(isActive_in) {}
    /** \brief . */
    KeyObjectPair(const KeyObjectPair<ObjType> &kop)
      : first(key), second(kop.second), key(kop.key), isActive_(kop.isActive_) {}
    /** \brief . */
    KeyObjectPair<ObjType>& operator=(const KeyObjectPair<ObjType> &kop)
      {
        second = kop.second;
        key = kop.key;
        isActive_ = kop.isActive_;
        return *this;
      }
    /** \brief . */
    static KeyObjectPair<ObjType> makeInvalid()
      { return KeyObjectPair("", ObjType(), false); }
    /** \brief . */
    bool isActive() const { return isActive_; }
  private:
    bool isActive_;
  };

  /** \brief Predicate for selecting active object entries in filtered iterator. */
  template<class ObjType>
  class SelectActive {
  public:
    bool operator()(const KeyObjectPair<ObjType> &key_and_obj) const
      { return key_and_obj.isActive(); }
  };

  /** \brief Thrown if an invalid ordinal index is passed in. */
  class InvalidOrdinalIndexError : public ExceptionBase
  {public:InvalidOrdinalIndexError(const std::string& what_arg) : ExceptionBase(what_arg) {}};

  /** \brief Thrown if an invalid string is passed in. */
  class InvalidKeyError : public ExceptionBase
  {public:InvalidKeyError(const std::string& what_arg) : ExceptionBase(what_arg) {}};

};


/** \brief String indexed ordered value-type object container class.
 *
 * This class is a simple utility class for managing the storage and retrievel
 * of value-type objects which the following features/properties:
 * <ul>
 * <li> Objects are held by-value (i.e., like std::vector, std::map, etc.) </li>
 * <li> Objects are indexed by a unique string key value. </li>
 * <li> Objects can also be indexed by a unique ordering ordinal. </li>
 * <li> Iterators are orderded according to insertion order. </li>
 * </ul>
 *
 * The design of this class comes with a few important limitations:
 * <ul>
 * <li> Storage grows accoding to the number of insertions and removing
 *    objects does not decrease it.  This is a result of the simple
 *    implementation needed to preserve ordering <li>
 * </ul>
 *
 * \todo Implement compression of unused entries.  This will invalidate the
 * indexes but will allow handling of lots of inserts and deletes of elements.
 */
template<class ObjType>
class StringIndexedOrderedValueObjectContainer
  : private StringIndexedOrderedValueObjectContainerBase
{
private:

  /** \brief . */
  typedef KeyObjectPair<ObjType> key_and_obj_t;
  /** \brief . */
  typedef std::deque<key_and_obj_t> key_and_obj_array_t; // See notes below!
  /** \brief . */
  typedef std::map<std::string, OrdinalIndex> key_to_idx_map_t;

public:

  /** \name Public types. */
  //@{

  /** \brief Ordinal used for the index. */
  typedef StringIndexedOrderedValueObjectContainerBase::Ordinal Ordinal;

  /** \brief The non-const iterator type. */
  typedef FilteredIterator<typename key_and_obj_array_t::iterator, SelectActive<ObjType> >
  Iterator;

  /** \brief The const iterator type. */
  typedef FilteredIterator<typename key_and_obj_array_t::const_iterator, SelectActive<ObjType> >
  ConstIterator;

  //@}

  /** \name Constructors/Destructors/Info */
  //@{

  /** \brief . */
  StringIndexedOrderedValueObjectContainer();

  /** \brief . */
  Ordinal numObjects() const;

  /** \brief . */
  Ordinal numStorage() const;

  //@}

  /** \name Set, get, and remove functions */
  //@{

  /** \brief Set (or reset) object by value and return its ordinal index.
   *
   * If the object with the given key index does not exist, it will be added.
   * If an object with the given key does not exist, it will be created.
   *
   * \return Returns the ordinal index by which the object can be looked up
   * with.
   */
  Ordinal setObj(const std::string &key, const ObjType &obj);

  /** \brief Get the ordinal index given the string key.
   *
   * If the key does not exist, then getInvalidOrdinal() is returned.
   */
  inline Ordinal getObjOrdinalIndex(const std::string &key) const;

  /** \brief Get a nonconst semi-persisting association with the stored object
   * indexed by ordinal.
   */
  inline Ptr<ObjType> getNonconstObjPtr(const Ordinal &idx);

  /** \brief Get a const semi-persisting association with the stored object
   * indexed by ordinal.
   */
  inline Ptr<const ObjType> getObjPtr(const Ordinal &idx) const;

  /** \brief Get a nonconst semi-persisting association with the stored object
   * indexed by string key.
   */
  inline Ptr<ObjType> getNonconstObjPtr(const std::string &key);

  /** \brief Get a const semi-persisting association with the stored object
   * indexed by string key.
   */
  inline Ptr<const ObjType> getObjPtr(const std::string &key) const;

  /** \brief Remove an object given its ordinal index.
   *
   * Each object is errased by assigning to a default-constructed ObjType().
   * This, for example, will wipe out the reference count for a smart pointer
   * class or will unsize an array, etc..
   */
  void removeObj(const Ordinal &idx);

  /** \brief Remove an object given its string key. */
  void removeObj(const std::string &key);

  //@}

  /** \name Iterator access */
  //@{

  /** \brief . */
  inline Iterator nonconstBegin();

  /** \brief . */
  inline Iterator nonconstEnd();

  /** \brief . */
  inline ConstIterator begin() const;

  /** \brief . */
  inline ConstIterator end() const;

  //@}

private: // data members

  /** \brief Stories objects contiguously along with key strings. */
  key_and_obj_array_t key_and_obj_array_;
  /** \brief Provides lookups of key -> ordinal index into above array. */
  key_to_idx_map_t key_to_idx_map_;

  // The above is a fairly simple data-structure.
  //
  // key_and_obj_array_: Array that stores the objects (and key names), by
  // value, in the order they are inserted.  Any removed objects will have the
  // index valuie of getInvalidOrdinal().  The key strings are also storied
  // with the objects so that a clean iterator can over the objects has access
  // to both the key and the object value.
  //
  // key_to_idx_map_: Maps the unique string key to the unigue ordinal index
  // for an object.
  //
  // NOTES:
  // 
  // A) This data-structure stores the key names twice in order to allow for
  // optimal iterator performance.  The array key_and_obj_array_ allows fast
  // ordered iterators through the data but in order to also provide the names
  // in a fast manner, they have to be stored twice.
  //
  // B) Deleting objects is done by just removing an entry from
  // key_to_idx_map_ but the corresponding entry in key_and_obj_array_ is just
  // abandoned with the object value set to ObjType().

private: // member functions

  /** \brief . */
  void assertOrdinalIndex(const Ordinal idx) const;

  /** \brief . */
  key_and_obj_t& getNonconstKeyAndObject(const Ordinal idx);

  /** \brief . */
  const key_and_obj_t& getKeyAndObject(const Ordinal idx) const;

  /** \brief . */
  void throwInvalidKeyError(const Ordinal idx, const std::string &key) const;

  /** \brief . */
  Ordinal assertKeyGetOrdinal(const std::string &key) const;

};


//
// StringIndexedOrderedValueObjectContainer: Inline Implementations
//


// Set, get, and remove functions


template<class ObjType>
inline
Ptr<ObjType>
StringIndexedOrderedValueObjectContainer<ObjType>::getNonconstObjPtr(const Ordinal &idx)
{
  return ptrFromRef(getNonconstKeyAndObject(idx).second);
}


template<class ObjType>
inline
Ptr<const ObjType>
StringIndexedOrderedValueObjectContainer<ObjType>::getObjPtr(const Ordinal &idx) const
{
  return ptrFromRef(getKeyAndObject(idx).second);
}


template<class ObjType>
inline
Ptr<ObjType>
StringIndexedOrderedValueObjectContainer<ObjType>::getNonconstObjPtr(const std::string &key)
{
  return getNonconstObjPtr(assertKeyGetOrdinal(key));
}


template<class ObjType>
inline
Ptr<const ObjType>
StringIndexedOrderedValueObjectContainer<ObjType>::getObjPtr(const std::string &key) const
{
  return getObjPtr(assertKeyGetOrdinal(key));
}


// Iterator access


template<class ObjType>
inline
typename StringIndexedOrderedValueObjectContainer<ObjType>::Iterator
StringIndexedOrderedValueObjectContainer<ObjType>::nonconstBegin()
{
  return Iterator(key_and_obj_array_.begin(), key_and_obj_array_.begin(),
    key_and_obj_array_.end());
}


template<class ObjType>
inline
typename StringIndexedOrderedValueObjectContainer<ObjType>::Iterator
StringIndexedOrderedValueObjectContainer<ObjType>::nonconstEnd()
{
  return Iterator(key_and_obj_array_.end(), key_and_obj_array_.begin(),
    key_and_obj_array_.end());
}


template<class ObjType>
inline
typename StringIndexedOrderedValueObjectContainer<ObjType>::ConstIterator
StringIndexedOrderedValueObjectContainer<ObjType>::begin() const
{
  return ConstIterator(key_and_obj_array_.begin(), key_and_obj_array_.begin(),
    key_and_obj_array_.end());
}


template<class ObjType>
inline
typename StringIndexedOrderedValueObjectContainer<ObjType>::ConstIterator
StringIndexedOrderedValueObjectContainer<ObjType>::end() const
{
  return ConstIterator(key_and_obj_array_.end(), key_and_obj_array_.begin(),
    key_and_obj_array_.end());
}


//
// StringIndexedOrderedValueObjectContainer: Template Implementations
//


// Constructors/Destructors/Info


template<class ObjType>
StringIndexedOrderedValueObjectContainer<ObjType>::StringIndexedOrderedValueObjectContainer()
{}


template<class ObjType>
typename StringIndexedOrderedValueObjectContainer<ObjType>::Ordinal
StringIndexedOrderedValueObjectContainer<ObjType>::numObjects() const
{
  return key_to_idx_map_.size();
}


template<class ObjType>
typename StringIndexedOrderedValueObjectContainer<ObjType>::Ordinal
StringIndexedOrderedValueObjectContainer<ObjType>::numStorage() const
{
  return key_and_obj_array_.size();
}


// Set, get, and remove functions


template<class ObjType>
inline
typename StringIndexedOrderedValueObjectContainer<ObjType>::Ordinal
StringIndexedOrderedValueObjectContainer<ObjType>::getObjOrdinalIndex(const std::string &key) const
{
  key_to_idx_map_t::const_iterator itr = key_to_idx_map_.find(key);
  if (itr != key_to_idx_map_.end()) {
    return itr->second.idx;
  }
  return getInvalidOrdinal();
}


template<class ObjType>
typename StringIndexedOrderedValueObjectContainer<ObjType>::Ordinal
StringIndexedOrderedValueObjectContainer<ObjType>::setObj(const std::string &key,
  const ObjType &obj)
{
  typename key_to_idx_map_t::iterator obj_idx_itr = key_to_idx_map_.find(key);
  if (obj_idx_itr != key_to_idx_map_.end()) {
    // Object with the given key already exists
    const Ordinal obj_idx = obj_idx_itr->second.idx;
    key_and_obj_array_[obj_idx].second = obj;
    return obj_idx;
  }
  // Object with the given key does not already exist so create a new one.
  key_and_obj_array_.push_back(key_and_obj_t(key, obj));
  const Ordinal new_idx = key_and_obj_array_.size()-1;
  key_to_idx_map_[key] = new_idx;
  return new_idx;
}


template<class ObjType>
void StringIndexedOrderedValueObjectContainer<ObjType>::removeObj(const Ordinal &idx)
{
  key_and_obj_t &key_and_obj = getNonconstKeyAndObject(idx); 
  key_to_idx_map_.erase(key_and_obj.first);
  key_and_obj = key_and_obj_t::makeInvalid();
}


template<class ObjType>
void StringIndexedOrderedValueObjectContainer<ObjType>::removeObj(const std::string &key)
{
  typename key_to_idx_map_t::iterator itr = key_to_idx_map_.find(key);
  if (itr == key_to_idx_map_.end()) {
    throwInvalidKeyError(getInvalidOrdinal(), key);
  } 
  const Ordinal idx = itr->second.idx;
  key_to_idx_map_.erase(itr);
  key_and_obj_array_[idx] = key_and_obj_t::makeInvalid(); 
}


// private


template<class ObjType>
void StringIndexedOrderedValueObjectContainer<ObjType>::assertOrdinalIndex(const Ordinal idx) const
{
  TEUCHOS_TEST_FOR_EXCEPTION( !(0 <= idx && idx < numStorage()),
    InvalidOrdinalIndexError,
    "Error, the ordinal index " << idx << " is invalid"
    << " because it falls outside of the range of valid objects"
    << " [0,"<<numStorage()-1<<"]!");
}


template<class ObjType>
typename StringIndexedOrderedValueObjectContainer<ObjType>::key_and_obj_t&
StringIndexedOrderedValueObjectContainer<ObjType>::getNonconstKeyAndObject(const Ordinal idx)
{
  assertOrdinalIndex(idx);
  key_and_obj_t &key_and_obj = key_and_obj_array_[idx];
  TEUCHOS_TEST_FOR_EXCEPTION( !key_and_obj.isActive(),
    InvalidOrdinalIndexError,
    "Error, the ordinal index " << idx << " is invalid"
    << " because the object has been deleted!");
  return key_and_obj;
}


template<class ObjType>
const typename StringIndexedOrderedValueObjectContainer<ObjType>::key_and_obj_t&
StringIndexedOrderedValueObjectContainer<ObjType>::getKeyAndObject(const Ordinal idx) const
{
  assertOrdinalIndex(idx);
  const key_and_obj_t &key_and_obj = key_and_obj_array_[idx];
  TEUCHOS_TEST_FOR_EXCEPTION( !key_and_obj.isActive(),
    InvalidOrdinalIndexError,
    "Error, the ordinal index " << idx << " is invalid"
    << " because the object has been deleted!");
  return key_and_obj;
}


template<class ObjType>
void
StringIndexedOrderedValueObjectContainer<ObjType>::throwInvalidKeyError(
  const Ordinal idx, const std::string &key) const
{
  TEUCHOS_TEST_FOR_EXCEPTION( idx == getInvalidOrdinal(), InvalidKeyError,
    "Error, the key '" << key << "' does not exist!");
}


template<class ObjType>
typename StringIndexedOrderedValueObjectContainer<ObjType>::Ordinal
StringIndexedOrderedValueObjectContainer<ObjType>::assertKeyGetOrdinal(const std::string &key) const
{
  const Ordinal idx = getObjOrdinalIndex(key);
  throwInvalidKeyError(idx, key);
  return idx;
}

  
} // end of Teuchos namespace

/* Notes:
 *
 * This class may have a bit of a fragile implemenation.  It uses std::deque
 * instead of std::vector to hold the stored objects.  This is so that once an
 * object is added, it will keep the exact same address no matter how many
 * other objects are added.  This is not the case with std::vector because
 * when the memory is resized it will copy the value objects making the old
 * objects invalid.  My guess with the std::deque class is that it would
 * allocate and store the chunks such that once an objects was added to a
 * chunk that it would not get reallocated and moved like in std::vector.  I
 * don't feel very good about relying on this behavior but my guess is that
 * this will be pretty portable.  If this turns out not to be portable, we can
 * always just use RCP<ObjType> to store the object and that will guarantee
 * that the object address will never be moved.  Going with an RCP<ObjType>
 * implementation would also allow the Ptr<ObjType> views to catch dangling
 * references in a debug-mode build.
 */


#endif // TEUCHOS_STRING_INDEXED_ORDERED_VALUE_OBJECT_CONTAINER_HPP


