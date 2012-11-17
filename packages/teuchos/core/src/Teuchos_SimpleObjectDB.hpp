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

#ifndef TEUCHOS_SIMPLE_OBJECT_DB_HPP
#define TEUCHOS_SIMPLE_OBJECT_DB_HPP


#include "Teuchos_Array.hpp"
#include "Teuchos_ConstNonconstObjectContainer.hpp"
#include "Teuchos_as.hpp"


/*! \file Teuchos_SimpleObjectDB.hpp
    \brief A simple object table class for Teuchos
*/

/*! \class Teuchos::SimpleObjectDB
    \brief This class provides a central place to store objects
*/

namespace Teuchos
{


/** \brief Simple object object database.
 *
 * This simple class stores and retrieves objects (using an simple integral
 * type index).
 *
 * The main advantages of using this class over a simple Teuchos::Array<RCP<T>
 * > object are:
 * <ul>
 *
 * <li>Object indexes are reused to allow large numbers of
 * store[Nonconst,Const]Obj()/removeObj() calls without growing memory.  A
 * simple Teuchos::Array implementation would not allow for that.</li>
 *
 * <li>Seemlessly handle both const and non-const objects using the runtime
 * protection of const.</li>
 *
 * </ul>
 *
 * NOTE: The type should always be a nonconst type (i.e. <tt>T=SomeType</tt>)
 * and not a const type (i.e. not <tt>T=const SomeType</tt>).  The class will
 * not compile if given a const type.
 */
template <class T>
class SimpleObjectDB
{
public:
  
  /** \brief Construct an empty DB. */
  SimpleObjectDB();

  /** \brief Return the current size of the table.
   *
   * This number will be greater or equal to the number of currently stored
   * objects.  This function is really only of interest to unit testing code.
   */
  int tableSize() const;

  /** \brief Return number of free indexes.
   *
   * This is the number of free table indexes that are ready to be recycled
   * for objects to be added.
   */
  int numFreeIndexes() const;

  /** \brief Return number of non-null stored objects.
   */
  int numObjects() const;

  /** \brief Store a non-const object.
   *
   * \return Index to the object.
   */
  int storeNonconstObj(const RCP<T> &obj);

  /** \brief Store a const object.
   *
   * \return Index to the object.
   */
  int storeConstObj(const RCP<const T> &obj);
  
  /** \brief Performs an rcp_dynamic_cast<>() to store the obejct.
   *
   * If the dynamic cast fails, then an exception is thrown.
   *
   * \return Index to the object.
   */
  template <class TOld>
  int storeCastedNonconstObj(const RCP<TOld> & robj_old);

  /** \brief Remove a stored object without returning it.
   *
   * This function avoids reference counting overhead when the user
   * just wants the object to go away and does not care to grab it.
   */
  void removeObj(const int index);

  /** \breif Remove a stored non-const object from the table and return it.
   *
   * This will throw NonconstAccessError if a const object was added with
   * <tt>storeConstObj()</tt>.
   */
  RCP<T> removeNonconstObj(const int index);

  /** \breif Remove a stored const object from the table and return it.
   */
  RCP<const T> removeConstObj(const int index);
  
  /** \brief Remove an indexed object from the table.
   *
   * If this object is solely owned by this table, then it will be destroyed
   * according to the embedded RCP destruction policy object.
   *
   * \return The reference count of the object after it removed from the
   * table.  If this count is 0, then the object was freed as a side-effect.
   */
  int removeRCP(int &index);
  
  /** \brief Get an object (nonconst persisting association).
   *
   * This will throw NonconstAccessError if a const object was added with
   * <tt>storeConstObj()</tt>.
   */
  RCP<T> getNonconstObjRCP(const int index);
  
  /** \brief Get an object (const persisting association).
   */
  RCP<const T> getConstObjRCP(const int index) const;
  
  /** \brief Get an object (nonconst semi-persisting association).
   *
   * This will throw NonconstAccessError if a const object was added with
   * <tt>storeConstObj()</tt>.
   */
  Ptr<T> getNonconstObjPtr(const int index);
  
  /** \brief Get an object (const semi-persisting association).
   */
  Ptr<const T> getConstObjPtr(const int index) const;
  
  /** \brief Clear out all storage.
   *
   * Clears out the table storage and all RCPs to stored objects.  Any objects
   * solely owned by this table will be freed.
   */
  void purge();
  
private:

  typedef Array<ConstNonconstObjectContainer<T> > tableOfObjects_t;
  typedef Array<int> freedIndices_t;
  
  tableOfObjects_t tableOfObjects_;
  freedIndices_t freedIndices_;

  void validateIndex(const int index) const;

  template <class T2>
  int storeObjectImpl(const RCP<T2> &robj);
  
  void removeObjImpl(const int index);

};


/** \brief Nonmember constructor */
template <class T>
RCP<SimpleObjectDB<T> > createSimpleObjectDB()
{
  return rcp(new SimpleObjectDB<T>);
}


//
// Template definitions
//


template <class T>
SimpleObjectDB<T>::SimpleObjectDB()
{}


template <class T>
int SimpleObjectDB<T>::tableSize() const
{
  return tableOfObjects_.size();
}


template <class T>
int SimpleObjectDB<T>::numFreeIndexes() const
{
  return freedIndices_.size();
}


template <class T>
int SimpleObjectDB<T>::numObjects() const
{
  return tableSize() - numFreeIndexes();
}


template <class T>
int SimpleObjectDB<T>::storeNonconstObj(const RCP<T> &obj)
{
  return storeObjectImpl(obj);
}


template <class T>
int SimpleObjectDB<T>::storeConstObj(const RCP<const T> &obj)
{
  return storeObjectImpl(obj);
}


template <class T>
template <class TOld>
int SimpleObjectDB<T>::storeCastedNonconstObj(const RCP<TOld> & robj_old)
{
  return storeNonconstObj(rcp_dynamic_cast<T>(robj_old, true));
}


template <class T>
void SimpleObjectDB<T>::removeObj(const int index)
{
  validateIndex(index);
  removeObjImpl(index);
}


template <class T>
RCP<T> SimpleObjectDB<T>::removeNonconstObj(const int index)
{
  validateIndex(index);
  const RCP<T> obj = tableOfObjects_[index].getNonconstObj();
  removeObjImpl(index);
  return obj;
}


template <class T>
RCP<const T> SimpleObjectDB<T>::removeConstObj(const int index)
{
  validateIndex(index);
  const RCP<const T> obj = tableOfObjects_[index].getConstObj();
  removeObjImpl(index);
  return obj;
}


template <class T>
int SimpleObjectDB<T>::removeRCP(int &index)
{
  const int index_in = index;
  validateIndex(index);
  const int cnt = tableOfObjects_[index_in].count();
  removeObjImpl(index_in);
  index = -1;
  return (cnt - 1);
}


template <class T>
RCP<T> SimpleObjectDB<T>::getNonconstObjRCP(const int index)
{
  validateIndex(index);
  return tableOfObjects_[index].getNonconstObj();
}


template <class T>
RCP<const T> SimpleObjectDB<T>::getConstObjRCP(const int index) const
{
  validateIndex(index);
  return tableOfObjects_[index].getConstObj();
}


template <class T>
Ptr<T> SimpleObjectDB<T>::getNonconstObjPtr(const int index)
{
  validateIndex(index);
  return tableOfObjects_[index].getNonconstObj().ptr();
}


template <class T>
Ptr<const T> SimpleObjectDB<T>::getConstObjPtr(const int index) const
{
  validateIndex(index);
  return tableOfObjects_[index].getConstObj().ptr();
}


template <class T>
void SimpleObjectDB<T>::purge()
{
  // Wipe out all memory (see Item 82 in "C++ Coding Standards")
  tableOfObjects_t().swap(tableOfObjects_);
  freedIndices_t().swap(freedIndices_);
}


// private


template <class T>
void SimpleObjectDB<T>::validateIndex(const int index) const
{
  using Teuchos::as;
  TEUCHOS_TEST_FOR_EXCEPTION(
    !(0 <= index && index < as<int>(tableOfObjects_.size())),
    RangeError,
    "Error, the object index = " << index << " falls outside of the range"
    << " of valid objects [0,"<<tableOfObjects_.size()<<"]");
  const RCP<const T> &obj = tableOfObjects_[index].getConstObj();
  TEUCHOS_TEST_FOR_EXCEPTION(is_null(obj), NullReferenceError,
    "Error, the object at index "<<index<<" of type "
    <<TypeNameTraits<T>::name()<<" has already been deleted!");
}


template <class T>
template <class T2>
int SimpleObjectDB<T>::storeObjectImpl(const RCP<T2> & robj)
{
  robj.assert_not_null();

  int index = -1;

  if (freedIndices_.size() != 0) {
    index = freedIndices_.back();
    freedIndices_.pop_back();
    tableOfObjects_[index].initialize(robj);
  } else {
    tableOfObjects_.push_back(robj);
    index = tableOfObjects_.size() - 1;
  }

  return index;
}


template <class T>
void SimpleObjectDB<T>::removeObjImpl(const int index)
{
  tableOfObjects_[index] = null;
  freedIndices_.push_back(index);
}


} // end namespace Teuchos


#endif // TEUCHOS_SIMPLE_OBJECT_DB_HPP

