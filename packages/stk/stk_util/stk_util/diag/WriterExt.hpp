/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef STK_UTIL_DIAG_WRITEREXT_HPP
#define STK_UTIL_DIAG_WRITEREXT_HPP

#include <bitset>
#include <iostream>
#include <list>
#include <map>
#include <set>
#include <memory>
#include <stack>
#include <string>
#include <typeinfo>
#include <vector>
#include <stk_util/stk_config.h>
#include <stk_util/util/VecMap.hpp>
#include <stk_util/util/VecSet.hpp>

#include <stk_util/diag/Writer.hpp>


#include <stk_util/diag/String.hpp>
#include <stk_util/diag/StringUtil.hpp>

#include <stk_util/util/FeatureTest.hpp>
#include <stk_util/util/Array.hpp>
#include <stk_util/util/FArrayPrint.hpp>
#include <stk_util/parallel/MPI.hpp>
#include <stk_util/diag/Mapv.hpp>
#include <stk_util/diag/Option.hpp>


namespace stk {
namespace diag {

///
/// @addtogroup diag_writer_detail
/// @{
///

/**
 * @brief Function <b>operator<<</b> wrties a std::type_info name to the diagnostic
 * writer.
 *
 * @param dout    a <b>Writer</b> reference to the diagnostic writer to
 *      write the std::type_info object to.
 *
 * @param t    a <b>std::type_info</b> const reference to the std::typeinfo
 *      object.
 *
 * @return    a <b>Writer</b> reference to this object
 */
Writer &operator<<(Writer &dout, const std::type_info &t);

/**
 * @brief Template function <b>operator<<</b> writes an std::auto_ptr object
 * address and content to the diagnostic writer.
 *
 * @param dout    a <b>Writer</b> reference to the diagnostic writer to
 *      write the std::auto_ptr object.
 *
 * @param t    a <b>std::auto_ptr</b> const reference to the object.
 *
 * @return    a <b>Writer</b> reference to this object
 */
template <class T>
Writer &operator<<(Writer &dout, const std::auto_ptr<T> &t) {
  if (t.get())
    dout << " " << typeid(t) << ", " << t.get() << ", " << *t;
  else
    dout << " " << typeid(t) << ", <not created or not owner>";

  return dout;
}

/**
 * @brief Template function <b>operator<<</b> writes the members of an arbitrary
 * std::pair object to the diagnostic writer.
 *
 * @param dout    a <b>Writer</b> reference to the diagnostic writer to
 *      write the std::pair members to.
 *
 * @param pair    a <b>std::pair</b> const reference to the pair of objects.
 *
 * @return    a <b>Writer</b> reference to this object
 */
template <class T, class U>
Writer &operator<<(Writer & dout, const std::pair<T, U> &pair) {
//  dout << typeid(pair) << "(" << pair.first << ":" << pair.second << ")";
  dout << "(" << pair.first << ":" << pair.second << ")";

  return dout;
}

/**
 * @brief Template <b>dump</b> prints the object contained within a std::vector
 * object to the diagnostic writer.
 *
 * @param dout    a <b>Writer</b> reference to the diagnostic writer to
 *      write the std::vector to.
 *
 * @param t    a <b>std::vector</b> of objects.
 *
 * @return    a <b>Writer</b> reference to this object
 */
template <class T>
Writer &
dump(
  Writer &      dout,
  const std::vector<T> &  t)
{
  if (dout.shouldPrint()) {
    dout << typeid(t) << ", size " << t.size() << push << dendl;

    if (t.size() <= 10) {
      for (typename std::vector<T>::const_iterator it = t.begin(); it != t.end(); ++it)
        dout << (*it) << " ";
      dout << dendl;
    }
    else {
      int i = 0;
      for (typename std::vector<T>::const_iterator it = t.begin(); it != t.end(); ++it, ++i)
        dout << "[" << i << "] " << (*it) << dendl;
    }

    dout << pop;
  }

  return dout;
}

/**
 * @brief Template function <b>dump</b> prints the object pointed to that are
 * contained within a std::vector object to the diagnostic writer.
 *
 * @param dout    a <b>Writer</b> reference to the diagnostic writer to
 *      write the std::vector to.
 *
 * @param t    a <b>std::vector</b> of objects.
 *
 * @return    a <b>Writer</b> reference to this object
 */
template <class T>
Writer &
dump(
  Writer &      dout,
  const std::vector<T *> &  t)
{
  if (dout.shouldPrint()) {
    dout << typeid(t) << ", size " << t.size() << push << dendl;

    int i = 0;
    for (typename std::vector<T *>::const_iterator it = t.begin(); it != t.end(); ++it, ++i)
      dout << "[" << i << "] " << c_ptr_<T>(*it) << dendl;

    dout << pop;
  }

  return dout;
}

/**
 * @brief Template function <b>dump</b> prints the object contained within a
 * std::list object to the diagnostic writer.
 *
 * @param dout    a <b>Writer</b> reference to the diagnostic writer to
 *      write the std::list to.
 *
 * @param t    a <b>std::list</b> of objects.
 *
 * @return    a <b>Writer</b> reference to this object
 */
template <class T>
Writer &
dump(
  Writer &      dout,
  const std::list<T> &    t)
{
  if (dout.shouldPrint()) {
    dout << typeid(t) << ", size " << t.size() << push << dendl;

    int i = 0;
    for (typename std::list<T>::const_iterator it = t.begin(); it != t.end(); ++it, ++i)
      dout << "[" << i << "] " << (*it) << dendl;

    dout << pop;
  }

  return dout;
}

/**
 * @brief Template function <b>dump</b> prints the object pointed to that are
 * contained within a std::list object to the diagnostic writer.
 *
 * @param dout    a <b>Writer</b> reference to the diagnostic writer to
 *      write the std::list to.
 *
 * @param t    a <b>std::list</b> of objects.
 *
 * @return    a <b>Writer</b> reference to this object
 */
template <class T>
Writer &
dump(
  Writer &      dout,
  const std::list<T *> &  t)
{
  if (dout.shouldPrint()) {
    dout << typeid(t) << ", size " << t.size() << push << dendl;

    int i = 0;
    for (typename std::list<T *>::const_iterator it = t.begin(); it != t.end(); ++it, ++i)
      dout << "[" << i << "] " << c_ptr_<T>(*it) << dendl;

    dout << pop;
  }

  return dout;
}

/**
 * @brief Template function <b>dump</b> prints the object contained within a
 * std::map object to the diagnostic writer.
 *
 * @param dout    a <b>Writer</b> reference to the diagnostic writer to
 *      write the std::map to.
 *
 * @param t    a <b>std::map</b> of objects.
 *
 * @return    a <b>Writer</b> reference to this object
 */
template <class Key, class T, class L>
Writer &
dump(
  Writer &      dout,
  const std::map<Key, T, L> &  t)
{
  if (dout.shouldPrint()) {
    dout << typeid(t) << ", size " << t.size() << push << dendl;

    for (typename std::map<Key, T, L>::const_iterator it = t.begin(); it != t.end(); ++it)
      dout << "[" << (*it).first << "] " << (*it).second << dendl;

    dout << pop;
  }

  return dout;
}

/**
 * @brief Template function <b>dump</b> prints the object pointed to that are
 * contained within a std::map to the diagnostic writer.
 *
 * @param dout    a <b>Writer</b> reference to the diagnostic writer to
 *      write the std::map to.
 *
 * @param t    a <b>std::map</b> of objects.
 *
 * @return    a <b>Writer</b> reference to this object
 */
template <class Key, class T, class L>
Writer &
dump(
  Writer &      dout,
  const std::map<Key, T *, L> &  t)
{
  if (dout.shouldPrint()) {
    dout << typeid(t) << ", size " << t.size() << push << dendl;

    for (typename std::map<Key, T *, L>::const_iterator it = t.begin(); it != t.end(); ++it)
      dout << "[" << (*it).first << "] " << c_ptr_<T>((*it).second) << std::endl;

    dout << pop;
  }

  return dout;
}

/**
 * @brief Template function <b>dump</b> prints the object contained within a
 * std::multimap object to the diagnostic writer.
 *
 * @param dout    a <b>Writer</b> reference to the diagnostic writer to
 *      write the std::multimap to.
 *
 * @param t    a <b>std::multimap</b> of objects.
 *
 * @return    a <b>Writer</b> reference to this object
 */
template <class Key, class T, class L>
Writer &
dump(
  Writer &      dout,
  const std::multimap<Key, T, L> &  t)
{
  if (dout.shouldPrint()) {
    dout << typeid(t) << ", size " << t.size() << push << dendl;

    for (typename std::multimap<Key, T, L>::const_iterator it = t.begin(); it != t.end(); ++it)
      dout << "[" << (*it).first << "] " << (*it).second << dendl;

    dout << pop;
  }

  return dout;
}

/**
 * @brief Template function <b>dump</b> prints the object pointed to that are
 * contained within a std::multimap to the diagnostic writer.
 *
 * @param dout    a <b>Writer</b> reference to the diagnostic writer to
 *      write the std::multimap to.
 *
 * @param t    a <b>std::multimap</b> of objects.
 *
 * @return    a <b>Writer</b> reference to this object
 */
template <class Key, class T, class L>
Writer &
dump(
  Writer &      dout,
  const std::multimap<Key, T *, L> &  t)
{
  if (dout.shouldPrint()) {
    dout << typeid(t) << ", size " << t.size() << push << dendl;

    for (typename std::multimap<Key, T *, L>::const_iterator it = t.begin(); it != t.end(); ++it)
      dout << "[" << (*it).first << "] " << c_ptr_<T>((*it).second) << std::endl;

    dout << pop;
  }

  return dout;
}

/**
 * @brief Template function <b>dump</b> prints the object contained within a
 * std::set object to the diagnostic writer.
 *
 * @param dout    a <b>Writer</b> reference to the diagnostic writer to
 *      write the std::set to.
 *
 * @param t    a <b>std::set</b> of objects.
 *
 * @return    a <b>Writer</b> reference to this object
 */
template <class Key, class L>
Writer &
dump(
  Writer &      dout,
  const std::set<Key, L> &  t)
{
  if (dout.shouldPrint()) {
    dout << typeid(t) << ", size " << t.size() << push << dendl;

    for (typename std::set<Key, L>::const_iterator it = t.begin(); it != t.end(); ++it)
      dout << (*it) << dendl;

    dout << pop;
  }

  return dout;
}

/**
 * @brief Template function <b>dump</b> prints the object contained within a
 * std::set object to the diagnostic writer.
 *
 * @param dout    a <b>Writer</b> reference to the diagnostic writer to
 *      write the std::set to.
 *
 * @param t    a <b>std::set</b> of objects.
 *
 * @return    a <b>Writer</b> reference to this object
 */
template <class Key, class L>
Writer &
dump(
  Writer &      dout,
  const std::set<Key *, L> &  t)
{
  if (dout.shouldPrint()) {
    dout << typeid(t) << ", size " << t.size() << push << dendl;

    for (typename std::set<Key *, L>::const_iterator it = t.begin(); it != t.end(); ++it)
      dout << c_ptr_<Key>((*it)) << dendl;

    dout << pop;
  }

  return dout;
}

/**
 * @brief Template function <b>dump</b> prints the object contained within a
 * std::multiset object to the diagnostic writer.
 *
 * @param dout    a <b>Writer</b> reference to the diagnostic writer to
 *      write the std::multiset to.
 *
 * @param t    a <b>std::multiset</b> of objects.
 *
 * @return    a <b>Writer</b> reference to this object
 */
template <class Key, class L>
Writer &
dump(
  Writer &      dout,
  const std::multiset<Key, L> &  t)
{
  if (dout.shouldPrint()) {
    dout << typeid(t) << ", size " << t.size() << push << dendl;

    for (typename std::multiset<Key, L>::const_iterator it = t.begin(); it != t.end(); ++it)
      dout << (*it) << dendl;

    dout << pop;
  }

  return dout;
}

/**
 * @brief Template function <b>dump</b> prints the object contained within a
 * std::multiset object to the diagnostic writer.
 *
 * @param dout    a <b>Writer</b> reference to the diagnostic writer to
 *      write the std::multiset to.
 *
 * @param t    a <b>std::multiset</b> of objects.
 *
 * @return    a <b>Writer</b> reference to this object
 */
template <class Key, class L>
Writer &
dump(
  Writer &      dout,
  const std::multiset<Key *, L> &  t)
{
  if (dout.shouldPrint()) {
    dout << typeid(t) << ", size " << t.size() << push << dendl;

    for (typename std::multiset<Key *, L>::const_iterator it = t.begin(); it != t.end(); ++it)
      dout << c_ptr_<Key>((*it)) << dendl;

    dout << pop;
  }

  return dout;
}

// /**
//  * @brief Template <b>dump</b> prints the object contained within a
//  * hash_map to the diagnostic writer.
//  *
//  * @param dout    a <b>Writer</b> reference to the diagnostic writer to
//  *      write the hash_map to.
//  *
//  * @param t    a <b>hash_map</b> of objects.
//  *
//  * @return    a <b>Writer</b> reference to this object
//  */
// template <class Key, class T>
// Writer &
// dump(
//   Writer &      dout,
//   const hash_map<Key, T> &  t)
// {
//   if (dout.shouldPrint()) {
//     dout << typeid(t) << ", size " << t.size() << push << dendl;

//     for (typename hash_map<Key, T>::const_iterator it = t.begin(); it != t.end(); ++it)
//       dout <<  "[" << (*it).first << "] " << (*it).second << dendl;

//     dout << pop;
//   }

//   return dout;
// }

// /**
//  * @brief Template <b>dump</b> prints the object pointed to that are
//  * contained within a hash_map to the diagnostic writer.
//  *
//  * @param dout    a <b>Writer</b> reference to the diagnostic writer to
//  *      write the hash_map to.
//  *
//  * @param t    a <b>hash_map</b> of objects.
//  *
//  * @return    a <b>Writer</b> reference to this object
//  */
// template <class Key, class T>
// Writer &
// dump(
//   Writer &      dout,

//   const hash_map<Key, T *> &  t)
// {
//   if (dout.shouldPrint()) {
//     dout << typeid(t) << ", size " << t.size() << push << dendl;

//     for (typename hash_map<Key, T *>::const_iterator it = t.begin(); it != t.end(); ++it)
//       dout << "[" << (*it).first << "] " << (*it)->second << dendl;

//     dout << pop;
//   }

//   return dout;
// }


template <size_t n>
Writer &operator<<(Writer &dout, const std::bitset<n> &t) {
  if (dout.shouldPrint())
    dout.getStream() << t;

  return dout;
}


/**
 * @brief Member function <b>operator<<</b> write the std::vector object to
 * the diagnostic writer.
 *
 * @param dout    a <b>Writer</b> reference to the diagnostic writer to
 *      write the std::list to.
 *
 * @param t    a <b>std::vector</b> const reference to the std::vector.
 *
 * @return    a <b>Writer</b> reference to this object
 */
template <class T>
Writer &operator<<(Writer &dout, const std::vector<T> &t) {
  return dump(dout, t);
}

/**
 * @brief Template function <b>operator<<</b> write the std::list object to the
 * diagnostic writer.
 *
 * @param dout    a <b>Writer</b> reference to the diagnostic writer to
 *      write the std::list to.
 *
 * @param t    a <b>std::list</b> const reference to the std::list.
 *
 * @return    a <b>Writer</b> reference to this object
 */
template <class T>
Writer &operator<<(Writer &dout, const std::list<T> &t) {
  return dump(dout, t);
}

/**
 * @brief Template function <b>operator<<</b> writes the std::map object to the
 * diagnostic writer.
 *
 * @param dout    a <b>Writer</b> reference to the diagnostic writer to
 *      write the std::map to.
 *
 * @param t    a <b>std::map</b> const reference to the std::map.
 *
 * @return    a <b>Writer</b> reference to this object
 */
template <class Key, class T, class L>
Writer &operator<<(Writer &dout, const std::map<Key, T, L> &t) {
  return dump(dout, t);
}

/**
 * @brief Template function <b>operator<<</b> writes the std::multimap object to the
 * diagnostic writer.
 *
 * @param dout    a <b>Writer</b> reference to the diagnostic writer to
 *      write the std::multimap to.
 *
 * @param t    a <b>std::multimap</b> const reference to the std::multimap.
 *
 * @return    a <b>Writer</b> reference to this object
 */
template <class Key, class T, class L>
Writer &operator<<(Writer &dout, const std::multimap<Key, T, L> &t) {
  return dump(dout, t);
}

/**
 * @brief Template function <b>operator<<</b> writes the std::set object to the
 * diagnostic writer.
 *
 * @param dout    a <b>Writer</b> reference to the diagnostic writer to
 *      write the std::set to.
 *
 * @param t    a <b>std::set</b> const reference to the std::set.
 *
 * @return    a <b>Writer</b> reference to this object
 */
template <class Key, class L>
Writer &operator<<(Writer &dout, const std::set<Key, L> &t) {
  return dump(dout, t);
}

/**
 * @brief Template function <b>operator<<</b> writes the std::multiset object to the
 * diagnostic writer.
 *
 * @param dout    a <b>Writer</b> reference to the diagnostic writer to
 *      write the std::multiset to.
 *
 * @param t    a <b>std::multiset</b> const reference to the std::multiset.
 *
 * @return    a <b>Writer</b> reference to this object
 */
template <class Key, class L>
Writer &operator<<(Writer &dout, const std::multiset<Key, L> &t) {
  return dump(dout, t);
}


/**
 * @brief Function <b>operator<<</b> writes a sierra String object to the diagnostic
 * writer.
 *
 * @param dout		a <b>Writer</b> reference to the diagnostic writer to
 *			write the sierra string to.
 *
 * @param s		a <b>sierra::String</b> const reference to the sierra string to write.
 *
 * @return		a <b>Writer</b> reference to this object
 */
Writer &operator<<(Writer &dout, const sierra::String &s);

/**
 * @brief Function <b>operator<<</b> writes a sierra Identifier object to the
 * diagnostic writer.
 *
 * @param dout		a <b>Writer</b> reference to the diagnostic writer to
 *			write the sierra identifier to.
 *
 * @param s		a <b>sierra::String</b> const reference to the sierra identifier to write.
 *
 * @return		a <b>Writer</b> reference to this object
 */
Writer &operator<<(Writer &dout, const sierra::Identifier &s);

#if defined( STK_HAS_MPI )
/**
 * @brief Function <b>operator<<</b> writes the MPI::Loc<int> type to the output stream.
 *
 * @param dout		a <b>Writer</b> reference to the diagnostic writer to write the
 *			c style string to.
 *
 * @param loc		a <b>MPI::Loc<int></b> const reference to the
 *			MPI::MaxLoc/MPI::MinLoc operator object.
 *
 * @return		a <b>Writer</b> reference to this object
 */
Writer &operator<<(Writer &dout, const sierra::MPI::Loc<int> &loc);

/**
 * @brief Function <b>operator<<</b> writes the sierra::MPI::Loc<double> type to the output stream.
 *
 * @param dout		a <b>Writer</b> reference to the diagnostic writer to write the
 *			c style string to.
 *
 * @param loc		a <b>sierra::MPI::Loc<int></b> const reference to the
 *			sierra::MPI::MaxLoc/sierra::MPI::MinLoc operator object.
 *
 * @return		a <b>Writer</b> reference to this object
 */
Writer &operator<<(Writer &dout, const sierra::MPI::Loc<double> &loc);

/**
 * @brief Function <b>operator<<</b> writes the sierra::MPI::Loc<float> type to the output stream.
 *
 * @param dout		a <b>Writer</b> reference to the diagnostic writer to write the
 *			c style string to.
 *
 * @param loc		a <b>sierra::MPI::Loc<int></b> const reference to the
 *			sierra::MPI::MaxLoc/sierra::MPI::MinLoc operator object.
 *
 * @return		a <b>Writer</b> reference to this object
 */
Writer &operator<<(Writer &dout, const sierra::MPI::Loc<float> &loc);

/**
 * @brief Function <b>operator<<</b> writes the TempLoc type to the output stream.
 *
 * @param dout		a <b>Writer</b> reference to the diagnostic writer to write the
 *			c style string to.
 *
 * @param loc		a <b>sierra::MPI::Loc<int></b> const reference to the
 *			sierra::MPI::MaxLoc/sierra::MPI::MinLoc operator object.
 *
 * @return		a <b>Writer</b> reference to this object
 */
Writer &operator<<(Writer &dout, const sierra::MPI::TempLoc &loc);

#endif // if defined( STK_HAS_MPI )

/**
 * @brief Template function <b>c_ptr_name</b> implements c_ptr_func with the function
 * 'name'.
 *
 * @param t		a <b>T</b> pointer to an object that is call the name
 *			member function.
 *
 * @return		a <b>c_ptr_func_</b> object which contains a member function
 *			pointer to name();
 */
template <class T>
c_ptr_func_<T, const sierra::String &> c_ptr_name(const T *t) {
  return c_ptr_func_<T, const sierra::String &>(t, &T::name);
}

template< class ElementType,
	  class Tag0,
	  class Tag1,
	  class Tag2,
	  class Tag3,
	  class Tag4,
	  class Tag5,
	  class Tag6,
	  class Tag7 >
Writer &
operator<<(Writer &dout, const sierra::Array<ElementType, Tag0, Tag1, Tag2, Tag3, Tag4, Tag5, Tag6, Tag7 > &array)
{
  typedef sierra::Array<ElementType, Tag0, Tag1, Tag2, Tag3, Tag4, Tag5, Tag6, Tag7 > X;

  if (dout.shouldPrint()) {
    std::string type = sierra::demangle(typeid(X).name());
    dout << type.substr(0, type.find(", sierra::TypeListEnd")) << ">" << push << dendl;
    dout.getStream() << array; // ArrayVerbosePrint<X::NumDim>::dump(dout, array.dimension(), array.ptr(), array.stride());
    dout << pop << dendl;
  }
  return dout;
}


template< class ElementType,
	  class Tag0,
	  class Tag1,
	  class Tag2,
	  class Tag3,
	  class Tag4,
	  class Tag5,
	  class Tag6,
	  class Tag7 >
Writer &
operator<<(Writer &dout, const sierra::ArrayContainer<ElementType, Tag0, Tag1, Tag2, Tag3, Tag4, Tag5, Tag6, Tag7 > &array)
{
  typedef sierra::ArrayContainer<ElementType, Tag0, Tag1, Tag2, Tag3, Tag4, Tag5, Tag6, Tag7 > X;

  if (dout.shouldPrint()) {
    std::string type = sierra::demangle(typeid(X).name());
    dout << type.substr(0, type.find(", sierra::TypeListEnd")) << ">" << push << dendl;
    dout.getStream() << array; // ArrayVerbosePrint<X::NumDim>::dump(dout, array.dimension(), array.ptr(), array.stride());
    dout << pop << dendl;
  }
  return dout;
}


template< class ElementType, int Dimension>
Writer &
operator<<(Writer &dout, const sierra::FArray<ElementType, Dimension> &array)
{
  typedef sierra::FArray<ElementType, Dimension> X;

  if (dout.shouldPrint()) {
    dout << sierra::demangle(typeid(X).name()) << push << dendl;
    dout.getStream() << array; // ArrayVerbosePrint<X::NumDim>::dump(dout, array.dimension(), array.ptr(), array.stride());
    dout << pop << dendl;
  }
  return dout;
}

template< class ElementType, int Dimension>
Writer &
operator<<(Writer &dout, const sierra::FArrayContainer<ElementType, Dimension> &array)
{
  typedef sierra::FArrayContainer<ElementType, Dimension> X;

  if (dout.shouldPrint()) {
    dout << sierra::demangle(typeid(X).name()) << push << dendl;
    dout.getStream() << array; // ArrayVerbosePrint<X::NumDim>::dump(dout, array.dimension(), array.ptr(), array.stride());
    dout << pop << dendl;
  }
  return dout;
}

/**
 * @brief Template function <b>dump</b> writes a Mapv_no_delete object to the
 * diagnostic writer.
 *
 * @param dout		a <b>Writer</b> reference to the diagnostic writer to
 *			write the Mapv_no_delete to.
 *
 * @param t		a <b>sierra::String</b> const reference to the Mapv_no_delete to write.
 *
 * @return		a <b>Writer</b> reference to this object
 */
template <class T>
Writer &
dump(
  Writer &			dout,
  const sierra::Mapv_no_delete<T> &	t)
{
  if (dout.shouldPrint()) {
    dout << typeid(t) << ", size " << t.size() << push << dendl;

    for (typename sierra::Mapv_no_delete<T>::const_iterator it = t.begin(); it != t.end(); ++it) {
      dout << "[" << (*it).mapv_key() << "] " << (*it) << dendl;
    }

    dout << pop;
  }

  return dout;
}

/**
 * @brief Template function <b>dump</b> writes the vecmap object to the
 * diagnostic writer.
 *
 * @param dout		a <b>Writer</b> reference to the diagnostic writer to
 *			write the vecmap to.
 *
 * @param t		a <b>vecmap</b> const reference to the vecmap.
 *
 * @return		a <b>Writer</b> reference to this object
 */
template <class Key, class T, class U>
Writer &
dump(
  Writer &			dout,
  const sierra::vecmap<Key, T, U> &	t)
{
  if (dout.shouldPrint()) {
    dout << typeid(t) << ", size " << t.size() << push << dendl;

    for (typename sierra::vecmap<Key, T, U>::const_iterator it = t.begin(); it != t.end(); ++it)
      dout << "[" << (*it).first << "] " << (*it).second << dendl;

    dout << pop;
  }

  return dout;
}

/**
 * @brief Template function <b>dump</b> writes a vecmap of pointers object to the
 * diagnostic writer.
 *
 * @param dout		a <b>Writer</b> reference to the diagnostic writer to
 *			write the vecmap of pointers to.
 *
 * @param t		a <b>vecmap</b> const reference to the vecmap of pointers.
 *
 * @return		a <b>Writer</b> reference to this object
 */
template <class Key, class T, class U>
Writer &
dump(
  Writer &			dout,
  const sierra::vecmap<Key, T *, U> &	t)
{
  if (dout.shouldPrint()) {
    dout << typeid(t) << ", size " << t.size() << push << dendl;

    for (typename sierra::vecmap<Key, T *, U>::const_iterator it = t.begin(); it != t.end(); ++it)
      dout << "[" << (*it).first << "] " << *(*it).second << dendl;

    dout << pop;
  }

  return dout;
}

/**
 * @brief Template function <b>dump</b> writes a vecmap of pointers object to the
 * dignostic writer.
 *
 * @param dout		a <b>Writer</b> reference to the diagnostic writer to
 *			write the vecmap of pointers to.
 *
 * @param t		a <b>vecmap</b> const reference to the vecmap of pointers.
 *
 * @return		a <b>Writer</b> reference to this object
 */
template <class Key, class T, class U>
Writer &
dump(
  Writer &				dout,
  const sierra::vecmap<Key *, T *, U> &	t)
{
  if (dout.shouldPrint()) {
    dout << typeid(t) << ", size " << t.size() << push << dendl;

    for (typename sierra::vecmap<Key *, T *, U>::const_iterator it = t.begin(); it != t.end(); ++it)
      dout << "[" << (*it).first << "] " << *(*it).second << dendl;

    dout << pop;
  }

  return dout;
}

/**
 * @brief Template function <b>dump</b> writes a vecset object to the diagnostic
 * writer.
 *
 * @param dout		a <b>Writer</b> reference to the diagnostic writer to
 *			write the vecset to.
 *
 * @param t		a <b>vecset</b> const reference to the vecset.
 *
 * @return		a <b>Writer</b> reference to this object
 */
template <class T, class U>
Writer &
dump(
  Writer &		dout,
  const sierra::vecset<T, U> &	t)
{
  if (dout.shouldPrint()) {
    dout << typeid(t) << ", size " << t.size() << push << dendl;

    int i = 0;
    for (typename sierra::vecset<T, U>::const_iterator it = t.begin(); it != t.end(); ++it, ++i)
      dout << "[" << i << "] " << (*it) << dendl;

    dout << pop;
  }

  return dout;
}

/**
 * @brief Template function <b>dump</b> writes a vecset of pointers object to the
 * diagnostic writer.
 *
 * @param dout		a <b>Writer</b> reference to the diagnostic writer to
 *			write the vecset of pointers to.
 *
 * @param t		a <b>vecset</b> const reference to the vecset of pointers.
 *
 * @return		a <b>Writer</b> reference to this object
 */
template <class T, class U>
Writer &
dump(
  Writer &			dout,
  const sierra::vecset<T *, U> &	t)
{
  if (dout.shouldPrint()) {
    dout << typeid(t) << ", size " << t.size() << push << dendl;

    int i = 0;
    for (typename sierra:: vecset<T *, U>::const_iterator it = t.begin(); it != t.end(); ++it, ++i)
      dout << "[" << i << "] " << *(*it) << dendl;

    dout << pop;
  }

  return dout;
}

/**
 * @brief Template <b>dump</b> writes a MapvNode object to the diagnostic writer.
 *
 * @param dout		a <b>Writer</b> reference to the diagnostic writer to write the
 *			Mapvnod to.
 *
 * @param t		a <b>MapvNode</b> const reference to the MapvNode to write.
 *
 * @return		a <b>Writer</b> reference to this object
 */
template <class T, class U>
Writer &
dump(
  Writer &			dout,
  const sierra::MapvNode<T, U> &	t)
{
  if (dout.shouldPrint()) {
    dout << typeid(t) << ", " << t.mapv_key();
  }

  return dout;
}

/**
 * @brief Template function <b>dump</b> writes a Mapv object to the diagnostic writer.
 *
 * @param dout		a <b>Writer</b> reference to the diagnostic writer to write the
 *			Mapv to.
 *
 * @param t		a <b>std::vector</b> const reference to the Mapv to write.
 *
 * @return		a <b>Writer</b> reference to this object
 */
template <class T, class U>
Writer &
dump(
  Writer &		dout,
  const sierra::Mapv<T, U> &	t)
{
  if (dout.shouldPrint()) {
    dout << typeid(t) << ", size " << t.size() << push << dendl;

    for (typename sierra::Mapv<T, U>::const_iterator it = t.begin(); it != t.end(); ++it) {
      dout << "[" << (*it).mapv_key() << "] " << (*it) << dendl;
    }

    dout << pop;
  }

  return dout;
}

/**
 * @brief Member function <b>operator<<</b> writer a vecset object the
 * diagnostic writer.
 *
 * @param dout		a <b>Writer</b> reference to the diagnostic writer to
 *			write the vecset to.
 *
 * @param t		a <b>vecset</b> const reference to the vecset.
 *
 * @return		a <b>Writer</b> reference to this object
 */
template <class T, class U>
Writer &operator<<(Writer &dout, const sierra::vecset<T, U> &t) {
  return dump(dout, t);
}

/**
 * @brief Template function <b>operator<<</b> writes the vecmap object to the
 * diagnostic writer.
 *
 * @param dout		a <b>Writer</b> reference to the diagnostic writer to
 *			write the vecmap to.
 *
 * @param t		a <b>vecmap</b> const reference to the vecmap.
 *
 * @return		a <b>Writer</b> reference to this object
 */
template <class Key, class T, class U>
Writer &operator<<(Writer &dout, const sierra::vecmap<Key, T, U> &t) {
  return dump(dout, t);
}

/**
 * @brief Member function <b>operator<<</b> writes a vecset of pointers object to
 * the diagnostic writer.
 *
 * @param dout		a <b>Writer</b> reference to the diagnostic writer to
 *			write the vecset of pointers to.
 *
 * @param t		a <b>vecset</b> const reference to the vecset of pointers.
 *
 * @return		a <b>Writer</b> reference to this object
 */
template <class T, class U>
Writer &operator<<(Writer &dout, const sierra::vecset<T *, U> &t) {
  return dump(dout, t);
}

/**
 * @brief Template function <b>operator<<</b> writea a vecmap of pointers with key
 * pointers object to the diagnostic writer.
 *
 * @param dout		a <b>Writer</b> reference to the diagnostic writer to
 *			write the vecmap of pointers to.
 *
 * @param t		a <b>vecmap</b> const reference to the vecmap of pointers.
 *
 * @return		a <b>Writer</b> reference to this object
 */
template <class Key, class T, class U>
Writer &operator<<(Writer &dout, const sierra::vecmap<Key *, T *, U> &t) {
  return dump(dout, t);
}

/**
 * @brief Template function <b>operator<<</b> writes a Mpav_no_delete object to the
 * diagnostic writer.
 *
 * @param dout		a <b>Writer</b> reference to the diagnostic writer to
 *			write the Mapv_no_delete to.
 *
 * @param t		a <b>vecmap</b> const reference to the Mapv_no_delete.
 *
 * @return		a <b>Writer</b> reference to this object
 */
template <class T>
Writer &operator<<(Writer &dout, const sierra::Mapv_no_delete<T> &t) {
  return dump(dout, t);
}

/**
 * @brief Member function <b>operator<<</b> writes a Mapv object to the diagnostic
 * writer.
 *
 * @param dout		a <b>Writer</b> reference to the diagnostic writer to
 *			write the Mapv to.
 *
 * @param t		a <b>vecmap</b> const reference to the Mapv.
 *
 * @return		a <b>Writer</b> reference to this object
 */
template <class T, class U>
Writer &operator<<(Writer &dout, const sierra::Mapv<T, U> &t) {
  return dump(dout, t);
}

/**
 * @brief Template function <b>operator<<</b> writes a MapvNode object to the
 * diagnostic writer.
 *
 * @param dout		a <b>Writer</b> reference to the diagnostic writer to
 *			write the MapvNode to.
 *
 * @param t		a <b>MapvNode</b> const reference to the MapvNode.
 *
 * @return		a <b>Writer</b> reference to this object
 */
template <class T, class U>
Writer &operator<<(Writer &dout, const sierra::MapvNode<T, U> &t) {
  return dump(dout, t);
}

/**
 * @brief Template function <b>operator<<</b> writes a vecmap of pointers object to
 * the dignostic writer.
 *
 * @param dout		a <b>Writer</b> reference to the diagnostic writer to
 *			write the vecmap of pointers to.
 *
 * @param t		a <b>vecmap</b> const reference to the vecmap of pointers.
 *
 * @return		a <b>Writer</b> reference to this object
 */
template <class Key, class T, class U>
Writer &operator<<(Writer &dout, const sierra::vecmap<Key, T *, U> &t) {
  return dump(dout, t);
}


///
/// @}
///

} // namespace diag
} // namespace stk

#endif // STK_UTIL_DIAG_WRITEREXT_HPP
