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

#include <stk_util/diag/Writer.hpp>

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
 * @param dout		a <b>Writer</b> reference to the diagnostic writer to
 *			write the std::type_info object to.
 *
 * @param t		a <b>std::type_info</b> const reference to the std::typeinfo
 *			object.
 *
 * @return		a <b>Writer</b> reference to this object
 */
Writer &operator<<(Writer &dout, const std::type_info &t);

/**
 * @brief Template function <b>operator<<</b> writes an std::auto_ptr object
 * address and content to the diagnostic writer.
 *
 * @param dout		a <b>Writer</b> reference to the diagnostic writer to
 *			write the std::auto_ptr object.
 *
 * @param t		a <b>std::auto_ptr</b> const reference to the object.
 *
 * @return		a <b>Writer</b> reference to this object
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
 * @param dout		a <b>Writer</b> reference to the diagnostic writer to
 *			write the std::pair members to.
 *
 * @param pair		a <b>std::pair</b> const reference to the pair of objects.
 *
 * @return		a <b>Writer</b> reference to this object
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
 * @param dout		a <b>Writer</b> reference to the diagnostic writer to
 *			write the std::vector to.
 *
 * @param t		a <b>std::vector</b> of objects.
 *
 * @return		a <b>Writer</b> reference to this object
 */
template <class T>
Writer &
dump(
  Writer &			dout,
  const std::vector<T> &	t)
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
 * @param dout		a <b>Writer</b> reference to the diagnostic writer to
 *			write the std::vector to.
 *
 * @param t		a <b>std::vector</b> of objects.
 *
 * @return		a <b>Writer</b> reference to this object
 */
template <class T>
Writer &
dump(
  Writer &			dout,
  const std::vector<T *> &	t)
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
 * @param dout		a <b>Writer</b> reference to the diagnostic writer to
 *			write the std::list to.
 *
 * @param t		a <b>std::list</b> of objects.
 *
 * @return		a <b>Writer</b> reference to this object
 */
template <class T>
Writer &
dump(
  Writer &			dout,
  const std::list<T> &		t)
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
 * @param dout		a <b>Writer</b> reference to the diagnostic writer to
 *			write the std::list to.
 *
 * @param t		a <b>std::list</b> of objects.
 *
 * @return		a <b>Writer</b> reference to this object
 */
template <class T>
Writer &
dump(
  Writer &			dout,
  const std::list<T *> &	t)
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
 * @param dout		a <b>Writer</b> reference to the diagnostic writer to
 *			write the std::map to.
 *
 * @param t		a <b>std::map</b> of objects.
 *
 * @return		a <b>Writer</b> reference to this object
 */
template <class Key, class T, class L>
Writer &
dump(
  Writer &			dout,
  const std::map<Key, T, L> &	t)
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
 * @param dout		a <b>Writer</b> reference to the diagnostic writer to
 *			write the std::map to.
 *
 * @param t		a <b>std::map</b> of objects.
 *
 * @return		a <b>Writer</b> reference to this object
 */
template <class Key, class T, class L>
Writer &
dump(
  Writer &			dout,
  const std::map<Key, T *, L> &	t)
{
  if (dout.shouldPrint()) {
    dout << typeid(t) << ", size " << t.size() << push << dendl;

    for (typename std::map<Key, T *, L>::const_iterator it = t.begin(); it != t.end(); ++it)
      dout << "[" << (*it).first << "] " << c_ptr_<T>((*it)->second) << std::endl;

    dout << pop;
  }

  return dout;
}

/**
 * @brief Template function <b>dump</b> prints the object contained within a
 * std::multimap object to the diagnostic writer.
 *
 * @param dout		a <b>Writer</b> reference to the diagnostic writer to
 *			write the std::multimap to.
 *
 * @param t		a <b>std::multimap</b> of objects.
 *
 * @return		a <b>Writer</b> reference to this object
 */
template <class Key, class T, class L>
Writer &
dump(
  Writer &			dout,
  const std::multimap<Key, T, L> &	t)
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
 * @param dout		a <b>Writer</b> reference to the diagnostic writer to
 *			write the std::multimap to.
 *
 * @param t		a <b>std::multimap</b> of objects.
 *
 * @return		a <b>Writer</b> reference to this object
 */
template <class Key, class T, class L>
Writer &
dump(
  Writer &			dout,
  const std::multimap<Key, T *, L> &	t)
{
  if (dout.shouldPrint()) {
    dout << typeid(t) << ", size " << t.size() << push << dendl;

    for (typename std::multimap<Key, T *, L>::const_iterator it = t.begin(); it != t.end(); ++it)
      dout << "[" << (*it).first << "] " << c_ptr_<T>((*it)->second) << std::endl;

    dout << pop;
  }

  return dout;
}

/**
 * @brief Template function <b>dump</b> prints the object contained within a
 * std::set object to the diagnostic writer.
 *
 * @param dout		a <b>Writer</b> reference to the diagnostic writer to
 *			write the std::set to.
 *
 * @param t		a <b>std::set</b> of objects.
 *
 * @return		a <b>Writer</b> reference to this object
 */
template <class Key, class L>
Writer &
dump(
  Writer &			dout,
  const std::set<Key, L> &	t)
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
 * @param dout		a <b>Writer</b> reference to the diagnostic writer to
 *			write the std::set to.
 *
 * @param t		a <b>std::set</b> of objects.
 *
 * @return		a <b>Writer</b> reference to this object
 */
template <class Key, class L>
Writer &
dump(
  Writer &			dout,
  const std::set<Key *, L> &	t)
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
 * @param dout		a <b>Writer</b> reference to the diagnostic writer to
 *			write the std::multiset to.
 *
 * @param t		a <b>std::multiset</b> of objects.
 *
 * @return		a <b>Writer</b> reference to this object
 */
template <class Key, class L>
Writer &
dump(
  Writer &			dout,
  const std::multiset<Key, L> &	t)
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
 * @param dout		a <b>Writer</b> reference to the diagnostic writer to
 *			write the std::multiset to.
 *
 * @param t		a <b>std::multiset</b> of objects.
 *
 * @return		a <b>Writer</b> reference to this object
 */
template <class Key, class L>
Writer &
dump(
  Writer &			dout,
  const std::multiset<Key *, L> &	t)
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
//  * @param dout		a <b>Writer</b> reference to the diagnostic writer to
//  *			write the hash_map to.
//  *
//  * @param t		a <b>hash_map</b> of objects.
//  *
//  * @return		a <b>Writer</b> reference to this object
//  */
// template <class Key, class T>
// Writer &
// dump(
//   Writer &			dout,
//   const hash_map<Key, T> &	t)
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
//  * @param dout		a <b>Writer</b> reference to the diagnostic writer to
//  *			write the hash_map to.
//  *
//  * @param t		a <b>hash_map</b> of objects.
//  *
//  * @return		a <b>Writer</b> reference to this object
//  */
// template <class Key, class T>
// Writer &
// dump(
//   Writer &			dout,

//   const hash_map<Key, T *> &	t)
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
 * @param dout		a <b>Writer</b> reference to the diagnostic writer to
 *			write the std::list to.
 *
 * @param t		a <b>std::vector</b> const reference to the std::vector.
 *
 * @return		a <b>Writer</b> reference to this object
 */
template <class T>
Writer &operator<<(Writer &dout, const std::vector<T> &t) {
  return dump(dout, t);
}

/**
 * @brief Template function <b>operator<<</b> write the std::list object to the
 * diagnostic writer.
 *
 * @param dout		a <b>Writer</b> reference to the diagnostic writer to
 *			write the std::list to.
 *
 * @param t		a <b>std::list</b> const reference to the std::list.
 *
 * @return		a <b>Writer</b> reference to this object
 */
template <class T>
Writer &operator<<(Writer &dout, const std::list<T> &t) {
  return dump(dout, t);
}

/**
 * @brief Template function <b>operator<<</b> writes the std::map object to the
 * diagnostic writer.
 *
 * @param dout		a <b>Writer</b> reference to the diagnostic writer to
 *			write the std::map to.
 *
 * @param t		a <b>std::map</b> const reference to the std::map.
 *
 * @return		a <b>Writer</b> reference to this object
 */
template <class Key, class T, class L>
Writer &operator<<(Writer &dout, const std::map<Key, T, L> &t) {
  return dump<Key, T, L>(dout, t);
}

/**
 * @brief Template function <b>operator<<</b> writes the std::multimap object to the
 * diagnostic writer.
 *
 * @param dout		a <b>Writer</b> reference to the diagnostic writer to
 *			write the std::multimap to.
 *
 * @param t		a <b>std::multimap</b> const reference to the std::multimap.
 *
 * @return		a <b>Writer</b> reference to this object
 */
template <class Key, class T, class L>
Writer &operator<<(Writer &dout, const std::multimap<Key, T, L> &t) {
  return dump<Key, T, L>(dout, t);
}

/**
 * @brief Template function <b>operator<<</b> writes the std::set object to the
 * diagnostic writer.
 *
 * @param dout		a <b>Writer</b> reference to the diagnostic writer to
 *			write the std::set to.
 *
 * @param t		a <b>std::set</b> const reference to the std::set.
 *
 * @return		a <b>Writer</b> reference to this object
 */
template <class Key, class L>
Writer &operator<<(Writer &dout, const std::set<Key, L> &t) {
  return dump<Key, L>(dout, t);
}

/**
 * @brief Template function <b>operator<<</b> writes the std::multiset object to the
 * diagnostic writer.
 *
 * @param dout		a <b>Writer</b> reference to the diagnostic writer to
 *			write the std::multiset to.
 *
 * @param t		a <b>std::multiset</b> const reference to the std::multiset.
 *
 * @return		a <b>Writer</b> reference to this object
 */
template <class Key, class L>
Writer &operator<<(Writer &dout, const std::multiset<Key, L> &t) {
  return dump<Key, L>(dout, t);
}

///
/// @}
///

} // namespace diag
} // namespace stk

#endif // STK_UTIL_DIAG_WRITEREXT_HPP
