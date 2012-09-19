/*------------------------------------------------------------------------*/
/*                 Copyright 2010 - 2011 Sandia Corporation.              */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef STK_UTIL_UTIL_RANGE_HPP
#define STK_UTIL_UTIL_RANGE_HPP

#include <boost/range.hpp>

namespace stk {

template <class It>
inline
typename boost::iterator_range<It>::const_iterator const_begin(const boost::iterator_range<It> &range) {
  return boost::begin(range);
}

template <class It>
inline
typename boost::iterator_range<It>::const_iterator const_end(const boost::iterator_range<It> &range) {
  return boost::end(range);
}

template <class It>
inline
typename boost::iterator_range<It>::iterator begin(const boost::iterator_range<It> &range) {
  return boost::begin(range);
}

template <class It>
inline
typename boost::iterator_range<It>::iterator end(const boost::iterator_range<It> &range) {
  return boost::end(range);
}

template <class T>
inline
typename boost::iterator_range<T>::const_iterator const_begin(const boost::sub_range<T> &range) {
  return boost::begin(range);
}

template <class T>
inline
typename boost::iterator_range<T>::const_iterator const_end(const boost::sub_range<T> &range) {
  return boost::end(range);
}

template <class T>
inline
typename boost::iterator_range<T>::iterator begin(const boost::sub_range<T> &range) {
  return boost::begin(range);
}

template <class T>
inline
typename boost::iterator_range<T>::iterator end(const boost::sub_range<T> &range) {
  return boost::end(range);
}

} // namespace stk

#endif // STK_UTIL_UTIL_RANGE_HPP
