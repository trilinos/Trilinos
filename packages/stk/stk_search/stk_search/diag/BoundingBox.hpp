/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef stk_search_diag_BoundingBox_hpp
#define stk_search_diag_BoundingBox_hpp

#include <ostream>
#include <stk_search/BoundingBox.hpp>

namespace stk { namespace search { namespace box {

template <class Key, class Data, int Dim>
std::ostream & operator<<(std::ostream & out, PointBoundingBox<Key,Data,Dim> const& box)
{
  out << "min: ( ";
  for (int i=0; i<Dim; ++i) {
    out << box.lower(i) << ", ";
  }
  out << "\b\b )  max: ( ";
  for (int i=0; i<Dim; ++i) {
    out << box.upper(i) << ", ";
  }
  out << "\b\b )";

  return out;
}

template <class Key, class Data, int Dim>
std::ostream & operator<<(std::ostream & out, SphereBoundingBox<Key,Data,Dim> const& box)
{
  out << "min: ( ";
  for (int i=0; i<Dim; ++i) {
    out << box.lower(i) << ", ";
  }
  out << "\b\b )  max: ( ";
  for (int i=0; i<Dim; ++i) {
    out << box.upper(i) << ", ";
  }
  out << "\b\b )";

  return out;
}

template <class Key, class Data, int Dim>
std::ostream & operator<<(std::ostream & out, AxisAlignedBoundingBox<Key,Data,Dim> const& box)
{
  out << "min: ( ";
  for (int i=0; i<Dim; ++i) {
    out << box.lower(i) << ", ";
  }
  out << "\b\b )  max: ( ";
  for (int i=0; i<Dim; ++i) {
    out << box.upper(i) << ", ";
  }
  out << "\b\b )";

  return out;
}

}}} //namespace stk::search::box

#endif //stk_search_diag_BoundingBox_hpp
