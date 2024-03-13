// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef SIERRA_Akri_Utility_h
#define SIERRA_Akri_Utility_h

#include <memory>
#include <vector>
#include <limits>
#include <cmath>

namespace krino {
namespace utility {

  inline bool sign_change( double f1, double f2 ) {
    return ( (f1 < 0.) ? (f2 >= 0.) : (f2 < 0.) ); // GOMA sign convention
    //return ( (f1 > 0.) ? (f2 <= 0.) : (f2 > 0.) ); // Marching cubes sign convention
  }

  inline int sign( double f ) {
    return ( (f < 0.) ? -1 : 1 ); // GOMA sign convention
    //return ( (f > 0.) ? 1 : -1 ); // Marching cubes sign convention
  }

  inline bool is_less(double f1, double f2, double tol)      { return (f2-f1 > tol*(std::fabs(f1)+std::fabs(f2))); }
  inline bool is_more(double f1, double f2, double tol)      { return (is_less(f2,f1,tol)); }
  inline bool is_not_equal(double f1, double f2, double tol) { return (is_less(f1,f2,tol) || is_less(f2,f1,tol)); }
  inline bool is_equal(double f1, double f2, double tol)     { return (!is_less(f1,f2,tol) && !is_less(f2,f1,tol)); }

  inline bool is_less(double f1, double f2)      { return is_less(f1,f2,100.0*std::numeric_limits<double>::epsilon()); }
  inline bool is_more(double f1, double f2)      { return is_more(f1,f2,100.0*std::numeric_limits<double>::epsilon()); }
  inline bool is_not_equal(double f1, double f2) { return is_not_equal(f1,f2,100.0*std::numeric_limits<double>::epsilon()); }
  inline bool is_equal(double f1, double f2)     { return is_equal(f1,f2,100.0*std::numeric_limits<double>::epsilon()); }

  inline bool is_less(double f1, int i1, double f2, int i2, double tol) { return is_less(f1,f2,tol) ? true : (is_more(f1,f2,tol) ? false : (i1<i2)); }
  inline bool is_more(double f1, int i1, double f2, int i2, double tol) { return (is_less(f2,i2,f1,i1,tol)); }
  inline bool is_less(double f1, int i1, double f2, int i2) { return is_less(f1,i1,f2,i2,100.0*std::numeric_limits<double>::epsilon()); }
  inline bool is_more(double f1, int i1, double f2, int i2) { return is_more(f1,i1,f2,i2,100.0*std::numeric_limits<double>::epsilon()); }

  template<class T>
  inline void
  free_all(std::vector<T> & vec) {std::vector<T> empty_vec; vec.swap(empty_vec);}

} // namespace utility
} // namespace krino

#endif // SIERRA_Akri_Utility_h
