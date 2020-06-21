// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef PerceptBoostArray_hpp
#define PerceptBoostArray_hpp

#if defined(__INTEL_COMPILER) && (__INTEL_COMPILER == 1210)
#pragma warning disable 2196 2536 279
#endif

#include <iostream>
#include <array>

namespace percept {

  template<class T, std::size_t N>
  inline std::ostream& operator<< (std::ostream& out, const std::array<T,N>& coll)
  {
    typedef std::array<T,N> Array;
    typename Array::const_iterator pos;

    for (pos=coll.begin(); pos!=coll.end(); ++pos) {
      out << *pos << ' ';
    }
    out << std::endl;
    return out;
  }

}
#endif
