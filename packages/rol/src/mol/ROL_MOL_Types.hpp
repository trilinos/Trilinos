// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_MOL_TYPES_HPP
#define ROL_MOL_TYPES_HPP

#include "ROL_Ptr.hpp"
#include "ROL_Types.hpp"
#include "ROL_Vector.hpp"
#include <string>
#include <vector>

namespace ROL {

enum EMOType{
  MOTYPE_CC = 0,
  MOTYPE_NBI,
  MOTYPE_LAST
};

inline std::string EMOTypeToString(EMOType tr) {
  std::string retString;
  switch(tr) {
    case MOTYPE_CC:   retString = "Convex Combination";           break;
    case MOTYPE_NBI:  retString = "Normal Boundary Intersection"; break;
    case MOTYPE_LAST: retString = "Last Type (Dummy)";            break;
    default:          retString = "INVALID EMOType";
  }
  return retString;
}

inline int isValidMOType(EMOType ls) {
  return( (ls == MOTYPE_CC) ||
          (ls == MOTYPE_NBI) );
}

inline EMOType & operator++(EMOType &type) {
  return type = static_cast<EMOType>(type+1);
}

inline EMOType operator++(EMOType &type, int) {
  EMOType oldval = type;
  ++type;
  return oldval;
}

inline EMOType & operator--(EMOType &type) {
  return type = static_cast<EMOType>(type-1);
}

inline EMOType operator--(EMOType &type, int) {
  EMOType oldval = type;
  --type;
  return oldval;
}

inline EMOType StringToEMOType(std::string s) {
  s = removeStringFormat(s);
  for ( EMOType st = MOTYPE_CC; st < MOTYPE_LAST; ++st ) {
    if ( !s.compare(removeStringFormat(EMOTypeToString(st))) ) {
      return st;
    }
  }
  return MOTYPE_LAST;
}

template<typename Real>
struct ParetoData {
public:
  const Ptr<Vector<Real>> solution; // Solution vector
  const std::vector<Real> lambda;   // Weight vector
  const std::vector<Real> values;   // Objective function values
  const EExitStatus exitstatus;     // Reason algorithm stopped
  ParetoData(const Ptr<Vector<Real>>& x, const std::vector<Real>& l, const std::vector<Real>& v, EExitStatus s)
    : solution(x), lambda(l), values(v), exitstatus(s) {}
};

}

#endif
