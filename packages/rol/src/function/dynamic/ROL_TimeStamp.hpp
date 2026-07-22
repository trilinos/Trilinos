// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#pragma once
#ifndef ROL_TIMESTAMP_HPP
#define ROL_TIMESTAMP_HPP


#include <vector>

#include "ROL_Ptr.hpp"

/** @ingroup dynamic_group
    \class ROL::TimeStamp
    \brief Contains local time step information.
*/

namespace ROL {


template<typename> struct TimeStamp;

// Alias for pointer vector of TimeStamp objects
template<typename Real>
using TimeStampsPtr = Ptr<std::vector<TimeStamp<Real>>>;


template<typename Real>
struct TimeStamp {

  using size_type = typename std::vector<Real>::size_type;

  size_type         k;   // Time-step number  
  std::vector<Real> t;   // Time points for this time step

  TimeStamp& operator= ( const TimeStamp& ts ) {
    k = ts.k; t = ts.t;
    return *this;
  }

  /** \brief Create a vector of uniform TimeStamp objects for the interval [t_initial,t_final] 
             where each step has local time points based on t_ref defined on the standard 
             interval [0,1] */
  static TimeStampsPtr<Real> make_uniform( Real t_initial, 
                                           Real t_final, 
                                           const std::vector<Real>& t_ref, 
                                           size_type num_steps ) {

    auto timeStamp = ROL::makePtr<std::vector<ROL::TimeStamp<Real>>>(num_steps);
  
    Real dt = (t_final-t_initial)/num_steps; // size of each time step
    size_type nt = t_ref.size();             // number of time points per step
  
    for( size_type k=0; k<num_steps; ++k ) {
      (*timeStamp)[k].t.resize(nt);
      for( size_type l=0; l<nt; ++l ) (*timeStamp)[k].t[l] = dt*(k+t_ref[l]);
      (*timeStamp)[k].k = k;
    }
    return timeStamp;
  }
};

} // namespace ROL


#endif // ROL_TIMESTAMP_HPP

