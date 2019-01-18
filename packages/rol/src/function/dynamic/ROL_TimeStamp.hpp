// @HEADER
// ************************************************************************
//
//               Rapid Optimization Library (ROL) Package
//                 Copyright (2014) Sandia Corporation
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
// Questions? Contact lead developers:
//              Drew Kouri   (dpkouri@sandia.gov) and
//              Denis Ridzal (dridzal@sandia.gov)
//
// ************************************************************************
// @HEADER


#pragma once
#ifndef ROL_TIMESTAMP_HPP
#define ROL_TIMESTAMP_HPP


#include <vector>

/** @ingroup dynamic_group
    \class ROL::TimeStamp
    \brief Contains local time step information.
*/

namespace ROL {

template<typename Real>
struct TimeStamp {

  using size_type = typename std::vector<Real>::size_type;

  size_type         k;   // Time-step number  
  std::vector<Real> t;   // Time points for this time step

  /** \brief Create a vector of uniform TimeStamp objects for the interval [t_initial,t_final] where each step
             has local time points based on t_ref defined on the standard interval [0,1] */
  static ROL::Ptr<std::vector<TimeStamp<Real>>> make_uniform( Real t_initial, Real t_final, const std::vector<Real>& t_ref, size_type num_steps ) {

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

