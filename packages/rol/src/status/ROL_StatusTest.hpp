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

#ifndef ROL_STATUSTEST_H
#define ROL_STATUSTEST_H

#include "ROL_Types.hpp"

/** \class ROL::StatusTest
    \brief Provides an interface to check status of optimization algorithms.
*/


namespace ROL {


template <class Real>
class StatusTest {
private:

  Real gtol_;
  Real stol_;
  int  max_iter_;

public:

  virtual ~StatusTest() {}

  StatusTest( Real gtol = 1.e-6, Real stol = 1.e-12, int max_iter = 100 ) :  
    gtol_(gtol), stol_(stol), max_iter_(max_iter) {}

  /** \brief Check algorithm status.
  */
  virtual bool check( AlgorithmState<Real> &state ) {
     if ( (state.gnorm > this->gtol_) && 
          (state.snorm > this->stol_) && 
          (state.iter  < this->max_iter_) ) {
       return true;
     }
     else {
       return false;
     }
  }

}; // class StatusTest


template <class Real>
class StatusTestSQP : public StatusTest<Real> {
private:

  Real gtol_;
  Real ctol_;
  Real stol_;
  int  max_iter_;

public:

  virtual ~StatusTestSQP() {}

  StatusTestSQP( Real gtol = 1e-6, Real ctol = 1e-6, Real stol = 1e-12, int max_iter = 100 ) :  
    gtol_(gtol), ctol_(ctol), stol_(stol), max_iter_(max_iter) {}

  /** \brief Check algorithm status.
  */
  virtual bool check( AlgorithmState<Real> &state ) {
     if ( ((state.gnorm > this->gtol_) || (state.cnorm > this->ctol_)) && 
          (state.snorm > this->stol_) && 
          (state.iter  < this->max_iter_) ) {
       return true;
     }
     else {
       return false;
     }
  }

}; // class StatusTestSQP

} // namespace ROL

#endif
