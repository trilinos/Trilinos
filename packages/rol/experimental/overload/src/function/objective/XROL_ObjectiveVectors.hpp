
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

#include "XROL.hpp"

/** 
    \file  XROL_ObjectiveVectors.hpp
    \brief Simple containers for primal and dual
           vectors needed to store intermediate 
           values
*/

namespace XROL {


template<class X> 
class ObjectiveVectors {

  template<class T> using vector     = std::vector<T>;
  template<class T> using unique_ptr = std::unique_ptr<T>;

  using size_type = typename vector<X>::size_type;
 
public:
 
  void allocate_x( const X& x, size_type d=1 ) {
    if( !is_xprim_allocated_ ) {
      xprim_.resize(d);
      for( size_type i=0; i<d; ++i ) 
        xprim_.push_back(std::move(clone(x)));
      is_xprim_allocated_ = true;
    }
  }

  void allocate_g( const dual_t<X>& g, size_type d=1 ) {
    if( !is_xdual_allocated_ ) {
      xdual_.resize(d);
      for( size_type i=0; i<d; ++i ) 
        xdual_.push_back(std::move(clone(g)));
      is_xdual_allocated_ = true;
    }
  }

  X& x( size_type i = 0 ) { 
    if( !is_xprim_allocated_ ) 
      throw UnallocatedVector(": optimization space"); 
    return *(xprim_[i]); 
  }


  dual_t<X>& g( size_type i = 0 ) { 
    if( !is_xdual_allocated_ )
      throw UnallocatedVector(": optimization dual space");
    return *(xdual_[i]); 
  }

private:
 
  vector<unique_ptr<X>>         xprim_;
  vector<unique_ptr<dual_t<X>>> xdual_;

  bool is_xprim_allocated_; 
  bool is_xdual_allocated_; 
 
};

} // namespace XROL

