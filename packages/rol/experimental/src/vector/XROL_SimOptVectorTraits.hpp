
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


namespace XROL {

template<class A, class B> using pair = std::pair<A,B>;
template<class T>          using ptr  = std::shared_ptr<T>;


template<class Sim, class Opt>
using VSO = pair<ptr<Sim>,ptr<Opt>>;

template<class U, class Z>
struct VectorCreationTraits<VSO<U,Z>> {

  using V = VSO<U,Z>;

  

  static auto clone( const V& x ) {

    auto up = ST::clone(*std::get<0>(x));
    auto zp = OT::clone(*std::get<1>(x));
              
    return std::make_shared<V>(std::make_pair(up,zp));
  }

  // TODO: Implement basis
};

  
template<class U, class Z>
struct VectorSpaceTraits<VSO<U,Z>> {

  using ST = VectorTraits<U>;
  using OT = VectorTraits<Z>;

  

  static auto dimension( const V& x ) {
    auto up = std::get<0>(x);
    auto zp = std::get<1>(x);
    return ST::dimension(*up) + OT::dimension(*zp);
  }

  static void dual( V& xdual, const V& x ) {
    auto up  = std::get<0>(x);
    auto zp  = std::get<1>(x);
    auto udp = std::get<0>(xdual);
    auto zdp = std::get<1>(xdual);

    ST::dual(*udp,*up);
    OT::dual(*udp,*up);
    
  }

};




} // namespace XROL

