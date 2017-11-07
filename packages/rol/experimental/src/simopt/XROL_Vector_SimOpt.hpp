
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

#include "XROL_VectorTraits.hpp"


namespace XROL {

template<class S, class O> 
using Vector_SimOpt = std::pair<std::unique_ptr<S>,
                                std::unique_ptr<O>>;

// Traits specialization

template<class S, class O>
struct VectorIndex<Vector_SimOpt<S,O>> {
  using type = std::common_type_t<index_t<S>,index_t<O>>;
};

template<class S, class O>
struct VectorElement<Vector_SimOpt<S,O>> {
  using type = std::common_type_t<element_t<S>,element_t<O>>;
};

template<class S, class O>
struct VectorMagnitude<Vector_SimOpt<S,O>> {
  using type = std::common_type_t<magnitude_t<S>,magnitude_t<O>>;
};


template<class S, class O>
struct VectorDual<Vector_SimOpt<S,O>> {
  using type = Vector_SimOpt<dual_t<S>,dual_t<O>>;
};

template<class S, class O>
struct implements_elementwise<Vector_SimOpt<S,O>> : 
  std::integral_constant<bool,implements_elementwise<S>() 
                          and implements_elementwise<O>()>{}; 

template<class S, class O>
struct implements_core<Vector_SimOpt<S,O>> : 
  std::integral_constant<bool,implements_core<S>() 
                          and implements_core<O>()>{}; 




// Functions

template<class S, class O>
std::unique_ptr<Vector_SimOpt<S,O>> 
clone( const Vector_SimOpt<S,O>& v ) {
  using namespace std;
  auto ptr = make_unique<Vector_SimOpt<S,O>>( move( clone( *(v.first) ) ) ,
                                              move( clone( *(v.second) ) ) );
  return std::move( ptr );
}

template<class S, class O>
index_t<Vector_SimOpt<S,O>>
dimension( const Vector_SimOpt<S,O>& x ) {
  return dimension( *(x.first) ) + 
         dimension( *(x.second) ); 
}

template<class S, class O>
void set( Vector_SimOpt<S,O>& x, const Vector_SimOpt<S,O>& y ) {
  set( *(x.first),  *(y.first)  );
  set( *(x.second), *(y.second) );
}

template<class S, class O>
void dual( dual_t<Vector_SimOpt<S,O>>&  xdual, 
           const Vector_SimOpt<S,O>& xprim ) {
  set(xdual,xprim);
}

template<class S, class O>
void plus( Vector_SimOpt<S,O>& x, const Vector_SimOpt<S,O>& y ) {
  plus( *(x.first),  *(y.first)  );
  plus( *(x.second), *(y.second) );
}

template<class S, class O>
void scale( Vector_SimOpt<S,O>& x, const element_t<Vector_SimOpt<S,O>> alpha ) {
  scale( *(x.first),  alpha );
  scale( *(x.second), alpha );
}

template<class S, class O>
void fill( Vector_SimOpt<S,O>& x, const element_t<Vector_SimOpt<S,O>> alpha ) {
  fill( *(x.first),  alpha );
  fill( *(x.second), alpha );
}

template<class S, class O>
void scale( Vector_SimOpt<S,O>& x,  
            const element_t<Vector_SimOpt<S,O>> alpha,
            const Vector_SimOpt<S,O>& y ) {
  axpy( *(x.first),  alpha, *(y.first)  );
  axpy( *(x.second), alpha, *(y.second) );
}


template<class S, class O>
void basis( Vector_SimOpt<S,O>& b,
            index_t<Vector_SimOpt<S,O>> i ) {
  auto d = dimension(*(b.first));
  if(i < d ) {
    basis( *(b.first) );   
    fill( *(b.second), 0 );
  }
  else {
    fill( *(b.first), 0 );
    basis( *(b.second), i-d );
  }
}

template<class S, class O>
element_t<Vector_SimOpt<S,O>>
dot( const Vector_SimOpt<S,O>& x,
     const Vector_SimOpt<S,O>& y ) {
  return dot( *(x.first),  *(y.first)  ) + 
         dot( *(x.second), *(y.second) );
}

template<class S, class O>
magnitude_t<Vector_SimOpt<S,O>>
norm( const Vector_SimOpt<S,O>& x ) {
  return std::sqrt( dot( *(x.first),  *(x.first)  ) + 
                    dot( *(x.second), *(x.second) ) );
}

template<class S, class O>
void print( const Vector_SimOpt<S,O>& x, std::ostream& os ) {
  os << "Sim: "; print( *(x.first)  );
  os << "Opt: "; print( *(x.second) );
}

template<class R, class S, class O>
auto reduce( const R& r, const Vector_SimOpt<S,O>& x ) {
  auto rfirst  = reduce( r, *(x.first)  );
  auto rsecond = reduce( r, *(x.second) );
  return r(rfirst,rsecond);
}

template<class Generator, class Distribution, class S, class O>
void randomize( Generator& g, Distribution& d, Vector_SimOpt<S,O>& x ) {
  randomize( g, d, *(x.first)  );
  randomize( g, d, *(x.second) );
}

template<class S, class O, class Function, class... Vs>
void eval_function( Vector_SimOpt<S,O>& x, const Function &f, const Vs&... vs ) {
  eval_function( *(x.first),  f, std::make_tuple( *(vs.first)  ... ) );
  eval_function( *(x.second), f, std::make_tuple( *(vs.second) ... ) );
}


template<class R, class F, class S, class O, class... Vs>
auto eval_function_and_reduce( const R& r, const F& f, 
  const Vector_SimOpt<S,O>& x, const Vs&... vs ) {
  auto rfirst  = eval_function_and_reduce( r, f, *(x.first),  *(vs.first)  ... );
  auto rsecond = eval_function_and_reduce( r, f, *(x.second), *(vs.second) ... );
  return r( rfirst, rsecond );
}

} // namespace XROL
