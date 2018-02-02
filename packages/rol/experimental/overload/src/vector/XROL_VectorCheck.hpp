
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

template<class V> 
auto checkVector( const V& x, const V& y, const V& z, std::ostream &os = std::cout ) {

  using namespace std;

  element_t<V> zero{0.0};
  element_t<V> one{1.0};
  element_t<V> a{1.234};
  element_t<V> b{-0.4321};

  unique_ptr<V> v_ptr = clone(z);

  auto v    = *v_ptr;
  auto vtmp = *clone(z);
  auto xtmp = *clone(x);
  auto ytmp = *clone(y);

  vector<element_t<V>> vCheck;
  //vCheck.resize(dimension(x));

  auto store = [&vCheck,v]() { vCheck.push_back(norm(v)); };
  auto reset = [x,y,z,&xtmp,&ytmp,&v]() { set(v,z);  set(xtmp,x);  set(ytmp,y); };
  auto println = [&os,&vCheck]( auto text ) { 
    to_ostream(os, setw(76), left, text);
    os << ". Consistency error:  " << vCheck.back() << "\n";
  };



  os << "\n" << setw(96) << left << setfill('*') << 
    "********** Begin verification of linear algebra. " << "\n\n";
  
  reset();
  plus(v,x); plus(xtmp,z); axpy(v,-one,xtmp); 
  store();
  os << scientific << setprecision(12) << setfill('.');
  println("Commutativity of addition");

  reset(); 
  plus(ytmp,x); plus(v,ytmp); plus(xtmp,z); plus(xtmp,y); 
  axpy(v,-one,xtmp);
  store();
  println("Associativity of addition");

  reset();
  fill(v,0); plus(v,x); axpy(v, -one, x); 
  store();
  println("Identity element of addition");

  reset();
  scale(v,-one); plus(v,z);
  store();
  println("Inverse elements of addition");

  reset();
  scale(v,one); axpy(v,-one, z);
  store();
  println("Identity element of scalar multiplication");

  reset();
  set(v,z); set(vtmp,z);
  scale(v,b); scale(v,a); scale(vtmp,a*b); axpy(v,-one,vtmp); 
  store();
  println("Consistency of scalar multiplication with field multiplication");

  set(v,z); set(vtmp,z);
  scale(v,a+b); scale(vtmp,a); axpy(vtmp,b,z); axpy(v,-one,vtmp); 
  store();
  println("Distributivity of scalar multiplication with respect to field addition");

  reset();
  plus(v,x); scale(v,a); scale(xtmp,a); axpy(xtmp,a,z); axpy(v,-one,xtmp);
  store();
  println("Distributivity of scalar multiplication with respect to vector addition");

  vCheck.push_back(abs(dot(z,x) - dot(x,z)));
  println("Commutativity of dot (inner) product over the field of reals");

  set(xtmp,x);
  plus(xtmp,y); 
  auto zxt = dot(z,xtmp);
  auto zx  = dot(z,x);
  auto zy  = dot(z,y);
  {
    vCheck.push_back(abs(zxt - zx - zy)/max(abs(zxt), max(abs(zx), abs(zy))));
  }
  println("Additivity of dot (inner) product");

  set(v,z);
  auto vnorm = norm(v);
  if (vnorm == zero) {
    scale(v,a);
    vCheck.push_back(std::abs(norm(v) - zero));
  } else {
    scale(v,one/vnorm);
    vCheck.push_back(std::abs(norm(v) - one));
  }
  println("Consistency of scalar multiplication and norm");

  // Reflexivity.
  set(v,z);
  dual(xtmp,z);
  dual(ytmp,xtmp);
  axpy(v,-one, ytmp); 
  store();
  println("Reflexivity"); 
  os << "\n";
  os << setfill('*');
  // Restore format state of pStream used for the header info.
  os << std::setw(96) << std::left << "********** End verification of linear algebra." << "\n\n";

  return vCheck;

}



} // namespace XROL
