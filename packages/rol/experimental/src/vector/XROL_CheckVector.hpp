
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

#include "XROL_StdVectorTraits.hpp"

namespace XROL {

template<class V> 
auto checkVector( const V& x, const V& y, const V& z, std::ostream &os = std::cout ) {

  using VT = VectorTraits<V>;

  using IndexType     = typename ElementTraits<V>::IndexType;
  using ElementType   = typename ElementTraits<V>::ElementType;
  using MagnitudeType = typename ElementTraits<V>::MagnitudeType;

  ElementType zero{0.0};
  ElementType one{1.0};
  ElementType a{1.234};
  ElementType b{-0.4321};

  int width =  94;

  std::vector<ElementType> vCheck;

  auto v    = *VT::clone(z);
  auto vtmp = *VT::clone(z);
  auto xtmp = *VT::clone(x);
  auto ytmp = *VT::clone(y);

  os << "\n" << std::setw(width) << std::left << std::setfill('*') << "********** Begin verification of linear algebra. " << "\n\n";

  // Commutativity of addition.
  VT::set(v,z);  VT::set(xtmp,x);  VT::set(ytmp,y);
  VT::plus(v,x); VT::plus(xtmp,z); VT::axpy(v,-one,xtmp); vCheck.push_back(VT::norm(v));
  *pStream << std::scientific << std::setprecision(12) << std::setfill('>');
  *pStream << std::setw(width) << std::left << "Commutativity of addition. Consistency error: " << " " << vCheck.back() << "\n";

  // Associativity of addition.
  VT::set(v,z);     VT::set(xtmp,x);  VT::set(ytmp,y);
  VT::plus(ytmp,x); VT::plus(v,ytmp); VT::plus(xtmp,z); VT::plus(xtmp,y); 
  *pStream << std::setw(width) << std::left << "Associativity of addition. Consistency error: " << " " << vCheck.back() << "\n";
/*
  // Identity element of addition.
  v->set(*this); xtmp->set(x); ytmp->set(y);
  v->zero(); v->plus(x); v->axpy(-one, x); vCheck.push_back(v->norm());
  *pStream << std::setw(width) << std::left << "Identity element of addition. Consistency error: " << " " << vCheck.back() << "\n";

  // Inverse elements of addition.
  v->set(*this); xtmp->set(x); ytmp->set(y);
  v->scale(-one); v->plus(*this); vCheck.push_back(v->norm());
  *pStream << std::setw(width) << std::left << "Inverse elements of addition. Consistency error: " << " " << vCheck.back() << "\n";

  // Identity element of scalar multiplication.
  v->set(*this); xtmp->set(x); ytmp->set(y);
  v->scale(one); v->axpy(-one, *this); vCheck.push_back(v->norm());
  *pStream << std::setw(width) << std::left << "Identity element of scalar multiplication. Consistency error: " << " " << vCheck.back() << "\n";

  // Consistency of scalar multiplication with field multiplication.
  v->set(*this); vtmp->set(*this);
  v->scale(b); v->scale(a); vtmp->scale(a*b); v->axpy(-one, *vtmp); vCheck.push_back(v->norm());
  *pStream << std::setw(width) << std::left << "Consistency of scalar multiplication with field multiplication. Consistency error: " << " " << vCheck.back() << "\n";

  // Distributivity of scalar multiplication with respect to field addition.
  v->set(*this); vtmp->set(*this);
  v->scale(a+b); vtmp->scale(a); vtmp->axpy(b, *this); v->axpy(-one, *vtmp); vCheck.push_back(v->norm());
  *pStream << std::setw(width) << std::left << "Distributivity of scalar multiplication with respect to field addition. Consistency error: " << " " << vCheck.back() << "\n";

  // Distributivity of scalar multiplication with respect to vector addition.
  v->set(*this); xtmp->set(x); ytmp->set(y);
  v->plus(x); v->scale(a); xtmp->scale(a); xtmp->axpy(a, *this); v->axpy(-one, *xtmp); vCheck.push_back(v->norm());
  *pStream << std::setw(width) << std::left << "Distributivity of scalar multiplication with respect to vector addition. Consistency error: " << " " << vCheck.back() << "\n";

  // Commutativity of dot (inner) product over the field of reals.
  vCheck.push_back(std::abs(this->dot(x) - x.dot(*this)));
  *pStream << std::setw(width) << std::left << "Commutativity of dot (inner) product over the field of reals. Consistency error: " << " " << vCheck.back() << "\n";

  // Additivity of dot (inner) product.
  xtmp->set(x);
  xtmp->plus(y); vCheck.push_back(std::abs(this->dot(*xtmp) - this->dot(x) - this->dot(y))/std::max(std::abs(this->dot(*xtmp)), std::max(std::abs(this->dot(x)), std::abs(this->dot(y)))));
  *pStream << std::setw(width) << std::left << "Additivity of dot (inner) product. Consistency error: " << " " << vCheck.back() << "\n";

  // Consistency of scalar multiplication and norm.
  v->set(*this);
  Real vnorm = v->norm();
  if (vnorm == zero) {
    v->scale(a);
    vCheck.push_back(std::abs(v->norm() - zero));
  } else {
    v->scale(one/vnorm);
    vCheck.push_back(std::abs(v->norm() - one));
  }
  *pStream << std::setw(width) << std::left << "Consistency of scalar multiplication and norm. Consistency error: " << " " << vCheck.back() << "\n";

  // Reflexivity.
  v->set(*this);
  xtmp = Teuchos::rcp_const_cast<Vector>(Teuchos::rcpFromRef(this->dual()));
  ytmp = Teuchos::rcp_const_cast<Vector>(Teuchos::rcpFromRef(xtmp->dual()));
  v->axpy(-one, *ytmp); vCheck.push_back(v->norm());
  *pStream << std::setw(width) << std::left << "Reflexivity. Consistency error: " << " " << vCheck.back() << "\n\n";


  // Restore format state of pStream used for the header info.
  pStream->copyfmt(headerFormatState);
  *pStream << std::setw(width) << std::left << "********** End verification of linear algebra. " << "\n\n";

  // Restore format state of the original pStream.
  pStream->copyfmt(oldFormatState);
*/
  return vCheck;

  

}



} // namespace XROL
