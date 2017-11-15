
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

#include "XROL_Core.hpp"

namespace XROL {

namespace details {

using namespace std;

template<class V>
class Vector : public CRTP<V> {
public:

  using IndexT   = index_t<V>;
  using ElementT = element_t<V>;
  using NormT    = norm_t<V>;
  using DualT    = dual_t<V>;
  
  void plus( const Vector& x )                       { this->impl().plus(x.impl());        }
  void set( const Vector& x )                        { this->impl().set(x.impl());         }
  NormT dot( const Vector& x ) const                 { return this->impl().dot(x.impl());  }
  NormT norm() const                                 { return this->impl().norm();         }
  unique_ptr<Vector> clone() const                   { return this->impl().clone();        }
  void axpy( const ElementT alpha, const Vector& x ) { this->impl().axpy(alpha,x.impl());  }
  void fill( const ElementT alpha )                  { this->impl().fill(alpha);           }
  void scale( const ElementT alpha )                 { this->impl().scale(alpha);          }
  unique_ptr<Vector> basis( IndexT i ) const         { return this->impl().basis(i);       }
  IndexT dimension() const                           { return this->impl().dimension();    }
  void dual(DualT& x) const                          { return this->impl().dual(x.dual()); }

  void print( ostream& os, const string& delimiter=" " ) const { 
    this->impl().print(os);            
  } 

  // Elementwise functions

  // y = f(x1,x2,...)
  template<class F, class... Vs>
  void applyFunction( F&& f, Vs&&... vs ) {
    this->impl().applyFunction( forward<F>(f), forward<Vs>(vs)... );
  }
  
  // result = r(result,y_i) for all i
  template<class R>
  NormT reduce( R&& r ) const {
    this->impl().reduce(forward<R>(r));
  }

  // result = r(result, f(y_i,x1_i,x2_i,...) for all i
  template<class F, class R, class... Vs>
  NormT applyFunctionAndReduce( F&& f, R&& r, Vs&&... vs ) const {
    return this->impl().applyFunctionAndReduce( forward<F>(f), forward<R>(r), forward<Vs>(vs)... );  
  }

   
 vector<NormT> checkVector( const Vector &x,
                            const Vector &y,
                            const bool printToStream = true,
                            ostream & outStream = cout ) const {
    NormT zero =  0.0;
    NormT one  =  1.0;
    NormT a    =  1.234;
    NormT b    = -0.4321;
    int width =  94;
    vector<NormT> vCheck;

    Teuchos::oblackholestream bhs; // outputs nothing

    Teuchos::RCP<ostream> pStream;
    if (printToStream) {
      pStream = Teuchos::rcp(&outStream, false);
    } else {
      pStream = Teuchos::rcp(&bhs, false);
    }

    // Save the format state of the original pStream.
    Teuchos::oblackholestream oldFormatState, headerFormatState;
    oldFormatState.copyfmt(*pStream);

    Teuchos::RCP<Vector> v    = this->clone();
    Teuchos::RCP<Vector> vtmp = this->clone();
    Teuchos::RCP<Vector> xtmp = x.clone();
    Teuchos::RCP<Vector> ytmp = y.clone();

    //*pStream << "\n************ Begin verification of linear algebra.\n\n";
    *pStream << "\n" << setw(width) << left << setfill('*') << "********** Begin verification of linear algebra. " << "\n\n";
    headerFormatState.copyfmt(*pStream);

    // Commutativity of addition.
    v->set(*this); xtmp->set(x); ytmp->set(y);
    v->plus(x); xtmp->plus(*this); v->axpy(-one, *xtmp); vCheck.push_back(v->norm());
    *pStream << scientific << setprecision(12) << setfill('>');
    *pStream << setw(width) << left << "Commutativity of addition. Consistency error: " << " " << vCheck.back() << "\n";

    // Associativity of addition.
    v->set(*this); xtmp->set(x); ytmp->set(y);
    ytmp->plus(x); v->plus(*ytmp); xtmp->plus(*this); xtmp->plus(y); v->axpy(-one, *xtmp); vCheck.push_back(v->norm());
    *pStream << setw(width) << left << "Associativity of addition. Consistency error: " << " " << vCheck.back() << "\n";

    // Identity element of addition.
    v->set(*this); xtmp->set(x); ytmp->set(y);
    v->zero(); v->plus(x); v->axpy(-one, x); vCheck.push_back(v->norm());
    *pStream << setw(width) << left << "Identity element of addition. Consistency error: " << " " << vCheck.back() << "\n";

    // Inverse elements of addition.
    v->set(*this); xtmp->set(x); ytmp->set(y);
    v->scale(-one); v->plus(*this); vCheck.push_back(v->norm());
    *pStream << setw(width) << left << "Inverse elements of addition. Consistency error: " << " " << vCheck.back() << "\n";

    // Identity element of scalar multiplication.
    v->set(*this); xtmp->set(x); ytmp->set(y);
    v->scale(one); v->axpy(-one, *this); vCheck.push_back(v->norm());
    *pStream << setw(width) << left << "Identity element of scalar multiplication. Consistency error: " << " " << vCheck.back() << "\n";

    // Consistency of scalar multiplication with field multiplication.
    v->set(*this); vtmp->set(*this);
    v->scale(b); v->scale(a); vtmp->scale(a*b); v->axpy(-one, *vtmp); vCheck.push_back(v->norm());
    *pStream << setw(width) << left << "Consistency of scalar multiplication with field multiplication. Consistency error: " << " " << vCheck.back() << "\n";

    // Distributivity of scalar multiplication with respect to field addition.
    v->set(*this); vtmp->set(*this);
    v->scale(a+b); vtmp->scale(a); vtmp->axpy(b, *this); v->axpy(-one, *vtmp); vCheck.push_back(v->norm());
    *pStream << setw(width) << left << "Distributivity of scalar multiplication with respect to field addition. Consistency error: " << " " << vCheck.back() << "\n";

    // Distributivity of scalar multiplication with respect to vector addition.
    v->set(*this); xtmp->set(x); ytmp->set(y);
    v->plus(x); v->scale(a); xtmp->scale(a); xtmp->axpy(a, *this); v->axpy(-one, *xtmp); vCheck.push_back(v->norm());
    *pStream << setw(width) << left << "Distributivity of scalar multiplication with respect to vector addition. Consistency error: " << " " << vCheck.back() << "\n";

    // Commutativity of dot (inner) product over the field of reals.
    vCheck.push_back(abs(this->dot(x) - x.dot(*this)));
    *pStream << setw(width) << left << "Commutativity of dot (inner) product over the field of reals. Consistency error: " << " " << vCheck.back() << "\n";

    // Additivity of dot (inner) product.
    xtmp->set(x);
    xtmp->plus(y); vCheck.push_back(abs(this->dot(*xtmp) - this->dot(x) - this->dot(y))/max(abs(this->dot(*xtmp)), max(abs(this->dot(x)), abs(this->dot(y)))));
    *pStream << setw(width) << left << "Additivity of dot (inner) product. Consistency error: " << " " << vCheck.back() << "\n";

    // Consistency of scalar multiplication and norm.
    v->set(*this);
    NormT vnorm = v->norm();
    if (vnorm == zero) {
      v->scale(a);
      vCheck.push_back(abs(v->norm() - zero));
    } else {
      v->scale(one/vnorm);
      vCheck.push_back(abs(v->norm() - one));
    }
    *pStream << setw(width) << left << "Consistency of scalar multiplication and norm. Consistency error: " << " " << vCheck.back() << "\n";

    // Reflexivity.
    v->set(*this);
    xtmp = Teuchos::rcp_const_cast<Vector>(Teuchos::rcpFromRef(this->dual()));
    ytmp = Teuchos::rcp_const_cast<Vector>(Teuchos::rcpFromRef(xtmp->dual()));
    v->axpy(-one, *ytmp); vCheck.push_back(v->norm());
    *pStream << setw(width) << left << "Reflexivity. Consistency error: " << " " << vCheck.back() << "\n\n";

    //*pStream << "************   End verification of linear algebra.\n\n";

    // Restore format state of pStream used for the header info.
    pStream->copyfmt(headerFormatState);
    *pStream << setw(width) << left << "********** End verification of linear algebra. " << "\n\n";
    // Restore format state of the original pStream.
    pStream->copyfmt(oldFormatState);

    return vCheck;
  }

}; // class Vector

} // namespace details

template<class V> using Vector = details::Vector<V>;

//template<class V> 
//struct VectorFactory {
//  static unique_ptr<V> create() {}
//}; 


} // namespace XROL
