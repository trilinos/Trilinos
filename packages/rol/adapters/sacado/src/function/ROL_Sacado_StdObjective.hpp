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

#ifndef ROL_SACADO_STDOBJECTIVE_HPP
#define ROL_SACADO_STDOBJECTIVE_HPP

#include "Sacado.hpp"
#include "ROL_Objective.hpp"
#include "ROL_StdVector.hpp"

namespace ROL {

/** \brief Generic StdObjective wrapper class for class that uses Sacado */
template<class Real, template<class> class Obj>
class Sacado_StdObjective : public Objective<Real> {

  template <typename T> using vector = std::vector<T>;
  
  typedef Vector<Real>    V;
  typedef StdVector<Real> SV;

protected:

  // Template Object should inherit from ROL::StdObjective
  Obj<Real> obj_;

  /* Evaluate the gradient at x */
  template<class ScalarT> 
  void gradientAD( vector<ScalarT> &g, const vector<ScalarT> &x, Real &tol );

  /* Compute the action of the Hessian evaluated at x on a vector v */
  template<class ScalarT> 
  void hessVecAD( vector<ScalarT> &hv, const vector<ScalarT> &v, 
                  const vector<ScalarT> &x, Real &tol ); 

public:

  /* Evaluate the objective function at x */
  using Objective<Real>::value;
  Real value( const V &x, Real &tol ) {
    const SV xs  = dynamic_cast<const SV&>(x);
    return value(*(xs.getVector()), tol);            
  }

  Real value( const vector<Real> &x, Real &tol ) {
    return obj_.value(x,tol);
  }

  /* Evaluate the gradient at x */
  using Objective<Real>::gradient; 
  void gradient( V &g, const V &x, Real &tol ) {
    SV gs  = dynamic_cast<SV&>(g);
    const SV xs  = dynamic_cast<const SV&>(x);
    this->gradient(*(gs.getVector()),*(xs.getVector()),tol);
  }

  void gradient( vector<Real> &g, const vector<Real> &x, Real &tol ) {
      this->gradientAD(g,x,tol); 
  }

  /* Compute the action of the Hessian evaluated at x on a vector v */
  using Objective<Real>::hessVec;
  void hessVec( V &hv, const V &v, const V &x, Real &tol ) {
    SV hvs = dynamic_cast<SV&>(hv);
    const SV vs = dynamic_cast<const SV&>(v);
    const SV xs = dynamic_cast<const SV&>(x);
    hessVec(*(hvs.getVector()),*(vs.getVector()),*(xs.getVector()),tol);
  }

  void hessVec( vector<Real> &hv, const vector<Real> &v, const vector<Real> &x, Real &tol ) {
      this->hessVecAD(hv,v,x,tol);
  }
};



template<class Real, template<class> class Obj>
template<class ScalarT>
void Sacado_StdObjective<Real,Obj>::gradientAD(vector<ScalarT> &g, const vector<ScalarT> &x, Real &tol) { 

    // Data type which supports automatic differentiation 
    typedef Sacado::Fad::DFad<ScalarT> FadType;
    typedef vector<FadType>            Fadvector;

    // Get a pointer to the gradient vector
    int n = x.size();

    // Create a vector of independent variables
    ROL::Ptr<Fadvector> x_fad = ROL::makePtr<Fadvector>();

    x_fad->reserve(n);   

    // Initialize constructor for each element
    for(int i=0; i<n; ++i) {
        x_fad->push_back(FadType(n,i,x[i])); 
    }

    // AD access to objective function
    FadType J_fad = obj_.value(*x_fad,tol);

    // Evaluate gradient
    for(int i=0; i<n; ++i) {
        g[i] = J_fad.dx(i);
    }
}



template <class Real, template<class> class Obj>
template <class ScalarT>
void Sacado_StdObjective<Real,Obj>::hessVecAD( vector<ScalarT> &hv, const vector<ScalarT> &v, 
                                               const vector<ScalarT> &x, Real &tol ) {

    // Data type which supports automatic differentiation 
    typedef Sacado::Fad::SFad<ScalarT,1> FadType;
    typedef vector<FadType>              Fadvector;

    int n = x.size();
   
    // Create a vector of independent variables
    ROL::Ptr<Fadvector> x_fad = ROL::makePtr<Fadvector>();

    x_fad->reserve(n); 

    // Allocate for gradient   
    Fadvector g_fad(n);

    for(int i=0; i<n; ++i) {
        x_fad->push_back( FadType(1,x[i]) );
    }

    // Set directional derivative    
    for(int i=0; i<n; ++i) {
        (*x_fad)[i].fastAccessDx(0) = v[i];
    }
    
    this->gradientAD(g_fad,*x_fad,tol);

    for(int i=0; i<n; ++i) {
        hv[i] = g_fad[i].dx(0);            
    }

} // class Sacado_StdObjective
 
} // namespace ROL

#endif // ROL_SACADO_STDOBJECTIVE_HPP
