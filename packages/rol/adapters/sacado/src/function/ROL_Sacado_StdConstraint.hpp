// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_SACADO_STDEQUALITYCONSTRAINT_HPP
#define ROL_SACADO_STDEQUALITYCONSTRAINT_HPP

#include "Sacado.hpp"
#include "ROL_StdConstraint.hpp"

namespace ROL {

//! \brief ROL interface wrapper for Sacado Constraint
template<class Real, template<class> class Constr>
class Sacado_StdConstraint : public Constraint<Real> {

  template <typename T> using vector = std::vector<T>;

  typedef Vector<Real>     V;
  typedef StdVector<Real>  SV; 

protected:

  // Template Object should inherit from ROL::StdConstraint
  Constr<Real> constr_;

  template<class ScalarT>
  void applyJacobianAD( vector<ScalarT> &jv, const vector<ScalarT> &v,
                        const vector<ScalarT> &x, Real &tol );   
 
  template<class ScalarT>
  void applyAdjointJacobianAD( vector<ScalarT> &aju, const vector<ScalarT> &u,
                               const vector<ScalarT> &x, Real &tol);

  template<class ScalarT>
  void applyAdjointHessianAD( vector<ScalarT> &ahuv, const vector<ScalarT> &u,
                              const vector<ScalarT> &v, const vector<ScalarT> &x, Real &tol);


public: 
 
  using Constraint<Real>::value;
  void value(V &c, const V &x, Real &tol ) {
    SV cs = dynamic_cast<SV&>(c);
    const SV xs = dynamic_cast<const SV&>(x);
    value(*(cs.getVector()),*(xs.getVector()),tol); 
  }

  void value(vector<Real> &c, const vector<Real>  &x, Real &tol) {
    constr_.value(c,x,tol);
  }
  
  using Constraint<Real>::applyJacobian;
  void applyJacobian( V &jv, const V &v, const V &x, Real &tol ) {
    SV jvs = dynamic_cast<SV&>(jv);
    const SV vs = dynamic_cast<const SV&>(v);
    const SV xs = dynamic_cast<const SV&>(x);
    applyJacobian(*(jvs.getVector()), *(vs.getVector()), *(xs.getVector()), tol);   
  }  
  
  void applyJacobian(vector<Real>  &jv, const vector<Real>  &v, 
                             const vector<Real>  &x, Real &tol) {
    this->applyJacobianAD(jv,v,x,tol);
  }

  using Constraint<Real>::applyAdjointJacobian;
  void applyAdjointJacobian( V &aju, const V &u, const V &x, Real &tol ) {
    SV ajus = dynamic_cast<SV&>(aju);
    const SV us = dynamic_cast<const SV&>(u);
    const SV xs = dynamic_cast<const SV&>(x);
    applyAdjointJacobian(*(ajus.getVector()),*(us.getVector()),*(xs.getVector()),tol);  
  }

  void applyAdjointJacobian(vector<Real>  &aju, const vector<Real>  &u,
                                    const vector<Real>  &x, Real &tol) {
    this->applyAdjointJacobianAD(aju,u,x,tol);
  } 


  using Constraint<Real>::applyAdjointHessian;
  void applyAdjointHessian( V &ahuv, const V &u, const V& v, const V &x, Real &tol ) {
    SV ahuvs = dynamic_cast<SV&>(ahuv);
    const SV us = dynamic_cast<const SV&>(u);
    const SV vs = dynamic_cast<const SV&>(v);
    const SV xs = dynamic_cast<const SV&>(x);
    applyAdjointHessian(*(ahuvs.getVector()),*(us.getVector()),
                        *(vs.getVector()),*(xs.getVector()),tol);  
  }

  void applyAdjointHessian(vector<Real>  &ahuv, const vector<Real>  &u,
                           const vector<Real>  &v, const vector<Real>  &x, Real &tol){
    this->applyAdjointHessianAD(ahuv,u,v,x,tol);
  }

}; // class Sacado_StdConstraint


template<class Real, template<class> class Constr>
template<class ScalarT>
void Sacado_StdConstraint<Real,Constr>::applyJacobianAD(vector<ScalarT> &jv, const vector<ScalarT> &v, 
                                                                const vector<ScalarT> &x, Real &tol) {

    // Data type which supports automatic differentiation 
    typedef Sacado::Fad::DFad<ScalarT> FadType;
    typedef vector<FadType>       Fadvector;

    int n = x.size();
    int m = jv.size();

    // Create a vector of independent variables
    ROL::Ptr<Fadvector> x_fad = ROL::makePtr<Fadvector>();

    x_fad->reserve(n);

    // Initialize constructor for each element
    for(int i=0; i<n; ++i) {
        x_fad->push_back(FadType(n,i,x[i])); 
    }

    // Create a vector of independent variables
    Fadvector c_fad(m);

    constr_.value(c_fad,*x_fad,tol);

    for(int i=0; i<m; ++i) {
        jv[i] = 0;
        for(int j=0; j<n; ++j) {
            jv[i] += v[j]*c_fad[i].dx(j); 
        }   
    }       
}



template<class Real, template<class> class Constr>
template<class ScalarT>
void Sacado_StdConstraint<Real,Constr>::applyAdjointJacobianAD(vector<ScalarT> &aju, const vector<ScalarT> &u, 
                                                                    const vector<ScalarT> &x, Real &tol) {

   // Data type which supports automatic differentiation 
    typedef Sacado::Fad::DFad<ScalarT> FadType;
    typedef vector<FadType>            Fadvector;
  
    int n = x.size();
    int m = u.size();

    // Create a vector of independent variables
    ROL::Ptr<Fadvector> x_fad = ROL::makePtr<Fadvector>();

    x_fad->reserve(n);

    // Initialize constructor for each element
    for(int i=0; i<n; ++i) {
        x_fad->push_back(FadType(n,i,x[i])); 
    }

    Fadvector c_fad(m);
 
    // Evaluate constraint
    constr_.value(c_fad,*x_fad,tol);
    
    FadType udotc = 0;
    
    for(int j=0;j<m;++j){ 
        udotc += c_fad[j]*u[j];
    } 

    // Take gradient of the dot product
    for(int i=0;i<n;++i){
        aju[i] = udotc.dx(i);
    } 
}


template<class Real, template<class> class Constr>
template<class ScalarT>
void Sacado_StdConstraint<Real,Constr>::applyAdjointHessianAD(vector<ScalarT> &ahuv, const vector<ScalarT> &u,
                                                                   const vector<ScalarT> &v, const vector<ScalarT> &x, 
                                                                   Real &tol){
  // Data type which supports automatic differentiation 
  typedef Sacado::Fad::SFad<ScalarT,1> FadType;
  typedef vector<FadType>              Fadvector;

  int m = u.size();

  // Number of optimization variables
  int n = x.size();
  
  // Create a vector of independent variables
  ROL::Ptr<Fadvector> x_fad = ROL::makePtr<Fadvector>();

  x_fad->reserve(n);

  // Allocate for directional adjoint Jacobian
  Fadvector aju_fad(n);

  for(int i=0; i<n; ++i) {
    x_fad->push_back( FadType(1,x[i]) );

    // Set derivative direction
    (*x_fad)[i].fastAccessDx(0) = v[i];     
  }

  // Allocate for constraint vector
  Fadvector c_fad(m);

  // Allocate for dual constraint vector
  Fadvector u_fad(m);

  for(int j=0; j<m; ++j) {
    u_fad[j] = u[j];
  }

  // Evaluate constraint adjoint Jacobian direction
  this->applyAdjointJacobianAD( aju_fad, u_fad, *x_fad, tol);

  for(int i=0; i<n; ++i) {
    ahuv[i] = aju_fad[i].dx(0);             
  }
}

} // namespace ROL

#endif // ROL_SACADO_STDEQUALITYCONSTRAINT_HPP
