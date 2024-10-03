// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   example_01.cpp
    \brief  Minimize the Gross-Pitaevskii functional and demonstrate 
            the effect of choice of function space of the Gradient on
            convergence.

    \details Minimize the one-dimensional Gross-Pitaevskii (GP) energy 
             functional
             \f[ J[\psi] = \int \frac{1}{2} |\nabla\psi|^2 + V(x)|\psi|^2 
                           +g|\psi|^4 \,\mathrm{d}x \f]
             Subject to the equality constraint that the particle density be
             normalized. 
             \f[ e(\psi) = \int |\psi|^2\,\mathrm{d}x - 1 = 0 \f]
             For simplicity, we will assume the wavefunction \f$\psi\f$ to 
             be real-valued, the potential function \f$ V(x)\geq 0\f$,
             the computational domain is the interval \f$[0,1]\f$, and that
             \f$\psi(0)=\psi(1)=0\f$. We also discretize the problem using
             second-order centered finite differences on a uniform grid. 

             \f[
             \psi''(x_i) \approx = \frac{\psi(x_{i-1})-2\psi(x_i)+\psi(x_{i+1})}{\Delta x^2}
             \f]

    \author Greg von Winckel
    \date   Mon Dec  1 12:55:12 MST 2014
*/

#include <iostream>

#include "ROL_Stream.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include "ROL_ParameterList.hpp"
#include "ROL_StdVector.hpp"
#include "ROL_Objective.hpp"
#include "ROL_Constraint.hpp"
#include "ROL_Algorithm.hpp"
#include "ROL_CompositeStep.hpp"
#include "ROL_ConstraintStatusTest.hpp"


using namespace ROL;

template<class Real>
class Objective_GrossPitaevskii : public Objective<Real> {

  typedef std::vector<Real>  vector;
  typedef Vector<Real>       V;
  typedef StdVector<Real>    SV;

  typedef typename vector::size_type uint;

private:


  /** \var Real g_ appearing before quartic term in GP functional    */ 
  Real g_;    

  /*! \var int nx_ Number of interior nodes  */ 
  uint  nx_;     

  /*! \var int nx_ Mesh spacing \f$ \Delta x = \frac{1}{n_x+1} \f$  */ 
  Real dx_;     
        
  /*! \var ptr Vp_ Pointer to potential vector  */ 
  ROL::Ptr<const vector> Vp_;    

  ROL::Ptr<const vector> getVector( const V& x ) {
    
    return dynamic_cast<const SV&>(x).getVector();
  }

  ROL::Ptr<vector> getVector( V& x ) {
    
    return dynamic_cast<SV&>(x).getVector();
  }



  //! Apply finite difference operator 
  /*! Compute \f$K\psi\f$, where \f$K\f$ is the finite difference approximation 
      of \f$-D_x^2\f$ */
  void applyK(const Vector<Real> &v, Vector<Real> &kv) {

    using namespace Teuchos;

    // Pointer to direction vector 
    ROL::Ptr<const vector> vp = getVector(v);

    // Pointer to action of Hessian on direction vector 
    ROL::Ptr<vector> kvp = getVector(kv);

    Real dx2 = dx_*dx_;

    (*kvp)[0] = (2.0*(*vp)[0]-(*vp)[1])/dx2;
  
    for(uint i=1;i<nx_-1;++i) {
      (*kvp)[i] = (2.0*(*vp)[i]-(*vp)[i-1]-(*vp)[i+1])/dx2;
    } 

    (*kvp)[nx_-1] = (2.0*(*vp)[nx_-1]-(*vp)[nx_-2])/dx2;
  } 

  public: 

  Objective_GrossPitaevskii(const Real &g, const Vector<Real> &V) : g_(g),  
    Vp_(getVector(V))  {
    nx_ = Vp_->size(); 
    dx_ = (1.0/(1.0+nx_));
  }
           
  //! Evaluate \f$J[\psi]\f$
  /*! \f[ J[\psi]=\frac{1}{2} \int\limits_0^1 |\psi'|^2 + 
          V(x)|\psi|^2+g|\psi|^4\,\mathrm{d}x \f] 
          where the integral is approximated with the trapezoidal rule and
          the derivative is approximated using finite differences */
  Real value( const Vector<Real> &psi, Real &tol ) {

        
    
    // Pointer to opt vector 
    ROL::Ptr<const vector> psip = getVector(psi);

    // Pointer to K applied to opt vector 
    ROL::Ptr<V> kpsi = psi.clone();
    ROL::Ptr<vector> kpsip = getVector(*kpsi);

    Real J = 0;

    this->applyK(psi,*kpsi);

    for(uint i=0;i<nx_;++i) {
      J += (*psip)[i]*(*kpsip)[i] + (*Vp_)[i]*pow((*psip)[i],2) + g_*pow((*psip)[i],4);
    } 
      
    // Rescale for numerical integration by trapezoidal rule
    J *= 0.5*dx_;

    return J;
  }

  //! Evaluate \f$\nabla J[\psi]\f$
  /*! \f[ \nabla J[\psi] = -\psi'' + V(x)\psi+2g|\psi|^3 \f] */
  void gradient( Vector<Real> &g, const Vector<Real> &psi, Real &tol ) {

    

    // Pointer to opt vector 
    ROL::Ptr<const vector> psip = getVector(psi);

    // Pointer to gradient vector 
    ROL::Ptr<vector> gp = getVector(g);

    // Pointer to K applied to opt vector 
    ROL::Ptr<V> kpsi = psi.clone();
    ROL::Ptr<vector> kpsip = getVector(*kpsi);

    applyK(psi,*kpsi);

    for(uint i=0;i<nx_;++i) {
      (*gp)[i] = ((*kpsip)[i] + (*Vp_)[i]*(*psip)[i] + 2.0*g_*pow((*psip)[i],3))*dx_;
    } 
      
  }



  //! Evaluate \f$\nabla^2 J[\psi] v\f$
  /*! \f[ \nabla^2 J[\psi]v = -v'' + V(x)v+6g|\psi|^2 v \f] */
  void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &psi, Real &tol ) {

    

    // Pointer to opt vector 
    ROL::Ptr<const vector> psip = getVector(psi);

    // Pointer to direction vector 
    ROL::Ptr<const vector> vp = getVector(v);

    // Pointer to action of Hessian on direction vector 
    ROL::Ptr<vector> hvp = getVector(hv);

    applyK(v,hv);
 
    for(uint i=0;i<nx_;++i) {
      (*hvp)[i] *= dx_;
      (*hvp)[i] += ( (*Vp_)[i] + 6.0*g_*pow((*psip)[i],2) )*(*vp)[i]*dx_;
    } 
  }

};


template<class Real>
class Normalization_Constraint : public Constraint<Real> {

  typedef std::vector<Real>  vector;
  typedef Vector<Real>       V;
  typedef StdVector<Real>    SV;

  typedef typename vector::size_type uint;


private:     
  uint nx_;
  Real dx_;

  ROL::Ptr<const vector> getVector( const V& x ) {
    
    return dynamic_cast<const SV&>(x).getVector();
  }

  ROL::Ptr<vector> getVector( V& x ) {
    
    return dynamic_cast<SV&>(x).getVector();
  }

public:
  Normalization_Constraint(int n, Real dx) : nx_(n), dx_(dx) {}          

  //! Evaluate \f$c[\psi]\f$
  /*! \f[ c[\psi]= \int\limits_0^1 |\psi|^2\,\mathrm{d}x - 1 \f] 
      where the integral is approximated with the trapezoidal rule and
      the derivative is approximated using finite differences. 
      This constraint is a scalar */
  void value(Vector<Real> &c, const Vector<Real> &psi, Real &tol){

    

    // Pointer to constraint vector (only one element)
    ROL::Ptr<vector> cp = getVector(c);

    // Pointer to optimization vector     
    ROL::Ptr<const vector> psip = getVector(psi);

    (*cp)[0] = -1.0;
    for(uint i=0;i<nx_;++i) {
      (*cp)[0] += dx_*pow((*psip)[i],2);
    } 
  }

  //! Evaluate \f$c'[\psi]v\f$
  /*! \f[ c'[\psi]v= 2 \int\limits_0^1 \psi v\,\mathrm{d}x  \f]
      The action of the Jacobian on a vector produces a scalar */
  void applyJacobian(Vector<Real> &jv, const Vector<Real> &v, const Vector<Real> &psi, Real &tol){

    

    // Pointer to action of Jacobian of constraint on direction vector (yields scalar)
    ROL::Ptr<vector> jvp = getVector(jv);

    // Pointer to direction vector     
    ROL::Ptr<const vector> vp = getVector(v);

    // Pointer to optimization vector     
    ROL::Ptr<const vector> psip = getVector(psi);
      
        (*jvp)[0] = 0;
        for(uint i=0;i<nx_;++i) {
            (*jvp)[0] += 2.0*dx_*(*psip)[i]*(*vp)[i];
        }
    }

  //! Evaluate \f$(c'[\psi])^\ast v\f$
  /*! \f[ (c'[\psi])^\ast v = 2 \int\limits_0^1 \psi v\,\mathrm{d}x  \f] 
       The action of the Jacobian adjoint on a scalar produces a vector */
  void applyAdjointJacobian(Vector<Real> &ajv, const Vector<Real> &v, const Vector<Real> &psi, Real &tol){

     

    // Pointer to action of adjoint of Jacobian of constraint on direction vector (yields vector)
    ROL::Ptr<vector> ajvp = getVector(ajv);

    // Pointer to direction vector     
    ROL::Ptr<const vector> vp = getVector(v);

    // Pointer to optimization vector     
    ROL::Ptr<const vector> psip = getVector(psi);

    for(uint i=0;i<nx_;++i) {
      (*ajvp)[i] = 2.0*dx_*(*psip)[i]*(*vp)[0];
    }
  }

  //! Evaluate \f$((c''[\psi])^\ast v)u\f$
  /*! \f[ ((c''[\psi])^\ast v)u = 2 v u   \f] 
      The action of the Hessian adjoint on a on a vector v in a direction u produces a vector of
      the same size as \f$\psi\f$ */
  void applyAdjointHessian(Vector<Real> &ahuv, const Vector<Real> &u, const Vector<Real> &v, 
                           const Vector<Real> &psi, Real &tol){

    

    // The pointer to action of constraint Hessian in u,v inner product
    ROL::Ptr<vector> ahuvp = getVector(ahuv);

    // Pointer to direction vector u     
    ROL::Ptr<const vector> up = getVector(u);

    // Pointer to direction vector v     
    ROL::Ptr<const vector> vp = getVector(v);

    // Pointer to optimization vector     
    ROL::Ptr<const vector> psip = getVector(psi);

    for(uint i=0;i<nx_;++i) {
      (*ahuvp)[i] = 2.0*dx_*(*vp)[i]*(*up)[0];        
    }  
  }
};


