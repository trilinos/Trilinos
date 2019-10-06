#pragma once
#ifndef RAYLEIGHQUOTIENT_HPP
#define RAYLEIGHQUOTIENT_HPP

#include "ROL_ComplexStdVector.hpp"
#include "HermitianMatrix.hpp"
#include "ROL_Objective.hpp"

/** \class RayleighQuotient
    \brief Objective function for the Rayleigh quotient of a 
           Hermitian Matrix
    
    The Rayleigh quotient is defined for a Hermitian matrix \f$A\f$ as

    \f[ \rho(z) = \frac{ z^\ast A z } { z^\ast z } \f]

    where \f$z\in\mathbb{C}^n\f$, \f$A\in\mathbb{C}^{n\times n}\f$, 
    \f$ A^\ast = A \f$, and \f$ z^\ast = \bar z^\top \f$

    Using the decomposition \f$ z = x + i y \f$, where \f$x,y\in\mathbb{R}^n\f$
    and the definition of the Wirtinger gradients

    \f[ \nabla_z        = \frac{1}{2}\left( \nabla_x - i \nabla_y \right),\quad
        \nabla_{\bar z} = \frac{1}{2}\left( \nabla_x + i \nabla_y \right) \f]
   
    Let \f$ \langle u, v\lange = 2\Re[u^\ast v] \f$
 
    For convenience, introduce the mappings $p,q:\mathbb{C}^n\rightarrow\mathbb{R}$ as 

    \f[ p = \langle z, A z\rangle \f$ and \f$ q = \langle z, z\rangle \f] 

    So that we may express the Rayleigh quotient simply as 

    \f[ \rho(z) = \frac{p}{q} \f]

    The gradient of the Rayleigh quotient with respect to \f$z\f$ is

    \f[ \nabla_z \rho = \frac{ q A z^\ast - p z^\ast }{q^2} \f]
     
    This is understood to be a covector or one-form. We can also write a 
    gradient with respect to \f$ \bar z \f$ as 

    \f[ \nabla_{\bar z} \rho = \frac{ q A z^\top - p z^\top }{q^2} \f]
    
    Directional derivatives can be evaluated given a direction vector \f$ v \in \mathbb{C}^n \f$ 
    as follows:

    \f[ \rho_v = (\nabla_z \rho)^\ast v + (\nabla_{\bar z} \rho)^\ast \bar v 
               = 2\Re[(\nabla_z \rho)^\ast v] \f] 

    The Hessian can be thought of as a \f$2\times 2\f$ block system 

    \f[ \nabla^2 \rho = 
        \begin{pmatrix}
        \nabla_z (\nabla_z \rho)^\ast & \nabla_{\bar z} (\nabla_z \rho)^\ast \\
        \nabla_{\bar z} (\nabla_z \rho)^\ast & \nabla_{\bar z} (\nabla_z \rho)^\ast \\
        \end{pmatrix} 
    \f]
    
     Hessian-vector products can be expressed through block multiplication

     \f[ Hv = [\nabla^2 \rho]v = 
         \begin{pmatrix}
         \nabla_z (\nabla_z \rho)^\ast & \nabla_{\bar z} (\nabla_z \rho)^\ast \\
         \nabla_{\bar z} (\nabla_z \rho)^\ast & \nabla_{\bar z} (\nabla_z \rho)^\ast \\
         \end{pmatrix} \begin{pmatrix} v \\  \bar v \end{pmatrix} 
     \f]

     This can be expressed as

     \f[ Hv = \frac{(4p\alpha-2q\beta)z + q^2 Av - qp v - 2\alpha q y }{q^3} \f]
     
     where \f$ \alpha = \langle z, v\rangle \f$, and \f$ \beta = \langle z, A v \rangle \f$ 
      
 */

template<typename Real>
class RayleighQuotient : public ROL::Objective<Real> {
public:

  using V = ROL::StdVector<Real,std::complex<Real>>;

  RayleighQuotient( HermitianMatrix<Real> A_in ) : 
    A_(A_in), y_(A_in.size()) {}

  Real value( const ROL::Vector<Real>& z, Real& tol ) override {

    A_.apply(y_,z,tol);   
    auto p = z.dot(y_);
    auto q = z.dot(z);
    return p/q;    
  }

  void gradient(       ROL::Vector<Real>& g, 
                 const ROL::Vector<Real>& z, 
                 Real& tol ) override {
    
    A_.apply(g,z,tol);      // y = Az
    auto p = z.dot(g);
    auto q = z.dot(z);
    g.axpy(-p/q,z);         // g = Az - rho * z
    g.scale(2/q);           // g = (Az - rho * z)/q
  }

  void hessVec(       ROL::Vector<Real>& Hv, 
                const ROL::Vector<Real>& v, 
                const ROL::Vector<Real>& z, 
                Real& tol ) override {

    A_.apply(y_,z,tol);     // y = Az
    A_.apply(Hv,v,tol);

    auto p     = z.dot(y_);    
    auto q     = z.dot(z);  
    auto alpha = z.dot(v);  
    auto beta  = z.dot(Hv);  // beta  = Re[(Az,v)]
    auto rho   = p/q;           

    Hv.scale(q);
    Hv.axpy(-2*alpha,y_);
    Hv.axpy(-p,v);
    Hv.axpy(4*rho*alpha-2*beta,z);
    Hv.scale(2/(q*q));
  }

private:
  HermitianMatrix<Real>                   A_;
  ROL::StdVector<Real,std::complex<Real>> y_;

};

template<typename Real>
inline ROL::Ptr<RayleighQuotient<Real>> 
make_RayleighQuotient( HermitianMatrix<Real> A ) {
  return ROL::makePtr<RayleighQuotient<Real>>(A);
}

#endif //RAYLEIGHQUOTIENT_HPP

