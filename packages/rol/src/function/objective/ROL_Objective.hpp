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

#ifndef ROL_OBJECTIVE_H
#define ROL_OBJECTIVE_H

#include "ROL_Vector.hpp"
#include "ROL_Types.hpp"
#include <iostream>

/** @ingroup func_group
    \class ROL::Objective
    \brief Provides the interface to evaluate objective functions.

    ROL's objective function interface is designed for Fr$eacute;chet differentiable 
    functionals \f$f:\mathcal{X}\to\mathbb{R}\f$, where \f$\mathcal{X}\f$ is a Banach
    space.  The basic operator interace, to be implemented by the user, requires:
    \li #value -- objective function evaluation.

    It is strongly recommended that the user additionally overload:
    \li #gradient -- the objective function gradient -- the default is a finite-difference approximation;
    \li #hessVec  -- the action of the Hessian -- the default is a finite-difference approximation.

    The user may also overload:
    \li #update     -- update the objective function at each new iteration;
    \li #dirDeriv   -- compute the directional derivative -- the default is a finite-difference approximation;
    \li #invHessVec -- the action of the inverse Hessian;
    \li #precond    -- the action of a preconditioner for the Hessian.

    ---
*/


namespace ROL {

template <class Real>
class Objective {
public:

  virtual ~Objective() {}

  /** \brief Update objective function. 

      This function updates the objective function at new iterations. 
      @param[in]          x      is the new iterate. 
      @param[in]          flag   is true if the iterate has changed.
      @param[in]          iter   is the outer algorithm iterations count.
  */
  virtual void update( const Vector<Real> &x, bool flag = true, int iter = -1 ) {
    ROL_UNUSED(x);
    ROL_UNUSED(flag);
    ROL_UNUSED(iter);
  }

  /** \brief Compute value.

      This function returns the objective function value.
      @param[in]          x   is the current iterate.
      @param[in]          tol is a tolerance for inexact objective function computation.
  */
  virtual Real value( const Vector<Real> &x, Real &tol ) = 0;

  /** \brief Compute gradient.

      This function returns the objective function gradient.
      @param[out]         g   is the gradient.
      @param[in]          x   is the current iterate.
      @param[in]          tol is a tolerance for inexact objective function computation.

      The default implementation is a finite-difference approximation based on the function value.
      This requires the definition of a basis \f$\{\phi_i\}\f$ for the optimization vectors x and
      the definition of a basis \f$\{\psi_j\}\f$ for the dual optimization vectors (gradient vectors g).
      The bases must be related through the Riesz map, i.e., \f$ R \{\phi_i\} = \{\psi_j\}\f$,
      and this must be reflected in the implementation of the ROL::Vector::dual() method.
  */
  virtual void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) ;

  /** \brief Compute directional derivative.

      This function returns the directional derivative of the objective function in the \f$d\f$ direction.
      @param[in]          x   is the current iterate.
      @param[in]          d   is the direction.
      @param[in]          tol is a tolerance for inexact objective function computation.
  */
  virtual Real dirDeriv( const Vector<Real> &x, const Vector<Real> &d, Real &tol ) ;

  /** \brief Apply Hessian approximation to vector.

      This function applies the Hessian of the objective function to the vector \f$v\f$.
      @param[out]         hv  is the the action of the Hessian on \f$v\f$.
      @param[in]          v   is the direction vector.
      @param[in]          x   is the current iterate.
      @param[in]          tol is a tolerance for inexact objective function computation.
  */
  virtual void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol );

  /** \brief Apply inverse Hessian approximation to vector.

      This function applies the inverse Hessian of the objective function to the vector \f$v\f$.
      @param[out]         hv  is the action of the inverse Hessian on \f$v\f$.
      @param[in]          v   is the direction vector.
      @param[in]          x   is the current iterate.
      @param[in]          tol is a tolerance for inexact objective function computation.
  */
  virtual void invHessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
    ROL_UNUSED(hv);
    ROL_UNUSED(v);
    ROL_UNUSED(x);
    ROL_UNUSED(tol);
    ROL_TEST_FOR_EXCEPTION(true, std::invalid_argument,
      ">>> ERROR (ROL::Objective): invHessVec not implemented!"); 
    //hv.set(v.dual());
  }

  /** \brief Apply preconditioner to vector.

      This function applies a preconditioner for the Hessian of the objective function to the vector \f$v\f$.
      @param[out]         Pv  is the action of the Hessian preconditioner on \f$v\f$.
      @param[in]          v   is the direction vector.
      @param[in]          x   is the current iterate.
      @param[in]          tol is a tolerance for inexact objective function computation.
  */
  virtual void precond( Vector<Real> &Pv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
    ROL_UNUSED(x);
    ROL_UNUSED(tol);
    Pv.set(v.dual());
  }

  /** \brief Finite-difference gradient check.

      This function computes a sequence of one-sided finite-difference checks for the gradient.  
      At each step of the sequence, the finite difference step size is decreased.  The output 
      compares the error 
      \f[
          \left| \frac{f(x+td) - f(x)}{t} - \langle \nabla f(x),d\rangle_{\mathcal{X}^*,\mathcal{X}}\right|.
      \f]
      if the approximation is first order. More generally, difference approximation is
      \f[
          \frac{1}{t} \sum\limits_{i=1}^m w_i f(x+t c_i d)     
      \f]
      where m = order+1, \f$w_i\f$ are the difference weights and \f$c_i\f$ are the difference steps
      @param[in]      x             is an optimization variable.
      @param[in]      d             is a direction vector.
      @param[in]      printToStream is a flag that turns on/off output.
      @param[out]     outStream     is the output stream.
      @param[in]      numSteps      is a parameter which dictates the number of finite difference steps.
      @param[in]      order         is the order of the finite difference approximation (1,2,3,4)
  */
  virtual std::vector<std::vector<Real> > checkGradient( const Vector<Real> &x,
                                                         const Vector<Real> &d,
                                                         const bool printToStream = true,
                                                         std::ostream & outStream = std::cout,
                                                         const int numSteps = ROL_NUM_CHECKDERIV_STEPS,
                                                         const int order = 1 ) {
    return checkGradient(x, x.dual(), d, printToStream, outStream, numSteps, order);
  }

  /** \brief Finite-difference gradient check.

      This function computes a sequence of one-sided finite-difference checks for the gradient.  
      At each step of the sequence, the finite difference step size is decreased.  The output 
      compares the error 
      \f[
          \left| \frac{f(x+td) - f(x)}{t} - \langle \nabla f(x),d\rangle_{\mathcal{X}^*,\mathcal{X}}\right|.
      \f]
      if the approximation is first order. More generally, difference approximation is
      \f[
          \frac{1}{t} \sum\limits_{i=1}^m w_i f(x+t c_i d)     
      \f]
      where m = order+1, \f$w_i\f$ are the difference weights and \f$c_i\f$ are the difference steps
 
      @param[in]      x             is an optimization variable.
      @param[in]      g             is used to create a temporary gradient vector.
      @param[in]      d             is a direction vector.
      @param[in]      printToStream is a flag that turns on/off output.
      @param[out]     outStream     is the output stream.
      @param[in]      numSteps      is a parameter which dictates the number of finite difference steps.
      @param[in]      order         is the order of the finite difference approximation (1,2,3,4)
  */
  virtual std::vector<std::vector<Real> > checkGradient( const Vector<Real> &x,
                                                         const Vector<Real> &g,
                                                         const Vector<Real> &d,
                                                         const bool printToStream = true,
                                                         std::ostream & outStream = std::cout,
                                                         const int numSteps = ROL_NUM_CHECKDERIV_STEPS,
                                                         const int order = 1 );


  /** \brief Finite-difference gradient check with specified step sizes.

      This function computes a sequence of one-sided finite-difference checks for the gradient.  
      At each step of the sequence, the finite difference step size is decreased.  The output 
      compares the error 
      \f[
          \left| \frac{f(x+td) - f(x)}{t} - \langle \nabla f(x),d\rangle_{\mathcal{X}^*,\mathcal{X}}\right|.
      \f]
      if the approximation is first order. More generally, difference approximation is
      \f[
          \frac{1}{t} \sum\limits_{i=1}^m w_i f(x+t c_i d)     
      \f]
      where m = order+1, \f$w_i\f$ are the difference weights and \f$c_i\f$ are the difference steps
      @param[in]      x             is an optimization variable.
      @param[in]      d             is a direction vector.
      @param[in]      steps         is vector of steps of user-specified size.
      @param[in]      printToStream is a flag that turns on/off output.
      @param[out]     outStream     is the output stream.
      @param[in]      order         is the order of the finite difference approximation (1,2,3,4)
  */
  virtual std::vector<std::vector<Real> > checkGradient( const Vector<Real> &x,
                                                         const Vector<Real> &d,
                                                         const std::vector<Real> &steps,
                                                         const bool printToStream = true,
                                                         std::ostream & outStream = std::cout,
                                                         const int order = 1 ) {

    return checkGradient(x, x.dual(), d, steps, printToStream, outStream, order);

  }


  /** \brief Finite-difference gradient check with specified step sizes.

      This function computes a sequence of one-sided finite-difference checks for the gradient.  
      At each step of the sequence, the finite difference step size is decreased.  The output 
      compares the error 
      \f[
          \left| \frac{f(x+td) - f(x)}{t} - \langle \nabla f(x),d\rangle_{\mathcal{X}^*,\mathcal{X}}\right|.
      \f]
      if the approximation is first order. More generally, difference approximation is
      \f[
          \frac{1}{t} \sum\limits_{i=1}^m w_i f(x+t c_i d)     
      \f]
      where m = order+1, \f$w_i\f$ are the difference weights and \f$c_i\f$ are the difference steps
 
      @param[in]      x             is an optimization variable.
      @param[in]      g             is used to create a temporary gradient vector.
      @param[in]      d             is a direction vector.
      @param[in]      steps         is vector of steps of user-specified size.
      @param[in]      printToStream is a flag that turns on/off output.
      @param[out]     outStream     is the output stream.
      @param[in]      order         is the order of the finite difference approximation (1,2,3,4)
  */
  virtual std::vector<std::vector<Real> > checkGradient( const Vector<Real> &x,
                                                         const Vector<Real> &g,
                                                         const Vector<Real> &d,
                                                         const std::vector<Real> &steps,
                                                         const bool printToStream = true,
                                                         std::ostream & outStream = std::cout,
                                                         const int order = 1 );

  /** \brief Finite-difference Hessian-applied-to-vector check.

      This function computes a sequence of one-sided finite-difference checks for the Hessian.  
      At each step of the sequence, the finite difference step size is decreased.  The output 
      compares the error 
      \f[
          \left\| \frac{\nabla f(x+tv) - \nabla f(x)}{t} - \nabla^2 f(x)v\right\|_{\mathcal{X}^*},
      \f]
      if the approximation is first order. More generally, difference approximation is
      \f[
          \frac{1}{t} \sum\limits_{i=1}^m w_i \nabla f(x+t c_i v),     
      \f]
      where m = order+1, \f$w_i\f$ are the difference weights and \f$c_i\f$ are the difference steps
      @param[in]      x             is an optimization variable.
      @param[in]      v             is a direction vector.
      @param[in]      printToStream is a flag that turns on/off output.
      @param[out]     outStream     is the output stream.
      @param[in]      numSteps      is a parameter which dictates the number of finite difference steps.
      @param[in]      order         is the order of the finite difference approximation (1,2,3,4)
  */
  virtual std::vector<std::vector<Real> > checkHessVec( const Vector<Real> &x,
                                                        const Vector<Real> &v,
                                                        const bool printToStream = true,
                                                        std::ostream & outStream = std::cout,
                                                        const int numSteps = ROL_NUM_CHECKDERIV_STEPS,
                                                        const int order = 1 ) {

    return checkHessVec(x, x.dual(), v, printToStream, outStream, numSteps, order);

  }

  /** \brief Finite-difference Hessian-applied-to-vector check.

      This function computes a sequence of one-sided finite-difference checks for the Hessian.  
      At each step of the sequence, the finite difference step size is decreased.  The output 
      compares the error 
      \f[
          \left\| \frac{\nabla f(x+tv) - \nabla f(x)}{t} - \nabla^2 f(x)v\right\|_{\mathcal{X}^*},
      \f]
      if the approximation is first order. More generally, difference approximation is
      \f[
          \frac{1}{t} \sum\limits_{i=1}^m w_i \nabla f(x+t c_i v),     
      \f]
      where m = order+1, \f$w_i\f$ are the difference weights and \f$c_i\f$ are the difference steps
      @param[in]      x             is an optimization variable.
      @param[in]      hv            is used to create temporary gradient and Hessian-times-vector vectors.
      @param[in]      v             is a direction vector.
      @param[in]      printToStream is a flag that turns on/off output.
      @param[out]     outStream     is the output stream.
      @param[in]      numSteps      is a parameter which dictates the number of finite difference steps.
      @param[in]      order         is the order of the finite difference approximation (1,2,3,4)
  */
  virtual std::vector<std::vector<Real> > checkHessVec( const Vector<Real> &x,
                                                        const Vector<Real> &hv,
                                                        const Vector<Real> &v,
                                                        const bool printToStream = true,
                                                        std::ostream & outStream = std::cout,
                                                        const int numSteps = ROL_NUM_CHECKDERIV_STEPS,
                                                        const int order = 1) ;


  /** \brief Finite-difference Hessian-applied-to-vector check with specified step sizes.

      This function computes a sequence of one-sided finite-difference checks for the Hessian.  
      At each step of the sequence, the finite difference step size is decreased.  The output 
      compares the error 
      \f[
          \left\| \frac{\nabla f(x+tv) - \nabla f(x)}{t} - \nabla^2 f(x)v\right\|_{\mathcal{X}^*},
      \f]
      if the approximation is first order. More generally, difference approximation is
      \f[
          \frac{1}{t} \sum\limits_{i=1}^m w_i \nabla f(x+t c_i v),     
      \f]
      where m = order+1, \f$w_i\f$ are the difference weights and \f$c_i\f$ are the difference steps
      @param[in]      x             is an optimization variable.
      @param[in]      v             is a direction vector.
      @param[in]      steps         is vector of steps of user-specified size.
      @param[in]      printToStream is a flag that turns on/off output.
      @param[out]     outStream     is the output stream.
      @param[in]      order         is the order of the finite difference approximation (1,2,3,4)
  */
  virtual std::vector<std::vector<Real> > checkHessVec( const Vector<Real> &x,
                                                        const Vector<Real> &v,
                                                        const std::vector<Real> &steps,
                                                        const bool printToStream = true,
                                                        std::ostream & outStream = std::cout,
                                                        const int order = 1 ) {

    return checkHessVec(x, x.dual(), v, steps, printToStream, outStream, order);

  }

  /** \brief Finite-difference Hessian-applied-to-vector check with specified step sizes.

      This function computes a sequence of one-sided finite-difference checks for the Hessian.  
      At each step of the sequence, the finite difference step size is decreased.  The output 
      compares the error 
      \f[
          \left\| \frac{\nabla f(x+tv) - \nabla f(x)}{t} - \nabla^2 f(x)v\right\|_{\mathcal{X}^*},
      \f]
      if the approximation is first order. More generally, difference approximation is
      \f[
          \frac{1}{t} \sum\limits_{i=1}^m w_i \nabla f(x+t c_i v),     
      \f]
      where m = order+1, \f$w_i\f$ are the difference weights and \f$c_i\f$ are the difference steps
      @param[in]      x             is an optimization variable.
      @param[in]      hv            is used to create temporary gradient and Hessian-times-vector vectors.
      @param[in]      v             is a direction vector.
      @param[in]      steps         is vector of steps of user-specified size.
      @param[in]      printToStream is a flag that turns on/off output.
      @param[out]     outStream     is the output stream.
      @param[in]      order         is the order of the finite difference approximation (1,2,3,4)
  */
  virtual std::vector<std::vector<Real> > checkHessVec( const Vector<Real> &x,
                                                        const Vector<Real> &hv,
                                                        const Vector<Real> &v,
                                                        const std::vector<Real> &steps,
                                                        const bool printToStream = true,
                                                        std::ostream & outStream = std::cout,
                                                        const int order = 1) ;


  /** \brief Hessian symmetry check.

      This function checks the symmetry of the Hessian by comparing 
      \f[
         \langle \nabla^2f(x)v,w\rangle_{\mathcal{X}^*,\mathcal{X}}
         \quad\text{and}\quad
         \langle \nabla^2f(x)w,v\rangle_{\mathcal{X}^*,\mathcal{X}}.
      \f]
      @param[in]      x             is an optimization variable.
      @param[in]      v             is a direction vector.
      @param[in]      w             is a direction vector.
      @param[in]      printToStream is a flag that turns on/off output.
      @param[out]     outStream     is the output stream.
  */
  virtual std::vector<Real> checkHessSym( const Vector<Real> &x,
                                          const Vector<Real> &v,
                                          const Vector<Real> &w,
                                          const bool printToStream = true,
                                          std::ostream & outStream = std::cout ) {

    return checkHessSym(x, x.dual(), v, w, printToStream, outStream);

  }

  /** \brief Hessian symmetry check.

      This function checks the symmetry of the Hessian by comparing 
      \f[
         \langle \nabla^2f(x)v,w\rangle_{\mathcal{X}^*,\mathcal{X}}
         \quad\text{and}\quad
         \langle \nabla^2f(x)w,v\rangle_{\mathcal{X}^*,\mathcal{X}}.
      \f]
      @param[in]      x             is an optimization variable.
      @param[in]      hv            is used to create temporary Hessian-times-vector vectors.
      @param[in]      v             is a direction vector.
      @param[in]      w             is a direction vector.
      @param[in]      printToStream is a flag that turns on/off output.
      @param[out]     outStream     is the output stream.
  */
  virtual std::vector<Real> checkHessSym( const Vector<Real> &x,
                                          const Vector<Real> &hv,
                                          const Vector<Real> &v,
                                          const Vector<Real> &w,
                                          const bool printToStream = true,
                                          std::ostream & outStream = std::cout );

// Definitions for parametrized (stochastic) objective functions
private:
  std::vector<Real> param_;

protected:
  const std::vector<Real> getParameter(void) const {
    return param_;
  }

public:
  virtual void setParameter(const std::vector<Real> &param) {
    param_.assign(param.begin(),param.end());
  }

}; // class Objective

} // namespace ROL

// include templated definitions
#include <ROL_ObjectiveDef.hpp>

#endif
