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

#ifndef ROL_EQUALITY_CONSTRAINT_H
#define ROL_EQUALITY_CONSTRAINT_H

#include "ROL_Vector.hpp"
#include "ROL_Types.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseVector.hpp"
#include "Teuchos_LAPACK.hpp"
#include <iostream>

/** @ingroup func_group
    \class ROL::EqualityConstraint
    \brief Defines the equality constraint operator interface.

    ROL's equality constraint interface is designed for Fr&eacute;chet differentiable
    operators \f$c:\mathcal{X} \rightarrow \mathcal{C}\f$, where \f$\mathcal{X}\f$
    and \f$\mathcal{C}\f$ are Banach spaces.  The constraints are of the form
    \f[
      c(x) = 0 \,.
    \f]
    The basic operator interface, to be
    implemented by the user, requires:
    \li #value -- constraint evaluation.

    It is strongly recommended that the user additionally overload:
    \li #applyJacobian        -- action of the constraint Jacobian --the default is
                                 a finite-difference approximation;
    \li #applyAdjointJacobian -- action of the adjoint of the constraint Jacobian --the default is
                                 a finite-difference approximation.

    The user may also overload:
    \li #applyAdjointHessian  -- action of the adjoint of the constraint Hessian --the default
                                 is a finite-difference approximation based on the adjoint Jacobian;
    \li #solveAugmentedSystem -- solution of the augmented system --the default is an iterative
                                 scheme based on the action of the Jacobian and its adjoint.
    \li #applyPreconditioner  -- action of a constraint preconditioner --the default is null-op.

    ---
*/


namespace ROL {

template <class Real>
class EqualityConstraint {
private:
  bool activated_;

public:

  virtual ~EqualityConstraint() {}

  /** \brief Evaluate the constraint operator \f$c:\mathcal{X} \rightarrow \mathcal{C}\f$
             at \f$x\f$.

             @param[out]      c   is the result of evaluating the constraint operator at @b x; a constraint-space vector
             @param[in]       x   is the constraint argument; an optimization-space vector
             @param[in,out]   tol is a tolerance for inexact evaluations; currently unused

             On return, \f$\mathsf{c} = c(x)\f$,
             where \f$\mathsf{c} \in \mathcal{C}\f$, \f$\mathsf{x} \in \mathcal{X}\f$.

             ---
  */
  virtual void value(Vector<Real> &c,
                     const Vector<Real> &x,
                     Real &tol) = 0;

 
  /** \brief Apply the constraint Jacobian at \f$x\f$, \f$c'(x) \in L(\mathcal{X}, \mathcal{C})\f$,
             to vector \f$v\f$.

             @param[out]      jv  is the result of applying the constraint Jacobian to @b v at @b x; a constraint-space vector
             @param[in]       v   is an optimization-space vector
             @param[in]       x   is the constraint argument; an optimization-space vector
             @param[in,out]   tol is a tolerance for inexact evaluations; currently unused

             On return, \f$\mathsf{jv} = c'(x)v\f$, where
             \f$v \in \mathcal{X}\f$, \f$\mathsf{jv} \in \mathcal{C}\f$. \n\n
             The default implementation is a finite-difference approximation.

             ---
  */
  virtual void applyJacobian(Vector<Real> &jv,
                             const Vector<Real> &v,
                             const Vector<Real> &x,
                             Real &tol);


  /** \brief Apply the adjoint of the the constraint Jacobian at \f$x\f$, \f$c'(x)^* \in L(\mathcal{C}^*, \mathcal{X}^*)\f$,
             to vector \f$v\f$.

             @param[out]      ajv is the result of applying the adjoint of the constraint Jacobian to @b v at @b x; a dual optimization-space vector
             @param[in]       v   is a dual constraint-space vector
             @param[in]       x   is the constraint argument; an optimization-space vector
             @param[in,out]   tol is a tolerance for inexact evaluations; currently unused

             On return, \f$\mathsf{ajv} = c'(x)^*v\f$, where
             \f$v \in \mathcal{C}^*\f$, \f$\mathsf{ajv} \in \mathcal{X}^*\f$. \n\n
             The default implementation is a finite-difference approximation.

             ---
  */
  virtual void applyAdjointJacobian(Vector<Real> &ajv,
                                    const Vector<Real> &v,
                                    const Vector<Real> &x,
                                    Real &tol);


  /** \brief Apply the adjoint of the constraint Hessian at \f$x\f$,
             \f$c''(x)^* \in L(L(\mathcal{C}^*, \mathcal{X}^*), \mathcal{X}^*)\f$,
             to vector \f$v\f$ in direction \f$u\f$.

             @param[out]      ahuv is the result of applying the adjoint of the constraint Hessian to @b v at @b x in direction @b u; a dual optimization-space vector
             @param[in]       u    is the direction vector; a dual constraint-space vector
             @param[in]       v    is a dual optimization-space vector
             @param[in]       x    is the constraint argument; an optimization-space vector
             @param[in,out]   tol  is a tolerance for inexact evaluations; currently unused

             On return, \f$ \mathsf{ahuv} = c''(x)^*(u,v) \f$, where
             \f$u \in \mathcal{C}^*\f$, \f$v \in \mathcal{X}^*\f$, and \f$\mathsf{ahuv} \in \mathcal{X}^*\f$. \n\n
             The default implementation is a finite-difference approximation based on the adjoint Jacobian.

             ---
  */
  virtual void applyAdjointHessian(Vector<Real> &ahuv,
                                   const Vector<Real> &u,
                                   const Vector<Real> &v,
                                   const Vector<Real> &x,
                                   Real &tol);


  /** \brief Approximately solves the <em> augmented system </em>
             \f[
                 \begin{pmatrix}
                   I     & c'(x)^* \\
                   c'(x) & 0
                 \end{pmatrix}
                 \begin{pmatrix}
                   v_{1} \\
                   v_{2}
                 \end{pmatrix}
                 =
                 \begin{pmatrix}
                   b_{1} \\
                   b_{2}
                 \end{pmatrix}
             \f]
             where \f$v_{1} \in \mathcal{X}\f$, \f$v_{2} \in \mathcal{C}^*\f$,
             \f$b_{1} \in \mathcal{X}^*\f$, \f$b_{2} \in \mathcal{C}\f$,
             \f$I : \mathcal{X} \rightarrow \mathcal{X}^*\f$ is an identity
             operator, and \f$0 : \mathcal{C}^* \rightarrow \mathcal{C}\f$
             is a zero operator.

             @param[out]      v1  is the optimization-space component of the result
             @param[out]      v2  is the dual constraint-space component of the result
             @param[in]       b1  is the dual optimization-space component of the right-hand side
             @param[in]       b2  is the constraint-space component of the right-hand side
             @param[in]       x   is the constraint argument; an optimization-space vector
             @param[in,out]   tol is the nominal relative residual tolerance

             On return, \f$ [\mathsf{v1} \,\, \mathsf{v2}] \f$ approximately
             solves the augmented system, where the size of the residual is
             governed by special stopping conditions. \n\n
             The default implementation is the preconditioned generalized
             minimal residual (GMRES) method, which enables the use of
             nonsymmetric preconditioners.

             ---
  */
  virtual std::vector<Real> solveAugmentedSystem(Vector<Real> &v1,
                                                 Vector<Real> &v2,
                                                 const Vector<Real> &b1,
                                                 const Vector<Real> &b2,
                                                 const Vector<Real> &x,
                                                 Real &tol);


  /** \brief Apply a constraint preconditioner at \f$x\f$, \f$P(x) \in L(\mathcal{C}, \mathcal{C})\f$,
             to vector \f$v\f$.  In general, this preconditioner satisfies the following relationship:
             \f[
               c'(x) c'(x)^* P(x) v \approx v \,.
             \f]
             It is used by the #solveAugmentedSystem method.

             @param[out]      pv  is the result of applying the constraint preconditioner to @b v at @b x; a constraint-space vector
             @param[in]       v   is a constraint-space vector
             @param[in]       x   is the preconditioner argument; an optimization-space vector
             @param[in,out]   tol is a tolerance for inexact evaluations

             On return, \f$\mathsf{pv} = P(x)v\f$, where
             \f$v \in \mathcal{C}\f$, \f$\mathsf{pv} \in \mathcal{C}\f$. \n\n
             The default implementation is the identity operator.

             ---
  */
  virtual void applyPreconditioner(Vector<Real> &pv,
                                   const Vector<Real> &v,
                                   const Vector<Real> &x,
                                   Real &tol) {
    pv.set(v);
  }


  EqualityConstraint(void) : activated_(true) {}

  /** \brief Update constraint functions.  
                x is the optimization variable, 
                flag = true if optimization variable is changed,
                iter is the outer algorithm iterations count.
  */
  virtual void update( const Vector<Real> &x, bool flag = true, int iter = -1 ) {}

  /** \brief Check if the vector, v, is feasible
  */
  virtual bool isFeasible( const Vector<Real> &v ) { return true; }

  /** \brief Turn on constraints 
  */
  void activate(void)    { this->activated_ = true;  }

  /** \brief Turn off constraints
  */
  void deactivate(void)  { this->activated_ = false; }

  /** \brief Check if constraints are on
  */
  bool isActivated(void) { return this->activated_;  }

  /** \brief Finite-difference check for the constraint Jacobian application.
  */
  virtual std::vector<std::vector<Real> > checkApplyJacobian( const Vector<Real> &x,
                                                              const Vector<Real> &v,
                                                              const Vector<Real> &jv,
                                                              const bool printToScreen = true,
                                                              const int numSteps = ROL_NUM_CHECKDERIV_STEPS ) ;

  /** \brief Finite-difference check for the application of the adjoint of constraint Jacobian.
  */
  virtual std::vector<std::vector<Real> > checkApplyAdjointJacobian(const Vector<Real> &x,
                                                                    const Vector<Real> &v,
                                                                    const bool printToScreen = true,
                                                                    const int numSteps = ROL_NUM_CHECKDERIV_STEPS ) ;

  /** \brief Finite-difference check for the application of the adjoint of constraint Hessian.
  */
  virtual std::vector<std::vector<Real> > checkApplyAdjointHessian(const Vector<Real> &x,
                                                                   const Vector<Real> &u,
                                                                   const Vector<Real> &v,
                                                                   const bool printToScreen = true,
                                                                   const int numSteps = ROL_NUM_CHECKDERIV_STEPS ) ;

}; // class EqualityConstraint

} // namespace ROL

#include "ROL_EqualityConstraintDef.hpp"

#endif
