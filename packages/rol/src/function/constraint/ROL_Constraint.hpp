// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_CONSTRAINT_H
#define ROL_CONSTRAINT_H

#include "ROL_Vector.hpp"
#include "ROL_UpdateType.hpp"
#include "ROL_Types.hpp"
#include <iostream>

/** @ingroup func_group
    \class ROL::Constraint
    \brief Defines the general constraint operator interface.

    ROL's constraint interface is designed for Fr&eacute;chet differentiable
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
class Constraint {
private:
  bool activated_;

public:
  virtual ~Constraint(void) {}

  Constraint(void) : activated_(true) {}

  /** \brief Update constraint function. 

      This function updates the constraint function at new iterations. 
      @param[in]          x      is the new iterate. 
      @param[in]          type   is the type of update requested.
      @param[in]          iter   is the outer algorithm iterations count.
  */
  virtual void update( const Vector<Real> &x, UpdateType type, int iter = -1 ) {
    ROL_UNUSED(x);
    ROL_UNUSED(type);
    ROL_UNUSED(iter);
  }

  /** \brief Update constraint functions.  
                x is the optimization variable, 
                flag = true if optimization variable is changed,
                iter is the outer algorithm iterations count.
  */
  virtual void update( const Vector<Real> &x, bool flag = true, int iter = -1 ) {}

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


  /** \brief Apply the adjoint of the the constraint Jacobian at \f$x\f$, \f$c'(x)^* \in L(\mathcal{C}^*, \mathcal{X}^*)\f$,
             to vector \f$v\f$.

             @param[out]      ajv is the result of applying the adjoint of the constraint Jacobian to @b v at @b x; a dual optimization-space vector
             @param[in]       v   is a dual constraint-space vector
             @param[in]       x   is the constraint argument; an optimization-space vector
             @param[in]       dualv  is a vector used for temporary variables; a constraint-space vector
             @param[in,out]   tol is a tolerance for inexact evaluations; currently unused

             On return, \f$\mathsf{ajv} = c'(x)^*v\f$, where
             \f$v \in \mathcal{C}^*\f$, \f$\mathsf{ajv} \in \mathcal{X}^*\f$. \n\n
             The default implementation is a finite-difference approximation.

             ---
  */

  virtual void applyAdjointJacobian(Vector<Real> &ajv,
                                    const Vector<Real> &v,
                                    const Vector<Real> &x,
                                    const Vector<Real> &dualv,
                                    Real &tol);


  /** \brief Apply the derivative of the adjoint of the constraint Jacobian at \f$x\f$
             to vector \f$u\f$ in direction \f$v\f$,
             according to \f$ v \mapsto c''(x)(v,\cdot)^*u \f$.

             @param[out]      ahuv is the result of applying the derivative of the adjoint of the constraint Jacobian at @b x to vector @b u in direction @b v; a dual optimization-space vector
             @param[in]       u    is the direction vector; a dual constraint-space vector
             @param[in]       v    is an optimization-space vector
             @param[in]       x    is the constraint argument; an optimization-space vector
             @param[in,out]   tol  is a tolerance for inexact evaluations; currently unused

             On return, \f$ \mathsf{ahuv} = c''(x)(v,\cdot)^*u \f$, where
             \f$u \in \mathcal{C}^*\f$, \f$v \in \mathcal{X}\f$, and \f$\mathsf{ahuv} \in \mathcal{X}^*\f$. \n\n
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
             \f$I : \mathcal{X} \rightarrow \mathcal{X}^*\f$ is an identity or Riesz
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


  /** \brief Apply a constraint preconditioner at \f$x\f$, \f$P(x) \in L(\mathcal{C}, \mathcal{C}^*)\f$,
             to vector \f$v\f$.  Ideally, this preconditioner satisfies the following relationship:
             \f[
               \left[c'(x) \circ R \circ c'(x)^* \circ P(x)\right] v = v \,,
             \f]
             where R is the appropriate Riesz map in \f$L(\mathcal{X}^*, \mathcal{X})\f$.  It is used by the #solveAugmentedSystem method.

             @param[out]      pv  is the result of applying the constraint preconditioner to @b v at @b x; a dual constraint-space vector
             @param[in]       v   is a constraint-space vector
             @param[in]       x   is the preconditioner argument; an optimization-space vector
             @param[in]       g   is the preconditioner argument; a dual optimization-space vector, unused
             @param[in,out]   tol is a tolerance for inexact evaluations

             On return, \f$\mathsf{pv} = P(x)v\f$, where
             \f$v \in \mathcal{C}\f$, \f$\mathsf{pv} \in \mathcal{C}^*\f$. \n\n
             The default implementation is the Riesz map in \f$L(\mathcal{C}, \mathcal{C}^*)\f$.

             ---
  */
  virtual void applyPreconditioner(Vector<Real> &pv,
                                   const Vector<Real> &v,
                                   const Vector<Real> &x,
                                   const Vector<Real> &g,
                                   Real &tol) {
    pv.set(v.dual());
  }

  /** \brief Turn on constraints 
  */
  void activate(void)    { activated_ = true;  }

  /** \brief Turn off constraints
  */
  void deactivate(void)  { activated_ = false; }

  /** \brief Check if constraints are on
  */
  bool isActivated(void) { return activated_;  }

  /** \brief Finite-difference check for the constraint Jacobian application.

      Details here.
  */
  virtual std::vector<std::vector<Real> > checkApplyJacobian( const Vector<Real> &x,
                                                              const Vector<Real> &v,
                                                              const Vector<Real> &jv,
                                                              const std::vector<Real> &steps,
                                                              const bool printToStream = true,
                                                              std::ostream & outStream = std::cout,
                                                              const int order = 1 ) ;


  /** \brief Finite-difference check for the constraint Jacobian application.

      Details here.

  */
  virtual std::vector<std::vector<Real> > checkApplyJacobian( const Vector<Real> &x,
                                                              const Vector<Real> &v,
                                                              const Vector<Real> &jv,
                                                              const bool printToStream = true,
                                                              std::ostream & outStream = std::cout,
                                                              const int numSteps = ROL_NUM_CHECKDERIV_STEPS,
                                                              const int order = 1 ) ;

  /** \brief Finite-difference check for the application of the adjoint of constraint Jacobian.

      Details here. (This function should be deprecated)

  */
  virtual std::vector<std::vector<Real> > checkApplyAdjointJacobian(const Vector<Real> &x,
                                                                    const Vector<Real> &v,
                                                                    const Vector<Real> &c,
                                                                    const Vector<Real> &ajv,
                                                                    const bool printToStream = true,
                                                                    std::ostream & outStream = std::cout,
                                                                    const int numSteps = ROL_NUM_CHECKDERIV_STEPS ) ;

  /* \brief Check the consistency of the Jacobian and its adjoint. Verify that the deviation 
     \f$|\langle w^\top,Jv\rangle-\langle adj(J)w,v|\f$ is sufficiently small. 

     @param[in]      w              is a dual constraint-space vector \f$w\in \mathcal{C}^\ast\f$
     @param[in]      v              is an optimization space vector \f$v\in \mathcal{X}\f$
     @param[in]      x              is the constraint argument \f$x\in\mathcal{X}\f$
     @param[in]      printToStream  is is a flag that turns on/off output
     @param[in]      outStream      is the output stream

     Returns the deviation.
 */

  virtual Real checkAdjointConsistencyJacobian(const Vector<Real> &w,
                                               const Vector<Real> &v,
                                               const Vector<Real> &x,
                                               const bool printToStream = true,
                                               std::ostream & outStream = std::cout) {
    return checkAdjointConsistencyJacobian(w, v, x, w.dual(), v.dual(), printToStream, outStream);
  }

  virtual Real checkAdjointConsistencyJacobian(const Vector<Real> &w,
                                               const Vector<Real> &v,
                                               const Vector<Real> &x,
                                               const Vector<Real> &dualw, 
                                               const Vector<Real> &dualv, 
                                               const bool printToStream = true,
                                               std::ostream & outStream = std::cout);


  /** \brief Finite-difference check for the application of the adjoint of constraint Hessian.

      Details here.

  */
  virtual std::vector<std::vector<Real> > checkApplyAdjointHessian(const Vector<Real> &x,
                                                                   const Vector<Real> &u,
                                                                   const Vector<Real> &v,
                                                                   const Vector<Real> &hv,
                                                                   const std::vector<Real> &step,
                                                                   const bool printToScreen = true,
                                                                   std::ostream & outStream = std::cout,
                                                                   const int order = 1 ) ;
  /** \brief Finite-difference check for the application of the adjoint of constraint Hessian.

      Details here.

  */
  virtual std::vector<std::vector<Real> > checkApplyAdjointHessian(const Vector<Real> &x,
                                                                   const Vector<Real> &u,
                                                                   const Vector<Real> &v,
                                                                   const Vector<Real> &hv,
                                                                   const bool printToScreen = true,
                                                                   std::ostream & outStream = std::cout,
                                                                   const int numSteps = ROL_NUM_CHECKDERIV_STEPS,
                                                                   const int order = 1 ) ;

// Definitions for parametrized (stochastic) constraints
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

}; // class Constraint

} // namespace ROL

#include "ROL_ConstraintDef.hpp"

#endif
