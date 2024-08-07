// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef NONLINEAR_PROGRAM_HPP
#define NONLINEAR_PROGRAM_HPP

#include "Sacado.hpp"

#include "ROL_StdVector.hpp"
#include "ROL_Bounds.hpp"
#include "ROL_PartitionedVector.hpp"
#include "ROL_OptimizationProblem.hpp"
#include "ROL_Sacado_StdObjective.hpp"
#include "ROL_Sacado_StdConstraint.hpp"
#include <cmath>

/* \class ROL::NonlinearProgram 
 * \brief Provides a simple interface to formulating and solving
 * modest sized nonlinear programs consisting of objective functions with 
 * optional equality, inequality, and bound constraints. Creates 
 * OptimizationProblem objects which encapsulate the quantities needed
 * to represent a problem as Type-U, Type-E, Type-B, or Type-EB.
 *
 * Automatic differentiation to obtain gradients, Hessians, and Jacobians
 * is provided by Sacado, so the user only needs to implement value functions
 * for the objective and constraints and provide upper and lower bounds 
 * when they exist.     
 *
 * NonlinearProgram assumes the user's Objective and Constraints are 
 * implemented using ROL::StdVector
 */

namespace ROL {

template<class Real>
class NonlinearProgram {

  template <typename T> using vector = std::vector<T>;  
 
protected:

  typedef Vector<Real>                V;
  typedef StdVector<Real>             SV;
  typedef PartitionedVector<Real>     PV;

  typedef Objective<Real>             OBJ;
  typedef Constraint<Real>            CON;
  typedef Bounds<Real>                BND;
  typedef OptimizationProblem<Real>   OPT;

  typedef typename vector<Real>::size_type size_type;

private:

  ROL::Ptr<vector<Real> > lp_;
  ROL::Ptr<vector<Real> > up_;
  ROL::Ptr<BND> bnd_;

protected:

  virtual ~NonlinearProgram(){}

  /* Set the ith element of the of lower bound equal to val \f$\ell_i=\text{val}\f$ 
     @param[in]   i   is the element to set
     @param[in]   val is the value to set it to 
  */
  void setLower(int i, Real val) {
    (*lp_)[i] = val;  
  }

  /* Set the ith element of the of upper bound equal to val \f$u_=\text{val}\f$
     @param[in]   i   is the element to set
     @param[in]   val is the value to set it to  
  */
  void setUpper(int i, Real val) {
    (*up_)[i] = val;
  }

  /* Set the lower bound to l */
  void setLower(Real l[]) {
    lp_ = ROL::makePtr<vector<Real>>(l, l+dimension_x() );
  }

  /* Set the upper bound to u */
  void setUpper(Real u[]) {
    up_ = ROL::makePtr<vector<Real>>(u, u+dimension_x() );
  }

  /* \brief Call this function in your problem's constructor to turn 
      off bound constraints */
  void noBound() {
    bnd_->deactivate();   
  }
  
  /* \brief Return the objective function. */
  virtual const ROL::Ptr<OBJ> getObjective() = 0;

  /* \brief Return the equality constraint. Do not overload if you 
            have no equality constraint. */
  virtual const ROL::Ptr<CON> getEqualityConstraint() {
    return ROL::nullPtr;
  }

  /* \brief Return the equality constraint vector (used for cloning).
            Returns a ROL::nullPtr if there is no equality constraint */
  const ROL::Ptr<V> getEqualityMultiplier() { 
    int n = dimension_ce();
    if( n > 0 ) {
      ROL::Ptr<vector<Real> > l = ROL::makePtr<vector<Real>>(n,1.0);
      return ROL::makePtr<SV>(l);
    }
    else {
      return ROL::nullPtr;
    }
  }

  /* \brief Return the inequality constraint. Do not overload if you 
            have no inequality constraint. */
  virtual const ROL::Ptr<CON> getInequalityConstraint() {
    return ROL::nullPtr;
  }

  /* \brief Return the inequality constraint vector (used for cloning). 
            Returns a ROL::nullPtr if there is no inequality constraint */
  const ROL::Ptr<V> getInequalityMultiplier() {
    int n = dimension_ci();
    if( n > 0 ) {
      ROL::Ptr<vector<Real> > l = ROL::makePtr<vector<Real>>(n,1.0);
      return ROL::makePtr<SV>(l);
    }
    else {
      return ROL::nullPtr;
    }
  }

  /* \brief Return the bounds on the inequality constraint.
            Returns a ROL::nullPtr if there is no inequality constraint */
  const ROL::Ptr<BND> getInequalityBoundConstraint() {
    int n = dimension_ci();
    if( n > 0 ) {
      const Real lval(0), uval(ROL_INF<Real>());
      ROL::Ptr<V> l = ROL::makePtr<SV>(
                          ROL::makePtr<vector<Real>>(n,lval) );
      ROL::Ptr<V> u = ROL::makePtr<SV>(
                          ROL::makePtr<vector<Real>>(n,uval) );
      return ROL::makePtr<BND>(l,u);
    }
    else {
      return ROL::nullPtr;
    }
  }

  /* \brief Create vector */
  ROL::Ptr<V> createOptVector( const Real * const array ) {
    ROL::Ptr<vector<Real> > x = 
      ROL::makePtr<vector<Real>>(array,array+dimension_x());
    return ROL::makePtr<SV>( x ); 
  }

  // Create an ROL::Ptr to a std::vector with dimensionality of the optimization space
//  ROL::Ptr<vector<Real> > createOptVector() {
//    int n = dimension_x();
//    ROL::Ptr<V> x = ROL::makePtr<vector<Real>>(n);
//    return x;
//  } 

  // Default deactivated bound
  ROL::Ptr<BND> getBoundConstraint() {
    return bnd_;
  }

public:
  /** \brief Return the dimension of the optimization (and bound) vectors. */
  virtual int dimension_x() { return 1; }

  /** \brief Return the dimension of the equality constraint vector  */
  virtual int dimension_ce() { return 0; }

  /** \brief Return the dimension of the inequality constrant vector */
  virtual int dimension_ci() { return 0; }


  NonlinearProgram(int n) {
     // Create lower and upper bounds of negative and positive infinity respectively
     lp_ = ROL::makePtr<vector<Real>>(n,ROL::ROL_NINF<Real>());
     up_ = ROL::makePtr<vector<Real>>(n,ROL::ROL_INF<Real>());
     ROL::Ptr<V> l = ROL::makePtr<SV>( lp_ );
     ROL::Ptr<V> u = ROL::makePtr<SV>( up_ );
     bnd_ = ROL::makePtr<BND>( l, u ); 
  }

   /* \brief Create the OptimizationProblem from the supplied components and 
             return it */
  ROL::Ptr<OPT> getOptimizationProblem() {
    ROL::Ptr<V> x = getInitialGuess()->clone();
    x->set(*getInitialGuess());
    return ROL::makePtr<OPT>( getObjective(),
                                   x,
                                   getBoundConstraint(),
                                   getEqualityConstraint(),
                                   getEqualityMultiplier(),
                                   getInequalityConstraint(),
                                   getInequalityMultiplier(),
                                   getInequalityBoundConstraint() );
  }

  /* \brief Return the initial guess for the optimization vector */
  virtual const ROL::Ptr<const V> getInitialGuess() = 0;

  /* \brief Return whether or not the initial guess is feasible */
  virtual bool initialGuessIsFeasible() = 0;

  /* \brief Return the value of the Objective function for the initial guess */
  virtual Real getInitialObjectiveValue() = 0;

  /* \brief Return the set of vectors that solve the nonlinear program if
            they are known */
  virtual ROL::Ptr<const V> getSolutionSet() { return ROL::nullPtr; } 
  
  /* \brief Return the value of the objective function for a solution vector. 
            If not known, return infinity */
  virtual Real getSolutionObjectiveValue() { 
     return ROL_INF<Real>(); 
  };
   
   /* \brief If the problem has known solutions, return whether ROL
             has acceptibly solved problem */
  bool foundAcceptableSolution( const V &x, const Real &tolerance=std::sqrt(ROL_EPSILON<Real>()) ) {
    ROL::Ptr<const PV> sol = ROL::dynamicPtrCast<const PV>( getSolutionSet() );
    ROL::Ptr<V> error;
    ROL::Ptr<const V> xv;
    if ( dimension_ci() > 0 ) {
      xv = dynamic_cast<const PV&>(x).get(0);
      error = xv->clone();
    }
    else {
      xv = ROL::makePtrFromRef(x);
      error = x.clone();
    }
 
    Real minerror = ROL::ROL_INF<Real>();
    for( size_type k=0; k<sol->numVectors(); ++k) {
      error->set(*xv);
      error->axpy(-1.0,*(sol->get(k)));
      minerror = std::min(minerror,error->norm());
    }
    return (minerror < tolerance);
  }
  
   

}; // class NonlinearProgram

} // namespace ROL

#endif // NONLINEAR_PROGRAM_HPP
