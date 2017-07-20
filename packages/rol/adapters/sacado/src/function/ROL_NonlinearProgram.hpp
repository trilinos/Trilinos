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

  Teuchos::RCP<vector<Real> > lp_;
  Teuchos::RCP<vector<Real> > up_;
  Teuchos::RCP<BND> bnd_;

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
    lp_ = Teuchos::rcp( new vector<Real>(l, l+dimension_x() ) );
  }

  /* Set the upper bound to u */
  void setUpper(Real u[]) {
    up_ = Teuchos::rcp( new vector<Real>(u, u+dimension_x() ) );
  }

  /* \brief Call this function in your problem's constructor to turn 
      off bound constraints */
  void noBound() {
    bnd_->deactivate();   
  }
  
  /* \brief Return the objective function. */
  virtual const Teuchos::RCP<OBJ> getObjective() = 0;

  /* \brief Return the equality constraint. Do not overload if you 
            have no equality constraint. */
  virtual const Teuchos::RCP<CON> getEqualityConstraint() {
    return Teuchos::null;
  }

  /* \brief Return the equality constraint vector (used for cloning).
            Returns a Teuchos::null if there is no equality constraint */
  const Teuchos::RCP<V> getEqualityMultiplier() { 
    int n = dimension_ce();
    if( n > 0 ) {
      Teuchos::RCP<vector<Real> > l = Teuchos::rcp( new vector<Real>(n,1.0) );
      return Teuchos::rcp( new SV(l) );
    }
    else {
      return Teuchos::null;
    }
  }

  /* \brief Return the inequality constraint. Do not overload if you 
            have no inequality constraint. */
  virtual const Teuchos::RCP<CON> getInequalityConstraint() {
    return Teuchos::null;
  }

  /* \brief Return the inequality constraint vector (used for cloning). 
            Returns a Teuchos::null if there is no inequality constraint */
  const Teuchos::RCP<V> getInequalityMultiplier() {
    int n = dimension_ci();
    if( n > 0 ) {
      Teuchos::RCP<vector<Real> > l = Teuchos::rcp( new vector<Real>(n,1.0) );
      return Teuchos::rcp( new SV(l) );
    }
    else {
      return Teuchos::null;
    }
  }

  /* \brief Return the bounds on the inequality constraint.
            Returns a Teuchos::null if there is no inequality constraint */
  const Teuchos::RCP<BND> getInequalityBoundConstraint() {
    int n = dimension_ci();
    if( n > 0 ) {
      const Real lval(0), uval(ROL_INF<Real>());
      Teuchos::RCP<V> l = Teuchos::rcp( new SV(
                          Teuchos::rcp( new vector<Real>(n,lval) ) ) );
      Teuchos::RCP<V> u = Teuchos::rcp( new SV(
                          Teuchos::rcp( new vector<Real>(n,uval) ) ) );
      return Teuchos::rcp( new BND(l,u) );
    }
    else {
      return Teuchos::null;
    }
  }

  /* \brief Create vector */
  Teuchos::RCP<V> createOptVector( const Real * const array ) {
    Teuchos::RCP<vector<Real> > x = 
      Teuchos::rcp( new vector<Real>(array,array+dimension_x()) );
    return Teuchos::rcp( new SV( x ) ); 
  }

  // Create an RCP to a std::vector with dimensionality of the optimization space
//  Teuchos::RCP<vector<Real> > createOptVector() {
//    int n = dimension_x();
//    Teuchos::RCP<V> x = Teuchos::rcp( new vector<Real>(n) );
//    return x;
//  } 

  // Default deactivated bound
  Teuchos::RCP<BND> getBoundConstraint() {
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
     lp_ = Teuchos::rcp( new vector<Real>(n,ROL::ROL_NINF<Real>()) );
     up_ = Teuchos::rcp( new vector<Real>(n,ROL::ROL_INF<Real>()) );
     Teuchos::RCP<V> l = Teuchos::rcp( new SV( lp_ ) );
     Teuchos::RCP<V> u = Teuchos::rcp( new SV( up_ ) );
     bnd_ = Teuchos::rcp( new BND( l, u ) ); 
  }

   /* \brief Create the OptimizationProblem from the supplied components and 
             return it */
  Teuchos::RCP<OPT> getOptimizationProblem() {
    Teuchos::RCP<V> x = getInitialGuess()->clone();
    x->set(*getInitialGuess());
    return Teuchos::rcp( new OPT ( getObjective(),
                                   x,
                                   getBoundConstraint(),
                                   getEqualityConstraint(),
                                   getEqualityMultiplier(),
                                   getInequalityConstraint(),
                                   getInequalityMultiplier(),
                                   getInequalityBoundConstraint() ) );
  }

  /* \brief Return the initial guess for the optimization vector */
  virtual const Teuchos::RCP<const V> getInitialGuess() = 0;

  /* \brief Return whether or not the initial guess is feasible */
  virtual bool initialGuessIsFeasible() = 0;

  /* \brief Return the value of the Objective function for the initial guess */
  virtual Real getInitialObjectiveValue() = 0;

  /* \brief Return the set of vectors that solve the nonlinear program if
            they are known */
  virtual Teuchos::RCP<const V> getSolutionSet() { return Teuchos::null; } 
  
  /* \brief Return the value of the objective function for a solution vector. 
            If not known, return infinity */
  virtual Real getSolutionObjectiveValue() { 
     return ROL_INF<Real>(); 
  };
   
   /* \brief If the problem has known solutions, return whether ROL
             has acceptibly solved problem */
  bool foundAcceptableSolution( const V &x, const Real &tolerance=std::sqrt(ROL_EPSILON<Real>()) ) {
    Teuchos::RCP<const PV> sol = Teuchos::rcp_dynamic_cast<const PV>( getSolutionSet() );
    Teuchos::RCP<V> error;
    Teuchos::RCP<const V> xv;
    if ( dimension_ci() > 0 ) {
      xv = Teuchos::dyn_cast<const PV>(x).get(0);
      error = xv->clone();
    }
    else {
      xv = Teuchos::rcp(&x, false);
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
