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

#ifndef ROL_COMPOSITESTEPSQP_H
#define ROL_COMPOSITESTEPSQP_H

#include "ROL_Types.hpp"
#include "ROL_Step.hpp"
#include <sstream>
#include <iomanip>

/** \class ROL::CompositeStepSQP
    \brief Implements the computation of optimization steps
           with composite-step trust-region SQP methods.
*/


namespace ROL {

template <class Real>
class CompositeStepSQP : public Step<Real> {
private:

  int TRflag_;
  int CGflag_;
  int CGiter_;

  Real lmhtol_;
  Real qntol_;
  Real pgtol_;
  Real projtol_;
  Real tangtol_;
  Real tntmax_;

public:

  virtual ~CompositeStepSQP() {}

  CompositeStepSQP( Teuchos::ParameterList & parlist ) : Step<Real>() {
    //Teuchos::RCP<StepState<Real> > step_state = Step<Real>::getState();
    TRflag_ = 0;
    CGflag_ = 0;
    CGiter_ = 0;
    Real nominal_tol = parlist.get("Nominal SQP Optimality Solver Tolerance", 1e-3);
    lmhtol_  = nominal_tol;
    qntol_   = nominal_tol;
    pgtol_   = nominal_tol;       
    projtol_ = nominal_tol;     
    tangtol_ = nominal_tol;     
  }

  /** \brief Initialize step.
  */
  void initialize( Vector<Real> &x, Vector<Real> &l,
                   Objective<Real> &obj, EqualityConstraint<Real> &con, 
                   AlgorithmState<Real> &algo_state ) {
    //Teuchos::RCP<StepState<Real> > step_state = Step<Real>::getState();

    algo_state.nfval = 0;
    algo_state.ncval = 0;
    algo_state.ngrad = 0;

    Real zerotol = 0.0;

    // Update objective and constraint.
    obj.update(x,true,algo_state.iter);    
    algo_state.value = obj.value(x, zerotol); 
    algo_state.nfval++;
    con.update(x,true,algo_state.iter);
    Teuchos::RCP<Vector<Real> > c = l.clone();
    con.value(*c, x, zerotol);    
    algo_state.cnorm = c->norm(); 
    algo_state.ncval++;
    Teuchos::RCP<Vector<Real> > g = x.clone();
    obj.gradient(*g, x, zerotol);

    // Compute gradient of Lagrangian at new multiplier guess.
    computeLagrangeMultiplier(l, x, *g, con);
    Teuchos::RCP<Vector<Real> > ajl = x.clone();
    con.applyAdjointJacobian(*ajl, l, x, zerotol);
    Teuchos::RCP<Vector<Real> > gl = x.clone();
    gl->set(*g); gl->plus(*ajl);
    algo_state.ngrad++;
    algo_state.gnorm = gl->norm();
  }

  /** \brief Compute step.
  */
  void compute( Vector<Real> &s, const Vector<Real> &x, const Vector<Real> &l,
                Objective<Real> &obj, EqualityConstraint<Real> &con, 
                AlgorithmState<Real> &algo_state ) {
    //Teuchos::RCP<StepState<Real> > step_state = Step<Real>::getState();
  }

  /** \brief Update step, if successful.
  */
  void update( Vector<Real> &x, Vector<Real> &l, const Vector<Real> &s,
               Objective<Real> &obj, EqualityConstraint<Real> &con, 
               AlgorithmState<Real> &algo_state ) {
    //Teuchos::RCP<StepState<Real> > state = Step<Real>::getState();

    //Real tol = std::sqrt(ROL_EPSILON);

    //Real eps = 0.0;
  
    // Update algorithm state
    //(algo_state.iterateVec)->set(x);
  }

  /** \brief Compute step for bound constraints; here only to satisfy the
             interface requirements, does nothing, needs refactoring.
  */
  void compute( Vector<Real> &s, const Vector<Real> &x, Objective<Real> &obj,
                        BoundConstraint<Real> &con,
                        AlgorithmState<Real> &algo_state ) {}

  /** \brief Update step, for bound constraints; here only to satisfy the 
             interface requirements, does nothing, needs refactoring.
  */
  void update( Vector<Real> &x, const Vector<Real> &s, Objective<Real> &obj,
                       BoundConstraint<Real> &con,
                       AlgorithmState<Real> &algo_state ) {} 

  /** \brief Print iterate header.
  */
  std::string printHeader( void ) const  {
    std::stringstream hist;
    hist << "  ";
    hist << std::setw(6)  << std::left << "iter";
    hist << std::setw(15) << std::left << "fval";
    hist << std::setw(15) << std::left << "cnorm";
    hist << std::setw(15) << std::left << "gLnorm";
    hist << std::setw(15) << std::left << "snorm";
    hist << std::setw(15) << std::left << "delta";
    hist << std::setw(10) << std::left << "#fval";
    hist << std::setw(10) << std::left << "#grad";
    hist << std::setw(10) << std::left << "flagTR";
    hist << std::setw(10) << std::left << "iterCG";
    hist << std::setw(10) << std::left << "flagCG";
    hist << "\n";
    return hist.str();
  }

  std::string printName( void ) const {
    std::stringstream hist;
    hist << "\n" << " Composite-step trust-region SQP solver";
    hist << "\n";
    return hist.str();
  }

  /** \brief Print iterate status.
  */
  std::string print( AlgorithmState<Real> & algo_state, bool pHeader = false ) const  {
    //const Teuchos::RCP<const StepState<Real> >& step_state = Step<Real>::getStepState();

    std::stringstream hist;
    hist << std::scientific << std::setprecision(6);
    if ( algo_state.iter == 0 ) {
      hist << printName();
    }
    if ( pHeader ) {
      hist << printHeader();
    }
    if ( algo_state.iter == 0 ) {
      hist << "  ";
      hist << std::setw(6)  << std::left << algo_state.iter;
      hist << std::setw(15) << std::left << algo_state.value;
      hist << std::setw(15) << std::left << algo_state.cnorm;
      hist << std::setw(15) << std::left << algo_state.gnorm;
      hist << "\n";
    }
    else {
      hist << "  "; 
      hist << std::setw(6)  << std::left << algo_state.iter;  
      hist << std::setw(15) << std::left << algo_state.value; 
      hist << std::setw(15) << std::left << algo_state.cnorm; 
      hist << std::setw(15) << std::left << algo_state.gnorm; 
      hist << std::setw(15) << std::left << algo_state.snorm; 
      hist << std::setw(10) << std::left << algo_state.nfval;              
      hist << std::setw(10) << std::left << algo_state.ngrad;              
      hist << std::setw(10) << std::left << TRflag_;              
      hist << std::setw(10) << std::left << CGiter_;
      hist << std::setw(10) << std::left << CGflag_;
      hist << "\n";
    }
    return hist.str();
  }

  /** \brief Compute Lagrange multipliers by solving the least-squares
             problem minimizing the gradient of the Lagrangian, via the
             augmented system formulation.

             @param[out]      l   is the updated Lagrange multiplier; a dual constraint-space vector
             @param[in]       x   is the current iterate; an optimization-space vector
             @param[in]       gf  is the gradient of the objective function; a dual optimization-space vector
             @param[in]       con is the equality constraint object

             On return ... fill the blanks.
  */
  void computeLagrangeMultiplier(Vector<Real> &l, const Vector<Real> &x, const Vector<Real> &gf, EqualityConstraint<Real> &con) {

    Real zerotol = 0.0;

    /* Apply adjoint of constraint Jacobian to current multiplier. */
    Teuchos::RCP<Vector<Real> > ajl = x.clone();
    con.applyAdjointJacobian(*ajl, l, x, zerotol);

    /* Form right-hand side of the augmented system. */
    Teuchos::RCP<Vector<Real> > b1 = x.clone();
    Teuchos::RCP<Vector<Real> > b2 = l.clone();
    // b1 is the negative gradient of the Lagrangian
    b1->set(gf); b1->plus(*ajl); b1->scale(-1.0);
    // b2 is zero
    b2->zero();

    /* Declare left-hand side of augmented system. */
    Teuchos::RCP<Vector<Real> > v1 = x.clone();
    Teuchos::RCP<Vector<Real> > v2 = l.clone();
    
    /* Compute linear solver tolerance. */
    Real b1norm  = b1->norm();
    Real tol = lmhtol_*b1norm;

    /* Solve augmented system. */
    con.solveAugmentedSystem(*v1, *v2, *b1, *b2, x, tol);

    /* Return updated Lagrange multiplier. */
    // v2 is the multiplier update
    l.plus(*v2);
  }

  /** \brief Compute quasi-normal step by minimizing the norm of
             the linearized constraint.

             Compute an approximate solution of the problem
             \f[
               \begin{array}{rl}
                 \min_{n} & \|c'(x_k)n + c(x_k)\|^2_{\mathcal{X}} \\
                 \mbox{subject to} & \|n\|_{\mathcal{X}} \le \delta
               \end{array}
             \f]
             The approximate solution is computed using the Powell dogleg
             method. The dogleg path is computed using the Cauchy point and
             a full Newton step.  Its intersection with the trust-region
             constraint gives the quasi-normal step.

             @param[out]      n     is the quasi-normal step; an optimization-space vector
             @param[in]       c     is the value of equality constraints; a constraint-space vector
             @param[in]       x     is the current iterate; an optimization-space vector
             @param[in]       delta is the trust-region radius for the quasi-normal step
             @param[in]       con   is the equality constraint object

  */
  void computeQuasinormalStep(Vector<Real> &n, const Vector<Real> &c, const Vector<Real> &x, Real delta, EqualityConstraint<Real> &con) {
    
  }

}; // class CompositeStepSQP

} // namespace ROL

#endif
