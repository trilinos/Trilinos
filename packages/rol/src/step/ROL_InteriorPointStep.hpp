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

#ifndef ROL_INTERIORPOINTSTEP_H
#define ROL_INTERIORPOINTSTEP_H

#include "ROL_CompositeStep.hpp"
#include "ROL_ConstraintStatusTest.hpp"
#include "ROL_InteriorPoint.hpp"
#include "ROL_Types.hpp"

namespace ROL {

template <class Real>
class InteriorPointStep : public Step<Real> {

typedef InteriorPoint::PenalizedObjective<Real>   IPOBJ;
typedef InteriorPoint::CompositeConstraint<Real>  IPCON;

private:

  Teuchos::RCP<StatusTest<Real> >          status_;
  Teuchos::RCP<Step<Real> >                step_;  
  Teuchos::RCP<Objective<Real> >           ipobj_;
  Teuchos::RCP<EqualityConstraint<Real> >  ipcon_;
  Teuchos::RCP<Algorithm<Real> >           algo_;
  Teuchos::RCP<Teuchos::ParameterList>     parlist_;

  // Storage
  Teuchos::RCP<Vector<Real> > x_;
  Teuchos::RCP<Vector<Real> > g_;
  Teuchos::RCP<Vector<Real> > l_;
  Teuchos::RCP<Vector<Real> > c_;

  Real mu_;      // Barrier parameter
  Real eps_;     // Minimal value of barrier parameter
  Real rho_;     // Barrier parameter reduction factor
  int  maxit_;   // Maximum number of interior point subproblem solves

  // For the subproblem
  Real gtol_;           // Status test gradient tolerance
  Real ctol_;           // Status test constraint tolerance
  Real stol_;           // Status test step tolerance
  int subproblemIter_;  // Status test maximum number of iterations

public:
 
  ~InteriorPointStep() {}

  InteriorPointStep(Teuchos::ParameterList &parlist) :
    Step<Real>(), 
    status_(Teuchos::null), 
    step_(Teuchos::null),
    ipobj_(Teuchos::null),
    ipcon_(Teuchos::null),
    algo_(Teuchos::null), 
    x_(Teuchos::null),
    g_(Teuchos::null),
    l_(Teuchos::null),
    c_(Teuchos::null) {

    using Teuchos::ParameterList;
    
    // List of general Interior Point parameters
    ParameterList& iplist  = parlist.sublist("Step").sublist("Interior Point");

    mu_             = iplist.get("Initial Barrier Penalty",1.0);
    eps_            = iplist.get("Minimum Barrier Penalty",1.e-4);
    rho_            = iplist.get("Barrier Penalty Reduction Factor",0.5);
    subproblemIter_ = iplist.get("Subproblem Iteration Limit",10);


    // List of Status Test parameters
    ParameterList& stlist  = parlist.sublist("Status Test");

    gtol_  = stlist.get("Gradient Tolerance", 1.e-8);
    ctol_  = stlist.get("Constraint Tolerance", 1.e-8);
    stol_  = stlist.get("Step Tolerance", 1.e-8);
    maxit_ = stlist.get("Iteration Limit", 100);
 
    // List of Composite Step SQP parameters
    ParameterList& cslist  = parlist.sublist("Step").sublist("Composite Step");  
    
     
    parlist_ = Teuchos::rcp(&parlist, false);

    // Create a Composite Step SQP subproblem solver
    step_ = Teuchos::rcp(new CompositeStep<Real>(cslist) );

  }

  /** \brief Initialize step with equality constraint 
   */
  virtual void initialize( Vector<Real> &x, const Vector<Real> &g, 
                           Vector<Real> &l, const Vector<Real> &c,
                           Objective<Real> &obj, EqualityConstraint<Real> &con, 
                           AlgorithmState<Real> &algo_state ) {

//    std::cout << "InteriorPointStep::initialize()" << std::endl;

    Teuchos::RCP<StepState<Real> > state = Step<Real>::getState();
    state->descentVec    = x.clone();
    state->gradientVec   = g.clone();
    state->constraintVec = c.clone();

    // Initialize storage
    x_ = x.clone();
    g_ = g.clone();
    l_ = l.clone();
    c_ = c.clone();

    x_->set(x);

    // Downcast Objective -> InteriorPointObjective
    IPOBJ &ipobj = Teuchos::dyn_cast<IPOBJ>(obj);
    IPCON &ipcon = Teuchos::dyn_cast<IPCON>(con);

    // Set initial penalty
    ipobj.updatePenalty(mu_);

    algo_state.nfval = 0;
    algo_state.ncval = 0;
    algo_state.ngrad = 0;

    Real zerotol = 0.0;
    obj.update(*x_,true,algo_state.iter);
    algo_state.value = obj.value(*x_,zerotol);

    obj.gradient(*g_,*x_,zerotol);
    algo_state.gnorm = g_->norm();

    con.value(*c_,*x_,zerotol);
    algo_state.cnorm = c_->norm();

    algo_state.nfval += ipobj.getNumberFunctionEvaluations();
    algo_state.ngrad += ipobj.getNumberGradientEvaluations();
    algo_state.ncval += ipcon.getNumberConstraintEvaluations(); 

  }


  /** \brief Compute step (equality constraints).
  */
  void compute( Vector<Real> &s, const Vector<Real> &x, const Vector<Real> &l,
                Objective<Real> &obj, EqualityConstraint<Real> &con, 
                AlgorithmState<Real> &algo_state ) {

//    std::cout << "InteriorPointStep::compute()" << std::endl;

    // Reset the status test
    status_ = Teuchos::rcp( new ConstraintStatusTest<Real>(gtol_,ctol_,stol_,maxit_) );

    // Create the algorithm 
    algo_ = Teuchos::rcp( new Algorithm<Real>(step_,status_,false) );

    x_->set(x);

    //  Run the algorithm
    algo_->run(*x_,*g_,*l_,*c_,obj,con,false);

    s.set(*x_); s.axpy(-1.0,x);

    // Get number of iterations from the subproblem solve
    subproblemIter_ = (algo_->getState())->iter;
    
  }


  /** \brief Update step, if successful (equality constraints).
  */
  void update( Vector<Real> &x, Vector<Real> &l, const Vector<Real> &s, Objective<Real> &obj, 
               EqualityConstraint<Real> &con,  AlgorithmState<Real> &algo_state ) {


//    std::cout << "InteriorPointStep::update()" << std::endl;

    Teuchos::RCP<StepState<Real> > state = Step<Real>::getState();
 
    // Update optimization vector
    x.plus(s);

    algo_state.iterateVec->set(x);
    state->descentVec->set(s);
    algo_state.snorm = s.norm();
    algo_state.iter++;

    // Downcast Objective -> InteriorPointObjective
    IPOBJ &ipobj = Teuchos::dyn_cast<IPOBJ>(obj);
    IPCON &ipcon = Teuchos::dyn_cast<IPCON>(con);

    Real zerotol = 0.0;


    algo_state.value = obj.value(x,zerotol);
    obj.gradient(*g_,x,zerotol);
    con.value(*c_,x,zerotol);

    algo_state.gnorm = g_->norm();
    algo_state.cnorm = c_->norm();
    algo_state.snorm = s.norm();

    algo_state.nfval += ipobj.getNumberFunctionEvaluations();
    algo_state.ngrad += ipobj.getNumberGradientEvaluations();
    algo_state.ncval += ipcon.getNumberConstraintEvaluations();

    // If we can reduce the barrier parameter, do so
    if(mu_ > eps_) {
      mu_ *= rho_;
      ipobj.updatePenalty(mu_);
    }
    
  }

  /** \brief Compute step for bound constraints; here only to satisfy the
             interface requirements, does nothing, needs refactoring.
  */
  virtual void compute( Vector<Real> &s, const Vector<Real> &x, Objective<Real> &obj, 
                        BoundConstraint<Real> &con, 
                        AlgorithmState<Real> &algo_state ) {}

  /** \brief Update step, for bound constraints; here only to satisfy the
             interface requirements, does nothing, needs refactoring.
  */
  virtual void update( Vector<Real> &x, const Vector<Real> &s, Objective<Real> &obj, 
                       BoundConstraint<Real> &con,
                       AlgorithmState<Real> &algo_state ) {}

  /** \brief Print iterate header.
  */
  std::string printHeader( void ) const {
    std::stringstream hist;
    hist << "  ";
    hist << std::setw(9)  << std::left  << "IPiter";
    hist << std::setw(9)  << std::left  << "CSiter";
    hist << std::setw(15) << std::left  << "penalty";
    hist << std::setw(15) << std::left  << "fval";
    hist << std::setw(15) << std::left  << "cnorm";
    hist << std::setw(15) << std::left  << "gLnorm";
    hist << std::setw(15) << std::left  << "snorm";
    hist << std::setw(8)  << std::left  << "#fval";
    hist << std::setw(8)  << std::left  << "#grad";
    hist << std::setw(8)  << std::left  << "#cval";

    hist << "\n";
    return hist.str();
  }

  /** \brief Print step name.
  */
  std::string printName( void ) const {
    std::stringstream hist;
    hist << "\n" << "Composite Step Interior Point Solver\n";
    return hist.str();
  }

  /** \brief Print iterate status.
  */
  std::string print( AlgorithmState<Real> &algo_state, bool pHeader = false ) const {
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
      hist << std::setw(9)  << std::left << algo_state.iter;
      hist << std::setw(9)  << std::left << subproblemIter_;
      hist << std::setw(15) << std::left << mu_;
      hist << std::setw(15) << std::left << algo_state.value;
      hist << std::setw(15) << std::left << algo_state.cnorm;
      hist << std::setw(15) << std::left << algo_state.gnorm;
      hist << "\n";
    }
    else {
      hist << "  ";
      hist << std::setw(9)  << std::left << algo_state.iter;
      hist << std::setw(9)  << std::left << subproblemIter_;
      hist << std::setw(15) << std::left << mu_;
      hist << std::setw(15) << std::left << algo_state.value;
      hist << std::setw(15) << std::left << algo_state.cnorm;
      hist << std::setw(15) << std::left << algo_state.gnorm;
      hist << std::setw(15) << std::left << algo_state.snorm;
//      hist << std::scientific << std::setprecision(6);
      hist << std::setw(8) << std::left << algo_state.nfval;
      hist << std::setw(8) << std::left << algo_state.ngrad;
      hist << std::setw(8) << std::left << algo_state.ncval;
      hist << "\n";
    }

    return hist.str(); 
  }





}; // class InteriorPointStep

} // namespace ROL

#endif // ROL_INTERIORPOINTSTEP_H
