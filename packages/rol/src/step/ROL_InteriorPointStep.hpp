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
#include "ROL_ObjectiveFromBoundConstraint.hpp"
#include "ROL_Types.hpp"


namespace ROL {

template <class Real>
class InteriorPointStep : public Step<Real> {

typedef InteriorPoint::PenalizedObjective<Real>   IPOBJ;
typedef InteriorPoint::CompositeConstraint<Real>  IPCON;

typedef PartitionedVector<Real> PV;
typedef typename PV::size_type  size_type; 

const static size_type OPT   = 0;
const static size_type SLACK = 1;

private:

  Teuchos::RCP<StatusTest<Real> >       status_;
  Teuchos::RCP<Step<Real> >             step_;  
  Teuchos::RCP<IPOBJ>                   ipobj_;
  Teuchos::RCP<IPCON>                   ipcon_;
  Teuchos::RCP<Algorithm<Real> >        algo_;
  Teuchos::RCP<Teuchos::ParameterList>  parlist_;

  // Storage
  Teuchos::RCP<PV> x_;
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

  int verbosity_;       // Adjust level of detail in printing step information

public:
 
  using Step<Real>::initialize;
  using Step<Real>::compute;
  using Step<Real>::update;

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
    
    verbosity_ = parlist.sublist("General").get("Print Verbosity",0);

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
 
    parlist_ = Teuchos::rcp(&parlist, false);

  }

  /** \brief Initialize step with equality constraint 
   */
  void initialize( Vector<Real> &x, const Vector<Real> &g, 
                   Vector<Real> &l, const Vector<Real> &c,
                   Objective<Real> &obj, EqualityConstraint<Real> &con, 
                   AlgorithmState<Real> &algo_state ) {

    Teuchos::RCP<StepState<Real> > state = Step<Real>::getState();
    state->descentVec    = x.clone();
    state->gradientVec   = g.clone();
    state->constraintVec = c.clone();

    // Initialize storage
    x_ = Teuchos::rcp_static_cast<PV>(x.clone());
    g_ = g.clone();
    l_ = l.clone();
    c_ = c.clone();

    x_->set(x);

    ipobj_ = Teuchos::rcp(&Teuchos::dyn_cast<IPOBJ>(obj),false);
    ipcon_ = Teuchos::rcp(&Teuchos::dyn_cast<IPCON>(con),false);

    // Set initial penalty
    ipobj_->updatePenalty(mu_);

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

    algo_state.nfval += ipobj_->getNumberFunctionEvaluations();
    algo_state.ngrad += ipobj_->getNumberGradientEvaluations();
    algo_state.ncval += ipcon_->getNumberConstraintEvaluations(); 

  }


  
  void initialize( Vector<Real> &x, const Vector<Real> &g, Vector<Real> &l, const Vector<Real> &c,
                   Objective<Real> &obj, EqualityConstraint<Real> &con, BoundConstraint<Real> &bnd, 
                   AlgorithmState<Real> &algo_state ) {
    initialize(x,g,l,c,obj,con,algo_state);
  }




  /** \brief Compute step (equality constraints).
  */
  void compute( Vector<Real> &s, const Vector<Real> &x, const Vector<Real> &l,
                Objective<Real> &obj, EqualityConstraint<Real> &con, 
                AlgorithmState<Real> &algo_state ) {

    // Create the algorithm 
    algo_ = Teuchos::rcp( new Algorithm<Real>("Composite Step",*parlist_,false) );

    x_->set(x);

    //  Run the algorithm
    algo_->run(*x_,*g_,*l_,*c_,*ipobj_,*ipcon_,false);

    s.set(*x_); s.axpy(-1.0,x);

    // Get number of iterations from the subproblem solve
    subproblemIter_ = (algo_->getState())->iter;
    
  }

  virtual void compute( Vector<Real> &s, const Vector<Real> &x, const Vector<Real> &l,
                        Objective<Real> &obj, EqualityConstraint<Real> &con, 
                        BoundConstraint<Real> &bnd,
                        AlgorithmState<Real> &algo_state ) {
    compute(s,x,l,obj,con,algo_state); 
  }



  /** \brief Update step, if successful (equality constraints).
  */
  void update( Vector<Real> &x, Vector<Real> &l, const Vector<Real> &s, Objective<Real> &obj, 
               EqualityConstraint<Real> &con,  AlgorithmState<Real> &algo_state ) {

    // If we can reduce the barrier parameter, do so
    if(mu_ > eps_) {
      mu_ *= rho_;
      ipobj_->updatePenalty(mu_);
    }

    Teuchos::RCP<StepState<Real> > state = Step<Real>::getState();
 
    // Update optimization vector
    x.plus(s);

    algo_state.iterateVec->set(x);
    state->descentVec->set(s);
    algo_state.snorm = s.norm();
    algo_state.iter++;

    Real zerotol = 0.0;

    algo_state.value = ipobj_->value(x,zerotol);
    algo_state.value = ipobj_->getObjectiveValue();

    ipcon_->value(*c_,x,zerotol);
    state->constraintVec->set(*c_);

    ipobj_->gradient(*g_,x,zerotol);
    state->gradientVec->set(*g_);

    ipcon_->applyAdjointJacobian(*g_,*l_,x,zerotol);
    state->gradientVec->plus(*g_);    

    x_->set(x);
    x_->axpy(-1.0,state->gradientVec->dual());

    Elementwise::ThresholdUpper<Real> threshold(0.0);

    //PartitionedVector<Real> &xpv = Teuchos::dyn_cast<PartitionedVector<Real> >(*x_);

    Teuchos::RCP<Vector<Real> > slack = x_->get(SLACK);
   
    slack->applyUnary(threshold);

    x_->axpy(-1.0,x);

    algo_state.gnorm = x_->norm();
    algo_state.cnorm = state->constraintVec->norm();
    algo_state.snorm = s.norm();

    algo_state.nfval += ipobj_->getNumberFunctionEvaluations();
    algo_state.ngrad += ipobj_->getNumberGradientEvaluations();
    algo_state.ncval += ipcon_->getNumberConstraintEvaluations();
    
  }

  void update( Vector<Real> &x, Vector<Real> &l, const Vector<Real> &s,
               Objective<Real> &obj, EqualityConstraint<Real> &con,
               BoundConstraint<Real> &bnd,
               AlgorithmState<Real> &algo_state ) {
    update(x,l,s,obj,con,algo_state); 
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

    if( verbosity_ > 0 ) {

      hist << std::string(116,'-') << "\n";
      hist << "Interior Point status output definitions\n\n";
   
      hist << "  IPiter  - Number of interior point steps taken\n";
      hist << "  CSiter  - Number of Composite Steps taken in each subproblem\n";
      hist << "  penalty - Penalty parameter multiplying the barrier objective\n";
      hist << "  fval    - Number of objective evaluations\n";
      hist << "  cnorm   - Norm of the composite constraint\n";
      hist << "  gLnorm  - Norm of the Lagrangian's gradient\n";
      hist << "  snorm   - Norm of step (update to optimzation and slack vector)\n";
      hist << "  #fval   - Number of objective function evaluations\n";
      hist << "  #grad   - Number of gradient evaluations\n";
      hist << "  #cval   - Number of composite constraint evaluations\n"; 
      hist << std::string(116,'-') << "\n";
      
     
    }

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
