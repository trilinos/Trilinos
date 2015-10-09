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

#include "ROL_CompositeStepSQP.hpp"
#include "ROL_InteriorPoint.hpp"
#include "ROL_Types.hpp"

namespace ROL {

template <class Real>
class InteriorPointStep : public Step<Real> {

typedef InteriorPointObjective<Real>          IPOBJ;
typedef InteriorPointEqualityConstraint<Real> IPCON;

private:

  Teuchos::RCP<Vector<Real> >           xvec_;
  Teuchos::RCP<Vector<Real> >           gvec_;
  Teuchos::RCP<Vector<Real> >           lvec_;
  Teuchos::RCP<Vector<Real> >           cvec_;

  Teuchos::RCP<StatusTest<Real> >       status_;
  Teuchos::RCP<Step<Real> >             step_;  
  Teuchos::RCP<DefaultAlgorithm<Real> > algo_;
  Teuchos::RCP<Teuchos::ParameterList>  parlist_;

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
    Step<Real>(), step_(Teuchos::null), status_(Teuchos::null) {

    using Teuchos::ParameterList;
    
    ParameterList& iplist  = parlist.sublist("Step").sublist("Interior Point");
    ParameterList& stlist  = parlist.sublist("Status Test");
    ParameterList& sqplist = parlist.sublist("Step").sublist("Composite Step SQP");  


    // Interior Point parameters
    mu_             = iplist.get("Initial Barrier Penalty",1.0);
    eps_            = iplist.get("Minimum Barrier Penalty",1.e-4);
    rho_            = iplist.get("Barrier Penalty Reduction Factor",0.5);
    subproblemIter_ = iplist.get("Subproblem Iteration Limit",10);

 
    // Status test parameters
    gtol_  = stlist.get("Gradient Tolerance", 1.e-8);
    ctol_  = stlist.get("Constraint Tolerance", 1.e-8);
    stol_  = stlist.get("Step Tolerance", 1.e-8);
    maxit_ = stlist.get("Iteration Limit", 100);
  
    
     
    parlist_ = Teuchos::rcp(&parlist, false);

    step_ = Teuchos::rcp(new CompositeStepSQP<Real>(sqplist) );

  }

  /** \brief Initialize step with equality constraint 
   */
  void initialize( Vector<Real> &x, const Vector<Real> &g, Vector<Real> &l, const Vector<Real> &c,
                   Objective<Real> &obj, EqualityConstraint<Real> &con, AlgorithmState<Real> &algo_state ) {

    Teuchos::RCP<StepState<Real> > state = Step<Real>::getState();
    state->descentVec    = x.clone();
    state->gradientVec   = g.clone();
    state->constraintVec = c.clone();

    // Downcast Objective -> InteriorPointObjective
    IPOBJ &ipobj = Teuchos::dyn_cast<IPOBJ>(obj);
    IPCON &ipcon = Teuchos::dyn_cast<IPCON>(con);

    // Set initial penalty
    ipobj.updatePenalty(mu_);

    xvec_ = x.clone();
    gvec_ = g.clone();
    lvec_ = l.clone();
    cvec_ = c.clone();

    algo_state.nfval = 0;
    algo_state.ncval = 0;
    algo_state.ngrad = 0;

    Real zerotol = 0.0;
    obj.update(x,true,algo_state.iter);
    algo_state.value = obj.value(x,zerotol);

    obj.gradient(g,x,zerotol);
    algo_state.gnorm = g.norm();

    con.value(c,x,zerotol);
    algo_state.cnorm = c.norm();

    algo_state.nfval += ipobj.getNumberFunctionEvaluations();
    algo_state.ngval += ipobj.getNumberGradientEvaluations();
    algo_state.ncval += ipcon.getNumberConstraintEvaluations(); 

  }

  void compute( Vector<Real> &s, const Vector<Real> &x, 
                Objective<Real> &obj, EqualityConstraint<Real> &con, 
                AlgorithmState<Real> &algo_state ) {

    // Reset the status test
    status_ = Teuchos::rcp( new StatusTestSQP<Real>(gtol_,ctol_,stol_,maxit_) );

    // Create the algorithm 
    algo_ = Teuchos::rcp( new DefaultAlgorithm<Real>(step_,status_,false) );

    xvec_->set(x);

    //  Run the algorithm
    algo_->run(*xvec_,*gvec_,*lvec_,*cvec_,obj,con,false);

    s.set(*xvec_); s.axpy(-1.0,x);

    // Get number of iterations from the subproblem solve
    subproblemIter_ = (algo_->getState())->iter;
    
  }

  void update( Vector<Real> &x, const Vector<Real> &s, Objective<Real> &obj, 
               EqualityConstraint<Real> &con,  AlgorithmState<Real> &algo_state ) {

    Teuchos::RCP<StepState<Real> > state = Step<Real>::getState();


    x.plus(s);
    algo_state.iter++;

    // Downcast Objective -> InteriorPointObjective
    IPOBJ &ipobj = Teuchos::dyn_cast<IPOBJ>(obj);
    IPCON &ipcon = Teuchos::dyn_cast<IPCON>(con);

    algo_state_->nfval += ipobj.getNumberFunctionEvaluations();
    algo_state_->ngval += ipobj.getNumberGradientEvaluations();
    algo_state_->ncval += ipcon.getNumberConstraintEvaluations();


    // If we can reduce the barrier parameter, do so
    if(mu_ > eps_) {
      mu_ *= rho_;
      ipobj.updatePenalty(mu_);
    }
    
    

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
  std::string printHeader( void ) const {
    std::stringstream hist;
    hist << "  ";
    hist << std::setw(6)  << std::left << "iter";
    hist << std::setw(15) << std::left << "fval";
    hist << std::setw(15) << std::left << "cnorm";
    hist << std::setw(15) << std::left << "gnorm";
    hist << std::setw(15) << std::left << "snorm";
    hist << std::setw(15) << std::left << "penalty";
    hist << std::setw(8) << std::left << "#fval";
    hist << std::setw(8) << std::left << "#grad";
    hist << std::setw(8) << std::left << "#cval";
    hist << std::setw(8) << std::left << "subIter";
    hist << "\n";
    return hist.str();
  }

  /** \brief Print step name.
  */
  std::string printName( void ) const {
    std::stringstream hist;
    hist << "\n" << " Interior Point solver\n";
    return hist.str();
  }

  /** \brief Print iterate status.
  */
  std::string print( AlgorithmState<Real> &algo_state, bool printHeader = false ) const {
    std::stringstream hist;
    hist << std::scientific << std::setprecision(6);
    if ( algo_state.iter == 0 ) {
      hist << "  ";
      hist << std::setw(6)  << std::left << algo_state.iter;
      hist << std::setw(15) << std::left << algo_state.value;
      hist << std::setw(15) << std::left << algo_state.cnorm;
      hist << std::setw(15) << std::left << algo_state.gnorm;     
      hist << std::setw(15) << std::left << algo_state.snorm;     
      hist << std::setw(15) << std::left << mu_;     
      hist << std::setw(8)  << std::left << algo_state.nfval;
      hist << std::setw(8)  << std::left << algo_state.ngrad;
      hist << std::setw(8)  << std::left << algo_state.ncval;
      hist << std::setw(8)  << std::left << subproblemIter_;
      hist << "\n";
    }
    if ( print_header ) {
      hist << printHeader();
    }
    return hist.str(); 
  }





}; // class InteriorPointStep

} // namespace ROL

#endif // ROL_INTERIORPOINTSTEP_H
