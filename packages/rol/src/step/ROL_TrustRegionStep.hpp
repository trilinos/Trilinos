//@HEADER
// ***********************************************************************
//
//                     Rapid Optimization Library
//
// Questions? Contact:    Drew Kouri (dpkouri@sandia.gov)
//                      Denis Ridzal (dridzal@sandia.gov)
//
// ***********************************************************************
//@HEADER

#ifndef ROL_TRUSTREGIONSTEP_H
#define ROL_TRUSTREGIONSTEP_H

#include "ROL_Step.hpp"
#include "ROL_Secant.hpp"
#include "ROL_TrustRegion.hpp"
#include <sstream>
#include <iomanip>

/** \class ROL::TrustRegionStep
    \brief Provides the interface to compute optimization steps
           with trust regions.
*/


namespace ROL {

template <class Real>
class TrustRegionStep : public Step<Real> {
private:

  TrustRegionType TRtype_;
  TrustRegionStepType TRStype_;

  Teuchos::RCP<TrustRegion<Real> > trustRegion_;
  Teuchos::RCP<Secant<Real> > secant_;

  int CGiter_;
  int CGflag_;

  int TRflag_;

  int TR_nfval_;
  int TR_ngrad_;

  std::string step_;
  std::string TR_;

  Real del_;

public:

  virtual ~TrustRegionStep() {}

  TrustRegionStep( TrustRegionType TRtype = TrustRegionType_TruncatedCG, 
                   TrustRegionStepType TRStype = TrustRegionStep_NewtonKrylov,
                   int maxit = 20, Real tol1 = 1.e-4, Real tol2 = 1.e-2,
                   Real del = 100.0, Real delmin = 1.e-8, Real delmax = 5000.0,
                   Real eta0 = 0.05, Real eta1 = 0.05, Real eta2 = 0.9,
                   Real gamma0 = 0.0625, Real gamma1 = 0.25, Real gamma2 = 2.5, Real TRsafe = 1.0,
                   SecantType Stype = Secant_lBFGS, int L = 10, int BBtype = 1 ) 
    : TRtype_(TRtype), TRStype_(TRStype), del_(del) {
     
    trustRegion_ = Teuchos::rcp( new TrustRegion<Real>(TRtype,TRStype,maxit,tol1,tol2,delmin,delmax,eta0,eta1,eta2,
                                                       gamma0,gamma1,gamma2,TRsafe) );
    Teuchos::RCP<Secant<Real> > secant;
    if ( TRStype == TrustRegionStep_Newton ) {
      step_   = "Newton's Method";
    }
    else if ( TRStype == TrustRegionStep_NewtonKrylov ) {
      step_   = "Newton-CG";
    }
    else if ( TRStype == TrustRegionStep_NewtonKrylovSecantPreconditioning ) {
      if ( Stype == Secant_lBFGS ) {
        secant_ = Teuchos::rcp( new lBFGS<Real>(L) );
        step_   = "Newton-CG with Limited-Memory BFGS Preconditioning";
      }
      else if ( Stype == Secant_lDFP ) {
        secant_ = Teuchos::rcp( new lDFP<Real>(L) );
        step_   = "Newton-CG with Limited-Memory DFP Preconditioning";
      }
      else if ( Stype == Secant_BarzilaiBorwein ) {
        secant_ = Teuchos::rcp( new BarzilaiBorwein<Real>(BBtype) );
        step_   = "Newton-CG with Barzilai-Borwein Preconditioning";
      }
    }
    else if ( TRStype == TrustRegionStep_Secant ) {
      if ( Stype == Secant_lBFGS ) {
        secant_ = Teuchos::rcp( new lBFGS<Real>(L) );
        step_   = "Limited-Memory BFGS";
      }
      else if ( Stype == Secant_lDFP ) {
        secant_ = Teuchos::rcp( new lDFP<Real>(L) );
        step_   = "Limited-Memory DFP";
      }
      else if ( Stype == Secant_BarzilaiBorwein ) {
        secant_ = Teuchos::rcp( new BarzilaiBorwein<Real>(BBtype) );
        step_   = "Barzilai-Borwein";
      }
    }

    if ( TRtype == TrustRegionType_Dogleg ) {
      TR_ = "Dogleg";
    }
    else if ( TRtype == TrustRegionType_DoubleDogleg ) {
      TR_ = "Double Dogleg";
    }
    else if ( TRtype == TrustRegionType_TruncatedCG ) {
      TR_ = "Truncated CG";
    }
  }

  /** \brief Compute step.
  */
  void compute( Vector<Real> &s, const Vector<Real> &x, Objective<Real> &obj, AlgorithmState<Real> &algo_state ) {
    CGflag_ = 0;
    CGiter_ = 0;
    trustRegion_->run(s,algo_state.snorm,del_,CGflag_,CGiter_,
                      x,*(Step<Real>::state_->gradientVec),algo_state.gnorm,obj,secant_);
  }

  /** \brief Update step, if successful.
  */
  void update( Vector<Real> &x, const Vector<Real> &s, Objective<Real> &obj, AlgorithmState<Real> &algo_state ) {
    TRflag_   = 0;
    TR_nfval_ = 0;
    TR_ngrad_ = 0;
    Real fold = algo_state.value;
    Real fnew = 0.0;
    trustRegion_->update(x,fnew,del_,TR_nfval_,TR_ngrad_,TRflag_,
                         s,algo_state.snorm,fold,*(Step<Real>::state_->gradientVec),obj,secant_);
    algo_state.value = fnew;
    algo_state.nfval += TR_nfval_;
    algo_state.ngrad += TR_ngrad_;

    // Compute new gradient
    Teuchos::RCP<Vector<Real> > gp;
    if ( secant_ != Teuchos::null ) {
      gp = x.clone();
      gp->set(*(Step<Real>::state_->gradientVec));
    }
    obj.gradient(*(Step<Real>::state_->gradientVec),x);
    algo_state.ngrad++;
    
    // Update Secant Information
    if ( secant_ != Teuchos::null ) {
      secant_->update(*(Step<Real>::state_->gradientVec),*gp,s,algo_state.snorm,algo_state.iter+1);
    }
  
    // Update algorithm state
    (algo_state.iterateVec)->set(x);
    algo_state.gnorm = (Step<Real>::state_->gradientVec)->norm();
    algo_state.iter++;
  }

  /** \brief Print iterate status.
  */
  std::string print( AlgorithmState<Real> & algo_state, bool printHeader = false ) const  {
    std::stringstream hist;
    if ( algo_state.iter == 0 ) {
      hist << "\n" << step_ << " with " << TR_ << " Trust-Region Solver\n"; 
      hist << "  ";
      hist << std::setw(6)  << std::left << "iter";  
      hist << std::setw(15) << std::left << "value";
      hist << std::setw(15) << std::left << "gnorm"; 
      hist << std::setw(15) << std::left << "snorm";
      hist << std::setw(15) << std::left << "delta";
      hist << std::setw(10) << std::left << "#fval";
      hist << std::setw(10) << std::left << "#grad";
      hist << std::setw(10) << std::left << "tr_flag";
      if (    TRStype_ == TrustRegionStep_NewtonKrylov 
           || TRStype_ == TrustRegionStep_NewtonKrylovSecantPreconditioning ) {
        hist << std::setw(10) << std::left << "iterCG";
        hist << std::setw(10) << std::left << "flagCG";
      }
      hist << "\n";
      hist << std::setfill('-') << std::setw(120) << "-" << "\n";
      hist << std::setfill(' ') << "  ";
      hist << std::setw(6)  << std::left << algo_state.iter;
      hist << std::setw(15) << std::left << algo_state.value;
      hist << std::setw(15) << std::left << algo_state.gnorm;
      hist << "\n";
    }
    else {
      hist << "  "; 
      hist << std::setw(6)  << std::left << algo_state.iter;  
      hist << std::setw(15) << std::left << algo_state.value; 
      hist << std::setw(15) << std::left << algo_state.gnorm; 
      hist << std::setw(15) << std::left << algo_state.snorm; 
      hist << std::setw(15) << std::left << del_; 
      hist << std::setw(10) << std::left << algo_state.nfval;              
      hist << std::setw(10) << std::left << algo_state.ngrad;              
      hist << std::setw(10) << std::left << TRflag_;              
      if (    TRStype_ == TrustRegionStep_NewtonKrylov 
           || TRStype_ == TrustRegionStep_NewtonKrylovSecantPreconditioning ) {
        hist << std::setw(10) << std::left << CGiter_;
        hist << std::setw(10) << std::left << CGflag_;
      }
      hist << "\n";
    }
    return hist.str();
  }

}; // class Step

} // namespace ROL

#endif
