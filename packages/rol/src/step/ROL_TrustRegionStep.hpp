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

#include "ROL_Types.hpp"
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

  Teuchos::RCP<TrustRegion<Real> > trustRegion_;
  Teuchos::RCP<Secant<Real> >      secant_;

  int CGiter_;
  int CGflag_;

  int TRflag_;

  int TR_nfval_;
  int TR_ngrad_;

  ETrustRegion etr_;
  EDescent     edesc_;
  ESecant      esec_;

  std::string ssec_;

  Real del_;

public:

  virtual ~TrustRegionStep() {}

  TrustRegionStep( ETrustRegion etr   = TRUSTREGION_TRUNCATEDCG, 
                   EDescent     edesc = DESCENT_NEWTONKRYLOV,
                   int maxit = 20, Real tol1 = 1.e-4, Real tol2 = 1.e-2,
                   Real del = -1.0, Real delmin = 1.e-8, Real delmax = 5000.0,
                   Real eta0 = 0.05, Real eta1 = 0.05, Real eta2 = 0.9,
                   Real gamma0 = 0.0625, Real gamma1 = 0.25, Real gamma2 = 2.5, Real TRsafe = 1.0,
                   ESecant esec = SECANT_LBFGS, int L = 10, int BBtype = 1 ) 
    : etr_(etr), edesc_(edesc), esec_(esec), del_(del) {
     
    trustRegion_ = Teuchos::rcp( new TrustRegion<Real>(etr_,edesc_,maxit,tol1,tol2,delmin,delmax,eta0,eta1,eta2,
                                                       gamma0,gamma1,gamma2,TRsafe) );
    Teuchos::RCP<Secant<Real> > secant;
    if ( edesc_ == DESCENT_SECANTPRECOND || edesc_ == DESCENT_SECANT ) {
      if ( esec_ == SECANT_LBFGS ) {
        secant_ = Teuchos::rcp( new lBFGS<Real>(L) );
      }
      else if ( esec_ == SECANT_LDFP ) {
        secant_ = Teuchos::rcp( new lDFP<Real>(L) );
      }
      else if ( esec_ == SECANT_LSR1 ) {
        secant_ = Teuchos::rcp( new lSR1<Real>(L) );
      }
      else if ( esec_ == SECANT_BARZILAIBORWEIN ) {
        secant_ = Teuchos::rcp( new BarzilaiBorwein<Real>(BBtype) );
      }
    }
  }

  /** \brief Initialize step.
  */
  void initialize( const Vector<Real> &x, Objective<Real> &obj, AlgorithmState<Real> &algo_state ) {
    Teuchos::RCP<StepState<Real> >& state = Step<Real>::get_state();

    state->descentVec  = x.clone();
    state->gradientVec = x.clone();
    obj.gradient(*(state->gradientVec),x);
    algo_state.ngrad = 1;
    algo_state.gnorm = (state->gradientVec)->norm();
    algo_state.snorm = 1.e10;
    algo_state.value = obj.value(x);
    algo_state.nfval = 1;

    // Evaluate Objective Function at Cauchy Point
    if ( del_ <= ROL_EPSILON ) {
      Teuchos::RCP<Vector<Real> > Bg = x.clone();
      if ( edesc_ == DESCENT_SECANT ) {
        secant_->applyB(*Bg,*(state->gradientVec),x);
      }
      else {
        obj.hessVec(*Bg,*(state->gradientVec),x);
      }
      Real gBg   = Bg->dot(*(state->gradientVec));
      Real alpha = 1.0;
      if ( gBg > ROL_EPSILON ) {
        alpha = algo_state.gnorm*algo_state.gnorm/gBg;
      }
      Teuchos::RCP<Vector<Real> > cp = x.clone();
      cp->set(*(state->gradientVec)); 
      cp->scale(-alpha);
      Real fnew = obj.value(*cp);
      algo_state.nfval++;
      // Perform Quadratic Interpolation to Determine Initial Trust Region Radius
      Real gs = (state->gradientVec)->dot(*cp);
      del_    = -gs/(fnew - algo_state.value - gs)*alpha*algo_state.gnorm;
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
    if ( TRflag_ == 0 || TRflag_ == 1 ) {  
      if ( secant_ != Teuchos::null ) {
        gp = x.clone();
        gp->set(*(Step<Real>::state_->gradientVec));
      }
      obj.gradient(*(Step<Real>::state_->gradientVec),x);
      algo_state.ngrad++;
      if ( secant_ != Teuchos::null ) {
        secant_->update(*(Step<Real>::state_->gradientVec),*gp,s,algo_state.snorm,algo_state.iter+1);
      }
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
      hist << "\n" << EDescentToString(edesc_) 
           << " with " << ETrustRegionToString(etr_) 
           << " Trust-Region Solver\n"; 
      if ( edesc_ == DESCENT_SECANT || edesc_ == DESCENT_SECANTPRECOND ) { 
        hist << "Secant Type: " << ESecantToString(esec_) << "\n";
      }
      hist << "  ";
      hist << std::setw(6)  << std::left << "iter";  
      hist << std::setw(15) << std::left << "value";
      hist << std::setw(15) << std::left << "gnorm"; 
      hist << std::setw(15) << std::left << "snorm";
      hist << std::setw(15) << std::left << "delta";
      hist << std::setw(10) << std::left << "#fval";
      hist << std::setw(10) << std::left << "#grad";
      hist << std::setw(10) << std::left << "tr_flag";
      if (    edesc_ == DESCENT_NEWTONKRYLOV 
           || edesc_ == DESCENT_SECANTPRECOND
           || etr_   == TRUSTREGION_TRUNCATEDCG ) {
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
      if (    edesc_ == DESCENT_NEWTONKRYLOV 
           || edesc_ == DESCENT_SECANTPRECOND
           || etr_   == TRUSTREGION_TRUNCATEDCG ) {
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
