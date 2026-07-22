// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_PRIMALDUALINTERIORPOINTSTEP_H
#define ROL_PRIMALDUALINTERIORPOINTSTEP_H

#include "ROL_VectorsNorms.hpp"
#include "ROL_PartitionedVector.hpp"
#include "ROL_InteriorPointPenalty.hpp"
#include "ROL_PrimalDualInteriorPointResidual.hpp"
#include "ROL_Krylov.hpp"


/** @ingroup step_group
  
    Given the Type-EB problem
    
    \f[ \min_x f(x)\quad \text{s.t.}\quad c(x) = \f]
 
    with the bound constraint \f[ l \leq x \leq u \f]

    We introduce a barrier penalty to formuate a Type-E problem

    \f[ \min_x \varphi_\mu(x) = f(x) - \mu\sum\limits_i xscratch_i^l \ln(x_i-l_i) + xscratch_i^l \ln(u_i-x_i) \f]

    Where the weights are 
 
    \f[ l_i = \begin{cases} 0 & l_i=-\infty \\ 1 & l_i > -\infty\end{cases} ,\quad
        u_i = \begin{cases} 0 & u_i=\infty \\ 1 & u_i < \infty \end{cases} \f]


    Algorithm: 

    1) initialize parameters, counter, filter

    2) Evaluate \f$E_0\f$, check if (mu=0) problem converged

    3) Evaluate \f$E_\mu\f$, check convergence of (mu>0) problem
      3.1) Compute \f$\mu_{j+1}\f$ and \f$\tau_{j+1}\f$, increment j
      3.2) Reinitialize Filter
      3.3) If k=0, goto 3, else go to 4

    4) Compute search direction

    5) Do backtracking line search
      5.1) Initialize with max step length
      5.2) Compute trial point
      5.3) Check acceptability to filter. If good go to 5.4, else 5.5
      5.4) Check sufficient decrease. If accept go to 5.6, 5.5
      5.5) Initialize second-order correction
      5.6) Compute 2nd-order correction
      5.7) Check acceptability to filter. If reject, go to 5.10     
      5.8) Check sufficient decrease w.r.t current iterate. If accept, go to 5.6, else 5.9
      5.9) 2nd-order correction. If \f$p=p^\text{max}\f$, go to 5.6, else go to 5.10
      5.10) Choose new trial step size. If too small go to feasibility restoration phase 9, 
            else go to 5.6

   6) Accept trial point, set alpha, update muliplier estimates

   7) Augment filter if necesary

   8) Continue with next iteration, increase iteration counter and go to 2

   9) Feasibility restoration phase. Augmet filter and compute new iterate by 
      decreasing feasibility measure. Go to 8. 

   Jump to 5.10 if we encounter a NaN or Inf at a trial point  

*/

namespace ROL {
namespace InteriorPoint {
template <class Real> 
class PrimalDualInteriorPointStep : public Step<Real> {

  typedef Vector<Real>                                V;
  typedef PartitionedVector<Real>                     PV;
  typedef Objective<Real>                             OBJ;
  typedef BoundConstraint<Real>                       BND;
  typedef Krylov<Real>                                KRYLOV;
  typedef LinearOperator<Real>                        LINOP;
  typedef LinearOperatorFromConstraint<Real>          LOPEC;
  typedef Constraint<Real>                            EQCON;
  typedef StepState<Real>                             STATE;
  typedef InteriorPointPenalty<Real>                  PENALTY;
  typedef PrimalDualInteriorPointResidual<Real>       RESIDUAL;

private:

  ROL::Ptr<KRYLOV> krylov_;      // Krylov solver for the Primal Dual system 
  ROL::Ptr<LINOP>  precond_;     // Preconditioner for the Primal Dual system

  ROL::Ptr<BND> pbnd_;           // bound constraint for projecting x sufficiently away from the given bounds

  ROL::Ptr<V> x_;                // Optimization vector
  ROL::Ptr<V> g_;                // Gradient of the Lagrangian
  ROL::Ptr<V> l_;                // Lagrange multiplier

  ROL::Ptr<V> xl_;               // Lower bound vector
  ROL::Ptr<V> xu_;               // Upper bound vector

  ROL::Ptr<V> zl_;               // Lagrange multiplier for lower bound
  ROL::Ptr<V> zu_;               // Lagrange multiplier for upper bound

  ROL::Ptr<V> xscratch_;         // Scratch vector (size of x)
  ROL::Ptr<V> lscratch_;         // Scratch vector (size of l)

  ROL::Ptr<V> zlscratch_;        // Scratch vector (size of x)
  ROL::Ptr<V> zuscratch_;        // Scratch vector (size of x)

  ROL::Ptr<V> maskL_;            // Elements are 1 when xl>-INF, zero for xl = -INF
  ROL::Ptr<V> maskU_;            // Elements are 1 when xu< INF, zero for xu =  INF

  int iterKrylov_;
  int flagKrylov_;       

  bool symmetrize_;        // Symmetrize the Primal Dual system if true

  Elementwise::Multiply<Real> mult_; 

  // Parameters used by the Primal-Dual Interior Point Method

  Real mu_;                // Barrier penalty parameter
 
  Real eps_tol_;           // Error function tolerance   
  Real tau_;               // Current fraction-to-the-boundary parameter
  Real tau_min_;           // Minimum fraction-to-the-boundary parameter
  Real kappa_eps_;
  Real kappa_mu_;
  Real kappa1_;            // Feasibility projection parameter
  Real kappa2_;            // Feasibility projection parameter
  Real kappa_eps_;         // 
  Real lambda_max_;        //  multiplier maximum value
  Real theta_mu_;
  Real gamma_theta_;
  Real gamma_phi_
  Real sd_;                // Lagragian gradient scaling parameter
  Real scl_;               // Lower bound complementarity scaling parameter
  Real scu_;               // Upper bound complementarity scaling parameter
  Real smax_;              // Maximum scaling parameter

  Real diml_;              // Dimension of constaint
  Real dimx_;              // Dimension of optimization vector



  void updateState( const V& x, const V &l, OBJ &obj, 
                    EQCON &con, BND &bnd, ALGO &algo_state ) {
  
    Real tol = std::sqrt(ROL_EPSILON<Real>());
    ROL::Ptr<STATE> state = Step<Real>::getState();

    obj.update(x,true,algo_state.iter);
    con.update(x,true,algo_state.iter);

    algo_state.value = obj.value(tol);
    con.value(*(state->constraintVec),x,tol);

    obj.gradient(*(state->gradientVec),x,tol);   
    con.applyAdjointJacobian(*g_,l,x,tol);

    state->gradientVec->plus(*g_); // \f$ \nabla f(x)-\nabla c(x)\lambda \f$

    state->gradientVec->axpy(-1.0,*zl_); 
    state->gradientVec->axpy(-1.0,*zu_);

    // Scaled Lagrangian gradient sup norm
    algo_state.gnorm = normLinf(*(state->gradientVec))/sd_;

    // Constraint sup norm 
    algo_state.cnorm = normLinf(*(state->constraintVec));

    Elementwise::Multiply<Real> mult;

    Real lowerViolation;
    Real upperViolation;
     
    // Deviation from complementarity
    xscratch_->set(x);
    xscratch_->applyBinary(mult,*zl_);
    
    exactLowerViolation = normLinf(*xscratch_)/scl_;
 
    xscratch_->set(x);
    xscratch_->applyBunary(mult,*zu_);
    
    exactUpperBound = normLinf(*xscratch_)/scu_;

    // Measure ||xz||  
    algo_state.aggregateModelError = std::max(exactLowerViolation,
                                              exactUpperViolation);

  }

  /* When the constraint Jacobians are ill-conditioned, we can compute 
     multiplier vectors with very large norms, making it difficult to 
     satisfy the primal-dual equations to a small tolerance. These 
     parameters allow us to rescale the constraint qualification component
     of the convergence criteria 
  */
  void updateScalingParameters(void) {

    Real nl  = normL1(*l_);
    Real nzl = normL1(*zl_);
    Real nzu = normL1(*zu_);

    sd_  = (nl+nzl+nzu)/(diml_+2*dimx);
    sd_  = std::max(smax_,sd_)

    sd_ /= smax_;

    scl_  = nzl/dimx_;
    scl_  = std::max(smax_,scl_);
    scl_ /= smax_;

    scu_  = nzu/dimx_;
    scu_  = std::max(smax_,scu_);
    scu_ /= smax_;
  }

public:

  using Step<Real>::initialize;
  using Step<Real>::compute;
  using Step<Real>::update;

  PrimalDualInteriorPointStep( ROL::ParameterList &parlist, 
                               const ROL::Ptr<Krylov<Real> > &krylov = ROL::nullPtr,
                               const ROL::Ptr<LinearOperator<Real> > &precond = ROL::nullPtr ) : 
    Step<Real>(), krylov_(krylov), precond_(precond), iterKrylov_(0), flagKrylov_(0) {

    typedef ROL::ParameterList PL;

    PL &iplist = parlist.sublist("Step").sublist("Primal Dual Interior Point");

    kappa1_     = iplist.get("Bound Perturbation Coefficient 1",        1.e-2);
    kappa2_     = iplist.get("Bound Perturbation Coefficient 2",        1.e-2);
    lambda_max_ = iplist.get(" Multiplier Maximum Value",       1.e3 ); 
    smax_       = iplist.get("Maximum Scaling Parameter",               1.e2 );
    tau_min_    = iplist.get("Minimum Fraction-to-Boundary Parameter",  0.99 );
    kappa_mu_   = iplist.get("Multiplicative Penalty Reduction Factor", 0.2  );
    theta_mu_   = iplist.get("Penalty Update Power",                    1.5  );
    eps_tol_    = iplist.get("Error Tolerance",                         1.e-8);
    symmetrize_ = iplist.get("Symmetrize Primal Dual System",           true );

    PL &filter  = iplist.sublist("Filter Parameters");

    if(krylov_ == ROL::nullPtr) {
      krylov_ = KrylovFactory<Real>(parlist);
    }

    if( precond_ == ROL::nullPtr) {
      class IdentityOperator : public LINOP {
      public: 
        apply( V& Hv, const V &v, Real tol ) const {
          Hv.set(v);
        }
      }; // class IdentityOperator

      precond_ = ROL::makePtr<IdentityOperator<Real>>();
    }

  } 



  ~PrimalDualInteriorPointStep() {


  void initialize( Vector<Real> &x, const Vector<Real> &g, Vector<Real> &l, const Vector<Real> &c,
                   Objective<Real> &obj, Constraint<Real> &con, BoundConstraint<Real> &bnd,
                   AlgorithmState<Real> &algo_state ) {

     
    using Elementwise::ValueSet;   

    ROL::Ptr<PENALTY> &ipPen = dynamic_cast<PENALTY&>(obj);

    // Initialize step state
    ROL::Ptr<STATE> state = Step<Real>::getState();    
    state->descentVec   = x.clone();
    state->gradientVec  = g.clone();
    state->constaintVec = c.clone();

    diml_ = l.dimension();
    dimx_ = x.dimension();

    Real one(1.0); 
    Real zero(0.0);
    Real tol = std::sqrt(ROL_EPSILON<Real>());

    x_ = x.clone();
    g_ = g.clone();
    l_ = l.clone();
    c_ = c.clone();

    xscratch_ = x.clone(); 
    lscratch_ = l.clone(); 

    zlscratch_ = x.clone(); 
    zuscratch_ = x.clone(); 

    // Multipliers for lower and upper bounds 
    zl_ = x.clone();
    zu_ = x.clone();

    /*******************************************************************************************/
    /* Identify (implicitly) the index sets of upper and lower bounds by creating mask vectors */
    /*******************************************************************************************/
  
    xl_ = bnd.getLowerBound();
    xu_ = bnd.getUpperBound();

    maskl_ = ipPen.getLowerMask();
    masku_ = ipPen.getUpperMask();
    
    // Initialize bound constraint multipliers to 1 one where the corresponding bounds are finite
    zl_->set(maskl_);
    zu_->set(masku_);

    /*******************************************************************************************/
    /* Create a new bound constraint with perturbed bounds                                     */
    /*******************************************************************************************/
    
    ROL::Ptr<V> xdiff = xu_->clone();
    xdiff->set(*xu_); 
    xdiff->axpy(-1.0,*xl_);
    xdiff->scale(kappa2_);    

    class Max1X : public Elementwise::UnaryFunction<Real> {
    public:
      Real apply( const Real &x ) const {
        return std::max(1.0,x);
      }
    };

    Max1X                            max1x;
    Elementwise::AbsoluteValue<Real> absval;
    Elementwise::Min                 min;

    // Lower perturbation vector
    ROL::Ptr<V> pl = xl_->clone();
    pl->applyUnary(absval);
    pl->applyUnary(max1x);               // pl_i = max(1,|xl_i|)
    pl->scale(kappa1_);
    pl->applyBinary(min,xdiff);          // pl_i = min(kappa1*max(1,|xl_i|),kappa2*(xu_i-xl_i))

    // Upper perturbation vector
    ROL::Ptr<V> pu = xu_->clone();
    pu->applyUnary(absval);
    pu->applyUnary(max1x);              // pu_i = max(1,|xu_i|)
    pu->scale(kappa_1); 
    pu->applyBinary(min,xdiff);         // pu_u = min(kappa1*max(1,|xu_i|,kappa2*(xu_i-xl_i)))

    // Modified lower and upper bounds so that x in [xl+pl,xu-pu] using the above perturbation vectors
    pl->plus(*xl_);
    pu->scale(-1.0);
    pu->plus(*xu_);

    pbnd_ = ROL::makePtr<BoundConstraint<Real>>(pl,pu)

    // Project the initial guess onto the perturbed bounds
    pbnd_->project(x);


    /*******************************************************************************************/
    /* Solve least-squares problem for initial equality multiplier                             */
    /*******************************************************************************************/
   
    // \f$-(\nabla f-z_l+z_u) \f$
    g_->set(*zl_);
    g_->axpy(-1.0,g);
    g_->axpy(-1.0,*zu_)

    // We only need the updated multiplier
    lscratch_->zero();
    con.solveAugmentedSystem(*xscratch_,*l_,*g_,*lscratch_,x,tol);

    // If the multiplier supremum is too large, zero the vector as this appears to avoid poor
    // initial guesses of the multiplier in the case where the constraint Jacobian is 
    // ill-conditioned 
    if( normInf(l_) > lambda_max_ ) {
      l_->zero();
    }
 
    // Initialize the algorithm state
    algo_state.nfval = 0;
    algo_state.ncval = 0;
    algo_state.ngrad = 0;

    updateState(x,l,obj,con,bnd,algo_state);

  } // initialize()


  void compute( Vector<Real> &s, const Vector<Real> &x, const Vector<Real> &l, 
                Objective<Real> &obj, Constraint<Real> &con, 
                BoundConstraint<Real> &bnd, AlgorithmState<Real> &algo_state ) {


      

    Elementwise::Fill<Real>     minus_mu(-mu_); 
    Elementwise::Divide<Real>   div;
    Elementwise::Multiply<Real> mult;

    ROL::Ptr<STATE> state = Step<Real>::getState();    

    ROL::Ptr<OBJ>   obj_ptr = ROL::makePtrFromRef(obj);
    ROL::Ptr<EQCON> con_ptr = ROL::makePtrFromRef(con); 
    ROL::Ptr<BND>   bnd_ptr = ROL::makePtrFromRef(bnd);


    /*******************************************************************************************/
    /* Form Primal-Dual system residual and operator then solve for direction vector           */
    /*******************************************************************************************/
    
          
    ROL::Ptr<V> rhs = CreatePartitionedVector(state->gradientVec,
                                         state->constraintVec,
                                         resL_,
                                         resU_);
    
    ROL::Ptr<V> sysvec = CreatePartitionedVector( x_, l_, zl_, zu_ );


    ROL::Ptr<RESIDUAL> residual = ROL::makePtr<RESIDUAL>(obj,con,bnd,*sol,maskL_,maskU_,w_,mu_,symmetrize_); 

    residual->value(*rhs,*sysvec,tol);

    ROL::Ptr<V> sol = CreatePartitionedVector( xscratch_, lscratch_, zlscratch_, zuscratch_ );

    LOPEC jacobian( sysvec, residual ); 

    
    krylov_->run(*sol,jacobian,*residual,*precond_,iterKrylov_,flagKrylov_); 
    
    /*******************************************************************************************/
    /* Perform line search                                                                     */
    /*******************************************************************************************/



  } // compute() 
  


  void update( Vector<Real> &x, Vector<Real> &l, const Vector<Real> &s,
               Objective<Real> &obj, Constraint<Real> &con, 
               BoundConstraint<Real> &bnd, AlgorithmState<Real> &algo_state ) {

    // Check deviation from shifted complementarity
    Elementwise::Shift<Real>    minus_mu(-mu_);

    xscratch_->set(x);
    xscratch_->applyBinary(mult,*zl_);
    xscratch_->applyUnary(minus_mu);

    lowerViolation = normLinf(*xscratch_)/scl_; // \f$ \max_i xz_l^i - \mu \f$
 
    xscratch_->set(x);
    xscratch_->applyBinary(mult,*zu_);
    xscratch_->applyUnary(minus_mu);            
    
    upperBound = normLinf(*xscratch_)/scu_;

    // Evaluate \f$E_\mu(x,\lambda,z_l,z_u)\f$
    Real Emu = algo_state.gnorm;
    Emu = std::max(Emu,algo_state.cnorm);
    Emu = std::max(Emu,upperBound);
    Emu = std::max(Emu,lowerBound);

    // If sufficiently converged for the current mu, update it
    if(Emu < (kappa_epsilon_*mu_) ) {
      Real mu_old = mu_;
 
      /* \mu_{j+1} = \max\left{ \frac{\epsilon_\text{tol}}{10},
                                \min\{ \kappa_{\mu} \mu_j, 
                                       \mu_j^{\theta_\mu}\} \right\} */
      mu_ = std::min(kappa_mu_*mu_old,std::pow(mu_old,theta_mu_));
      mu_ = std::max(eps_tol_/10.0,mu_);
 
     // Update fraction-to-boundary parameter
     tau_ = std::max(tau_min_,1.0-mu_);     

          
 
    } 

  } // update() 




  // TODO: Implement header print out
  std::string printHeader( void ) const {
    std::string head("");
    return head;
  }

  // TODO: Implement name print out
  std::string printName( void ) const {
    std::string name("");
    return name;
  }


}; // class PrimalDualInteriorPointStep

} // namespace InteriorPoint
} // namespace ROL


#endif //  ROL_PRIMALDUALINTERIORPOINTSTEP_H
