//
// @HEADER
// ***********************************************************************
// 
//                           Rythmos Package
//                 Copyright (2005) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef Rythmos_STEPPER_IMPLICITBDF_H
#define Rythmos_STEPPER_IMPLICITBDF_H

#include "Rythmos_Stepper.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "Thyra_VectorBase.hpp"
#include "Thyra_ModelEvaluator.hpp"
#include "Thyra_ModelEvaluatorHelpers.hpp"
#include "Thyra_NonlinearSolverBase.hpp"
#include "Thyra_SingleResidSSDAEModelEvaluator.hpp"

namespace Rythmos {

/** \brief . */
template<class Scalar>
class ImplicitBDFStepper : public Stepper<Scalar>
{
  public:

    /** \brief . */
    ImplicitBDFStepper(
      const Teuchos::RefCountPtr<const Thyra::ModelEvaluator<Scalar> >        &model
      ,const Teuchos::RefCountPtr<const Thyra::NonlinearSolverBase<Scalar> >  &solver
      );

    /** \brief . */
    void setModel(const Teuchos::RefCountPtr<const Thyra::ModelEvaluator<Scalar> > &model);

    /** \brief . */
    void setSolver(const Teuchos::RefCountPtr<const Thyra::NonlinearSolverBase<Scalar> > &solver);

    /** \brief . */
    Scalar TakeStep(Scalar dt);
   
    /** \brief . */
    Scalar TakeStep();

    /** \brief . */
    Teuchos::RefCountPtr<const Thyra::VectorBase<Scalar> > get_solution() const;

    /** \brief . */
    Teuchos::RefCountPtr<const Thyra::VectorBase<Scalar> > get_residual() const;

    /** \brief . */
    std::string description() const;

    /** \brief . */
    std::ostream& describe(
      std::ostream                &out
      ,const Teuchos::EVerbosityLevel      verbLevel
      ,const std::string          leadingIndent
      ,const std::string          indentSpacer
      ) const;

    bool ErrWtVecSet(
      Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > &w, 
      Teuchos::RefCountPtr<const Thyra::VectorBase<Scalar> > &y
      );


  private:

    void obtainPredictor();
    void obtainResidual();
    void obtainJacobian();
    bool interpolateSolution(Scalar time);
    void updateHistory();
    void restoreHistory();
    void updateCoeffs();
    void initialize();
    void checkReduceOrder();
    void rejectStep();
    void completeStep();

    void setDefaultMagicNumbers(Teuchos::ParameterList &magicNumberList);

    // 05/05/06 tscoffe:  I hate the underscores for private variables!
    Teuchos::RefCountPtr<const Thyra::ModelEvaluator<Scalar> > model;
    Teuchos::RefCountPtr<const Thyra::NonlinearSolverBase<Scalar> > solver;
    Thyra::SingleResidSSDAEModelEvaluator<Scalar>   neModel;

    Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > xn0;
    Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > xpn0;
    std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > > xHistory;
    Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > residual;
    Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > errWtVec;
    Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > ee;

    Scalar time;


    ScalarMag relErrTol; // relative error tolerance
    ScalarMag absErrTol; // absolute error tolerance
    Scalar hh;        // Current step-size
    int currentOrder; // Current order of integration
    int oldOrder;     // previous order of integration
    int maxOrder;     // maximum order = min(5,user option maxord) - see below.
    int usedOrder;    // order used in current step (used after currentOrder is updated)
    Scalar alphas;    // $\alpha_s$ fixed-leading coefficient of this BDF method
    vector<Scalar> alpha;    // $\alpha_j(n)=h_n/\psi_j(n)$ coefficient used in local error test
                      // note:   $h_n$ = current step size, n = current time step
    Scalar alpha0;     // $-\sum_{j=1}^k \alpha_j(n)$ coefficient used in local error test
    Scalar cj ;        // $-\alpha_s/h_n$ coefficient used in local error test
    Scalar ck ;        // local error coefficient gamma[0] = 0; // $\gamma_j(n)=\sum_{l=1}^{j-1}1/\psi_l(n)$ coefficient used to
    vector<Scalar> gamma;    // calculate time derivative of history array for predictor 
    vector<Scalar> beta;     // coefficients used to evaluate predictor from history array
    vector<Scalar> psi;      // $\psi_j(n) = t_n-t_{n-j}$ intermediary variable used to 
                      // compute $\beta_j(n)$
    vector<Scalar> sigma;    // $\sigma_j(n) = \frac{h_n^j(j-1)!}{\psi_1(n)*\cdots *\psi_j(n)}$
    int numberOfSteps;// number of total time integration steps taken
    int  nef;
    Scalar usedStep;
    int nscsco;
    Scalar Ek;
    Scalar Ekm1;
    Scalar Ekm2;
    Scalar Ekp1;
    Scalar Est;
    Scalar Tk;
    Scalar Tkm1;
    Scalar Tkm2;
    Scalar Tkp1;
    int newOrder;
    bool initialPhase;
    Scalar stopTime;
    Scalar nextTimePt;
    Scalar minTimeStep;
    Scalar maxTimeStep;

    Teuchos::ParameterList magicNumber;

    enum actionFlag { ACTION_UNSET, ACTION_LOWER, ACTION_MAINTAIN, ACTION_RAISE };

};

// ////////////////////////////
// Defintions

template<class Scalar>
ImplicitBDFStepper<Scalar>::ImplicitBDFStepper(
  const Teuchos::RefCountPtr<const Thyra::ModelEvaluator<Scalar> > &model
  ,const Teuchos::RefCountPtr<const Thyra::NonlinearSolverBase<Scalar> > &solver
  )
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  // We should accept a parameter list to specify the maxord, reltol, abstol, and magic numbers
  // Initialize algorithm coefficients
  currentOrder=1; // Current order of integration
  oldOrder=1; // previous order of integration
  maxOrder=5;  // maximum order = min(5,user option maxord) - see below.
  usedOrder=1;  // order used in current step (used after currentOrder is updated)
  alphas=Scalar(-ST::one());  // $\alpha_s$ fixed-leading coefficient of this BDF method
  alpha.reserve(maxOrder+1);  // $\alpha_j(n)=h_n/\psi_j(n)$ coefficient used in local error test
                  // note:   $h_n$ = current step size, n = current time step
  gamma.reserve(maxOrder+1);  // calculate time derivative of history array for predictor 
  beta.reserve(maxOrder+1);   // coefficients used to evaluate predictor from history array
  psi.reserve(maxOrder+1);    // $\psi_j(n) = t_n-t_{n-j}$ intermediary variable used to 
                  // compute $\beta_j(n;$
  sigma.reserve(maxOrder+1);  // $\sigma_j(n) = \frac{h_n^j(j-1)!}{\psi_1(n)*\cdots *\psi_j(n)}$
  for (int i=0 ; i<maxOrder ; ++i)
  {
    alpha.push_back(ST::zero());
    beta.push_back(ST::zero());
    gamma.push_back(ST::zero());
    psi.push_back(ST::zero());
    sigma.push_back(ST::zero());
  }
  alpha0=ST::zero();   // $-\sum_{j=1}^k \alpha_j(n)$ coefficient used in local error test
  cj=ST::zero();      // $-\alpha_s/h_n$ coefficient used in local error test
  ck=ST::zero();      // local error coefficient gamma_[0] = 0; // $\gamma_j(n)=\sum_{l=1}^{j-1}1/\psi_l(n)$ coefficient used to
  numberOfSteps=0;   // number of total time integration steps taken
  nef=0;
  usedStep=ST::zero();
  nscsco=0;
  Ek=ST::zero();
  Ekm1=ST::zero();
  Ekm2=ST::zero();
  Ekp1=ST::zero();
  Est=ST::zero();
  Tk=ST::zero();
  Tkm1=ST::zero();
  Tkm2=ST::zero();
  Tkp1=ST::zero();
  newOrder=1;
  initialPhase=true;

  Teuchos::ParamterList magicNumber;
  setDefaultMagicNumbers(magicNumber);
  
  // Now we instantiate the model and the solver
  setModel(model);
  setSolver(solver);
}

template<class Scalar>
ImplicitBDFStepper<Scalar>::setDefaultMagicNumbers(
    Teuchos::ParameterList &magicNumberList)
{
  // Magic Number Defaults:
  magicNumberList.set( "h0_safety",      Scalar(2.0)     );
  magicNumberList.set( "h0_max_factor",  Scalar(0.0001)  );
  magicNumberList.set( "h_phase0_incr",  Scalar(2.0)     );
  magicNumberList.set( "h_max_inv",      ST::zero()      );
  magicNumberList.set( "Tkm1_Tk_safety", Scalar(2.0)     );
  magicNumberList.set( "Tkp_Tk_safety",  Scalar(0.5)     );
  magicNumberList.set( "r_factor",       Scalar(0.9)     );
  magicNumberList.set( "r_safety",       Scalar(2.0)     );
  magicNumberList.set( "r_fudge",        Scalar(0.0001)  );
  magicNumberList.set( "r_min",          Scalar(0.125)   );
  magicNumberList.set( "r_max",          Scalar(0.9)     );
  magicNumberList.set( "r_hincr_test",   Scalar(2.0)     );
  magicNumberList.set( "r_hincr",        Scalar(2.0)     );
  magicNumberList.set( "max_LET_fail",   Scalar(15)      );
}

template<class Scalar>
void ImplicitBDFStepper<Scalar>::setModel(const Teuchos::RefCountPtr<const Thyra::ModelEvaluator<Scalar> > &model_)
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  model = model_;
  t = ST::zero();
  xn0 = model->getNominalValues().get_x()->clone_v();
  xpn0 = model->getNominalValues().get_x_dot()->clone_v();
  residual = Thyra::createMember(model->get_f_space());
  errWtVec = xn0->clone_v();
  xHistory.push_back(xn0->clone_v());
  xHistory.push_back(xpn0->clone_v());
  for (int i=2 ; i<maxOrder ; ++i)
  {
    xHistory.push_back(xn0->clone_v()); 
    V_S(&*xHistory[i],ST::zero());
  }
  ee = Thyra::createMember(model->get_x_space());

}

template<class Scalar>
void ImplicitBDFStepper<Scalar>::setSolver(const Teuchos::RefCountPtr<const Thyra::NonlinearSolverBase<Scalar> > &solver_)
{
  solver = solver_;
}

template<class Scalar>
Scalar ImplicitBDFStepper<Scalar>::TakeStep()
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  // Set up problem coefficients (and handle first step)
  updateCoeffs();
  // compute predictor
  obtainPredictor();
  // solve nonlinear problem (as follows)
  
  //
  // Setup the nonlinear equations:
  //
  //   f_bar( x_dot_coeff * x_bar + x_dot_base, x_coeff * x_bar + x_base, t_base ) = 0
  //   x_dot_coeff = -alpha_s/dt
  //   x_dot_base = x_prime_pred + (alpha_s/dt) * x_pred
  //   x_coeff = 1
  //   x_base = 0
  //   t_base = tn+dt
  //
  Scalar coeff_x_dot = Scalar(-ST::one())*alpha_s/dt;
  V_StVpStV( &*x_dot_base, ST::one(), *y_dot_pred, alpha_s/dt, *y_pred );
  neModel.initialize(model,coeff_x_dot,x_dot_base,ST::one(),Teuchos::null,t_old+dt,y_pred);
  //
  // Solve the implicit nonlinear system to a tolerance of ???
  // 
  // 05/08/06 tscoffe:  I really need to get the update, not the solution from
  // the nonlinear solver.
  solver->solve( neModel, &*x, NULL, &*ee ); 
  
  // check error and evaluate LTE
  checkReduceOrder();
  
  // reject step and try again or complete step and finish
  rejectStep(); // this restores history and computes next step
  completeStep();  // this should only occur after the LTE passes
  
  return(dt);
}

template<class Scalar>
Scalar ImplicitBDFStepper<Scalar>::TakeStep(Scalar dt)
{
  //Not supported at this time.
  return(dt);
}

template<class Scalar>
Teuchos::RefCountPtr<const Thyra::VectorBase<Scalar> > ImplicitBDFStepper<Scalar>::get_solution() const
{
  return(x);
}

template<class Scalar>
Teuchos::RefCountPtr<const Thyra::VectorBase<Scalar> > ImplicitBDFStepper<Scalar>::get_residual() const
{
  return(f);
}

template<class Scalar>
std::string ImplicitBDFStepper<Scalar>::description() const
{
  std::string name = "Rythmos::ImplicitBDFStepper";
  return(name);
}

template<class Scalar>
std::ostream& ImplicitBDFStepper<Scalar>::describe(
      std::ostream                &out
      ,const Teuchos::EVerbosityLevel      verbLevel
      ,const std::string          leadingIndent
      ,const std::string          indentSpacer
      ) const
{
  if (verbLevel == Teuchos::VERB_EXTREME)
  {
    out << description() << "::describe" << std::endl;
    out << "model = " << std::endl;
    out << model->describe(out,verbLevel,leadingIndent,indentSpacer) << std::endl;
    out << "solver = " << std::endl;
    out << solver->describe(out,verbLevel,leadingIndent,indentSpacer) << std::endl;
    out << "x = " << std::endl;
    out << x->describe(out,verbLevel,leadingIndent,indentSpacer) << std::endl;
    out << "scaled_x_old = " << std::endl;
    out << scaled_x_old->describe(out,verbLevel,leadingIndent,indentSpacer) << std::endl;
    out << "f = " << std::endl;
    out << f->describe(out,verbLevel,leadingIndent,indentSpacer) << std::endl;
    out << "t = " << t << std::endl;
    out << "t_old = " << t_old << std::endl;
//    out << "neModel = " << std::endl;
//    out << neModel->describe(out,verbLevel,leadingIndent,indentSpacer) << std::endl;
  }
  return(out);
}

template<class Scalar>
void ImplicitBDFStepper<Scalar>
template<class Scalar>
void ImplicitBDFStepper<Scalar>::obtainPredictor()
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  
  // prepare history array for prediction
  for (int i=nscsco;i<=currentOrder;++i)
  {
    Vt_S(&*xHistory[i],beta[i]);
  }
  
  // evaluate predictor
  V_V(&*xn0,*xHistory[0]);
  V_S(&*xpn0,ST::zero());
  for (int i=1;i<=currentOrder;++i)
  {
    Vp_V(&*xn0,*xHistory[i]);
    Vp_StV(&*xpn0,gamma[i],*xhistory[i]);
  }
}

template<class Scalar>
void ImplicitBDFStepper<Scalar>::obtainResidual()
{
  // output: dsDae.RHSVectorPtr

  // This function returns the following residual:
  // $qpn0 - (alphas/hn)(Q(x)-qn0)+F(x)-B(t)$

  // Note:  dsDae.nextSolutionPtr is used to get Q,F,B in N_TIA_ControlAlgorithm::loadRHS.
  dsDae.RHSVectorPtr->linearCombo(1.0,*dsDae.daeQVectorPtr,-1.0,*dsDae.qn0Ptr);
  dsDae.RHSVectorPtr->linearCombo(1.0,*dsDae.qpn0Ptr,-alphas/hh,*dsDae.RHSVectorPtr);

  dsDae.RHSVectorPtr->linearCombo(1.0,*dsDae.RHSVectorPtr,+1.0,*dsDae.daeFVectorPtr);
  dsDae.RHSVectorPtr->linearCombo(1.0,*dsDae.RHSVectorPtr,-1.0,*dsDae.daeBVectorPtr);


}

template<class Scalar>
void ImplicitBDFStepper<Scalar>::obtainJacobian()
{

  // output: dsDae.JMatrixPtr

  // This function returns the following matrix:
  // $-(alphas/hn)dQdx(x)+dFdx$

  // Note:  dsDae.nextSolutionPtr is used to get dQdx,dFdx in N_TIA_ControlAlgorithm::loadJacobian.
//  dsDae.JmatrixPtr->linearCombo(-alphas/hh,*dsDae.dQdxMatrixPtr,+1.0,*dsDae.dFdxMatrixPtr);

  N_LAS_Matrix & dQdx = *(dsDae.dQdxMatrixPtr);
  N_LAS_Matrix & dFdx = *(dsDae.dFdxMatrixPtr);
  N_LAS_Matrix & Jac = *(dsDae.JMatrixPtr);

  Jac.add( dQdx );

  Jac.scale( -alphas/hh );

  Jac.add( dFdx );

}

template<class Scalar>
void ImplicitBDFStepper<Scalar>::interpolateSolution(Scalar timepoint, 
    	N_LAS_Vector * tmpSolVectorPtr)
{

// 03/23/04 tscoffe:  Currently this code is nearly identical to the IDA code
// for interpolating to an output time.  Either we acknowledge the copyright,
// the list of conditions in the license and the disclaimer or we rewrite this
// function.  The IDA license is included after this routine.
  Scalar tfuzz;   // fuzz factor to check for valid output time
  Scalar tp;      // approximately t{n-1}
  Scalar delt;    // distance between timepoint and currentTime
  Scalar c = 1.0; // coefficient for interpolation
  Scalar gam;     // coefficient for interpolation
  int kord;       // order of interpolation
  Scalar tn = currentTime;
  Scalar hh = hh;
  Scalar hused = usedStep;
  int kused = usedOrder;
  Scalar uround = ST::zero();  // unit round-off (set to zero for now)

  tfuzz = 100 * uround * (tn + hh);
  tp = tn - hused - tfuzz;
  if ( (timepoint - tp)*hh < ST::zero() ) 
    return false;

  *tmpSolVectorPtr = dsDae.xHistory[0];
  kord = kused;
  if ( (kused == 0) || (timepoint == tn) ) 
    kord = 1;

  delt = timepoint - tn;
  gam = delt/psi[0];
  for (int j=1 ; j <= kord ; ++j)
  {
    c = c*gam;
    gam = (delt + psi[j-1])/psi[j];
    tmpSolVectorPtr->linearCombo(1.0,*tmpSolVectorPtr,c,dsDae.xHistory[j]);
  }
  return true;
}

template<class Scalar>
void ImplicitBDFStepper<Scalar>::updateHistory()
{

  // Save Newton correction for potential order increase on next step.
  if (usedOrder < maxOrder)  
  {
    assign( &*xHistory[usedOrder+1], *ee );
  }
  // Update history arrays
  Vp_V( &*xHistory[usedOrder], *ee );
  for (int j=usedOrder-1;j>=0;j--) 
  {
    Vp_V( &*xHistory[j], *xHistory[j+1] );
  }

}

template<class Scalar>
void ImplicitBDFStepper<Scalar>::restoreHistory()
{

  // undo preparation of history array for prediction
  for (int i=nscsco;i<=currentOrder;++i)
  {
    Vt_S( &*xHistory[i], ST::one()/beta[i] );
  }
  for (int i=1;i<=currentOrder;++i)
  {
    psi[i-1] = psi[i] - hh;
  }
} 

template<class Scalar>
void ImplicitBDFStepper<Scalar>::updateCoeffs()
{
  // synchronize with Step Error Control
//  psi[0] = hh;
  // If the number of steps taken with constant order and constant stepsize is
  // more than the current order + 1 then we don't bother to update the
  // coefficients because we've reached a constant step-size formula.  When
  // this is is not true, then we update the coefficients for the variable
  // step-sizes. 
  if ((hh != usedStep) || (currentOrder != usedOrder))
    nscsco = 0;
  nscsco = min(nscsco+1,usedOrder+2);
  if (currentOrder+1 >= nscsco)
  {
    beta[0] = ST::one();
    alpha[0] = ST::one();
    Scalar temp1 = hh;
    sigma[0] = ST::one();
    gamma[0] = ST::zero();
    for (int i=1;i<=currentOrder;++i)
    {
      Scalar temp2 = psi[i-1];
      psi[i-1] = temp1;
      beta[i] = beta[i-1]*psi[i-1]/temp2;
      temp1 = temp2 + hh;
      alpha[i] = hh/temp1;
      sigma[i] = Scalar(i+1)*sigma[i-1]*alpha[i];
      gamma[i] = gamma[i-1]+alpha[i-1]/hh;
    }
    psi[currentOrder] = temp1;
    alphas = ST::zero();
    alpha0 = ST::zero();
    for (int i=0;i<currentOrder;++i)
    {
      alphas = alphas - Scalar(ST::one()/(i+ST::one()));
      alpha0 = alpha0 - alpha[i];
    }
    cj = -alphas/hh;
    ck = abs(alpha[currentOrder]+alphas-alpha0);
    ck = max(ck,alpha[currentOrder]);
  }
}

template<class Scalar>
void ImplicitBDFStepper<Scalar>::initialize()
{

  // we assume the solution vector is available here 
  // Note that I'm using currSolutionPtr instead of
  // nextSolutionPtr because this is the first step.

  // Update next stop time from StepErrorControl:
  updateStopTime();

  // Choose initial step-size
  Scalar time_to_stop = stopTime - currentTime;
  Scalar currentTimeStep;
  if (constantStepSize)
  {
    currentTimeStep = 0.1 * time_to_stop;
    currentTimeStep = Xycemin(startingTimeStep, currentTimeStep);
  }
  else
  {
    // compute an initial step-size based on rate of change in the solution initially
    Scalar dnorm_q = *(dsDae.qHistory[1].wRMSNorm(*dsDae.qErrWtVecPtr));
    if (dnorm_q > ST::zero())  // time-dependent DAE
    {
      currentTimeStep = Xycemin(h0_max_factor*abs(time_to_stop),sqrt(2.0)/(h0_safety*dnorm_q));
    } 
    else  // non-time-dependent DAE
    {
      currentTimeStep = h0_max_factor*abs(time_to_stop);
    }
    // choose min of user specified value and our value:
    if (startingTimeStep > ST::zero())
      currentTimeStep = Xycemin(startingTimeStep, currentTimeStep);
    // check for maximum step-size:
    Scalar rh = abs(currentTimeStep)*h_max_inv; 
    if (rh>1.0) currentTimeStep = currentTimeStep/rh;
  }
  hh = currentTimeStep;

  currentTimeStepRatio = 1.0;
  currentTimeStepSum   = 2.0*hh;

  lastTimeStep      = hh;
  lastTimeStepRatio = currentTimeStepRatio;
  lastTimeStepSum   = currentTimeStepSum;

  numberSuccessiveFailures = 0;
  stepAttemptStatus        = true;

//  tolAimFac = 0.5;

  nextTime = currentTime + hh;

  // x history
  dsDae.xHistory[0] = *(dsDae.currSolutionPtr);
  dsDae.xHistory[1].putScalar(ST::zero()); // no need to multiply by dt here

  // q history
  dsDae.qHistory[0] = *(dsDae.daeQVectorPtr);
  dsDae.qHistory[1].linearCombo(-1.0,*dsDae.daeFVectorPtr,1.0,*dsDae.daeBVectorPtr);
  dsDae.qHistory[1].scale(hh);

  // state history
  dsDae.sHistory[0] = *(dsDae.nextStatePtr);
  dsDae.sHistory[1].putScalar(ST::zero()); 

  // Coefficient initialization 
  numberOfSteps = 0;    // number of total time integration steps taken
  currentOrder = 1;
  usedOrder = 1;
  psi[0] = hh;
  cj = 1/psi[0];
  nscsco = 0;
}

template<class Scalar>
void ImplicitBDFStepper<Scalar>::checkReduceOrder()
{

// This routine puts its output in newOrder_

// This routine changes the following variables:
//    Ek, Tk, Est, newOrder, dsDae.delta_x, dsDae.delta_q,
//    Ekm1, Tkm1, Ekm2, Tkm2 

// This routine reads but does not change the following variables:
//    currentOrder, sigma, dsDae.newtonCorrectionPtr, dsDae.qNewtonCorrectionPtr,
//    dsDae.errWtVecPtr, dsDae.qErrWtVecPtr, dsDae.xHistory, dsDae.qHistory

  Scalar dnorm = norm_2(*errWtVec,*delta);
  Ek = sigma[currentOrder]*dnorm;
  Tk = Scalar(currentOrder+1)*Ek;
  Est = Ek;
  newOrder = currentOrder;
  if (currentOrder>1)
  {
    V_VpV(&*delta,*xHistory[currentOrder],ee);
    dnorm = norm_2(*errWtVec,*delta);
    Ekm1 = sigma[currentOrder-1]*dnorm;
    Tkm1 = currentOrder*Ekm1;
    if (currentOrder>2)
    {
      Vp_V(&*delta,*xHistory[currentOrder-1]);
      dnorm = norm_2(*errWtVec,*delta);
      Ekm2 = sigma[currentOrder-2]*dnorm;
      Tkm2 = (currentOrder-1)*Ekm2;
      if (max(Tkm1,Tkm2)<=Tk)
      {
        newOrder--;
        Est = Ekm1;
      }
    }
    else if (Tkm1 <= Tkm1_Tk_safety * Tk)
    {
      newOrder--;
      Est = Ekm1;
    }
  }
  // Check LTE here:
  if ((ck*dnorm) > ST::one())
  {
    // failed LTE
  } else {
    // passed LTE
  }
}

template<class Scalar>
void ImplicitBDFStepper<Scalar>::rejectStep()
{

// This routine puts its output in newTimeStep and newOrder

// This routine changes the following variables:
//    lastAttemptedTimeStep, initialPhase, nef, psi, newTimeStep,
//    newOrder, currentOrder, currentTimeStep, dsDae.xHistory,
//    dsDae.qHistory, nextTimePt, nextTime, currentTimeStepRatio,
//    currentTimeStepSum, nextTimePt

// This routine reads but does not change the following variables:
//    stepAttemptStatus, r_factor, r_safety, Est, r_fudge, r_min, r_max,
//    minTimeStep, maxTimeStep, currentTime, stopTime, lastTimeStep


  // Only update the time step if we are NOT running constant stepsize.
  bool adjustStep = (!constantStepSize);

  Scalar newTimeStep = hh;
  Scalar rr = 1.0; // step size ratio = new step / old step
  if ((stepAttemptStatus == false) && (adjustStep))
  {
    initialPhase = false;
    nef++;
    restoreHistory();
    // restore psi_
//    for (int i=1;i<=currentOrder;++i)
//      psi[i-1] = psi[i] - hh;

    if (nef >= max_LET_fail)  
    {
      string msg = "N_TIA_BackwardDifferentiation15::rejectStep: ";
      msg += "  Maximum number of local error test failures.  ";
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, msg);
    }

    if ((newtonConvergenceStatus <= 0))
    {
      /// 11/11/05 erkeite:  If the Newton solver fails, don't 
      // rely on the error estimate - it may be full of Nan's.
      rr = r_min;
      newTimeStep = rr * hh;

      if (nef > 2) newOrder = 1;//consistent with block below.
    }
    else
    {
      // 03/11/04 tscoffe:  Here is the block for choosing order & 
      // step-size when the local error test FAILS (but Newton 
      // succeeded). 
      if (nef == 1) // first local error test failure
      {
        rr = r_factor*pow(r_safety*Est+r_fudge,-1.0/(newOrder+1.0));
        rr = max(r_min,Xycemin(r_max,rr));
        newTimeStep = rr * hh;
      }
      else if (nef == 2) // second failure
      {
        rr = r_min;
        newTimeStep = rr * hh;
      }
      else if (nef > 2) // third and later failures
      {
        newOrder = 1;
        rr = r_min;
        newTimeStep = rr * hh;
      }
    }
    currentOrder = newOrder;
    if (numberOfSteps == 0) // still first step
    {
      psi[0] = newTimeStep;
      Vt_S(&*xHistory[1],rr);
    }
  }
  else if ((stepAttemptStatus == false) & (!adjustStep))
  {
    string tmp = "  BackwardDifferentiation15:rejectStep: Warning: Local error test failed with constant step-size.\n";
    N_ERH_ErrorMgr::report (N_ERH_ErrorMgr::USR_INFO_0, tmp);
  }

  // If the step needs to be adjusted:
  if (adjustStep)
  {
    newTimeStep = max(newTimeStep, minTimeStep);
    newTimeStep = min(newTimeStep, maxTimeStep);

    Scalar nextTimePt = currentTime + newTimeStep;

    if (nextTimePt > stopTime)
    {
      nextTimePt  = stopTime;
      newTimeStep = stopTime - currentTime;
    }

    hh = newTimeStep;
  }
  else // if time step is constant for this step:
  {
    Scalar nextTimePt = currentTime + hh;

    if (nextTimePt > stopTime)
    {
      nextTimePt      = stopTime;
      hh = stopTime - currentTime;
    }
  }
}

template<class Scalar>
void ImplicitBDFStepper<Scalar>::completeStep()
{

  numberOfSteps ++;
  nef = 0;
  lastTime    = currentTime;
  currentTime = nextTime;
  //
  // Only update the time step if we are NOT running constant stepsize.
  bool adjustStep = (!constantStepSize);

  lastAttemptedTimeStep = hh;

  Scalar newTimeStep = hh;
  Scalar rr = ST::one(); // step size ratio = new step / old step
  // 03/11/04 tscoffe:  Here is the block for choosing order & step-size when
  // the local error test PASSES (and Newton succeeded). 
  lastTimeStep = hh;
  lastTimeStepRatio = currentTimeStepRatio; // copied from calcTStep1
  lastTimeStepSum   = currentTimeStepSum; // copied from calcTStep1
  int orderDiff = currentOrder - usedOrder;
  usedOrder = currentOrder;
  usedStep = hh;
  if ((newOrder == currentOrder-1) || (currentOrder == maxOrder))
  {
    // If we reduced our order or reached max order then move to the next phase
    // of integration where we don't automatically double the step-size and
    // increase the order.
    initialPhase = false;
  }
  if (initialPhase)
  {
    currentOrder++;
    newTimeStep = h_phase0_incr * hh;
  }
  else // not in the initial phase of integration
  {
    int action = ACTION_UNSET;
    if (newOrder == currentOrder-1)
      action = ACTION_LOWER;
    else if (newOrder == maxOrder)
      action = ACTION_MAINTAIN;
    else if ((currentOrder+1>=nscsco) || (orderDiff == 1))
    {
      // If we just raised the order last time then we won't raise it again
      // until we've taken currentOrder+1 steps at order currentOrder.
      action = ACTION_MAINTAIN;
    }
    else // consider changing the order 
    {
      V_StVpStV(&*delta,ST::one(),*ee,Scalar(-ST::one()),*xHistory[currentOrder+1]);
      Scalar dNorm = norm_2(*errWtVec,*delta);
      Tkp1 = dnorm;
      Ekp1 = Tkp1/(currentOrder+2);
      if (currentOrder == 1)
      {
        if (Tkp1 >= Tkp1_Tk_safety * Tk)
          action = ACTION_MAINTAIN;
        else
          action = ACTION_RAISE;
      }
      else
      {
        if (Tkm1 <= min(Tk,Tkp1))
          action = ACTION_LOWER;
        else if (Tkp1 >= Tk)
          action = ACTION_MAINTAIN;
        else
          action = ACTION_RAISE;
      }
    }
    if (action == ACTION_RAISE)
    {
      currentOrder++;
      Est = Ekp1;
    }
    else if (action == ACTION_LOWER)
    {
      currentOrder--;
      Est = Ekm1;
    }
    newTimeStep = hh;
    rr = pow(r_safety*Est+r_fudge,-1.0/(currentOrder+1.0));
    if (rr >= r_hincr_test)
    {
      rr = r_hincr;
      newTimeStep = rr*hh;
    }
    else if (rr <= 1)
    {
      rr = max(r_min,Xycemin(r_max,rr));
      newTimeStep = rr*hh;
    }
  }
  // 03/22/04 tscoffe:  Note that updating the history has nothing to do with
  // the step-size and everything to do with the newton correction vectors.
  updateHistory();


  // 12/01/05 tscoffe:  This is necessary to avoid currentTimeStep == 0 right
  // before a breakpoint.  So I'm checking to see if currentTime is identically
  // equal to stopTime, in which case we are right before a breakpoint and we
  // should not adjust currentStepSize because that would result in
  // currentStepSize == 0.
  if (currentTime < stopTime)
  {
    // If the step needs to be adjusted:
    if (adjustStep)
    {
      newTimeStep = max(newTimeStep, minTimeStep);
      newTimeStep = min(newTimeStep, maxTimeStep);

      Scalar nextTimePt = currentTime + newTimeStep;

      if (nextTimePt > stopTime)
      {
        nextTimePt  = stopTime;
        newTimeStep = stopTime - currentTime;
      }

      nextTime = nextTimePt;

      hh = newTimeStep;
    }
    else // if time step is constant for this step:
    {
      Scalar nextTimePt = currentTime + hh;

      if (nextTimePt > stopTime)
      {
        nextTimePt      = stopTime;
        hh = stopTime - currentTime;
      }

      nextTime = nextTimePt;
    }
  }

// 11/02/04 tscoffe:  This should be done at the beginning of an integration
//                    stop rather than at the end, so its been moved into
//                    ControlAlgorithm::transientLoop_
  // Update next stop time in StepErrorControl:
  // updateStopTime();
}

template<class Scalar>
bool ImplicitBDFStepper<Scalar>::ErrWtVecSet(Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > &w, Teuchos::RefCountPtr<const Thyra::VectorBase<Scalar> > &y)
{
  abs(&*w,*y);
  V_StVpS(&*w,relErrTol,*w,absErrTol);
  reciprocal(&*w,*w);
  Vt_V(&*w,*w); // We square w because of how weighted norm_2 is computed.
  // divide by N to get RMS norm
  int N = y->size();
  Vt_S(&*w,Scalar(N));
  // Now you can compute WRMS norm as:
  // Scalar WRMSnorm = norm_2(w,y); // WRMS norm of y with respect to weights w.
}

} // namespace Rythmos

#endif //Rythmos_STEPPER_IMPLICITBDF_H
