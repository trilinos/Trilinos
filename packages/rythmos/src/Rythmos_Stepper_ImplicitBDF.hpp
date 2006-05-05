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

    // 05/05/06 tscoffe:  I hate the underscores for private variables!
    Teuchos::RefCountPtr<const Thyra::ModelEvaluator<Scalar> > model;
    Teuchos::RefCountPtr<const Thyra::NonlinearSolverBase<Scalar> > solver;
    Thyra::SingleResidSSDAEModelEvaluator<Scalar>   neModel;

    Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > xn0;
    Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > xpn0;
    std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > > xHistory;
    Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > residual;
    Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > errWtVec;

    Scalar time;


    ScalarMag relErrTol; // relative error tolerance
    ScalarMag absErrTol; // absolute error tolerance
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
    vector>Scalar> beta;     // coefficients used to evaluate predictor from history array
    vector>Scalar> psi;      // $\psi_j(n) = t_n-t_{n-j}$ intermediary variable used to 
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
    Scalar h0_safety;
    Scalar h0_max_factor;
    Scalar h_phase0_incr;
    Scalar h_max_inv;
    Scalar Tkm1_Tk_safety;
    Scalar Tkp1_Tk_safety;
    Scalar r_factor;
    Scalar r_safety;
    Scalar r_fudge;
    Scalar r_min;
    Scalar r_max;
    Scalar r_hincr_test;
    Scalar r_hincr;
    int max_LET_fail;

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
  // Magic Number Defaults:
  h0_safety=Scalar(2.0);
  h0_max_factor=Scalar(0.0001);
  h_phase0_incr=Scalar(2.0);
  h_max_inv=ST::zero();
  Tkm1_Tk_safety=Scalar(2.0);
  Tkp_Tk_safety=Scalar(0.5);
  r_factor=Scalar(0.9);
  r_safety=Scalar(2.0);
  r_fudge=Scalar(0.0001);
  r_min=Scalar(0.125);
  r_max=Scalar(0.9);
  r_hincr_test=Scalar(2.0);
  r_hincr=Scalar(2.0);
  max_LET_fail=Scalar(15);
  // Now we instantiate the model and the solver
  setModel(model);
  setSolver(solver);
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
  solver->solve( neModel, &*x ); // Note that x in input is x_old!
  
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

#ifdef RYTHMOS_STEPPER_DEBUG
  const string dashedline =
    "---------------------------------------------------------------"
    "-------------";
  std::cout << dashedline << std::endl;
  std::cout << "  Rythmos::ImplicitBDFStepper::obtainPredictor" << std::endl;
  std::cout << "nscsco: " << nscsco << std::endl;
  for (int i=0; i<=maxOrder ; ++i)
    std::cout << "beta[" << i << "] = " << beta[i] << std::endl;
  for (int i=0; i<=maxOrder ; ++i)
  {
    std::cout << "xHistory["<< i << "]: " << std::endl;
    dsDae.xHistory[i].printPetraObject();
    std::cout << std::endl;
  }
  for (int i=0; i<=maxOrder ; ++i)
  {
    std::cout << "qHistory["<< i << "]: " << std::endl;
    dsDae.qHistory[i].printPetraObject();
    std::cout << std::endl;
  }
  for (int i=0; i<=maxOrder ; ++i)
  {
    std::cout << " sHistory["<< i << "]: " << std::endl;
    dsDae.sHistory[i].printPetraObject();
    std::cout << std::endl;
  }
  std::cout << " xn0: " << std::endl;
  std::out << xn0Ptr->describe(out,verbLevel,leadingIndent,indentSpacer) << std::endl;
  std::cout << std::endl;
  std::cout << " qn0: " << std::endl;
  std::out << qn0Ptr->describe(out,verbLevel,leadingIndent,indentSpacer) << std::endl;
  std::cout << std::endl;
  std::cout << " qpn0: " << std::endl;
  std::out << qpn0Ptr->describe(out,verbLevel,leadingIndent,indentSpacer) << std::endl;
  std::cout << std::endl;
  std::cout << " sn0: " << std::endl;
  std::out << sn0Ptr->describe(out,verbLevel,leadingIndent,indentSpacer) << std::endl;
  std::cout << std::endl;
  std::cout << " spn0: " << std::endl;
  std::out << spn0Ptr->describe(out,verbLevel,leadingIndent,indentSpacer) << std::endl;
  std::cout << std::endl;
  std::cout << dashedline << std::endl;
#endif // RYTHMOS_STEPPER_DEBUG
  // copy the prediction into the next solution:
  *(dsDae.nextSolutionPtr) = *(dsDae.xn0Ptr);

  return;
}

template<class Scalar>
void ImplicitBDFStepper<Scalar>::obtainResidual()
{
  // output: dsDae.RHSVectorPtr

  // This function returns the following residual:
  // $qpn0 - (alphas/hn)(Q(x)-qn0)+F(x)-B(t)$

  // Note:  dsDae.nextSolutionPtr is used to get Q,F,B in N_TIA_ControlAlgorithm::loadRHS.
  dsDae.RHSVectorPtr->linearCombo(1.0,*dsDae.daeQVectorPtr,-1.0,*dsDae.qn0Ptr);
  dsDae.RHSVectorPtr->linearCombo(1.0,*dsDae.qpn0Ptr,-alphas/sec.currentTimeStep,*dsDae.RHSVectorPtr);

#ifdef RYTHMOS_STEPPER_DEBUG
  const string dashedline =
    "---------------------------------------------------------------"
    "-------------";
  std::cout << dashedline << std::endl;
  std::cout << "  Rythmos::ImplicitBDFStepper::obtainResidual");
  std::cout << " t = " << sec.nextTime << "" << std::endl;
  std::cout << " solution: " << std::endl;
  std::out << nextSolutionPtr->describe(out,verbLevel,leadingIndent,indentSpacer) << std::endl;
  std::cout << " daeQVector: " << std::endl;
  std::out << daeQVectorPtr->describe(out,verbLevel,leadingIndent,indentSpacer) << std::endl;
  std::cout << " qn0: " << std::endl;
  std::out << qn0Ptr->describe(out,verbLevel,leadingIndent,indentSpacer) << std::endl;
  std::cout << " qpn0: " << std::endl;
  std::out << qpn0Ptr->describe(out,verbLevel,leadingIndent,indentSpacer) << std::endl;
  std::cout << " alphas/hn: " << alphas/sec.currentTimeStep << "" << std::endl;
  std::cout << " daeFVector: " << std::endl;
  std::out << daeFVectorPtr->describe(out,verbLevel,leadingIndent,indentSpacer) << std::endl;
  std::cout << " daeBVector: " << std::endl;
  std::out << daeBVectorPtr->describe(out,verbLevel,leadingIndent,indentSpacer) << std::endl;
  std::cout << " dQdt-vector: " << std::endl;
  std::out << RHSVectorPtr->describe(out,verbLevel,leadingIndent,indentSpacer) << std::endl;
  std::cout << std::endl;
#endif

  dsDae.RHSVectorPtr->linearCombo(1.0,*dsDae.RHSVectorPtr,+1.0,*dsDae.daeFVectorPtr);
  dsDae.RHSVectorPtr->linearCombo(1.0,*dsDae.RHSVectorPtr,-1.0,*dsDae.daeBVectorPtr);

#ifdef RYTHMOS_STEPPER_DEBUG
  std::cout << " Residual-vector: " << std::endl;
  std::cout << "-(qpn0-(alpha_s/h)*(Q-qn0)+F-B) " << std::endl;
  std::out << RHSVectorPtr->describe(out,verbLevel,leadingIndent,indentSpacer) << std::endl;
  std::cout << dashedline << std::endl;
#endif

}

template<class Scalar>
void ImplicitBDFStepper<Scalar>::obtainJacobian()
{

#ifdef RYTHMOS_STEPPER_DEBUG
  const string dashedline =
    "---------------------------------------------------------------"
    "-------------";
  std::cout << dashedline << std::endl;
  std::cout << "  Rythmos::ImplicitBDFStepper::obtainJacobian" << std::endl;;
#endif
  // output: dsDae.JMatrixPtr

  // This function returns the following matrix:
  // $-(alphas/hn)dQdx(x)+dFdx$

  // Note:  dsDae.nextSolutionPtr is used to get dQdx,dFdx in N_TIA_ControlAlgorithm::loadJacobian.
//  dsDae.JmatrixPtr->linearCombo(-alphas/sec.currentTimeStep,*dsDae.dQdxMatrixPtr,+1.0,*dsDae.dFdxMatrixPtr);

  N_LAS_Matrix & dQdx = *(dsDae.dQdxMatrixPtr);
  N_LAS_Matrix & dFdx = *(dsDae.dFdxMatrixPtr);
  N_LAS_Matrix & Jac = *(dsDae.JMatrixPtr);

  Jac.add( dQdx );

#ifdef RYTHMOS_STEPPER_DEBUG
  if (tiaParams.debugLevel > 0)
  {
    std::cout << "\n dQdx: \n" <<std::endl;
    dQdx.printPetraObject();
  }
#endif

  Jac.scale( -alphas/sec.currentTimeStep );

#ifdef RYTHMOS_STEPPER_DEBUG
  if (tiaParams.debugLevel > 0)
  {
    double tmp_scalar =  -alphas/sec.currentTimeStep;
    std::cout << "\n scaled dQdx by " << tmp_scalar << " :" <<std::endl;
    Jac.printPetraObject();
  }
#endif

  Jac.add( dFdx );
#ifdef RYTHMOS_STEPPER_DEBUG
  if (tiaParams.debugLevel > 0)
  {
    std::cout << "\n dFdx:" <<std::endl;
    dFdx.printPetraObject();
    std::cout << "\n Total Jacobian:" <<std::endl;
    Jac.printPetraObject();
//    for (int i=0;i<3;++i)
//    {
//      printf("[ %25.20g\t%25.20g\t%25.20g ]\n",Jac[i][0],Jac[i][1],Jac[i][2]);
//    }

    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, dashedline);
    std::cout << std::endl;
  }
#endif

}

template<class Scalar>
void ImplicitBDFStepper<Scalar>::interpolateSolution(double timepoint, 
    	N_LAS_Vector * tmpSolVectorPtr)
{

// 03/23/04 tscoffe:  Currently this code is nearly identical to the IDA code
// for interpolating to an output time.  Either we acknowledge the copyright,
// the list of conditions in the license and the disclaimer or we rewrite this
// function.  The IDA license is included after this routine.
  double tfuzz;   // fuzz factor to check for valid output time
  double tp;      // approximately t{n-1}
  double delt;    // distance between timepoint and currentTime
  double c = 1.0; // coefficient for interpolation
  double gam;     // coefficient for interpolation
  int kord;       // order of interpolation
  double tn = sec.currentTime;
  double hh = sec.currentTimeStep;
  double hused = usedStep;
  int kused = usedOrder;
  double uround = ST::zero();  // unit round-off (set to zero for now)

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

#ifdef RYTHMOS_STEPPER_DEBUG
  const string dashedline =
    "---------------------------------------------------------------"
    "-------------";
  if (tiaParams.debugLevel > 0)
  {
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, "");
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, dashedline);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                           "  N_TIA_BackwardDifferentiation15::updateHistory");
    std::cout << "\n Before updates \n" << std::endl;
    for (int i=0; i<=maxOrder ; ++i)
    {
      std::cout << "\n xHistory["<< i << "]: \n" << std::endl;
      dsDae.xHistory[i].printPetraObject();
      std::cout << std::endl;
    }
    for (int i=0; i<=maxOrder ; ++i)
    {
      std::cout << "\n qHistory["<< i << "]: \n" << std::endl;
      dsDae.qHistory[i].printPetraObject();
      std::cout << std::endl;
    }
    for (int i=0; i<=maxOrder ; ++i)
    {
      std::cout << "\n sHistory["<< i << "]: \n" << std::endl;
      dsDae.sHistory[i].printPetraObject();
      std::cout << std::endl;
    }
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, dashedline);
  }
#endif // RYTHMOS_STEPPER_DEBUG

  // Save Newton correction for potential order increase on next step.
  if (usedOrder < maxOrder)  
  {
    dsDae.xHistory[usedOrder+1] = *dsDae.newtonCorrectionPtr;
    dsDae.qHistory[usedOrder+1] = *dsDae.qNewtonCorrectionPtr;
  }
  // Update history arrays
  dsDae.xHistory[usedOrder].linearCombo(1.0,dsDae.xHistory[usedOrder],1.0,*dsDae.newtonCorrectionPtr);
  dsDae.qHistory[usedOrder].linearCombo(1.0,dsDae.qHistory[usedOrder],1.0,*dsDae.qNewtonCorrectionPtr);
  for (int j=usedOrder-1;j>=0;j--) 
  {
    dsDae.xHistory[j].linearCombo(1.0,dsDae.xHistory[j],1.0,dsDae.xHistory[j+1]);
    dsDae.qHistory[j].linearCombo(1.0,dsDae.qHistory[j],1.0,dsDae.qHistory[j+1]);
  }

  // Update State History
  if (usedOrder < maxOrder)  
  {
    dsDae.sHistory[usedOrder+1] = *dsDae.sNewtonCorrectionPtr;
  }
  // Update history arrays
  dsDae.sHistory[usedOrder].linearCombo(1.0,dsDae.sHistory[usedOrder],1.0,*dsDae.sNewtonCorrectionPtr);
  for (int j=usedOrder-1;j>=0;j--) 
  {
    dsDae.sHistory[j].linearCombo(1.0,dsDae.sHistory[j],1.0,dsDae.sHistory[j+1]);
  }

#ifdef RYTHMOS_STEPPER_DEBUG
  if (tiaParams.debugLevel > 0)
  {
    std::cout << "\n After updates \n" << std::endl;
    std::cout << "\n newtonCorrectionPtr: " << std::endl;
    dsDae.newtonCorrectionPtr->printPetraObject();
    std::cout << "\n qnewtonCorrectionPtr: " << std::endl;
    dsDae.qNewtonCorrectionPtr->printPetraObject();
    for (int i=0; i<=maxOrder ; ++i)
    {
      std::cout << "\n xHistory["<< i << "]: \n" << std::endl;
      dsDae.xHistory[i].printPetraObject();
      std::cout << std::endl;
    }
    for (int i=0; i<=maxOrder ; ++i)
    {
      std::cout << "\n qHistory["<< i << "]: \n" << std::endl;
      dsDae.qHistory[i].printPetraObject();
      std::cout << std::endl;
    }
    std::cout << "\n sNewtonCorrectionPtr: " << std::endl;
    dsDae.sNewtonCorrectionPtr->printPetraObject();
    std::cout << std::endl;
    std::cout << "\n nextStatePtr: " << std::endl;
    dsDae.nextStatePtr->printPetraObject();
    std::cout << std::endl;
    for (int i=0; i<=maxOrder ; ++i)
    {
      std::cout << "\n sHistory["<< i << "]: \n" << std::endl;
      dsDae.sHistory[i].printPetraObject();
      std::cout << std::endl;
    }
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, dashedline);
  }
#endif // RYTHMOS_STEPPER_DEBUG
}

template<class Scalar>
void ImplicitBDFStepper<Scalar>::restoreHistory()
{

  // undo preparation of history array for prediction
  for (int i=nscsco;i<=currentOrder;++i)
  {
    dsDae.xHistory[i].scale(1/beta[i]);
    dsDae.qHistory[i].scale(1/beta[i]);
    dsDae.sHistory[i].scale(1/beta[i]);
  }
  for (int i=1;i<=currentOrder;++i)
  {
    psi[i-1] = psi[i] - (sec.currentTimeStep);
  }
#ifdef RYTHMOS_STEPPER_DEBUG
  const string dashedline =
    "---------------------------------------------------------------"
    "-------------";
  if (tiaParams.debugLevel > 0)
  {
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, "");
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, dashedline);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                           "  N_TIA_BackwardDifferentiation15::restoreHistory");
    for (int i=1;i<=currentOrder;++i)
      std::cout << "\n psi[i] = " << psi[i] << std::endl;
    for (int i=0; i<=maxOrder ; ++i)
    {
      std::cout << "\n xHistory["<< i << "]: \n" << std::endl;
      dsDae.xHistory[i].printPetraObject();
      std::cout << std::endl;
    }
    for (int i=0; i<=maxOrder ; ++i)
    {
      std::cout << "\n qHistory["<< i << "]: \n" << std::endl;
      dsDae.qHistory[i].printPetraObject();
      std::cout << std::endl;
    }
    for (int i=0; i<=maxOrder ; ++i)
    {
    std::cout << "\n sHistory["<< i << "]: \n" << std::endl;
    dsDae.sHistory[i].printPetraObject();
    std::cout << std::endl;
    }
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, dashedline);
  }
#endif // RYTHMOS_STEPPER_DEBUG
} 

template<class Scalar>
void ImplicitBDFStepper<Scalar>::updateCoeffs()
{
  // synchronize with Step Error Control
//  psi[0] = sec.currentTimeStep;
#ifdef RYTHMOS_STEPPER_DEBUG
  const string dashedline =
    "---------------------------------------------------------------"
    "-------------";
  if (tiaParams.debugLevel > 0)
  {
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, "");
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, dashedline);

    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                           "  N_TIA_BackwardDifferentiation15::updateCoeffs");
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
            "  currentTimeStep = ", sec.currentTimeStep);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
            "  numberOfSteps = ", numberOfSteps);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
            "  currentOrder = ", currentOrder);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
            "  nscsco = ", nscsco);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
            "  psi[0] = ", psi[0]);
  }
#endif
  // If the number of steps taken with constant order and constant stepsize is
  // more than the current order + 1 then we don't bother to update the
  // coefficients because we've reached a constant step-size formula.  When
  // this is is not true, then we update the coefficients for the variable
  // step-sizes. 
  if ((sec.currentTimeStep != usedStep) || (currentOrder != usedOrder))
    nscsco = 0;
  nscsco = min(nscsco+1,usedOrder+2);
  if (currentOrder+1 >= nscsco)
  {
    beta[0] = 1.0;
    alpha[0] = 1.0;
    double temp1 = sec.currentTimeStep;
    sigma[0] = 1.0;
    gamma[0] = ST::zero();
    for (int i=1;i<=currentOrder;++i)
    {
      double temp2 = psi[i-1];
      psi[i-1] = temp1;
      beta[i] = beta[i-1]*psi[i-1]/temp2;
      temp1 = temp2 + sec.currentTimeStep;
      alpha[i] = (sec.currentTimeStep)/temp1;
      sigma[i] = (i+1)*sigma[i-1]*alpha[i];
      gamma[i] = gamma[i-1]+alpha[i-1]/(sec.currentTimeStep);
    }
    psi[currentOrder] = temp1;
    alphas = ST::zero();
    alpha0 = ST::zero();
    for (int i=0;i<currentOrder;++i)
    {
      alphas = alphas - 1.0/(i+1.0);
      alpha0 = alpha0 - alpha[i];
    }
    cj = -alphas/(sec.currentTimeStep);
    ck = abs(alpha[currentOrder]+alphas-alpha0);
    ck = max(ck,alpha[currentOrder]);
  }
#ifdef RYTHMOS_STEPPER_DEBUG
  if (tiaParams.debugLevel > 0)
  {
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
            "  nscsco = ", nscsco);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
            "  beta[0] = ", beta[0]);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
            "  beta[1] = ", beta[1]);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
            "  beta[2] = ", beta[2]);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
            "  beta[3] = ", beta[3]);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
            "  beta[4] = ", beta[4]);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
            "  alpha[0] = ", alpha[0]);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
            "  alpha[1] = ", alpha[1]);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
            "  alpha[2] = ", alpha[2]);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
            "  alpha[3] = ", alpha[3]);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
            "  alpha[4] = ", alpha[4]);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
            "  alphas = ", alphas);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
            "  alpha0 = ", alpha0);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
            "  gamma[0] = ", gamma[0]);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
            "  gamma[1] = ", gamma[1]);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
            "  gamma[2] = ", gamma[2]);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
            "  gamma[3] = ", gamma[3]);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
            "  gamma[4] = ", gamma[4]);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
            "  psi[0] = ", psi[0]);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
            "  psi[1] = ", psi[1]);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
            "  psi[2] = ", psi[2]);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
            "  psi[3] = ", psi[3]);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
            "  psi[4] = ", psi[4]);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
            "  sigma[0] = ", sigma[0]);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
            "  sigma[1] = ", sigma[1]);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
            "  sigma[2] = ", sigma[2]);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
            "  sigma[3] = ", sigma[3]);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
            "  sigma[4] = ", sigma[4]);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
            "  ck = ", ck);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, dashedline);
  }
#endif // RYTHMOS_STEPPER_DEBUG
}

template<class Scalar>
void ImplicitBDFStepper<Scalar>::initialize()
{

  // we assume the solution vector is available here 
  // Note that I'm using currSolutionPtr instead of
  // nextSolutionPtr because this is the first step.

  // Update next stop time from StepErrorControl:
  sec.updateStopTime();

  // Choose initial step-size
  double time_to_stop = sec.stopTime - sec.currentTime;
  double currentTimeStep;
  if (tiaParams.constantStepSize)
  {
    currentTimeStep = 0.1 * time_to_stop;
    currentTimeStep = Xycemin(tiaParams.startingTimeStep, currentTimeStep);
  }
  else
  {
    // compute an initial step-size based on rate of change in the solution initially
    double dnorm_q = *(dsDae.qHistory[1].wRMSNorm(*dsDae.qErrWtVecPtr));
    if (dnorm_q > ST::zero())  // time-dependent DAE
    {
      currentTimeStep = Xycemin(h0_max_factor*abs(time_to_stop),sqrt(2.0)/(h0_safety*dnorm_q));
    } 
    else  // non-time-dependent DAE
    {
      currentTimeStep = h0_max_factor*abs(time_to_stop);
    }
    // choose min of user specified value and our value:
    if (tiaParams.startingTimeStep > ST::zero())
      currentTimeStep = Xycemin(tiaParams.startingTimeStep, currentTimeStep);
    // check for maximum step-size:
    double rh = abs(currentTimeStep)*h_max_inv; 
    if (rh>1.0) currentTimeStep = currentTimeStep/rh;
  }
  sec.currentTimeStep = currentTimeStep;

  sec.currentTimeStepRatio = 1.0;
  sec.currentTimeStepSum   = 2.0*sec.currentTimeStep;

  sec.lastTimeStep      = sec.currentTimeStep;
  sec.lastTimeStepRatio = sec.currentTimeStepRatio;
  sec.lastTimeStepSum   = sec.currentTimeStepSum;

  sec.numberSuccessiveFailures = 0;
  sec.stepAttemptStatus        = true;

//  sec.tolAimFac = 0.5;

  sec.nextTime = sec.currentTime + sec.currentTimeStep;

  // x history
  dsDae.xHistory[0] = *(dsDae.currSolutionPtr);
  dsDae.xHistory[1].putScalar(ST::zero()); // no need to multiply by dt here

  // q history
  dsDae.qHistory[0] = *(dsDae.daeQVectorPtr);
  dsDae.qHistory[1].linearCombo(-1.0,*dsDae.daeFVectorPtr,1.0,*dsDae.daeBVectorPtr);
  dsDae.qHistory[1].scale(sec.currentTimeStep);

  // state history
  dsDae.sHistory[0] = *(dsDae.nextStatePtr);
  dsDae.sHistory[1].putScalar(ST::zero()); 

  // Coefficient initialization 
  numberOfSteps = 0;    // number of total time integration steps taken
  currentOrder = 1;
  usedOrder = 1;
  psi[0] = sec.currentTimeStep;
  cj = 1/psi[0];
  nscsco = 0;
#ifdef RYTHMOS_STEPPER_DEBUG
  const string dashedline =
    "---------------------------------------------------------------"
    "-------------";
  if (tiaParams.debugLevel > 0)
  {
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, "");
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, dashedline);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                           "  N_TIA_BackwardDifferentiation15::initialize");
    std::cout << "\n xHistory: \n" << std::endl;
    dsDae.xHistory[0].printPetraObject();
    std::cout << std::endl;
    dsDae.xHistory[1].printPetraObject();
    std::cout << std::endl;
    std::cout << "\n qHistory: \n" << std::endl;
    dsDae.qHistory[0].printPetraObject();
    std::cout << std::endl;
    dsDae.qHistory[1].printPetraObject();
    std::cout << std::endl;
    std::cout << "\n sHistory: \n" << std::endl;
    dsDae.sHistory[0].printPetraObject();
    std::cout << std::endl;
    dsDae.sHistory[1].printPetraObject();
    std::cout << std::endl;
    std::cout << "\n" << "currentTimeStep = " << currentTimeStep << "\n" << std::endl;
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, dashedline);
  }
#endif // RYTHMOS_STEPPER_DEBUG
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

#ifdef RYTHMOS_STEPPER_DEBUG
  const string dashedline =
    "---------------------------------------------------------------"
    "-------------";
  if (tiaParams.debugLevel > 0)
  {
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, "");
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, dashedline);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                           "  N_TIA_BackwardDifferentiation15::checkReduceOrder");
  }
#endif // RYTHMOS_STEPPER_DEBUG

  // 03/11/04 tscoffe:  I only want to run this block after a step has been
  // attempted, but I want to do this regardless of the status from the local
  // error test.
  // 03/10/04 tscoffe:  Decide whether to reduce the order before considering
  // the local error test result.
  double dnorm_x = *(dsDae.newtonCorrectionPtr->wRMSNorm(*dsDae.errWtVecPtr)); // delta = newtonCorrection
  double dnorm_q = *(dsDae.qNewtonCorrectionPtr->wRMSNorm(*dsDae.qErrWtVecPtr)); // dnorm = norm of delta
  double dnorm = sqrt(0.5*dnorm_x*dnorm_x+0.5*dnorm_q*dnorm_q);
  Ek = sigma[currentOrder]*dnorm;
  Tk = (currentOrder+1)*Ek;
  Est = Ek;
  newOrder = currentOrder;
  if (currentOrder>1)
  {
    dsDae.delta_x->linearCombo(1.0,dsDae.xHistory[currentOrder],1.0,*dsDae.newtonCorrectionPtr);
    dnorm_x = *(dsDae.delta_x->wRMSNorm(*dsDae.errWtVecPtr));
    dsDae.delta_q->linearCombo(1.0,dsDae.qHistory[currentOrder],1.0,*dsDae.qNewtonCorrectionPtr);
    dnorm_q = *(dsDae.delta_q->wRMSNorm(*dsDae.qErrWtVecPtr));
    dnorm = sqrt(0.5*dnorm_x*dnorm_x+0.5*dnorm_q*dnorm_q);
    Ekm1 = sigma[currentOrder-1]*dnorm;
    Tkm1 = currentOrder*Ekm1;
    if (currentOrder>2)
    {
      dsDae.delta_x->linearCombo(1.0,dsDae.xHistory[currentOrder-1],1.0,*dsDae.delta_x);
      dnorm_x = *(dsDae.delta_x->wRMSNorm(*dsDae.errWtVecPtr));
      dsDae.delta_q->linearCombo(1.0,dsDae.qHistory[currentOrder-1],1.0,*dsDae.delta_q);
      dnorm_q = *(dsDae.delta_q->wRMSNorm(*dsDae.qErrWtVecPtr));
      dnorm = sqrt(0.5*dnorm_x*dnorm_x+0.5*dnorm_q*dnorm_q);
      Ekm2 = sigma[currentOrder-2]*dnorm;
      Tkm2 = (currentOrder-1)*Ekm2;
      if (Xycemax(Tkm1,Tkm2)<=Tk)
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
#ifdef RYTHMOS_STEPPER_DEBUG
  if (tiaParams.debugLevel > 0)
  {
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                          "  newOrder = ", newOrder);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, dashedline);
  }
#endif
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

// This routine reades but does not change the following variables:
//    stepAttemptStatus, r_factor, r_safety, Est, r_fudge, r_min, r_max,
//    minTimeStep, maxTimeStep, currentTime, stopTime, lastTimeStep


  // First we decide if we'll reduce the order independent of the local error test:
  checkReduceOrder();

#ifdef RYTHMOS_STEPPER_DEBUG
  const string dashedline =
    "---------------------------------------------------------------"
    "-------------";

  if (tiaParams.debugLevel > 0)
  {
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, "");
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, dashedline);

    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                           "  N_TIA_BackwardDifferentiation15::rejectStep");
  }
#endif

  // Only update the time step if we are NOT running constant stepsize.
  bool adjustStep = (!(tiaParams.constantStepSize) );

  sec.lastAttemptedTimeStep = sec.currentTimeStep;


  double newTimeStep = sec.currentTimeStep;
  double rr = 1.0; // step size ratio = new step / old step
  if ((sec.stepAttemptStatus == false) && (adjustStep))
  {
    initialPhase = false;
    nef++;
    restoreHistory();
    // restore psi_
//    for (int i=1;i<=currentOrder;++i)
//      psi[i-1] = psi[i] - sec.currentTimeStep;

    if (nef >= max_LET_fail)  
    {
      string msg = "N_TIA_BackwardDifferentiation15::rejectStep: ";
      msg += "  Maximum number of local error test failures.  ";
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, msg);
    }

    if ((sec.newtonConvergenceStatus <= 0))
    {
      /// 11/11/05 erkeite:  If the Newton solver fails, don't 
      // rely on the error estimate - it may be full of Nan's.
      rr = r_min;
      newTimeStep = rr * sec.currentTimeStep;

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
        rr = Xycemax(r_min,Xycemin(r_max,rr));
        newTimeStep = rr * sec.currentTimeStep;
      }
      else if (nef == 2) // second failure
      {
        rr = r_min;
        newTimeStep = rr * sec.currentTimeStep;
      }
      else if (nef > 2) // third and later failures
      {
        newOrder = 1;
        rr = r_min;
        newTimeStep = rr * sec.currentTimeStep;
      }
    }
    currentOrder = newOrder;
    if (numberOfSteps == 0) // still first step
    {
      psi[0] = newTimeStep;
      dsDae.xHistory[1].scale(rr);
      dsDae.qHistory[1].scale(rr);
    }
#ifdef RYTHMOS_STEPPER_DEBUG
    if (tiaParams.debugLevel > 0)
    {
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
              "  currentTimeStep = ", sec.currentTimeStep);
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
              "  numberOfSteps = ", numberOfSteps);
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
              "  currentOrder = ", currentOrder);
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
              "  nscsco = ", nscsco);
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
              "  alpha[0] = ", alpha[0]);
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
              "  alpha[1] = ", alpha[1]);
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
              "  alpha[2] = ", alpha[2]);
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
              "  alpha[3] = ", alpha[3]);
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
              "  alpha[4] = ", alpha[4]);
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
              "  psi[0] = ", psi[0]);
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
              "  psi[1] = ", psi[1]);
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
              "  psi[2] = ", psi[2]);
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
              "  psi[3] = ", psi[3]);
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
              "  psi[4] = ", psi[4]);
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
              "  sigma[0] = ", sigma[0]);
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
              "  sigma[1] = ", sigma[1]);
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
              "  sigma[2] = ", sigma[2]);
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
              "  sigma[3] = ", sigma[3]);
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
              "  sigma[4] = ", sigma[4]);
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                            "  rr = ", rr);
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                            "  r_factor = ", r_factor);
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                            "  r_safety = ", r_safety);
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                            "  Est = ", Est);
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                            "  r_fudge = ", r_fudge);
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                            "  newOrder = ", newOrder);

      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                            "  currentTimeStep = ", sec.currentTimeStep);
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                            "  newTimeStep = ", newTimeStep);
    }
#endif // RYTHMOS_STEPPER_DEBUG
  }
  else if ((sec.stepAttemptStatus == false) & (!adjustStep))
  {
    string tmp = "  BackwardDifferentiation15:rejectStep: Warning: Local error test failed with constant step-size.\n";
    N_ERH_ErrorMgr::report (N_ERH_ErrorMgr::USR_INFO_0, tmp);
  }

  // If the step needs to be adjusted:
  if (adjustStep)
  {
    newTimeStep = Xycemax(newTimeStep, sec.minTimeStep);
    newTimeStep = Xycemin(newTimeStep, sec.maxTimeStep);

    double nextTimePt = sec.currentTime + newTimeStep;

    if (nextTimePt > sec.stopTime)
    {
      nextTimePt  = sec.stopTime;
      newTimeStep = sec.stopTime - sec.currentTime;
    }

    sec.nextTime = nextTimePt;

    sec.currentTimeStepRatio = newTimeStep/sec.lastTimeStep;
    sec.currentTimeStepSum   = newTimeStep + sec.lastTimeStep;

#ifdef RYTHMOS_STEPPER_DEBUG
    if (tiaParams.debugLevel >0)
    {
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
              "  newTimeStep = ", newTimeStep);
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
              "  nextTime = ", sec.nextTime);
    }
#endif // RYTHMOS_STEPPER_DEBUG

    sec.currentTimeStep = newTimeStep;
  }
  else // if time step is constant for this step:
  {
    double nextTimePt = sec.currentTime + sec.currentTimeStep;

    if (nextTimePt > sec.stopTime)
    {
      nextTimePt      = sec.stopTime;
      sec.currentTimeStep = sec.stopTime - sec.currentTime;
    }

    sec.currentTimeStepRatio = sec.currentTimeStep / sec.lastTimeStep;
    sec.currentTimeStepSum   = sec.currentTimeStep + sec.lastTimeStep;

    sec.nextTime = nextTimePt;
  }
#ifdef RYTHMOS_STEPPER_DEBUG
  if (tiaParams.debugLevel >0)
  {
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, dashedline);
  }
#endif // RYTHMOS_STEPPER_DEBUG
}

template<class Scalar>
void ImplicitBDFStepper<Scalar>::completeStep()
{

  numberOfSteps ++;
  nef = 0;
  sec.lastTime    = sec.currentTime;
  sec.currentTime = sec.nextTime;
  // First we decide if we'll reduce the order independent of the local error test:
  checkReduceOrder();


#ifdef RYTHMOS_STEPPER_DEBUG
  const string dashedline =
    "---------------------------------------------------------------"
    "-------------";

  if (tiaParams.debugLevel > 0)
  {
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, "");
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, dashedline);

    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                           "  N_TIA_BackwardDifferentiation15::completeStep");
  }
#endif

  // Only update the time step if we are NOT running constant stepsize.
  bool adjustStep = (!(tiaParams.constantStepSize) );

  sec.lastAttemptedTimeStep = sec.currentTimeStep;

  double newTimeStep = sec.currentTimeStep;
  double rr = 1.0; // step size ratio = new step / old step
  // 03/11/04 tscoffe:  Here is the block for choosing order & step-size when
  // the local error test PASSES (and Newton succeeded). 
  sec.lastTimeStep = sec.currentTimeStep;
  sec.lastTimeStepRatio = sec.currentTimeStepRatio; // copied from calcTStep1
  sec.lastTimeStepSum   = sec.currentTimeStepSum; // copied from calcTStep1
  int orderDiff = currentOrder - usedOrder;
  usedOrder = currentOrder;
  usedStep = sec.currentTimeStep;
  if ((newOrder == currentOrder-1) || (currentOrder == tiaParams.maxOrder))
  {
    // If we reduced our order or reached max order then move to the next phase
    // of integration where we don't automatically double the step-size and
    // increase the order.
    initialPhase = false;
  }
  if (initialPhase)
  {
    currentOrder++;
    newTimeStep = h_phase0_incr * sec.currentTimeStep;
  }
  else // not in the initial phase of integration
  {
    int action = ACTION_UNSET;
    if (newOrder == currentOrder-1)
      action = ACTION_LOWER;
    else if (newOrder == tiaParams.maxOrder)
      action = ACTION_MAINTAIN;
    else if ((currentOrder+1>=nscsco) || (orderDiff == 1))
    {
      // If we just raised the order last time then we won't raise it again
      // until we've taken currentOrder+1 steps at order currentOrder.
      action = ACTION_MAINTAIN;
    }
    else // consider changing the order 
    {
      dsDae.delta_x->linearCombo(1.0,*dsDae.newtonCorrectionPtr,-1.0,dsDae.xHistory[currentOrder+1]);
      double dnorm_x = *(dsDae.delta_x->wRMSNorm(*dsDae.errWtVecPtr));
      dsDae.delta_q->linearCombo(1.0,*dsDae.qNewtonCorrectionPtr,-1.0,dsDae.qHistory[currentOrder+1]);
      double dnorm_q = *(dsDae.delta_q->wRMSNorm(*dsDae.qErrWtVecPtr));
      double dnorm = sqrt(0.5*dnorm_x*dnorm_x+0.5*dnorm_q*dnorm_q);
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
        if (Tkm1 <= Xycemin(Tk,Tkp1))
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
    newTimeStep = sec.currentTimeStep;
    rr = pow(r_safety*Est+r_fudge,-1.0/(currentOrder+1.0));
#ifdef RYTHMOS_STEPPER_DEBUG
    if (tiaParams.debugLevel > 0)
    {
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                            "  currentOrder = ", currentOrder);
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                            "  r_safety = ", r_safety);
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                            "  r_fudge = ", r_fudge);
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                            "  r_hincr = ", r_hincr);
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                            "  r_hincr_test = ", r_hincr_test);
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                            "  Est = ", Est);
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                            "  raw rr = ", rr);
    }
#endif
    if (rr >= r_hincr_test)
    {
      rr = r_hincr;
      newTimeStep = rr*sec.currentTimeStep;
    }
    else if (rr <= 1)
    {
      rr = Xycemax(r_min,Xycemin(r_max,rr));
      newTimeStep = rr*sec.currentTimeStep;
    }
  }
#ifdef RYTHMOS_STEPPER_DEBUG
    if (tiaParams.debugLevel > 0)
    {
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                            "  initialPhase = ", initialPhase);
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                            "  rr = ", rr);
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                            "  currentTimeStep = ", sec.currentTimeStep);
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                            "  currentTime = ", sec.currentTime);
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                            "  nextTime = ", sec.nextTime);
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                            "  newTimeStep = ", newTimeStep);
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                            "  minTimeStep = ", sec.minTimeStep);
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                            "  maxTimeStep = ", sec.maxTimeStep);
    }
#endif
  // 03/22/04 tscoffe:  Note that updating the history has nothing to do with
  // the step-size and everything to do with the newton correction vectors.
  updateHistory();


  // 12/01/05 tscoffe:  This is necessary to avoid currentTimeStep == 0 right
  // before a breakpoint.  So I'm checking to see if currentTime is identically
  // equal to stopTime, in which case we are right before a breakpoint and we
  // should not adjust currentStepSize because that would result in
  // currentStepSize == 0.
  if (sec.currentTime < sec.stopTime)
  {
    // If the step needs to be adjusted:
    if (adjustStep)
    {
      newTimeStep = Xycemax(newTimeStep, sec.minTimeStep);
      newTimeStep = Xycemin(newTimeStep, sec.maxTimeStep);

      double nextTimePt = sec.currentTime + newTimeStep;

      if (nextTimePt > sec.stopTime)
      {
        nextTimePt  = sec.stopTime;
        newTimeStep = sec.stopTime - sec.currentTime;
      }

      sec.nextTime = nextTimePt;

      sec.currentTimeStepRatio = newTimeStep/sec.lastTimeStep;
      sec.currentTimeStepSum   = newTimeStep + sec.lastTimeStep;

#ifdef RYTHMOS_STEPPER_DEBUG
      if (tiaParams.debugLevel >0)
      {
        N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                "  nextTime = ", sec.nextTime);
        N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                "  newTimeStep = ", newTimeStep);
      }
#endif

      sec.currentTimeStep = newTimeStep;
    }
    else // if time step is constant for this step:
    {
      double nextTimePt = sec.currentTime + sec.currentTimeStep;

      if (nextTimePt > sec.stopTime)
      {
        nextTimePt      = sec.stopTime;
        sec.currentTimeStep = sec.stopTime - sec.currentTime;
      }

      sec.currentTimeStepRatio = sec.currentTimeStep / sec.lastTimeStep;
      sec.currentTimeStepSum   = sec.currentTimeStep + sec.lastTimeStep;

      sec.nextTime = nextTimePt;
    }
  }
#ifdef RYTHMOS_STEPPER_DEBUG
  if (tiaParams.debugLevel >0)
  {
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, dashedline);
  }
#endif

// 11/02/04 tscoffe:  This should be done at the beginning of an integration
//                    stop rather than at the end, so its been moved into
//                    ControlAlgorithm::transientLoop_
  // Update next stop time in StepErrorControl:
  // sec.updateStopTime();
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
}

} // namespace Rythmos

#endif //Rythmos_STEPPER_IMPLICITBDF_H
