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


  private:

    void obtainPredictor();
    void obtainResidual();
    void obtainJacobian();
    void interpolateSolution(Scalar time);
    void updateHistory();
    void restoreHistory();
    int getMaxOrder();
    void updateCoeffs();
    void initialize();
    void checkReduceOrder();
    void rejectStep();
    void completeStep();

    Teuchos::RefCountPtr<const Thyra::ModelEvaluator<Scalar> > model_;
    Teuchos::RefCountPtr<const Thyra::NonlinearSolverBase<Scalar> > solver_;
    Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > x_;
    Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > scaled_x_old_;
    Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > f_;
    Scalar t_;
    Scalar t_old_;

    Thyra::SingleResidSSDAEModelEvaluator<Scalar>   neModel_;

    // Magic Numbers
    currentOrder_; // Current order of integration
    oldOrder_;     // previous order of integration
    maxOrder_;     // maximum order = min(5,user option maxord) - see below.
    usedOrder_;    // order used in current step (used after currentOrder_ is updated)
    alphas_;    // $\alpha_s$ fixed-leading coefficient of this BDF method
    alpha_;    // $\alpha_j(n)=h_n/\psi_j(n)$ coefficient used in local error test
                      // note:   $h_n$ = current step size, n = current time step
    alpha0_;     // $-\sum_{j=1}^k \alpha_j(n)$ coefficient used in local error test
    cj_ ;        // $-\alpha_s/h_n$ coefficient used in local error test
    ck_ ;        // local error coefficient gamma_[0] = 0; // $\gamma_j(n)=\sum_{l=1}^{j-1}1/\psi_l(n)$ coefficient used to
    gamma_;    // calculate time derivative of history array for predictor 
    beta_;     // coefficients used to evaluate predictor from history array
    psi_;      // $\psi_j(n) = t_n-t_{n-j}$ intermediary variable used to 
                      // compute $\beta_j(n)$
    sigma_;    // $\sigma_j(n) = \frac{h_n^j(j-1)!}{\psi_1(n)*\cdots *\psi_j(n)}$
    numberOfSteps_;// number of total time integration steps taken
    nef_;
    usedStep_;
    nscsco_;
    Ek_;
    Ekm1_;
    Ekm2_;
    Ekp1_;
    Est_;
    Tk_;
    Tkm1_;
    Tkm2_;
    Tkp1_;
    newOrder_;
    initialPhase_;
    h0_safety_;
    h0_max_factor_;
    h_phase0_incr_;
    h_max_inv_;
    Tkm1_Tk_safety_;
    Tkp1_Tk_safety_;
    r_factor_;
    r_safety_;
    r_fudge_;
    r_min_;
    r_max_;
    r_hincr_test_;
    r_hincr_;
    max_LET_fail_;

};

// ////////////////////////////
// Defintions

template<class Scalar>
ImplicitBDFStepper<Scalar>::ImplicitBDFStepper(
  const Teuchos::RefCountPtr<const Thyra::ModelEvaluator<Scalar> > &model
  ,const Teuchos::RefCountPtr<const Thyra::NonlinearSolverBase<Scalar> > &solver
  )
{
  setModel(model);
  setSolver(solver);
  // Magic Number Defaults:
  currentOrder=1; // Current order of integration
  oldOrder=1; // previous order of integration
  maxOrder=5;  // maximum order = min(5,user option maxord) - see below.
  usedOrder=1;  // order used in current step (used after currentOrder_ is updated)
  alphas=-1.0;  // $\alpha_s$ fixed-leading coefficient of this BDF method
  alpha=[0]=0.0;  // $\alpha_j(n)=h_n/\psi_j(n)$ coefficient used in local error test
                  // note:   $h_n$ = current step size, n = current time step
  alpha=[1]=0.0;
  alpha=[2]=0.0;
  alpha=[3]=0.0;
  alpha=[4]=0.0;
  alpha=[5]=0.0;
  alpha0=0.0;   // $-\sum_{j=1}^k \alpha_j(n)$ coefficient used in local error test
  cj=0.0;      // $-\alpha_s/h_n$ coefficient used in local error test
  ck=0.0;      // local error coefficient gamma_[0] = 0; // $\gamma_j(n)=\sum_{l=1}^{j-1}1/\psi_l(n)$ coefficient used to
  gamma=[0]=0.0;  // calculate time derivative of history array for predictor 
  gamma=[1]=0.0;
  gamma=[2]=0.0;
  gamma=[3]=0.0;
  gamma=[4]=0.0;
  gamma=[5]=0.0;
  beta=[0]=0.0;   // coefficients used to evaluate predictor from history array
  beta=[1]=0.0;
  beta=[2]=0.0;
  beta=[3]=0.0;
  beta=[4]=0.0;
  beta=[5]=0.0;
  psi=[0]=0.0;    // $\psi_j(n) = t_n-t_{n-j}$ intermediary variable used to 
                  // compute $\beta_j(n;$
  psi=[1]=0.0;
  psi=[2]=0.0;
  psi=[3]=0.0;
  psi=[4]=0.0;
  psi=[5]=0.0;
  sigma=[0]=0.0;  // $\sigma_j(n) = \frac{h_n^j(j-1)!}{\psi_1(n)*\cdots *\psi_j(n)}$
  sigma=[1]=0.0;
  sigma=[2]=0.0;
  sigma=[3]=0.0;
  sigma=[4]=0.0;
  sigma=[5]=0.0;
  numberOfSteps=0;   // number of total time integration steps taken
  nef=0;
  usedStep=0.0;
  nscsco=0;
  Ek=0.0;
  Ekm1=0.0;
  Ekm2=0.0;
  Ekp1=0.0;
  Est=0.0;
  Tk=0.0;
  Tkm1=0.0;
  Tkm2=0.0;
  Tkp1=0.0;
  newOrder=1;
  initialPhase=true;
  h0_safety=2.0;
  h0_max_factor=0.0001;
  h_phase0_incr=2.0;
  h_max_inv=0.0;
  Tkm1_Tk_safety=2.0;
  Tkp1_Tk_safety=0.5;
  r_factor=0.9;
  r_safety=2.0;
  r_fudge=0.0001;
  r_min=0.125;
  r_max=0.9;
  r_hincr_test=2.0;
  r_hincr=2.0;
  max_LET_fail=15;
}

template<class Scalar>
void ImplicitBDFStepper<Scalar>::setModel(const Teuchos::RefCountPtr<const Thyra::ModelEvaluator<Scalar> > &model)
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  model_ = model;
  t_ = ST::zero();
  x_ = model_->get_x_init()->clone_v();
  f_ = Thyra::createMember(model_->get_f_space());

  scaled_x_old_ = x_->clone_v();
  
}

template<class Scalar>
void ImplicitBDFStepper<Scalar>::setSolver(const Teuchos::RefCountPtr<const Thyra::NonlinearSolverBase<Scalar> > &solver)
{
  solver_ = solver;
}

template<class Scalar>
Scalar ImplicitBDFStepper<Scalar>::TakeStep()
{
  // print something out about this method not supporting automatic variable step-size
  typedef Teuchos::ScalarTraits<Scalar> ST;
  return(-ST::one());
}

template<class Scalar>
Scalar ImplicitBDFStepper<Scalar>::TakeStep(Scalar dt)
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  //
  // Setup the nonlinear equations:
  //
  //   f( (1/dt)* x + (-1/dt)*x_old), x, t ) = 0
  //
  V_StV( &*scaled_x_old_, Scalar(-ST::one()/dt), *x_ );
  t_old_ = t_;
  neModel_.initialize(model_,Scalar(ST::one()/dt),scaled_x_old_,ST::one(),Teuchos::null,t_old_+dt,Teuchos::null);
  //
  // Solve the implicit nonlinear system to a tolerance of ???
  //
  solver_->solve( neModel_, &*x_ ); // Note that x in input is x_old!
  //
  // Update the step
  //
  t_ += dt;

  return(dt);
}

template<class Scalar>
Teuchos::RefCountPtr<const Thyra::VectorBase<Scalar> > ImplicitBDFStepper<Scalar>::get_solution() const
{
  return(x_);
}

template<class Scalar>
Teuchos::RefCountPtr<const Thyra::VectorBase<Scalar> > ImplicitBDFStepper<Scalar>::get_residual() const
{
  return(f_);
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
    out << "model_ = " << std::endl;
    out << model_->describe(out,verbLevel,leadingIndent,indentSpacer) << std::endl;
    out << "solver_ = " << std::endl;
    out << solver_->describe(out,verbLevel,leadingIndent,indentSpacer) << std::endl;
    out << "x_ = " << std::endl;
    out << x_->describe(out,verbLevel,leadingIndent,indentSpacer) << std::endl;
    out << "scaled_x_old_ = " << std::endl;
    out << scaled_x_old_->describe(out,verbLevel,leadingIndent,indentSpacer) << std::endl;
    out << "f_ = " << std::endl;
    out << f_->describe(out,verbLevel,leadingIndent,indentSpacer) << std::endl;
    out << "t_ = " << t_ << std::endl;
    out << "t_old_ = " << t_old_ << std::endl;
//    out << "neModel_ = " << std::endl;
//    out << neModel_->describe(out,verbLevel,leadingIndent,indentSpacer) << std::endl;
  }
  return(out);
}

template<class Scalar>
void ImplicitBDFStepper<Scalar>
template<class Scalar>
void ImplicitBDFStepper<Scalar>::obtainPredictor()
{
  // prepare history array for prediction
  for (int i=nscsco_;i<=currentOrder_;++i)
  {
    dsDae.xHistory[i].scale(beta_[i]);
    dsDae.qHistory[i].scale(beta_[i]);
    dsDae.sHistory[i].scale(beta_[i]);
  }
  
  // evaluate predictor
  *dsDae.xn0Ptr = dsDae.xHistory[0];
  *dsDae.qn0Ptr = dsDae.qHistory[0];
  *dsDae.sn0Ptr = dsDae.sHistory[0];
  dsDae.qpn0Ptr->putScalar(0.0);
  dsDae.spn0Ptr->putScalar(0.0);
  for (int i=1;i<=currentOrder_;++i)
  {
    dsDae.xn0Ptr->linearCombo(1.0,dsDae.xHistory[i],1.0,*dsDae.xn0Ptr);
    dsDae.qn0Ptr->linearCombo(1.0,dsDae.qHistory[i],1.0,*dsDae.qn0Ptr);
    dsDae.sn0Ptr->linearCombo(1.0,dsDae.sHistory[i],1.0,*dsDae.sn0Ptr);
    dsDae.qpn0Ptr->linearCombo(gamma_[i],dsDae.qHistory[i],1.0,*dsDae.qpn0Ptr);
    dsDae.spn0Ptr->linearCombo(gamma_[i],dsDae.sHistory[i],1.0,*dsDae.spn0Ptr);
  }

#ifdef Xyce_DEBUG_TIME
  const string dashedline =
    "---------------------------------------------------------------"
    "-------------";
  if (tiaParams.debugLevel > 0)
  {
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, "");
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, dashedline);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                           "  N_TIA_BackwardDifferentiation15::obtainPredictor");
    cout << "\n nscsco_: " << nscsco_ << "\n" << endl;
    for (int i=0; i<=maxOrder_ ; ++i)
      cout << "\n beta_[" << i << "] = " << beta_[i] << "\n" << endl;
    for (int i=0; i<=maxOrder_ ; ++i)
    {
      cout << "\n xHistory["<< i << "]: \n" << endl;
      dsDae.xHistory[i].printPetraObject();
      cout << endl;
    }
    for (int i=0; i<=maxOrder_ ; ++i)
    {
      cout << "\n qHistory["<< i << "]: \n" << endl;
      dsDae.qHistory[i].printPetraObject();
      cout << endl;
    }
    for (int i=0; i<=maxOrder_ ; ++i)
    {
      cout << "\n sHistory["<< i << "]: \n" << endl;
      dsDae.sHistory[i].printPetraObject();
      cout << endl;
    }
    cout << "\n xn0: \n" << endl;
    dsDae.xn0Ptr->printPetraObject();
    cout << endl;
    cout << "\n qn0: \n" << endl;
    dsDae.qn0Ptr->printPetraObject();
    cout << endl;
    cout << "\n qpn0: \n" << endl;
    dsDae.qpn0Ptr->printPetraObject();
    cout << endl;
    cout << "\n sn0: \n" << endl;
    dsDae.sn0Ptr->printPetraObject();
    cout << endl;
    cout << "\n spn0: \n" << endl;
    dsDae.spn0Ptr->printPetraObject();
    cout << endl;
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, dashedline);
  }
#endif // Xyce_DEBUG_TIME
  // copy the prediction into the next solution:
  *(dsDae.nextSolutionPtr) = *(dsDae.xn0Ptr);

  return;
}

template<class Scalar>
void ImplicitBDFStepper<Scalar>::obtainResidual()
{
  // output: dsDae.RHSVectorPtr

  // This function returns the following residual:
  // $qpn0 - (alphas_/hn)(Q(x)-qn0)+F(x)-B(t)$

  // Note:  dsDae.nextSolutionPtr is used to get Q,F,B in N_TIA_ControlAlgorithm::loadRHS.
  dsDae.RHSVectorPtr->linearCombo(1.0,*dsDae.daeQVectorPtr,-1.0,*dsDae.qn0Ptr);
  dsDae.RHSVectorPtr->linearCombo(1.0,*dsDae.qpn0Ptr,-alphas_/sec.currentTimeStep,*dsDae.RHSVectorPtr);

#ifdef Xyce_DEBUG_TIME
  const string dashedline =
    "---------------------------------------------------------------"
    "-------------";
  if (tiaParams.debugLevel > 0)
  {
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, "");
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, dashedline);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                           "  N_TIA_BackwardDifferentiation15::obtainResidual");
    cout << "\n t = " << sec.nextTime << "\n" << endl;
    cout << "\n solution: \n" << endl;
    dsDae.nextSolutionPtr->printPetraObject();
    cout << "\n daeQVector: \n" << endl;
    dsDae.daeQVectorPtr->printPetraObject();
    cout << "\n qn0: \n" << endl;
    dsDae.qn0Ptr->printPetraObject();
    cout << "\n qpn0: \n" << endl;
    dsDae.qpn0Ptr->printPetraObject();
    cout << "\n alphas_/hn: " << alphas_/sec.currentTimeStep << "\n" << endl;
    cout << "\n daeFVector: \n" << endl;
    dsDae.daeFVectorPtr->printPetraObject();
    cout << "\n daeBVector: \n" << endl;
    dsDae.daeBVectorPtr->printPetraObject();
    cout << "\n dQdt-vector: \n" << endl;
    dsDae.RHSVectorPtr->printPetraObject();
    cout << endl;
  }
#endif

  dsDae.RHSVectorPtr->linearCombo(1.0,*dsDae.RHSVectorPtr,+1.0,*dsDae.daeFVectorPtr);
  dsDae.RHSVectorPtr->linearCombo(1.0,*dsDae.RHSVectorPtr,-1.0,*dsDae.daeBVectorPtr);

  // since the nonlinear solver is expecting a -f, scale by -1.0:
  dsDae.RHSVectorPtr->scale(-1.0);

  // if voltage limiting is on, add it in:
  if (dsDae.limiterFlag)
  {
    (dsDae.dQdxdVpVectorPtr)->scale( -alphas_/sec.currentTimeStep );

    (dsDae.RHSVectorPtr)->daxpy(
      *(dsDae.RHSVectorPtr), +1.0, *(dsDae.dQdxdVpVectorPtr));

    (dsDae.RHSVectorPtr)->daxpy(
      *(dsDae.RHSVectorPtr), +1.0, *(dsDae.dFdxdVpVectorPtr));
  }

#ifdef Xyce_DEBUG_TIME
  if (tiaParams.debugLevel > 0)
  {
    cout << "\n Residual-vector: \n" << endl;
    cout << "-(qpn0-(alpha_s/h)*(Q-qn0)+F-B) \n" << endl;
    dsDae.RHSVectorPtr->printPetraObject();
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, dashedline);
    cout << endl;
  }
#endif

}

template<class Scalar>
void ImplicitBDFStepper<Scalar>::obtainJacobian()
{

#ifdef Xyce_DEBUG_TIME
  const string dashedline =
    "---------------------------------------------------------------"
    "-------------";
  if (tiaParams.debugLevel > 0)
  {
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, "");
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, dashedline);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                           "  N_TIA_BackwardDifferentiation15::obtainJacobian");
  }
#endif
  // output: dsDae.JMatrixPtr

  // This function returns the following matrix:
  // $-(alphas_/hn)dQdx(x)+dFdx$

  // Note:  dsDae.nextSolutionPtr is used to get dQdx,dFdx in N_TIA_ControlAlgorithm::loadJacobian.
//  dsDae.JmatrixPtr->linearCombo(-alphas_/sec.currentTimeStep,*dsDae.dQdxMatrixPtr,+1.0,*dsDae.dFdxMatrixPtr);

  N_LAS_Matrix & dQdx = *(dsDae.dQdxMatrixPtr);
  N_LAS_Matrix & dFdx = *(dsDae.dFdxMatrixPtr);
  N_LAS_Matrix & Jac = *(dsDae.JMatrixPtr);

  Jac.add( dQdx );

#ifdef Xyce_DEBUG_TIME
  if (tiaParams.debugLevel > 0)
  {
    cout << "\n dQdx: \n" <<endl;
    dQdx.printPetraObject();
  }
#endif

  Jac.scale( -alphas_/sec.currentTimeStep );

#ifdef Xyce_DEBUG_TIME
  if (tiaParams.debugLevel > 0)
  {
    double tmp_scalar =  -alphas_/sec.currentTimeStep;
    cout << "\n scaled dQdx by " << tmp_scalar << " :" <<endl;
    Jac.printPetraObject();
  }
#endif

  Jac.add( dFdx );
#ifdef Xyce_DEBUG_TIME
  if (tiaParams.debugLevel > 0)
  {
    cout << "\n dFdx:" <<endl;
    dFdx.printPetraObject();
    cout << "\n Total Jacobian:" <<endl;
    Jac.printPetraObject();
//    for (int i=0;i<3;++i)
//    {
//      printf("[ %25.20g\t%25.20g\t%25.20g ]\n",Jac[i][0],Jac[i][1],Jac[i][2]);
//    }

    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, dashedline);
    cout << endl;
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
  double tp;      // approximately t_{n-1}
  double delt;    // distance between timepoint and currentTime
  double c = 1.0; // coefficient for interpolation
  double gam;     // coefficient for interpolation
  int kord;       // order of interpolation
  double tn = sec.currentTime;
  double hh = sec.currentTimeStep;
  double hused = usedStep_;
  int kused = usedOrder_;
  double uround = 0.0;  // unit round-off (set to zero for now)

  tfuzz = 100 * uround * (tn + hh);
  tp = tn - hused - tfuzz;
  if ( (timepoint - tp)*hh < 0.0 ) 
    return false;

  *tmpSolVectorPtr = dsDae.xHistory[0];
  kord = kused;
  if ( (kused == 0) || (timepoint == tn) ) 
    kord = 1;

  delt = timepoint - tn;
  gam = delt/psi_[0];
  for (int j=1 ; j <= kord ; ++j)
  {
    c = c*gam;
    gam = (delt + psi_[j-1])/psi_[j];
    tmpSolVectorPtr->linearCombo(1.0,*tmpSolVectorPtr,c,dsDae.xHistory[j]);
  }
  return true;
}

template<class Scalar>
void ImplicitBDFStepper<Scalar>::updateHistory()
{

#ifdef Xyce_DEBUG_TIME
  const string dashedline =
    "---------------------------------------------------------------"
    "-------------";
  if (tiaParams.debugLevel > 0)
  {
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, "");
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, dashedline);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                           "  N_TIA_BackwardDifferentiation15::updateHistory");
    cout << "\n Before updates \n" << endl;
    for (int i=0; i<=maxOrder_ ; ++i)
    {
      cout << "\n xHistory["<< i << "]: \n" << endl;
      dsDae.xHistory[i].printPetraObject();
      cout << endl;
    }
    for (int i=0; i<=maxOrder_ ; ++i)
    {
      cout << "\n qHistory["<< i << "]: \n" << endl;
      dsDae.qHistory[i].printPetraObject();
      cout << endl;
    }
    for (int i=0; i<=maxOrder_ ; ++i)
    {
      cout << "\n sHistory["<< i << "]: \n" << endl;
      dsDae.sHistory[i].printPetraObject();
      cout << endl;
    }
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, dashedline);
  }
#endif // Xyce_DEBUG_TIME

  // Save Newton correction for potential order increase on next step.
  if (usedOrder_ < maxOrder_)  
  {
    dsDae.xHistory[usedOrder_+1] = *dsDae.newtonCorrectionPtr;
    dsDae.qHistory[usedOrder_+1] = *dsDae.qNewtonCorrectionPtr;
  }
  // Update history arrays
  dsDae.xHistory[usedOrder_].linearCombo(1.0,dsDae.xHistory[usedOrder_],1.0,*dsDae.newtonCorrectionPtr);
  dsDae.qHistory[usedOrder_].linearCombo(1.0,dsDae.qHistory[usedOrder_],1.0,*dsDae.qNewtonCorrectionPtr);
  for (int j=usedOrder_-1;j>=0;j--) 
  {
    dsDae.xHistory[j].linearCombo(1.0,dsDae.xHistory[j],1.0,dsDae.xHistory[j+1]);
    dsDae.qHistory[j].linearCombo(1.0,dsDae.qHistory[j],1.0,dsDae.qHistory[j+1]);
  }

  // Update State History
  if (usedOrder_ < maxOrder_)  
  {
    dsDae.sHistory[usedOrder_+1] = *dsDae.sNewtonCorrectionPtr;
  }
  // Update history arrays
  dsDae.sHistory[usedOrder_].linearCombo(1.0,dsDae.sHistory[usedOrder_],1.0,*dsDae.sNewtonCorrectionPtr);
  for (int j=usedOrder_-1;j>=0;j--) 
  {
    dsDae.sHistory[j].linearCombo(1.0,dsDae.sHistory[j],1.0,dsDae.sHistory[j+1]);
  }

#ifdef Xyce_DEBUG_TIME
  if (tiaParams.debugLevel > 0)
  {
    cout << "\n After updates \n" << endl;
    cout << "\n newtonCorrectionPtr: " << endl;
    dsDae.newtonCorrectionPtr->printPetraObject();
    cout << "\n qnewtonCorrectionPtr: " << endl;
    dsDae.qNewtonCorrectionPtr->printPetraObject();
    for (int i=0; i<=maxOrder_ ; ++i)
    {
      cout << "\n xHistory["<< i << "]: \n" << endl;
      dsDae.xHistory[i].printPetraObject();
      cout << endl;
    }
    for (int i=0; i<=maxOrder_ ; ++i)
    {
      cout << "\n qHistory["<< i << "]: \n" << endl;
      dsDae.qHistory[i].printPetraObject();
      cout << endl;
    }
    cout << "\n sNewtonCorrectionPtr: " << endl;
    dsDae.sNewtonCorrectionPtr->printPetraObject();
    cout << endl;
    cout << "\n nextStatePtr: " << endl;
    dsDae.nextStatePtr->printPetraObject();
    cout << endl;
    for (int i=0; i<=maxOrder_ ; ++i)
    {
      cout << "\n sHistory["<< i << "]: \n" << endl;
      dsDae.sHistory[i].printPetraObject();
      cout << endl;
    }
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, dashedline);
  }
#endif // Xyce_DEBUG_TIME
}

template<class Scalar>
void ImplicitBDFStepper<Scalar>::restoreHistory()
{

  // undo preparation of history array for prediction
  for (int i=nscsco_;i<=currentOrder_;++i)
  {
    dsDae.xHistory[i].scale(1/beta_[i]);
    dsDae.qHistory[i].scale(1/beta_[i]);
    dsDae.sHistory[i].scale(1/beta_[i]);
  }
  for (int i=1;i<=currentOrder_;++i)
  {
    psi_[i-1] = psi_[i] - (sec.currentTimeStep);
  }
#ifdef Xyce_DEBUG_TIME
  const string dashedline =
    "---------------------------------------------------------------"
    "-------------";
  if (tiaParams.debugLevel > 0)
  {
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, "");
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, dashedline);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                           "  N_TIA_BackwardDifferentiation15::restoreHistory");
    for (int i=1;i<=currentOrder_;++i)
      cout << "\n psi_[i] = " << psi_[i] << endl;
    for (int i=0; i<=maxOrder_ ; ++i)
    {
      cout << "\n xHistory["<< i << "]: \n" << endl;
      dsDae.xHistory[i].printPetraObject();
      cout << endl;
    }
    for (int i=0; i<=maxOrder_ ; ++i)
    {
      cout << "\n qHistory["<< i << "]: \n" << endl;
      dsDae.qHistory[i].printPetraObject();
      cout << endl;
    }
    for (int i=0; i<=maxOrder_ ; ++i)
    {
    cout << "\n sHistory["<< i << "]: \n" << endl;
    dsDae.sHistory[i].printPetraObject();
    cout << endl;
    }
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, dashedline);
  }
#endif // Xyce_DEBUG_TIME
} 

template<class Scalar>
int ImplicitBDFStepper<Scalar>::getMaxOrder()
{
  return maxOrder_;
}

template<class Scalar>
void ImplicitBDFStepper<Scalar>::updateCoeffs()
{
  // synchronize with Step Error Control
//  psi_[0] = sec.currentTimeStep;
#ifdef Xyce_DEBUG_TIME
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
            "  numberOfSteps_ = ", numberOfSteps_);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
            "  currentOrder_ = ", currentOrder_);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
            "  nscsco_ = ", nscsco_);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
            "  psi_[0] = ", psi_[0]);
  }
#endif
  // If the number of steps taken with constant order and constant stepsize is
  // more than the current order + 1 then we don't bother to update the
  // coefficients because we've reached a constant step-size formula.  When
  // this is is not true, then we update the coefficients for the variable
  // step-sizes. 
  if ((sec.currentTimeStep != usedStep_) || (currentOrder_ != usedOrder_))
    nscsco_ = 0;
  nscsco_ = min(nscsco_+1,usedOrder_+2);
  if (currentOrder_+1 >= nscsco_)
  {
    beta_[0] = 1.0;
    alpha_[0] = 1.0;
    double temp1 = sec.currentTimeStep;
    sigma_[0] = 1.0;
    gamma_[0] = 0.0;
    for (int i=1;i<=currentOrder_;++i)
    {
      double temp2 = psi_[i-1];
      psi_[i-1] = temp1;
      beta_[i] = beta_[i-1]*psi_[i-1]/temp2;
      temp1 = temp2 + sec.currentTimeStep;
      alpha_[i] = (sec.currentTimeStep)/temp1;
      sigma_[i] = (i+1)*sigma_[i-1]*alpha_[i];
      gamma_[i] = gamma_[i-1]+alpha_[i-1]/(sec.currentTimeStep);
    }
    psi_[currentOrder_] = temp1;
    alphas_ = 0.0;
    alpha0_ = 0.0;
    for (int i=0;i<currentOrder_;++i)
    {
      alphas_ = alphas_ - 1.0/(i+1.0);
      alpha0_ = alpha0_ - alpha_[i];
    }
    cj_ = -alphas_/(sec.currentTimeStep);
    ck_ = abs(alpha_[currentOrder_]+alphas_-alpha0_);
    ck_ = max(ck_,alpha_[currentOrder_]);
  }
#ifdef Xyce_DEBUG_TIME
  if (tiaParams.debugLevel > 0)
  {
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
            "  nscsco_ = ", nscsco_);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
            "  beta_[0] = ", beta_[0]);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
            "  beta_[1] = ", beta_[1]);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
            "  beta_[2] = ", beta_[2]);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
            "  beta_[3] = ", beta_[3]);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
            "  beta_[4] = ", beta_[4]);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
            "  alpha_[0] = ", alpha_[0]);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
            "  alpha_[1] = ", alpha_[1]);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
            "  alpha_[2] = ", alpha_[2]);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
            "  alpha_[3] = ", alpha_[3]);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
            "  alpha_[4] = ", alpha_[4]);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
            "  alphas_ = ", alphas_);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
            "  alpha0_ = ", alpha0_);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
            "  gamma_[0] = ", gamma_[0]);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
            "  gamma_[1] = ", gamma_[1]);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
            "  gamma_[2] = ", gamma_[2]);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
            "  gamma_[3] = ", gamma_[3]);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
            "  gamma_[4] = ", gamma_[4]);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
            "  psi_[0] = ", psi_[0]);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
            "  psi_[1] = ", psi_[1]);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
            "  psi_[2] = ", psi_[2]);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
            "  psi_[3] = ", psi_[3]);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
            "  psi_[4] = ", psi_[4]);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
            "  sigma_[0] = ", sigma_[0]);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
            "  sigma_[1] = ", sigma_[1]);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
            "  sigma_[2] = ", sigma_[2]);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
            "  sigma_[3] = ", sigma_[3]);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
            "  sigma_[4] = ", sigma_[4]);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
            "  ck_ = ", ck_);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, dashedline);
  }
#endif // Xyce_DEBUG_TIME
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
    if (dnorm_q > 0.0)  // time-dependent DAE
    {
      currentTimeStep = Xycemin(h0_max_factor_*abs(time_to_stop),sqrt(2.0)/(h0_safety_*dnorm_q));
    } 
    else  // non-time-dependent DAE
    {
      currentTimeStep = h0_max_factor_*abs(time_to_stop);
    }
    // choose min of user specified value and our value:
    if (tiaParams.startingTimeStep > 0.0)
      currentTimeStep = Xycemin(tiaParams.startingTimeStep, currentTimeStep);
    // check for maximum step-size:
    double rh = abs(currentTimeStep)*h_max_inv_; 
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

//  sec.tolAimFac_ = 0.5;

  sec.nextTime = sec.currentTime + sec.currentTimeStep;

  // x history
  dsDae.xHistory[0] = *(dsDae.currSolutionPtr);
  dsDae.xHistory[1].putScalar(0.0); // no need to multiply by dt here

  // q history
  dsDae.qHistory[0] = *(dsDae.daeQVectorPtr);
  dsDae.qHistory[1].linearCombo(-1.0,*dsDae.daeFVectorPtr,1.0,*dsDae.daeBVectorPtr);
  dsDae.qHistory[1].scale(sec.currentTimeStep);

  // state history
  dsDae.sHistory[0] = *(dsDae.nextStatePtr);
  dsDae.sHistory[1].putScalar(0.0); 

  // Coefficient initialization 
  numberOfSteps_ = 0;    // number of total time integration steps taken
  currentOrder_ = 1;
  usedOrder_ = 1;
  psi_[0] = sec.currentTimeStep;
  cj_ = 1/psi_[0];
  nscsco_ = 0;
#ifdef Xyce_DEBUG_TIME
  const string dashedline =
    "---------------------------------------------------------------"
    "-------------";
  if (tiaParams.debugLevel > 0)
  {
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, "");
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, dashedline);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                           "  N_TIA_BackwardDifferentiation15::initialize");
    cout << "\n xHistory: \n" << endl;
    dsDae.xHistory[0].printPetraObject();
    cout << endl;
    dsDae.xHistory[1].printPetraObject();
    cout << endl;
    cout << "\n qHistory: \n" << endl;
    dsDae.qHistory[0].printPetraObject();
    cout << endl;
    dsDae.qHistory[1].printPetraObject();
    cout << endl;
    cout << "\n sHistory: \n" << endl;
    dsDae.sHistory[0].printPetraObject();
    cout << endl;
    dsDae.sHistory[1].printPetraObject();
    cout << endl;
    cout << "\n" << "currentTimeStep = " << currentTimeStep << "\n" << endl;
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, dashedline);
  }
#endif // Xyce_DEBUG_TIME
}

template<class Scalar>
void ImplicitBDFStepper<Scalar>::checkReduceOrder()
{

// This routine puts its output in newOrder_

// This routine changes the following variables:
//    Ek_, Tk_, Est_, newOrder_, dsDae.delta_x, dsDae.delta_q,
//    Ekm1_, Tkm1_, Ekm2_, Tkm2_ 

// This routine reads but does not change the following variables:
//    currentOrder_, sigma_, dsDae.newtonCorrectionPtr, dsDae.qNewtonCorrectionPtr,
//    dsDae.errWtVecPtr, dsDae.qErrWtVecPtr, dsDae.xHistory, dsDae.qHistory

#ifdef Xyce_DEBUG_TIME
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
#endif // Xyce_DEBUG_TIME

  // 03/11/04 tscoffe:  I only want to run this block after a step has been
  // attempted, but I want to do this regardless of the status from the local
  // error test.
  // 03/10/04 tscoffe:  Decide whether to reduce the order before considering
  // the local error test result.
  double dnorm_x = *(dsDae.newtonCorrectionPtr->wRMSNorm(*dsDae.errWtVecPtr)); // delta = newtonCorrection
  double dnorm_q = *(dsDae.qNewtonCorrectionPtr->wRMSNorm(*dsDae.qErrWtVecPtr)); // dnorm = norm of delta
  double dnorm = sqrt(0.5*dnorm_x*dnorm_x+0.5*dnorm_q*dnorm_q);
  Ek_ = sigma_[currentOrder_]*dnorm;
  Tk_ = (currentOrder_+1)*Ek_;
  Est_ = Ek_;
  newOrder_ = currentOrder_;
  if (currentOrder_>1)
  {
    dsDae.delta_x->linearCombo(1.0,dsDae.xHistory[currentOrder_],1.0,*dsDae.newtonCorrectionPtr);
    dnorm_x = *(dsDae.delta_x->wRMSNorm(*dsDae.errWtVecPtr));
    dsDae.delta_q->linearCombo(1.0,dsDae.qHistory[currentOrder_],1.0,*dsDae.qNewtonCorrectionPtr);
    dnorm_q = *(dsDae.delta_q->wRMSNorm(*dsDae.qErrWtVecPtr));
    dnorm = sqrt(0.5*dnorm_x*dnorm_x+0.5*dnorm_q*dnorm_q);
    Ekm1_ = sigma_[currentOrder_-1]*dnorm;
    Tkm1_ = currentOrder_*Ekm1_;
    if (currentOrder_>2)
    {
      dsDae.delta_x->linearCombo(1.0,dsDae.xHistory[currentOrder_-1],1.0,*dsDae.delta_x);
      dnorm_x = *(dsDae.delta_x->wRMSNorm(*dsDae.errWtVecPtr));
      dsDae.delta_q->linearCombo(1.0,dsDae.qHistory[currentOrder_-1],1.0,*dsDae.delta_q);
      dnorm_q = *(dsDae.delta_q->wRMSNorm(*dsDae.qErrWtVecPtr));
      dnorm = sqrt(0.5*dnorm_x*dnorm_x+0.5*dnorm_q*dnorm_q);
      Ekm2_ = sigma_[currentOrder_-2]*dnorm;
      Tkm2_ = (currentOrder_-1)*Ekm2_;
      if (Xycemax(Tkm1_,Tkm2_)<=Tk_)
      {
        newOrder_--;
        Est_ = Ekm1_;
      }
    }
    else if (Tkm1_ <= Tkm1_Tk_safety_ * Tk_)
    {
      newOrder_--;
      Est_ = Ekm1_;
    }
  }
#ifdef Xyce_DEBUG_TIME
  if (tiaParams.debugLevel > 0)
  {
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                          "  newOrder = ", newOrder_);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, dashedline);
  }
#endif
}

template<class Scalar>
void ImplicitBDFStepper<Scalar>::rejectStep()
{

// This routine puts its output in newTimeStep_ and newOrder_

// This routine changes the following variables:
//    lastAttemptedTimeStep, initialPhase_, nef_, psi_, newTimeStep_,
//    newOrder_, currentOrder_, currentTimeStep_, dsDae.xHistory,
//    dsDae.qHistory, nextTimePt, nextTime, currentTimeStepRatio,
//    currentTimeStepSum, nextTimePt

// This routine reades but does not change the following variables:
//    stepAttemptStatus, r_factor_, r_safety_, Est_, r_fudge_, r_min_, r_max_,
//    minTimeStep, maxTimeStep, currentTime, stopTime, lastTimeStep


  // First we decide if we'll reduce the order independent of the local error test:
  checkReduceOrder();

#ifdef Xyce_DEBUG_TIME
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


  double newTimeStep_ = sec.currentTimeStep;
  double rr = 1.0; // step size ratio = new step / old step
  if ((sec.stepAttemptStatus == false) && (adjustStep))
  {
    initialPhase_ = false;
    nef_++;
    restoreHistory();
    // restore psi_
//    for (int i=1;i<=currentOrder_;++i)
//      psi_[i-1] = psi_[i] - sec.currentTimeStep;

    if (nef_ >= max_LET_fail_)  
    {
      string msg = "N_TIA_BackwardDifferentiation15::rejectStep: ";
      msg += "  Maximum number of local error test failures.  ";
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, msg);
    }

    if ((sec.newtonConvergenceStatus <= 0))
    {
      /// 11/11/05 erkeite:  If the Newton solver fails, don't 
      // rely on the error estimate - it may be full of Nan's.
      rr = r_min_;
      newTimeStep_ = rr * sec.currentTimeStep;

      if (nef_ > 2) newOrder_ = 1;//consistent with block below.
    }
    else
    {
      // 03/11/04 tscoffe:  Here is the block for choosing order & 
      // step-size when the local error test FAILS (but Newton 
      // succeeded). 
      if (nef_ == 1) // first local error test failure
      {
        rr = r_factor_*pow(r_safety_*Est_+r_fudge_,-1.0/(newOrder_+1.0));
        rr = Xycemax(r_min_,Xycemin(r_max_,rr));
        newTimeStep_ = rr * sec.currentTimeStep;
      }
      else if (nef_ == 2) // second failure
      {
        rr = r_min_;
        newTimeStep_ = rr * sec.currentTimeStep;
      }
      else if (nef_ > 2) // third and later failures
      {
        newOrder_ = 1;
        rr = r_min_;
        newTimeStep_ = rr * sec.currentTimeStep;
      }
    }
    currentOrder_ = newOrder_;
    if (numberOfSteps_ == 0) // still first step
    {
      psi_[0] = newTimeStep_;
      dsDae.xHistory[1].scale(rr);
      dsDae.qHistory[1].scale(rr);
    }
#ifdef Xyce_DEBUG_TIME
    if (tiaParams.debugLevel > 0)
    {
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
              "  currentTimeStep = ", sec.currentTimeStep);
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
              "  numberOfSteps_ = ", numberOfSteps_);
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
              "  currentOrder_ = ", currentOrder_);
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
              "  nscsco_ = ", nscsco_);
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
              "  alpha_[0] = ", alpha_[0]);
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
              "  alpha_[1] = ", alpha_[1]);
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
              "  alpha_[2] = ", alpha_[2]);
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
              "  alpha_[3] = ", alpha_[3]);
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
              "  alpha_[4] = ", alpha_[4]);
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
              "  psi_[0] = ", psi_[0]);
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
              "  psi_[1] = ", psi_[1]);
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
              "  psi_[2] = ", psi_[2]);
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
              "  psi_[3] = ", psi_[3]);
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
              "  psi_[4] = ", psi_[4]);
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
              "  sigma_[0] = ", sigma_[0]);
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
              "  sigma_[1] = ", sigma_[1]);
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
              "  sigma_[2] = ", sigma_[2]);
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
              "  sigma_[3] = ", sigma_[3]);
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
              "  sigma_[4] = ", sigma_[4]);
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                            "  rr = ", rr);
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                            "  r_factor_ = ", r_factor_);
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                            "  r_safety_ = ", r_safety_);
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                            "  Est_ = ", Est_);
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                            "  r_fudge_ = ", r_fudge_);
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                            "  newOrder_ = ", newOrder_);

      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                            "  currentTimeStep = ", sec.currentTimeStep);
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                            "  newTimeStep_ = ", newTimeStep_);
    }
#endif // Xyce_DEBUG_TIME
  }
  else if ((sec.stepAttemptStatus == false) & (!adjustStep))
  {
    string tmp = "  BackwardDifferentiation15:rejectStep: Warning: Local error test failed with constant step-size.\n";
    N_ERH_ErrorMgr::report (N_ERH_ErrorMgr::USR_INFO_0, tmp);
  }

  // If the step needs to be adjusted:
  if (adjustStep)
  {
    newTimeStep_ = Xycemax(newTimeStep_, sec.minTimeStep);
    newTimeStep_ = Xycemin(newTimeStep_, sec.maxTimeStep);

    double nextTimePt = sec.currentTime + newTimeStep_;

    if (nextTimePt > sec.stopTime)
    {
      nextTimePt  = sec.stopTime;
      newTimeStep_ = sec.stopTime - sec.currentTime;
    }

    sec.nextTime = nextTimePt;

    sec.currentTimeStepRatio = newTimeStep_/sec.lastTimeStep;
    sec.currentTimeStepSum   = newTimeStep_ + sec.lastTimeStep;

#ifdef Xyce_DEBUG_TIME
    if (tiaParams.debugLevel >0)
    {
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
              "  newTimeStep_ = ", newTimeStep_);
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
              "  nextTime = ", sec.nextTime);
    }
#endif // Xyce_DEBUG_TIME

    sec.currentTimeStep = newTimeStep_;
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
#ifdef Xyce_DEBUG_TIME
  if (tiaParams.debugLevel >0)
  {
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, dashedline);
  }
#endif // Xyce_DEBUG_TIME
}

template<class Scalar>
void ImplicitBDFStepper<Scalar>::completeStep()
{

  numberOfSteps_ ++;
  nef_ = 0;
  sec.lastTime    = sec.currentTime;
  sec.currentTime = sec.nextTime;
  // First we decide if we'll reduce the order independent of the local error test:
  checkReduceOrder();


#ifdef Xyce_DEBUG_TIME
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

  double newTimeStep_ = sec.currentTimeStep;
  double rr = 1.0; // step size ratio = new step / old step
  // 03/11/04 tscoffe:  Here is the block for choosing order & step-size when
  // the local error test PASSES (and Newton succeeded). 
  sec.lastTimeStep = sec.currentTimeStep;
  sec.lastTimeStepRatio = sec.currentTimeStepRatio; // copied from calcTStep1
  sec.lastTimeStepSum   = sec.currentTimeStepSum; // copied from calcTStep1
  int orderDiff = currentOrder_ - usedOrder_;
  usedOrder_ = currentOrder_;
  usedStep_ = sec.currentTimeStep;
  if ((newOrder_ == currentOrder_-1) || (currentOrder_ == tiaParams.maxOrder))
  {
    // If we reduced our order or reached max order then move to the next phase
    // of integration where we don't automatically double the step-size and
    // increase the order.
    initialPhase_ = false;
  }
  if (initialPhase_)
  {
    currentOrder_++;
    newTimeStep_ = h_phase0_incr_ * sec.currentTimeStep;
  }
  else // not in the initial phase of integration
  {
    int action = TIAAction_UNSET;
    if (newOrder_ == currentOrder_-1)
      action = TIAAction_LOWER;
    else if (newOrder_ == tiaParams.maxOrder)
      action = TIAAction_MAINTAIN;
    else if ((currentOrder_+1>=nscsco_) || (orderDiff == 1))
    {
      // If we just raised the order last time then we won't raise it again
      // until we've taken currentOrder_+1 steps at order currentOrder_.
      action = TIAAction_MAINTAIN;
    }
    else // consider changing the order 
    {
      dsDae.delta_x->linearCombo(1.0,*dsDae.newtonCorrectionPtr,-1.0,dsDae.xHistory[currentOrder_+1]);
      double dnorm_x = *(dsDae.delta_x->wRMSNorm(*dsDae.errWtVecPtr));
      dsDae.delta_q->linearCombo(1.0,*dsDae.qNewtonCorrectionPtr,-1.0,dsDae.qHistory[currentOrder_+1]);
      double dnorm_q = *(dsDae.delta_q->wRMSNorm(*dsDae.qErrWtVecPtr));
      double dnorm = sqrt(0.5*dnorm_x*dnorm_x+0.5*dnorm_q*dnorm_q);
      Tkp1_ = dnorm;
      Ekp1_ = Tkp1_/(currentOrder_+2);
      if (currentOrder_ == 1)
      {
        if (Tkp1_ >= Tkp1_Tk_safety_ * Tk_)
          action = TIAAction_MAINTAIN;
        else
          action = TIAAction_RAISE;
      }
      else
      {
        if (Tkm1_ <= Xycemin(Tk_,Tkp1_))
          action = TIAAction_LOWER;
        else if (Tkp1_ >= Tk_)
          action = TIAAction_MAINTAIN;
        else
          action = TIAAction_RAISE;
      }
    }
    if (action == TIAAction_RAISE)
    {
      currentOrder_++;
      Est_ = Ekp1_;
    }
    else if (action == TIAAction_LOWER)
    {
      currentOrder_--;
      Est_ = Ekm1_;
    }
    newTimeStep_ = sec.currentTimeStep;
    rr = pow(r_safety_*Est_+r_fudge_,-1.0/(currentOrder_+1.0));
#ifdef Xyce_DEBUG_TIME
    if (tiaParams.debugLevel > 0)
    {
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                            "  currentOrder_ = ", currentOrder_);
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                            "  r_safety = ", r_safety_);
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                            "  r_fudge_ = ", r_fudge_);
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                            "  r_hincr_ = ", r_hincr_);
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                            "  r_hincr_test_ = ", r_hincr_test_);
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                            "  Est = ", Est_);
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                            "  raw rr = ", rr);
    }
#endif
    if (rr >= r_hincr_test_)
    {
      rr = r_hincr_;
      newTimeStep_ = rr*sec.currentTimeStep;
    }
    else if (rr <= 1)
    {
      rr = Xycemax(r_min_,Xycemin(r_max_,rr));
      newTimeStep_ = rr*sec.currentTimeStep;
    }
  }
#ifdef Xyce_DEBUG_TIME
    if (tiaParams.debugLevel > 0)
    {
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                            "  initialPhase_ = ", initialPhase_);
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                            "  rr = ", rr);
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                            "  currentTimeStep = ", sec.currentTimeStep);
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                            "  currentTime = ", sec.currentTime);
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                            "  nextTime = ", sec.nextTime);
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                            "  newTimeStep_ = ", newTimeStep_);
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
      newTimeStep_ = Xycemax(newTimeStep_, sec.minTimeStep);
      newTimeStep_ = Xycemin(newTimeStep_, sec.maxTimeStep);

      double nextTimePt = sec.currentTime + newTimeStep_;

      if (nextTimePt > sec.stopTime)
      {
        nextTimePt  = sec.stopTime;
        newTimeStep_ = sec.stopTime - sec.currentTime;
      }

      sec.nextTime = nextTimePt;

      sec.currentTimeStepRatio = newTimeStep_/sec.lastTimeStep;
      sec.currentTimeStepSum   = newTimeStep_ + sec.lastTimeStep;

#ifdef Xyce_DEBUG_TIME
      if (tiaParams.debugLevel >0)
      {
        N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                "  nextTime = ", sec.nextTime);
        N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0,
                "  newTimeStep_ = ", newTimeStep_);
      }
#endif

      sec.currentTimeStep = newTimeStep_;
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
#ifdef Xyce_DEBUG_TIME
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

} // namespace Rythmos

#endif //Rythmos_STEPPER_IMPLICITBDF_H
