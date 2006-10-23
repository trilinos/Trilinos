//@HEADER
// ***********************************************************************
//
//                           Rythmos Package
//                 Copyright (2006) Sandia Corporation
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
// Questions? Contact Todd S. Coffey (tscoffe@sandia.gov)
//
// ***********************************************************************
//@HEADER

#ifndef Rythmos_IMPLICITBDF_STEPPER_H
#define Rythmos_IMPLICITBDF_STEPPER_H

#include "Rythmos_Stepper.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "Thyra_VectorBase.hpp"
#include "Thyra_ModelEvaluator.hpp"
#include "Thyra_ModelEvaluatorHelpers.hpp"
#include "Thyra_NonlinearSolverBase.hpp"
#include "Thyra_SingleResidSSDAEModelEvaluator.hpp"
#include "Thyra_SolveSupportTypes.hpp"

namespace Rythmos {

enum BDFactionFlag { ACTION_UNSET, ACTION_LOWER, ACTION_MAINTAIN, ACTION_RAISE };
enum BDFstatusFlag { PREDICT_AGAIN, CONTINUE_ANYWAY, REP_ERR_FAIL, REP_CONV_FAIL };

/** \brief . */
template<class Scalar>
class ImplicitBDFStepper : virtual public Stepper<Scalar>
{
  public:

    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType ScalarMag;

    /** \brief . */
    ImplicitBDFStepper(
      const Teuchos::RefCountPtr<const Thyra::ModelEvaluator<Scalar> >  &model
      ,const Teuchos::RefCountPtr<Thyra::NonlinearSolverBase<Scalar> >  &solver
      );

    /** \brief . */
    ImplicitBDFStepper(
      const Teuchos::RefCountPtr<const Thyra::ModelEvaluator<Scalar> >  &model
      ,const Teuchos::RefCountPtr<Thyra::NonlinearSolverBase<Scalar> >  &solver
      ,Teuchos::RefCountPtr<Teuchos::ParameterList> &parameterList
      );

    /** \brief . */
    void InitializeStepper(Teuchos::RefCountPtr<Teuchos::ParameterList> &implicitBDFParameters);

    /** \brief . */
    void setModel(const Teuchos::RefCountPtr<const Thyra::ModelEvaluator<Scalar> > &model);

    /** \brief . */
    void setSolver(const Teuchos::RefCountPtr<Thyra::NonlinearSolverBase<Scalar> > &solver);

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
    void describe(
      Teuchos::FancyOStream       &out
      ,const Teuchos::EVerbosityLevel      verbLevel
      ) const;

    bool ErrWtVecSet(
      Thyra::VectorBase<Scalar> *w, 
      const Thyra::VectorBase<Scalar> &y
      );

    Scalar WRMSNorm(
      const Thyra::VectorBase<Scalar> &w
      , const Thyra::VectorBase<Scalar> &y
      ) const;
    
    /// Redefined from InterpolationBufferBase 
    /// Add points to buffer
    bool SetPoints(
      const std::vector<Scalar>& time_vec
      ,const std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > >& x_vec
      ,const std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > >& xdot_vec
      ,const std::vector<ScalarMag> & accuracy_vec 
      );
    
    /// Get values from buffer
    bool GetPoints(
      const std::vector<Scalar>& time_vec
      ,std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > >* x_vec
      ,std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > >* xdot_vec
      ,std::vector<ScalarMag>* accuracy_vec) const;

    /// Fill data in from another interpolation buffer
    bool SetRange(
      const Scalar& time_lower
      ,const Scalar& time_upper
      ,const InterpolationBufferBase<Scalar> & IB);

    /// Get interpolation nodes
    bool GetNodes(std::vector<Scalar>* time_vec) const;

    /// Remove interpolation nodes
    bool RemoveNodes(std::vector<Scalar>& time_vec);

    /// Get order of interpolation
    int GetOrder() const;

    /// Redefined from Teuchos::ParameterListAcceptor
    /** \brief . */
    void setParameterList(Teuchos::RefCountPtr<Teuchos::ParameterList> const& paramList);

    /** \brief . */
    Teuchos::RefCountPtr<Teuchos::ParameterList> getParameterList();

    /** \brief . */
    Teuchos::RefCountPtr<Teuchos::ParameterList> unsetParameterList();

    Teuchos::RefCountPtr<const Teuchos::ParameterList> getValidParameters() const;

  private:


    void obtainPredictor();
    void obtainResidual();
    void obtainJacobian();
    void updateHistory();
    void restoreHistory();
    void updateCoeffs();
    void initialize();
    Scalar checkReduceOrder();
    BDFstatusFlag rejectStep();
    void completeStep();

    void setDefaultMagicNumbers(Teuchos::ParameterList &magicNumberList);

    bool interpolateSolution(
        const Scalar& timepoint
        ,Thyra::VectorBase<Scalar>* x_ptr_
        ,ScalarMag* accuracy_ptr_
        ) const;

    // 05/05/06 tscoffe:  I hate the underscores for private variables!
    Teuchos::RefCountPtr<const Thyra::ModelEvaluator<Scalar> > model;
    Teuchos::RefCountPtr<Thyra::NonlinearSolverBase<Scalar> > solver;
    Thyra::SingleResidSSDAEModelEvaluator<Scalar>   neModel;

    Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > xn0;
    Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > xpn0;
    Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > x_dot_base;
    std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > > xHistory;
    Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > ee;
    Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > delta;
    Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > residual;
    Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > errWtVec;

    Scalar time;


    //typedef typename Thyra::ModelEvaluatorBase::InArgs<Scalar>::ScalarMag ScalarMag;
    ScalarMag relErrTol; // relative error tolerance
    ScalarMag absErrTol; // absolute error tolerance
    Scalar hh;        // Current step-size
    int currentOrder; // Current order of integration
    int oldOrder;     // previous order of integration
    int maxOrder;     // maximum order = min(5,user option maxord) - see below.
    int usedOrder;    // order used in current step (used after currentOrder is updated)
    Scalar alpha_s;    // $\alpha_s$ fixed-leading coefficient of this BDF method
    vector<Scalar> alpha;    // $\alpha_j(n)=h_n/\psi_j(n)$ coefficient used in local error test
                      // note:   $h_n$ = current step size, n = current time step
    Scalar alpha_0;     // $-\sum_{j=1}^k \alpha_j(n)$ coefficient used in local error test
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
    bool constantStepSize;

    // Magic Numbers:
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
    int    max_LET_fail;
    Scalar minTimeStep;
    Scalar maxTimeStep;

    int newtonConvergenceStatus;

    Teuchos::RefCountPtr<Teuchos::ParameterList> parameterList;

};

// ////////////////////////////
// Defintions

template<class Scalar>
void ImplicitBDFStepper<Scalar>::InitializeStepper(
    Teuchos::RefCountPtr<Teuchos::ParameterList> &implicitBDFParameters
    )
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  // Initialize algorithm coefficients
  setParameterList(implicitBDFParameters);
  maxOrder = parameterList->get("maxOrder",int(5)); // maximum order
  maxOrder = max(1,min(maxOrder,5)); // 1 <= maxOrder <= 5
  currentOrder=1; // Current order of integration
  oldOrder=1; // previous order of integration
  usedOrder=1;  // order used in current step (used after currentOrder is updated)
  alpha_s=Scalar(-ST::one());  // $\alpha_s$ fixed-leading coefficient of this BDF method
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
  alpha_0=ST::zero();   // $-\sum_{j=1}^k \alpha_j(n)$ coefficient used in local error test
  cj=ST::zero();      // $-\alpha_s/h_n$ coefficient used in local error test
  ck=ST::zero();      // local error coefficient gamma_[0] = 0; // $\gamma_j(n)=\sum_{l=1}^{j-1}1/\psi_l(n)$ coefficient used to
  hh=ST::zero();
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

  relErrTol = parameterList->get( "relErrTol", Scalar(1.0e-4) );
  absErrTol = parameterList->get( "absErrTol", Scalar(1.0e-6) );
  constantStepSize = parameterList->get( "constantStepSize", false );
  stopTime = parameterList->get( "stopTime", Scalar(10.0) );

  int outputLevel = parameterList->get( "outputLevel", int(-1) );
  outputLevel = min(max(outputLevel,-1),4);
  this->setVerbLevel(static_cast<Teuchos::EVerbosityLevel>(outputLevel));
  Teuchos::RefCountPtr<Teuchos::FancyOStream> out = this->getOStream();
  out->precision(15);
  out->setMaxLenLinePrefix(28);
  out->pushLinePrefix("Rythmos::ImplicitBDFStepper");
  out->setShowLinePrefix(true);
  out->setTabIndentStr("    ");

  setDefaultMagicNumbers(parameterList->sublist("magicNumbers"));

  Teuchos::OSTab ostab(out,1,"InitializeStepper");
  if ( static_cast<int>(this->getVerbLevel()) >= static_cast<int>(Teuchos::VERB_HIGH) )
  {
    *out << "maxOrder = " << maxOrder << endl;
    *out << "currentOrder = " << currentOrder << endl;
    *out << "oldOrder = " << oldOrder << endl;
    *out << "usedOrder = " << usedOrder << endl;
    *out << "alpha_s = " << alpha_s << endl;
    for (int i=0 ; i<maxOrder ; ++i)
    {
      *out << "alpha[" << i << "] = " << alpha[i] << endl;
      *out << "beta[" << i << "] = " << beta[i] << endl;
      *out << "gamma[" << i << "] = " << gamma[i] << endl;
      *out << "psi[" << i << "] = " << psi[i] << endl;
      *out << "sigma[" << i << "] = " << sigma[i] << endl;
    }
    *out << "alpha_0 = " << alpha_0 << endl;
    *out << "cj = " << cj << endl;
    *out << "ck = " << ck << endl;
    *out << "numberOfSteps = " << numberOfSteps << endl;
    *out << "nef = " << nef << endl;
    *out << "usedStep = " << usedStep << endl;
    *out << "nscsco = " << nscsco << endl;
    *out << "Ek = " << Ek << endl;
    *out << "Ekm1 = " << Ekm1 << endl;
    *out << "Ekm2 = " << Ekm2 << endl;
    *out << "Ekp1 = " << Ekp1 << endl;
    *out << "Est = " << Est << endl;
    *out << "Tk = " << Tk << endl;
    *out << "Tkm1 = " << Tkm1 << endl;
    *out << "Tkm2 = " << Tkm2 << endl;
    *out << "Tkp1 = " << Tkp1 << endl;
    *out << "newOrder = " << newOrder << endl;
    *out << "initialPhase = " << initialPhase << endl;
    *out << "relErrTol  = " << relErrTol  << endl;
    *out << "absErrTol  = " << absErrTol  << endl;
    *out << "constantStepSize  = " << constantStepSize  << endl;
    *out << "stopTime  = " << stopTime  << endl;
  }
}

template<class Scalar>
ImplicitBDFStepper<Scalar>::ImplicitBDFStepper(
  const Teuchos::RefCountPtr<const Thyra::ModelEvaluator<Scalar> > &model
  ,const Teuchos::RefCountPtr<Thyra::NonlinearSolverBase<Scalar> > &solver
  ,Teuchos::RefCountPtr<Teuchos::ParameterList> &parameterList
  )
{
  InitializeStepper(parameterList);

  // Now we instantiate the model and the solver
  setModel(model);
  setSolver(solver);
}

template<class Scalar>
ImplicitBDFStepper<Scalar>::ImplicitBDFStepper(
  const Teuchos::RefCountPtr<const Thyra::ModelEvaluator<Scalar> > &model
  ,const Teuchos::RefCountPtr<Thyra::NonlinearSolverBase<Scalar> > &solver
  )
{
  Teuchos::ParameterList emptyList;
  InitializeStepper(emptyList);

  // Now we instantiate the model and the solver
  setModel(model);
  setSolver(solver);
}

template<class Scalar>
void ImplicitBDFStepper<Scalar>::setDefaultMagicNumbers(
    Teuchos::ParameterList &magicNumberList)
{
  // Magic Number Defaults:
  h0_safety      = magicNumberList.get( "h0_safety",      Scalar(2.0)     );
  h0_max_factor  = magicNumberList.get( "h0_max_factor",  Scalar(0.001)   );
  h_phase0_incr  = magicNumberList.get( "h_phase0_incr",  Scalar(2.0)     );
  h_max_inv      = magicNumberList.get( "h_max_inv",      Scalar(0.0)     );
  Tkm1_Tk_safety = magicNumberList.get( "Tkm1_Tk_safety", Scalar(2.0)     );
  Tkp1_Tk_safety = magicNumberList.get( "Tkp1_Tk_safety", Scalar(0.5)     );
  r_factor       = magicNumberList.get( "r_factor",       Scalar(0.9)     );
  r_safety       = magicNumberList.get( "r_safety",       Scalar(2.0)     );
  r_fudge        = magicNumberList.get( "r_fudge",        Scalar(0.0001)  );
  r_min          = magicNumberList.get( "r_min",          Scalar(0.125)   );
  r_max          = magicNumberList.get( "r_max",          Scalar(0.9)     );
  r_hincr_test   = magicNumberList.get( "r_hincr_test",   Scalar(2.0)     );
  r_hincr        = magicNumberList.get( "r_hincr",        Scalar(2.0)     );
  max_LET_fail   = magicNumberList.get( "max_LET_fail",   int(15)         );
  minTimeStep    = magicNumberList.get( "minTimeStep",    Scalar(0.0)     );
  maxTimeStep    = magicNumberList.get( "maxTimeStep",    Scalar(10.0)    ); 

  Teuchos::RefCountPtr<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::OSTab ostab(out,1,"setDefaultMagicNumbers");
  if ( static_cast<int>(this->getVerbLevel()) >= static_cast<int>(Teuchos::VERB_HIGH) )
  {
    *out << "h0_safety = " << h0_safety << endl;
    *out << "h0_max_factor = " << h0_max_factor << endl;
    *out << "h_phase0_incr = " << h_phase0_incr << endl;
    *out << "h_max_inv = " << h_max_inv << endl;
    *out << "Tkm1_Tk_safety = " << Tkm1_Tk_safety  << endl;
    *out << "Tkp1_Tk_safety = " << Tkp1_Tk_safety << endl;
    *out << "r_factor = " << r_factor << endl;
    *out << "r_safety = " << r_safety << endl;
    *out << "r_fudge = " << r_fudge << endl;
    *out << "r_min = " << r_min << endl;
    *out << "r_max = " << r_max << endl;
    *out << "r_hincr_test = " << r_hincr_test << endl;
    *out << "r_hincr = " << r_hincr << endl;
    *out << "max_LET_fail = " << max_LET_fail << endl;
    *out << "minTimeStep = " << minTimeStep << endl;
    *out << "maxTimeStep = " << maxTimeStep << endl;
  }
}

template<class Scalar>
void ImplicitBDFStepper<Scalar>::setModel(const Teuchos::RefCountPtr<const Thyra::ModelEvaluator<Scalar> > &model_)
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  model = model_;
  time = ST::zero();
  xn0 = model->getNominalValues().get_x()->clone_v();
  xpn0 = model->getNominalValues().get_x_dot()->clone_v(); 
  x_dot_base = model->getNominalValues().get_x_dot()->clone_v();
  ee = model->getNominalValues().get_x()->clone_v();
  delta = model->getNominalValues().get_x()->clone_v();
  residual = Thyra::createMember(model->get_f_space());
  errWtVec = xn0->clone_v();
  ErrWtVecSet(&*errWtVec,*xn0);
  xHistory.push_back(xn0->clone_v());
  xHistory.push_back(xpn0->clone_v());
  for (int i=2 ; i<=maxOrder ; ++i) // Store maxOrder+1 vectors
  {
    xHistory.push_back(xn0->clone_v()); 
    V_S(&*xHistory[i],ST::zero());
  }
  initialize(); // Now that we have the model, we can do our initialization 
}

template<class Scalar>
void ImplicitBDFStepper<Scalar>::setSolver(const Teuchos::RefCountPtr<Thyra::NonlinearSolverBase<Scalar> > &solver_)
{
  solver = solver_;
}

template<class Scalar>
Scalar ImplicitBDFStepper<Scalar>::TakeStep()
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  typedef typename Thyra::ModelEvaluatorBase::InArgs<Scalar>::ScalarMag ScalarMag;
  Teuchos::RefCountPtr<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::OSTab ostab(out,1,"TakeStep");
  BDFstatusFlag status;
  while (1)
  {
    // Set up problem coefficients (and handle first step)
    updateCoeffs();
    // Set Error weight vector
    ErrWtVecSet(&*errWtVec,*xn0);
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
    Scalar coeff_x_dot = Scalar(-ST::one())*alpha_s/hh;
    V_StVpStV( &*x_dot_base, ST::one(), *xpn0, alpha_s/hh, *xn0 );
    neModel.initialize(model,coeff_x_dot,x_dot_base,ST::one(),Teuchos::null,time+hh,xn0);
    //
    // Solve the implicit nonlinear system to a tolerance of ???
    // 
    // 05/08/06 tscoffe:  I really need to get the update, not the solution from
    // the nonlinear solver.
    if(solver->getModel().get()!=&neModel)
      solver->setModel( Teuchos::rcp(&neModel,false) );
    /* // Thyra::TimeStepNewtonNonlinearSolver uses a built in solveCriteria, so you can't pass one in.
       // I believe this is the correct solveCriteria for IDA though.
    Thyra::SolveMeasureType nonlinear_solve_measure_type(Thyra::SOLVE_MEASURE_NORM_RESIDUAL,Thyra::SOLVE_MEASURE_ONE); 
    ScalarMag tolerance = relErrTol/ScalarMag(10.0); // This should be changed to match the condition in IDA
    Thyra::SolveCriteria<Scalar> nonlinearSolveCriteria(nonlinear_solve_measure_type, tolerance);
    Thyra::SolveStatus<Scalar> nonlinearSolveStatus = solver->solve( &*xn0, &nonlinearSolveCriteria, &*delta ); 
    */
    //Thyra::assign(&*xn0,ST::zero()); // 08/10/06 tscoffe:  what is this doing here?  It hoses the solve.
    Thyra::SolveStatus<Scalar> nonlinearSolveStatus = solver->solve( &*xn0, NULL, &*ee ); 
    if (nonlinearSolveStatus.solveStatus == Thyra::SOLVE_STATUS_CONVERGED) 
      newtonConvergenceStatus = 0;
    else 
      newtonConvergenceStatus = -1;

    // check error and evaluate LTE
    Scalar enorm = checkReduceOrder();
    
    if ( static_cast<int>(this->getVerbLevel()) >= static_cast<int>(Teuchos::VERB_HIGH) )
    {
      *out << "xn0 = " << std::endl;
      xn0->describe(*out,this->getVerbLevel());
      *out << "ee = " << std::endl;
      ee->describe(*out,this->getVerbLevel());
      *out << "delta = " << std::endl;
      delta->describe(*out,this->getVerbLevel());
      for (int i=0; i<max(2,maxOrder); ++i)
      {
        *out << "xHistory[" << i << "] = " << std::endl;
        xHistory[i]->describe(*out,this->getVerbLevel());
      }
      *out << "ck = " << ck << endl;
      *out << "enorm = " << enorm << endl;
      *out << "Local Truncation Error Check: (ck*enorm) < 1:  (" << ck*enorm << ") <?= 1" << endl;
    }
    // Check LTE here:
    if ((ck*enorm) > ST::one())
      status = rejectStep(); 
    else 
      break;
    if (status == CONTINUE_ANYWAY)
      break;
    if (status == REP_ERR_FAIL)
      return(Scalar(-ST::one()));
  }

  completeStep();  
  return(usedStep);
}

template<class Scalar>
Scalar ImplicitBDFStepper<Scalar>::TakeStep(Scalar dt)
{
  constantStepSize = true;
  hh = dt;
  dt = TakeStep();
  return(dt);
}

template<class Scalar>
Teuchos::RefCountPtr<const Thyra::VectorBase<Scalar> > ImplicitBDFStepper<Scalar>::get_solution() const
{
  return(xHistory[0]);
}

template<class Scalar>
Teuchos::RefCountPtr<const Thyra::VectorBase<Scalar> > ImplicitBDFStepper<Scalar>::get_residual() const
{
  return(residual);
}

template<class Scalar>
std::string ImplicitBDFStepper<Scalar>::description() const
{
  std::string name = "Rythmos::ImplicitBDFStepper";
  return(name);
}

template<class Scalar>
void ImplicitBDFStepper<Scalar>::describe(
      Teuchos::FancyOStream                &out
      ,const Teuchos::EVerbosityLevel      verbLevel
      ) const
{
  if (static_cast<int>(verbLevel) >= static_cast<int>(Teuchos::VERB_DEFAULT))
  {
    out << description() << "::describe" << std::endl;
    out << "model = " << model->description() << std::endl;
    out << "solver = " << solver->description() << std::endl;
    out << "neModel = " << neModel.description() << std::endl;
  }
  else if (static_cast<int>(verbLevel) >= static_cast<int>(Teuchos::VERB_LOW))
  {
    out << "time = " << time << std::endl;
    out << "hh = " << hh << std::endl;
    out << "currentOrder = " << currentOrder << std::endl;
  }
  else if (static_cast<int>(verbLevel) >= static_cast<int>(Teuchos::VERB_MEDIUM))
  {
  }
  else if (static_cast<int>(verbLevel) >= static_cast<int>(Teuchos::VERB_HIGH))
  {
    out << "model = " << std::endl;
    model->describe(out,verbLevel);
    out << "solver = " << std::endl;
    solver->describe(out,verbLevel);
    out << "neModel = " << std::endl;
    neModel.describe(out,verbLevel);
    out << "xn0 = " << std::endl;
    xn0->describe(out,verbLevel);
    out << "xpn0 = " << std::endl;
    xpn0->describe(out,verbLevel);
    out << "x_dot_base = " << std::endl;
    x_dot_base->describe(out,verbLevel);
    out << "xHistory = " << std::endl;
    for (int i=0 ; i < max(2,maxOrder) ; ++i)
    {
      out << "xHistory[" << i << "] = " << std::endl;
      xHistory[i]->describe(out,verbLevel);
    }
    out << "ee = " << std::endl;
    ee->describe(out,verbLevel);
    out << "delta = " << std::endl;
    delta->describe(out,verbLevel);
    out << "residual = " << std::endl;
    residual->describe(out,verbLevel);
    out << "errWtVec = " << std::endl;
    errWtVec->describe(out,verbLevel);
  }
}

template<class Scalar>
void ImplicitBDFStepper<Scalar>::obtainPredictor()
{
  typedef Teuchos::ScalarTraits<Scalar> ST;

  Teuchos::RefCountPtr<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::OSTab ostab(out,1,"obtainPredictor");
  if ( static_cast<int>(this->getVerbLevel()) >= static_cast<int>(Teuchos::VERB_HIGH) )
    *out << "currentOrder = " << currentOrder << std::endl;
  
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
    Vp_StV(&*xpn0,gamma[i],*xHistory[i]);
  }
  if ( static_cast<int>(this->getVerbLevel()) >= static_cast<int>(Teuchos::VERB_HIGH) )
  {
    *out << "xn0 = " << std::endl;
    xn0->describe(*out,this->getVerbLevel());
    *out << "xpn0 = " << std::endl;
    xpn0->describe(*out,this->getVerbLevel());
  }
}

template<class Scalar>
bool ImplicitBDFStepper<Scalar>::interpolateSolution(
        const Scalar& timepoint
        ,Thyra::VectorBase<Scalar>* x_ptr_
        ,ScalarMag* accuracy_ptr_
        ) const
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  Scalar tfuzz;   // fuzz factor to check for valid output time
  Scalar tp;      // approximately t{n-1}
  Scalar delt;    // distance between timepoint and time
  Scalar c = ST::one(); // coefficient for interpolation
  Scalar gam;     // coefficient for interpolation
  int kord;       // order of interpolation
  Scalar tn = time;
  Scalar hused = usedStep;
  int kused = usedOrder;
  Scalar uround = ST::zero();  // unit round-off (set to zero for now)

  tfuzz = 100 * uround * (tn + hh);
  tp = tn - hused - tfuzz;
  if ( (timepoint - tp)*hh < ST::zero() ) 
    return(false);
  if ( timepoint - (time-tfuzz) > ST::zero() )
    return(false);

  Thyra::V_V(x_ptr_,*xHistory[0]);
  kord = kused;
  if ( (kused == 0) || (timepoint == tn) ) 
    kord = 1;

  delt = timepoint - tn;
  gam = delt/psi[0];
  for (int j=1 ; j <= kord ; ++j)
  {
    c = c*gam;
    gam = (delt + psi[j-1])/psi[j];
    Thyra::Vp_StV(x_ptr_,c,*xHistory[j]);
  }
  *accuracy_ptr_ = pow(usedStep,kord);
  return(true);
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
  Teuchos::RefCountPtr<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::OSTab ostab(out,1,"updateHistory");
  if (static_cast<int>(this->getVerbLevel()) >= static_cast<int>(Teuchos::VERB_HIGH) )
  {
    for (int i=0;i<max(2,maxOrder);++i)
    {
      *out << "xHistory[" << i << "] = " << endl;
      xHistory[i]->describe(*out,this->getVerbLevel());
    }
  }

}

template<class Scalar>
void ImplicitBDFStepper<Scalar>::restoreHistory()
{
  typedef Teuchos::ScalarTraits<Scalar> ST;

  // undo preparation of history array for prediction
  for (int i=nscsco;i<=currentOrder;++i)
  {
    Vt_S( &*xHistory[i], ST::one()/beta[i] );
  }
  for (int i=1;i<=currentOrder;++i)
  {
    psi[i-1] = psi[i] - hh;
  }
  Teuchos::RefCountPtr<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::OSTab ostab(out,1,"restoreHistory");
  if (static_cast<int>(this->getVerbLevel()) >= static_cast<int>(Teuchos::VERB_HIGH) )
  {
    for (int i=0;i<maxOrder;++i)
      *out << "psi[" << i << "] = " << psi[i] << endl;
    for (int i=0;i<maxOrder;++i)
    {
      *out << "xHistory[" << i << "] = " << endl;
      xHistory[i]->describe(*out,this->getVerbLevel());
    }
  }
} 

template<class Scalar>
void ImplicitBDFStepper<Scalar>::updateCoeffs()
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
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
    alpha_s = ST::zero();
    alpha_0 = ST::zero();
    for (int i=0;i<currentOrder;++i)
    {
      alpha_s = alpha_s - Scalar(ST::one()/(i+ST::one()));
      alpha_0 = alpha_0 - alpha[i];
    }
    cj = -alpha_s/hh;
    ck = abs(alpha[currentOrder]+alpha_s-alpha_0);
    ck = max(ck,alpha[currentOrder]);
  }
  Teuchos::RefCountPtr<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::OSTab ostab(out,1,"updateCoeffs");
  if ( static_cast<int>(this->getVerbLevel()) >= static_cast<int>(Teuchos::VERB_HIGH) )
  {
    for (int i=0;i<=maxOrder;++i)
    {
      *out << "alpha[" << i << "] = " << alpha[i] << endl;
      *out << "beta[" << i << "] = " << beta[i] << endl;
      *out << "sigma[" << i << "] = " << sigma[i] << endl;
      *out << "gamma[" << i << "] = " << gamma[i] << endl;
      *out << "psi[" << i << "] = " << psi[i] << endl;
      *out << "alpha_s = " << alpha_s << endl;
      *out << "alpha_0 = " << alpha_0 << endl;
      *out << "cj = " << cj << endl;
      *out << "ck = " << ck << endl;
    }
  }

}

template<class Scalar>
void ImplicitBDFStepper<Scalar>::initialize()
{
  typedef Teuchos::ScalarTraits<Scalar> ST;

  // Choose initial step-size
  Scalar time_to_stop = stopTime - time;
  Scalar currentTimeStep;
  if (constantStepSize)
  {
    currentTimeStep = hh;
    //currentTimeStep = 0.1 * time_to_stop;
    //currentTimeStep = min(hh, currentTimeStep);
  }
  else
  {
    // compute an initial step-size based on rate of change in the solution initially
    Scalar ypnorm = WRMSNorm(*errWtVec,*xHistory[1]);
    if (ypnorm > ST::zero())  // time-dependent DAE
    {
      currentTimeStep = min(h0_max_factor*abs(time_to_stop),sqrt(2.0)/(h0_safety*ypnorm));
    } 
    else  // non-time-dependent DAE
    {
      currentTimeStep = h0_max_factor*abs(time_to_stop);
    }
    // choose min of user specified value and our value:
    if (hh > ST::zero())
      currentTimeStep = min(hh, currentTimeStep);
    // check for maximum step-size:
    Scalar rh = abs(currentTimeStep)*h_max_inv; 
    if (rh>1.0) currentTimeStep = currentTimeStep/rh;
  }
  hh = currentTimeStep;
  Teuchos::RefCountPtr<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::OSTab ostab(out,1,"initialize");
  if ( static_cast<int>(this->getVerbLevel()) >= static_cast<int>(Teuchos::VERB_HIGH) )
  {
    *out << "hh = " << hh << endl;
  }

  // x history
  assign(&*xHistory[0],*xn0);
  V_S(&*xHistory[1],ST::zero());

  // Coefficient initialization 
  numberOfSteps = 0;    // number of total time integration steps taken
  currentOrder = 1;
  usedOrder = 1;
  psi[0] = hh;
  cj = 1/psi[0];
  nscsco = 0;
}

template<class Scalar>
Scalar ImplicitBDFStepper<Scalar>::checkReduceOrder()
{
// This routine puts its output in newOrder_

// This routine changes the following variables:
//    Ek, Tk, Est, newOrder, dsDae.delta_x, dsDae.delta_q,
//    Ekm1, Tkm1, Ekm2, Tkm2 

// This routine reads but does not change the following variables:
//    currentOrder, sigma, dsDae.newtonCorrectionPtr, dsDae.qNewtonCorrectionPtr,
//    dsDae.errWtVecPtr, dsDae.qErrWtVecPtr, dsDae.xHistory, dsDae.qHistory

  Scalar enorm = WRMSNorm(*errWtVec,*ee);
  Ek = sigma[currentOrder]*enorm;
  Tk = Scalar(currentOrder+1)*Ek;
  Est = Ek;
  newOrder = currentOrder;
  Teuchos::RefCountPtr<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::OSTab ostab(out,1,"checkReduceOrder");
  if ( static_cast<int>(this->getVerbLevel()) >= static_cast<int>(Teuchos::VERB_HIGH) )
  {
    *out << "currentOrder = " << currentOrder << std::endl;
    *out << "Ek = " << Ek << std::endl;
    *out << "Tk = " << Tk << std::endl;
    *out << "enorm = " << enorm << std::endl;
  }
  if (currentOrder>1)
  {
    V_VpV(&*delta,*xHistory[currentOrder],*ee);
    Ekm1 = sigma[currentOrder-1]*WRMSNorm(*errWtVec,*delta);
    Tkm1 = currentOrder*Ekm1;
    if ( static_cast<int>(this->getVerbLevel()) >= static_cast<int>(Teuchos::VERB_HIGH) )
    {
      *out << "Ekm1 = " << Ekm1 << endl;
      *out << "Tkm1 = " << Tkm1 << endl;
    }
    if (currentOrder>2)
    {
      Vp_V(&*delta,*xHistory[currentOrder-1]);
      Ekm2 = sigma[currentOrder-2]*WRMSNorm(*errWtVec,*delta);
      Tkm2 = (currentOrder-1)*Ekm2;
      if ( static_cast<int>(this->getVerbLevel()) >= static_cast<int>(Teuchos::VERB_HIGH) )
      {
        *out << "Ekm2 = " << Ekm2 << endl;
        *out << "Tkm2 = " << Tkm2 << endl;
      }
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
  if ( static_cast<int>(this->getVerbLevel()) >= static_cast<int>(Teuchos::VERB_HIGH) )
  {
    *out << "Est = " << Est << endl;
    *out << "newOrder= " << newOrder << endl;
  }
  return(enorm);
}

template<class Scalar>
BDFstatusFlag ImplicitBDFStepper<Scalar>::rejectStep()
{

// This routine puts its output in newTimeStep and newOrder

// This routine changes the following variables:
//    initialPhase, nef, psi, newTimeStep,
//    newOrder, currentOrder, currentTimeStep, dsDae.xHistory,
//    dsDae.qHistory, nextTimePt, 
//    currentTimeStepSum, nextTimePt

// This routine reads but does not change the following variables:
//    r_factor, r_safety, Est, r_fudge, r_min, r_max,
//    minTimeStep, maxTimeStep, time, stopTime 

  // Only update the time step if we are NOT running constant stepsize.
  bool adjustStep = (!constantStepSize);

  Scalar newTimeStep = hh;
  Scalar rr = 1.0; // step size ratio = new step / old step
  nef++;
  Teuchos::RefCountPtr<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::OSTab ostab(out,1,"rejectStep");
  if ( static_cast<int>(this->getVerbLevel()) >= static_cast<int>(Teuchos::VERB_HIGH) )
  {
    *out << "adjustStep = " << adjustStep << endl;
    *out << "nef = " << nef << endl;
  }
  if (nef >= max_LET_fail)  
  {
    cerr << "Rythmos_Stepper_ImplicitBDF::rejectStep:  " 
          << "  Maximum number of local error test failures.  " << endl;
    return(REP_ERR_FAIL);
  }
  if (adjustStep)
  {
    initialPhase = false;
    if ( static_cast<int>(this->getVerbLevel()) >= static_cast<int>(Teuchos::VERB_HIGH) )
    {
      *out << "initialPhase = " << initialPhase << endl;
    }
    restoreHistory();
    // restore psi_
//    for (int i=1;i<=currentOrder;++i)
//      psi[i-1] = psi[i] - hh;

    if ((newtonConvergenceStatus < 0))
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
        rr = max(r_min,min(r_max,rr));
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
    if ( static_cast<int>(this->getVerbLevel()) >= static_cast<int>(Teuchos::VERB_HIGH) )
    {
      *out << "rr = " << rr << endl;
      *out << "newOrder = " << newOrder << endl;
    }
    currentOrder = newOrder;
    if (numberOfSteps == 0) // still first step
    {
      psi[0] = newTimeStep;
      Vt_S(&*xHistory[1],rr);
      if ( static_cast<int>(this->getVerbLevel()) >= static_cast<int>(Teuchos::VERB_HIGH) )
      {
        *out << "numberOfSteps == 0:" << endl;
        *out << "psi[0] = " << psi[0] << endl;
        *out << "xHistory[1] = " << std::endl;
        xHistory[1]->describe(*out,this->getVerbLevel());
      }
    }
  }
  else if (!adjustStep)
  {
    cerr << "Rythmos_Stepper_ImplicitBDF::rejectStep:  "
         << "Warning: Local error test failed with constant step-size." << endl;
  }

  BDFstatusFlag return_status = PREDICT_AGAIN;

  // If the step needs to be adjusted:
  if (adjustStep)
  {
    newTimeStep = max(newTimeStep, minTimeStep);
    newTimeStep = min(newTimeStep, maxTimeStep);

    Scalar nextTimePt = time + newTimeStep;

    if (nextTimePt > stopTime)
    {
      nextTimePt  = stopTime;
      newTimeStep = stopTime - time;
    }

    hh = newTimeStep;
  }
  else // if time step is constant for this step:
  {
    Scalar nextTimePt = time + hh;

    if (nextTimePt > stopTime)
    {
      nextTimePt      = stopTime;
      hh = stopTime - time;
    }
    return_status = CONTINUE_ANYWAY;
  }
  if ( static_cast<int>(this->getVerbLevel()) >= static_cast<int>(Teuchos::VERB_HIGH) )
  {
    *out << "hh = " << hh << endl;
  }
  return(return_status);
}

template<class Scalar>
void ImplicitBDFStepper<Scalar>::completeStep()
{
  typedef Teuchos::ScalarTraits<Scalar> ST;

  numberOfSteps ++;
  nef = 0;
  time += hh;
  Teuchos::RefCountPtr<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::OSTab ostab(out,1,"completeStep");
  if ( static_cast<int>(this->getVerbLevel()) >= static_cast<int>(Teuchos::VERB_HIGH) )
  {
    *out << "numberOfSteps = " << numberOfSteps << endl;
    *out << "nef = " << nef << endl;
    *out << "time = " << time << endl;
  }
  
  // Only update the time step if we are NOT running constant stepsize.
  bool adjustStep = (!constantStepSize);

  Scalar newTimeStep = hh;
  Scalar rr = ST::one(); // step size ratio = new step / old step
  // 03/11/04 tscoffe:  Here is the block for choosing order & step-size when
  // the local error test PASSES (and Newton succeeded). 
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
  if ( static_cast<int>(this->getVerbLevel()) >= static_cast<int>(Teuchos::VERB_HIGH) )
  {
    *out << "initialPhase = " << initialPhase << endl;
  }
  if (initialPhase)
  {
    currentOrder++;
    newTimeStep = h_phase0_incr * hh;
    if ( static_cast<int>(this->getVerbLevel()) >= static_cast<int>(Teuchos::VERB_HIGH) )
    {
      *out << "currentOrder = " << currentOrder << endl;
      *out << "newTimeStep = " << newTimeStep << endl;
    }
  }
  else // not in the initial phase of integration
  {
    BDFactionFlag action = ACTION_UNSET;
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
      Tkp1 = WRMSNorm(*errWtVec,*delta);
      Ekp1 = Tkp1/(currentOrder+2);
      if ( static_cast<int>(this->getVerbLevel()) >= static_cast<int>(Teuchos::VERB_HIGH) )
      {
        *out << "delta = " << endl;
        delta->describe(*out,this->getVerbLevel());
        *out << "Tkp1 = ||delta||_WRMS = " << Tkp1 << endl;
        *out << "Ekp1 = " << Ekp1 << endl;
      }
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
      rr = max(r_min,min(r_max,rr));
      newTimeStep = rr*hh;
    }
    if ( static_cast<int>(this->getVerbLevel()) >= static_cast<int>(Teuchos::VERB_HIGH) )
    {
      *out << "Est = " << Est << endl;
      *out << "currentOrder = " << currentOrder << endl;
      *out << "rr  = " << rr << endl;
      *out << "newTimeStep = " << newTimeStep << endl;
    }
  }
  // 03/22/04 tscoffe:  Note that updating the history has nothing to do with
  // the step-size and everything to do with the newton correction vectors.
  updateHistory();

  // 12/01/05 tscoffe:  This is necessary to avoid currentTimeStep == 0 right
  // before a breakpoint.  So I'm checking to see if time is identically
  // equal to stopTime, in which case we are right before a breakpoint and we
  // should not adjust currentStepSize because that would result in
  // currentStepSize == 0.
  if (time < stopTime)
  {
    // If the step needs to be adjusted:
    if (adjustStep)
    {
      newTimeStep = max(newTimeStep, minTimeStep);
      newTimeStep = min(newTimeStep, maxTimeStep);

      Scalar nextTimePt = time + newTimeStep;

      if (nextTimePt > stopTime)
      {
        nextTimePt  = stopTime;
        newTimeStep = stopTime - time;
      }

      hh = newTimeStep;
    }
    else // if time step is constant for this step:
    {
      Scalar nextTimePt = time + hh;

      if (nextTimePt > stopTime)
      {
        nextTimePt      = stopTime;
        hh = stopTime - time;
      }
    }
  }
  if ( static_cast<int>(this->getVerbLevel()) >= static_cast<int>(Teuchos::VERB_HIGH) )
  {
    *out << "hh = " << hh << endl;
  }
}

template<class Scalar>
bool ImplicitBDFStepper<Scalar>::ErrWtVecSet(Thyra::VectorBase<Scalar> *w_in, const Thyra::VectorBase<Scalar> &y)
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  Thyra::VectorBase<Scalar> &w = *w_in;
  abs(&w,y);
  Vt_S(&w,relErrTol);
  Vp_S(&w,absErrTol);
  reciprocal(&w,w);
  Vt_StV(&w,ST::one(),w); // We square w because of how weighted norm_2 is computed.
  // divide by N to get RMS norm
  int N = y.space()->dim();
  Vt_S(&w,Scalar(1.0/N));
  // Now you can compute WRMS norm as:
  // Scalar WRMSnorm = norm_2(w,y); // WRMS norm of y with respect to weights w.
  return true; // This should be updated to reflect success of reciprocal
}

template<class Scalar>
Scalar ImplicitBDFStepper<Scalar>::WRMSNorm(const Thyra::VectorBase<Scalar> &w, const Thyra::VectorBase<Scalar> &y) const
{
  return(norm_2(w,y));
}

template<class Scalar>
bool ImplicitBDFStepper<Scalar>::SetPoints(
    const std::vector<Scalar>& time_vec
    ,const std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > >& x_vec
    ,const std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > >& xdot_vec
    ,const std::vector<ScalarMag> & accuracy_vec 
    )
{
  return(false);
}

template<class Scalar>
bool ImplicitBDFStepper<Scalar>::GetPoints(
    const std::vector<Scalar>& time_vec
    ,std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > >* x_vec
    ,std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > >* xdot_vec
    ,std::vector<ScalarMag>* accuracy_vec) const
{
  bool status;
  Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > xdot_temp; // Teuchos::null
  for (int i=0 ; i<time_vec.size() ; ++i)
  {
    Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > x_temp = xn0->clone_v();
    ScalarMag accuracy;
    status = interpolateSolution(time_vec[i],&*x_temp,&accuracy);
    if (!status) return(status);
    x_vec->push_back(x_temp);
    xdot_vec->push_back(xdot_temp);
    accuracy_vec->push_back(accuracy);
  }
  if ( static_cast<int>(this->getVerbLevel()) >= static_cast<int>(Teuchos::VERB_HIGH) )
  {
    Teuchos::RefCountPtr<Teuchos::FancyOStream> out = this->getOStream();
    Teuchos::OSTab ostab(out,1,"BDFS::GetPoints");
    *out << "Passing out the interpolated values:" << std::endl;
    for (int i=0; i<time_vec.size() ; ++i)
    {
      *out << "time[" << i << "] = " << time_vec[i] << std::endl;
      *out << "x_vec[" << i << "] = " << std::endl;
      (*x_vec)[i]->describe(*out,this->getVerbLevel());
      if ( (*xdot_vec)[i] == Teuchos::null)
        *out << "xdot_vec[" << i << "] = Teuchos::null" << std::endl;
      else
      {
        *out << "xdot_vec[" << i << "] = " << std::endl;
        (*xdot_vec)[i]->describe(*out,this->getVerbLevel());
      }
      *out << "accuracy[" << i << "] = " << (*accuracy_vec)[i] << std::endl;
    }
  }
  return(status);
}

template<class Scalar>
bool ImplicitBDFStepper<Scalar>::SetRange(
    const Scalar& time_lower 
    ,const Scalar& time_upper
    ,const InterpolationBufferBase<Scalar>& IB)
{
  return(false);
}

template<class Scalar>
bool ImplicitBDFStepper<Scalar>::GetNodes(std::vector<Scalar>* time_vec) const
{
  if (numberOfSteps > 0)
    time_vec->push_back(time-usedStep);
  time_vec->push_back(time);
  return(true);
}

template<class Scalar>
bool ImplicitBDFStepper<Scalar>::RemoveNodes(std::vector<Scalar>& time_vec) 
{
  return(false);
}

template<class Scalar>
int ImplicitBDFStepper<Scalar>::GetOrder() const
{
  return(currentOrder);
}

template<class Scalar>
void ImplicitBDFStepper<Scalar>::setParameterList(Teuchos::RefCountPtr<Teuchos::ParameterList> const& paramList)
{
  parameterList = paramList;
}

template<class Scalar>
Teuchos::RefCountPtr<Teuchos::ParameterList> ImplicitBDFStepper<Scalar>::getParameterList()
{
  return(parameterList);
}

template<class Scalar>
Teuchos::RefCountPtr<Teuchos::ParameterList> ImplicitBDFStepper<Scalar>::unsetParameterList()
{
  Teuchos::RefCountPtr<Teuchos::ParameterList> temp_param_list = parameterList;
  parameterList = Teuchos::null;
  return(temp_param_list);
}

template<class Scalar>
Teuchos::RefCountPtr<const Teuchos::ParameterList> ImplicitBDFStepper<Scalar>::getValidParameters() const
{
  Teuchos::RefCountPtr<Teuchos::ParameterList> temp_param_list = Teuchos::rcp(new Teuchos::ParameterList);
  temp_param_list->set<int>   ( "maxOrder",         5              );
  temp_param_list->set<Scalar>( "relErrTol",        Scalar(1.0e-4) );
  temp_param_list->set<Scalar>( "absErrTol",        Scalar(1.0e-6) );
  temp_param_list->set<bool>  ( "constantStepSize", false          );
  temp_param_list->set<Scalar>( "stopTime",         Scalar(10.0)   );
  temp_param_list->set<int>   ( "outputLevel",       int(-1)        );

  Teuchos::ParameterList& magicNumberList = temp_param_list->sublist("magicNumbers");
  magicNumberList.set<Scalar>( "h0_safety",      Scalar(2.0)     );
  magicNumberList.set<Scalar>( "h0_max_factor",  Scalar(0.001)   );
  magicNumberList.set<Scalar>( "h_phase0_incr",  Scalar(2.0)     );
  magicNumberList.set<Scalar>( "h_max_inv",      Scalar(0.0)     );
  magicNumberList.set<Scalar>( "Tkm1_Tk_safety", Scalar(2.0)     );
  magicNumberList.set<Scalar>( "Tkp1_Tk_safety", Scalar(0.5)     );
  magicNumberList.set<Scalar>( "r_factor",       Scalar(0.9)     );
  magicNumberList.set<Scalar>( "r_safety",       Scalar(2.0)     );
  magicNumberList.set<Scalar>( "r_fudge",        Scalar(0.0001)  );
  magicNumberList.set<Scalar>( "r_min",          Scalar(0.125)   );
  magicNumberList.set<Scalar>( "r_max",          Scalar(0.9)     );
  magicNumberList.set<Scalar>( "r_hincr_test",   Scalar(2.0)     );
  magicNumberList.set<Scalar>( "r_hincr",        Scalar(2.0)     );
  magicNumberList.set<int>   ( "max_LET_fail",   int(15)         );
  magicNumberList.set<Scalar>( "minTimeStep",    Scalar(0.0)     );
  magicNumberList.set<Scalar>( "maxTimeStep",    Scalar(10.0)    ); 

  return(temp_param_list);
}

} // namespace Rythmos

#endif //Rythmos_IMPLICITBDF_STEPPER_H
