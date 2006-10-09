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
      ,Teuchos::ParameterList &parameterList
      );

    /** \brief . */
    void InitializeStepper(Teuchos::ParameterList &implicitBDFParameters);

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
      const std::vector<Scalar>& time_list
      ,const std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > >& x_list
      ,const std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > >& xdot_list);
    
    /// Get values from buffer
    bool GetPoints(
      const std::vector<Scalar>& time_list
      ,std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > >* x_list
      ,std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > >* xdot_list
      ,std::vector<ScalarMag>* accuracy_list) const;

    /// Fill data in from another interpolation buffer
    bool SetRange(
      const Scalar& time_lower
      ,const Scalar& time_upper
      ,const InterpolationBufferBase<Scalar> & IB);

    /// Get interpolation nodes
    bool GetNodes(std::vector<Scalar>* time_list) const;

    /// Remove interpolation nodes
    bool RemoveNodes(std::vector<Scalar>& time_list) const;

    /// Get order of interpolation
    int GetOrder() const;


  private:


    void obtainPredictor();
    void obtainResidual();
    void obtainJacobian();
    // bool interpolateSolution(Scalar time, Thyra::VectorBase<Scalar> &tmpSolVector);
    void updateHistory();
    void restoreHistory();
    void updateCoeffs();
    void initialize();
    Scalar checkReduceOrder();
    BDFstatusFlag rejectStep();
    void completeStep();

    void setDefaultMagicNumbers(Teuchos::ParameterList &magicNumberList);

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

#ifdef Rythmos_DEBUG
    int debugLevel;
    Teuchos::RefCountPtr<Teuchos::FancyOStream> debug_out;
#endif // Rythmos_DEBUG


};

// ////////////////////////////
// Defintions

template<class Scalar>
void ImplicitBDFStepper<Scalar>::InitializeStepper(
    Teuchos::ParameterList &implicitBDFParameters
    )
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  // Initialize algorithm coefficients
  maxOrder = implicitBDFParameters.get("maxOrder",int(5)); // maximum order
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

  relErrTol = implicitBDFParameters.get( "relErrTol", Scalar(1.0e-4) );
  absErrTol = implicitBDFParameters.get( "absErrTol", Scalar(1.0e-6) );
  constantStepSize = implicitBDFParameters.get( "constantStepSize", false );
  stopTime = implicitBDFParameters.get( "stopTime", Scalar(10.0) );


#ifdef Rythmos_DEBUG
  debugLevel = implicitBDFParameters.get( "debugLevel", int(1) );
  debug_out = Teuchos::VerboseObjectBase::getDefaultOStream();
  debug_out->precision(15);
  debug_out->setMaxLenLinePrefix(28);
  debug_out->pushLinePrefix("Rythmos::ImplicitBDFStepper");
  debug_out->setShowLinePrefix(true);
  debug_out->setTabIndentStr("    ");
#endif // Rythmos_DEBUG

  setDefaultMagicNumbers(implicitBDFParameters.sublist("magicNumbers"));

#ifdef Rythmos_DEBUG
  Teuchos::OSTab ostab(debug_out,1,"InitializeStepper");
  if (debugLevel > 1)
  {
    *debug_out << "maxOrder = " << maxOrder << endl;
    *debug_out << "currentOrder = " << currentOrder << endl;
    *debug_out << "oldOrder = " << oldOrder << endl;
    *debug_out << "usedOrder = " << usedOrder << endl;
    *debug_out << "alpha_s = " << alpha_s << endl;
    for (int i=0 ; i<maxOrder ; ++i)
    {
      *debug_out << "alpha[" << i << "] = " << alpha[i] << endl;
      *debug_out << "beta[" << i << "] = " << beta[i] << endl;
      *debug_out << "gamma[" << i << "] = " << gamma[i] << endl;
      *debug_out << "psi[" << i << "] = " << psi[i] << endl;
      *debug_out << "sigma[" << i << "] = " << sigma[i] << endl;
    }
    *debug_out << "alpha_0 = " << alpha_0 << endl;
    *debug_out << "cj = " << cj << endl;
    *debug_out << "ck = " << ck << endl;
    *debug_out << "numberOfSteps = " << numberOfSteps << endl;
    *debug_out << "nef = " << nef << endl;
    *debug_out << "usedStep = " << usedStep << endl;
    *debug_out << "nscsco = " << nscsco << endl;
    *debug_out << "Ek = " << Ek << endl;
    *debug_out << "Ekm1 = " << Ekm1 << endl;
    *debug_out << "Ekm2 = " << Ekm2 << endl;
    *debug_out << "Ekp1 = " << Ekp1 << endl;
    *debug_out << "Est = " << Est << endl;
    *debug_out << "Tk = " << Tk << endl;
    *debug_out << "Tkm1 = " << Tkm1 << endl;
    *debug_out << "Tkm2 = " << Tkm2 << endl;
    *debug_out << "Tkp1 = " << Tkp1 << endl;
    *debug_out << "newOrder = " << newOrder << endl;
    *debug_out << "initialPhase = " << initialPhase << endl;
    *debug_out << "relErrTol  = " << relErrTol  << endl;
    *debug_out << "absErrTol  = " << absErrTol  << endl;
    *debug_out << "constantStepSize  = " << constantStepSize  << endl;
    *debug_out << "stopTime  = " << stopTime  << endl;
  }
#endif // Rythmos_DEBUG
}

template<class Scalar>
ImplicitBDFStepper<Scalar>::ImplicitBDFStepper(
  const Teuchos::RefCountPtr<const Thyra::ModelEvaluator<Scalar> > &model
  ,const Teuchos::RefCountPtr<Thyra::NonlinearSolverBase<Scalar> > &solver
  ,Teuchos::ParameterList &parameterList
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

#ifdef Rythmos_DEBUG
  Teuchos::OSTab ostab(debug_out,1,"setDefaultMagicNumbers");
  if (debugLevel > 1)
  {
    *debug_out << "h0_safety = " << h0_safety << endl;
    *debug_out << "h0_max_factor = " << h0_max_factor << endl;
    *debug_out << "h_phase0_incr = " << h_phase0_incr << endl;
    *debug_out << "h_max_inv = " << h_max_inv << endl;
    *debug_out << "Tkm1_Tk_safety = " << Tkm1_Tk_safety  << endl;
    *debug_out << "Tkp1_Tk_safety = " << Tkp1_Tk_safety << endl;
    *debug_out << "r_factor = " << r_factor << endl;
    *debug_out << "r_safety = " << r_safety << endl;
    *debug_out << "r_fudge = " << r_fudge << endl;
    *debug_out << "r_min = " << r_min << endl;
    *debug_out << "r_max = " << r_max << endl;
    *debug_out << "r_hincr_test = " << r_hincr_test << endl;
    *debug_out << "r_hincr = " << r_hincr << endl;
    *debug_out << "max_LET_fail = " << max_LET_fail << endl;
    *debug_out << "minTimeStep = " << minTimeStep << endl;
    *debug_out << "maxTimeStep = " << maxTimeStep << endl;
  }
#endif // Rythmos_DEBUG
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
#ifdef Rythmos_DEBUG
  Teuchos::OSTab ostab(debug_out,1,"TakeStep");
#endif // Rythmos_DEBUG
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
    
#ifdef Rythmos_DEBUG
    if (debugLevel > 1)
    {
      *debug_out << "xn0 = " << std::endl;
      xn0->describe(*debug_out,Teuchos::VERB_EXTREME);
      *debug_out << "ee = " << std::endl;
      ee->describe(*debug_out,Teuchos::VERB_EXTREME);
      *debug_out << "delta = " << std::endl;
      delta->describe(*debug_out,Teuchos::VERB_EXTREME);
      for (int i=0; i<max(2,maxOrder); ++i)
      {
        *debug_out << "xHistory[" << i << "] = " << std::endl;
        xHistory[i]->describe(*debug_out,Teuchos::VERB_EXTREME);
      }
      *debug_out << "ck = " << ck << endl;
      *debug_out << "enorm = " << enorm << endl;
      *debug_out << "Local Truncation Error Check: (ck*enorm) < 1:  (" << ck*enorm << ") <?= 1" << endl;
    }
#endif // Rythmos_DEBUG
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
  if (verbLevel == Teuchos::VERB_EXTREME)
  {
    out << description() << "::describe" << std::endl;
    out << "model = " << std::endl;
    model->describe(out,verbLevel);
    out << "solver = " << std::endl;
    solver->describe(out,verbLevel);
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
    out << "time = " << time << std::endl;
    out << "hh = " << hh << std::endl;
    out << "currentOrder = " << currentOrder << std::endl;
    out << "neModel = " << std::endl;
    neModel.describe(out,verbLevel);
  }
}

template<class Scalar>
void ImplicitBDFStepper<Scalar>::obtainPredictor()
{
  typedef Teuchos::ScalarTraits<Scalar> ST;

#ifdef Rythmos_DEBUG
  Teuchos::OSTab ostab(debug_out,1,"obtainPredictor");
  if (debugLevel > 1)
  {
    *debug_out << "currentOrder = " << currentOrder << std::endl;
  }
#endif // Rythmos_DEBUG
  
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
#ifdef Rythmos_DEBUG
  if (debugLevel > 1)
  {
    *debug_out << "xn0 = " << std::endl;
    xn0->describe(*debug_out,Teuchos::VERB_EXTREME);
    *debug_out << "xpn0 = " << std::endl;
    xpn0->describe(*debug_out,Teuchos::VERB_EXTREME);
  }
#endif // Rythmos_DEBUG
}

/* // This function is not ready yet.
template<class Scalar>
void ImplicitBDFStepper<Scalar>::interpolateSolution(Scalar timepoint, 
    	Thyra::VectorBase<Scalar>  &tmpSolVector)
{
  typedef Teuchos::ScalarTraits<Scalar> ST;

// 03/23/04 tscoffe:  Currently this code is nearly identical to the IDA code
// for interpolating to an output time.  Either we acknowledge the copyright,
// the list of conditions in the license and the disclaimer or we rewrite this
// function.  The IDA license is included after this routine.
  Scalar tfuzz;   // fuzz factor to check for valid output time
  Scalar tp;      // approximately t{n-1}
  Scalar delt;    // distance between timepoint and time
  Scalar c = Scalar(1.0); // coefficient for interpolation
  Scalar gam;     // coefficient for interpolation
  int kord;       // order of interpolation
  Scalar tn = time;
  Scalar hh = hh;
  Scalar hused = usedStep;
  int kused = usedOrder;
  Scalar uround = ST::zero();  // unit round-off (set to zero for now)

  tfuzz = 100 * uround * (tn + hh);
  tp = tn - hused - tfuzz;
  if ( (timepoint - tp)*hh < ST::zero() ) 
    return false;

  assign(&*tmpSolVector,*xHistory[0]),
  kord = kused;
  if ( (kused == 0) || (timepoint == tn) ) 
    kord = 1;

  delt = timepoint - tn;
  gam = delt/psi[0];
  for (int j=1 ; j <= kord ; ++j)
  {
    c = c*gam;
    gam = (delt + psi[j-1])/psi[j];
    Vp_StV(&*tmpSolVector,c,*xHistory[j]);
  }
  return true;
}
*/

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
#ifdef Rythmos_DEBUG
  Teuchos::OSTab ostab(debug_out,1,"updateHistory");
  if (debugLevel > 1)
  {
    for (int i=0;i<max(2,maxOrder);++i)
    {
      *debug_out << "xHistory[" << i << "] = " << endl;
      xHistory[i]->describe(*debug_out,Teuchos::VERB_EXTREME);
    }
  }
#endif // Rythmos_DEBUG

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
#ifdef Rythmos_DEBUG
  Teuchos::OSTab ostab(debug_out,1,"restoreHistory");
  if (debugLevel > 1)
  {
    for (int i=0;i<maxOrder;++i)
      *debug_out << "psi[" << i << "] = " << psi[i] << endl;
    for (int i=0;i<maxOrder;++i)
    {
      *debug_out << "xHistory[" << i << "] = " << endl;
      xHistory[i]->describe(*debug_out,Teuchos::VERB_EXTREME);
    }
  }
#endif // Rythmos_DEBUG
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
#ifdef Rythmos_DEBUG
  Teuchos::OSTab ostab(debug_out,1,"updateCoeffs");
  if (debugLevel > 1)
  {
    for (int i=0;i<=maxOrder;++i)
    {
      *debug_out << "alpha[" << i << "] = " << alpha[i] << endl;
      *debug_out << "beta[" << i << "] = " << beta[i] << endl;
      *debug_out << "sigma[" << i << "] = " << sigma[i] << endl;
      *debug_out << "gamma[" << i << "] = " << gamma[i] << endl;
      *debug_out << "psi[" << i << "] = " << psi[i] << endl;
      *debug_out << "alpha_s = " << alpha_s << endl;
      *debug_out << "alpha_0 = " << alpha_0 << endl;
      *debug_out << "cj = " << cj << endl;
      *debug_out << "ck = " << ck << endl;
    }
  }
#endif // Rythmos_DEBUG

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
#ifdef Rythmos_DEBUG
  Teuchos::OSTab ostab(debug_out,1,"initialize");
  if (debugLevel > 1)
  {
    *debug_out << "hh = " << hh << endl;
  }
#endif // Rythmos_DEBUG

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
#ifdef Rythmos_DEBUG
  Teuchos::OSTab ostab(debug_out,1,"checkReduceOrder");
  if (debugLevel > 1)
  {
    *debug_out << "currentOrder = " << currentOrder << std::endl;
    *debug_out << "Ek = " << Ek << std::endl;
    *debug_out << "Tk = " << Tk << std::endl;
    *debug_out << "enorm = " << enorm << std::endl;
  }
#endif // Rythmos_DEBUG
  if (currentOrder>1)
  {
    V_VpV(&*delta,*xHistory[currentOrder],*ee);
    Ekm1 = sigma[currentOrder-1]*WRMSNorm(*errWtVec,*delta);
    Tkm1 = currentOrder*Ekm1;
#ifdef Rythmos_DEBUG
    if (debugLevel > 1)
    {
      *debug_out << "Ekm1 = " << Ekm1 << endl;
      *debug_out << "Tkm1 = " << Tkm1 << endl;
    }
#endif // Rythmos_DEBUG
    if (currentOrder>2)
    {
      Vp_V(&*delta,*xHistory[currentOrder-1]);
      Ekm2 = sigma[currentOrder-2]*WRMSNorm(*errWtVec,*delta);
      Tkm2 = (currentOrder-1)*Ekm2;
#ifdef Rythmos_DEBUG
      if (debugLevel > 1)
      {
        *debug_out << "Ekm2 = " << Ekm2 << endl;
        *debug_out << "Tkm2 = " << Tkm2 << endl;
      }
#endif // Rythmos_DEBUG
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
#ifdef Rythmos_DEBUG
  if (debugLevel > 1)
  {
    *debug_out << "Est = " << Est << endl;
    *debug_out << "newOrder= " << newOrder << endl;
  }
#endif // Rythmos_DEBUG
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
#ifdef Rythmos_DEBUG
  Teuchos::OSTab ostab(debug_out,1,"rejectStep");
  if (debugLevel > 1)
  {
    *debug_out << "adjustStep = " << adjustStep << endl;
    *debug_out << "nef = " << nef << endl;
  }
#endif // Rythmos_DEBUG
  if (nef >= max_LET_fail)  
  {
    cerr << "Rythmos_Stepper_ImplicitBDF::rejectStep:  " 
          << "  Maximum number of local error test failures.  " << endl;
    return(REP_ERR_FAIL);
  }
  if (adjustStep)
  {
    initialPhase = false;
#ifdef Rythmos_DEBUG
    if (debugLevel > 1)
    {
      *debug_out << "initialPhase = " << initialPhase << endl;
    }
#endif // Rythmos_DEBUG
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
#ifdef Rythmos_DEBUG
    if (debugLevel > 1)
    {
      *debug_out << "rr = " << rr << endl;
      *debug_out << "newOrder = " << newOrder << endl;
    }
#endif // Rythmos_DEBUG
    currentOrder = newOrder;
    if (numberOfSteps == 0) // still first step
    {
      psi[0] = newTimeStep;
      Vt_S(&*xHistory[1],rr);
#ifdef Rythmos_DEBUG
      if (debugLevel > 1)
      {
        *debug_out << "numberOfSteps == 0:" << endl;
        *debug_out << "psi[0] = " << psi[0] << endl;
        *debug_out << "xHistory[1] = " << std::endl;
        xHistory[1]->describe(*debug_out,Teuchos::VERB_EXTREME);
      }
#endif // Rythmos_DEBUG
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
#ifdef Rythmos_DEBUG
  if (debugLevel > 1)
  {
    *debug_out << "hh = " << hh << endl;
  }
#endif // Rythmos_DEBUG
  return(return_status);
}

template<class Scalar>
void ImplicitBDFStepper<Scalar>::completeStep()
{
  typedef Teuchos::ScalarTraits<Scalar> ST;

  numberOfSteps ++;
  nef = 0;
  time += hh;
#ifdef Rythmos_DEBUG
  Teuchos::OSTab ostab(debug_out,1,"completeStep");
  if (debugLevel > 1)
  {
    *debug_out << "numberOfSteps = " << numberOfSteps << endl;
    *debug_out << "nef = " << nef << endl;
    *debug_out << "time = " << time << endl;
  }
#endif // Rythmos_DEBUG
  
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
#ifdef Rythmos_DEBUG
  if (debugLevel > 1)
  {
    *debug_out << "initialPhase = " << initialPhase << endl;
  }
#endif // Rythmos_DEBUG
  if (initialPhase)
  {
    currentOrder++;
    newTimeStep = h_phase0_incr * hh;
#ifdef Rythmos_DEBUG
    if (debugLevel > 1)
    {
      *debug_out << "currentOrder = " << currentOrder << endl;
      *debug_out << "newTimeStep = " << newTimeStep << endl;
    }
#endif // Rythmos_DEBUG
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
#ifdef Rythmos_DEBUG
      if (debugLevel > 1)
      {
        *debug_out << "delta = " << endl;
        delta->describe(*debug_out,Teuchos::VERB_EXTREME);
        *debug_out << "Tkp1 = ||delta||_WRMS = " << Tkp1 << endl;
        *debug_out << "Ekp1 = " << Ekp1 << endl;
      }
#endif // Rythmos_DEBUG
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
#ifdef Rythmos_DEBUG
    if (debugLevel > 1)
    {
      *debug_out << "Est = " << Est << endl;
      *debug_out << "currentOrder = " << currentOrder << endl;
      *debug_out << "rr  = " << rr << endl;
      *debug_out << "newTimeStep = " << newTimeStep << endl;
    }
#endif // Rythmos_DEBUG
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
#ifdef Rythmos_DEBUG
  if (debugLevel > 1)
  {
    *debug_out << "hh = " << hh << endl;
  }
#endif // Rythmos_DEBUG
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
    const std::vector<Scalar>& time_list
    ,const std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > >& x_list
    ,const std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > >& xdot_list)
{
  return(false);
}

template<class Scalar>
bool ImplicitBDFStepper<Scalar>::GetPoints(
    const std::vector<Scalar>& time_list
    ,std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > >* x_list
    ,std::vector<Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > >* xdot_list
    ,std::vector<ScalarMag>* accuracy_list) const
{
  return(false);
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
bool ImplicitBDFStepper<Scalar>::GetNodes(std::vector<Scalar>* time_list) const
{
  return(false);
}

template<class Scalar>
bool ImplicitBDFStepper<Scalar>::RemoveNodes(std::vector<Scalar>& time_list) const
{
  return(false);
}

template<class Scalar>
int ImplicitBDFStepper<Scalar>::GetOrder() const
{
  return(currentOrder);
}


} // namespace Rythmos

#endif //Rythmos_IMPLICITBDF_STEPPER_H
