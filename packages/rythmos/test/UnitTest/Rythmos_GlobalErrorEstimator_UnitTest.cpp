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

#include "Teuchos_UnitTestHarness.hpp"

//#include "Rythmos_GlobalErrorEstimator.hpp"
#include "Rythmos_Types.hpp"
#include "Rythmos_IntegratorBuilder.hpp"
#include "../SinCos/SinCosModel.hpp"
#include "Rythmos_TimeStepNonlinearSolver.hpp"
#include "Rythmos_PointwiseInterpolationBufferAppender.hpp"
#include "Rythmos_CubicSplineInterpolator.hpp"
#include "Rythmos_InterpolationBuffer.hpp"
#include "Rythmos_TrailingInterpolationBufferAcceptingIntegratorBase.hpp"
#include "Rythmos_AdjointModelEvaluator.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_LinearNonlinearSolver.hpp"
#include "Thyra_DetachedVectorView.hpp"

namespace Rythmos {

//TEUCHOS_UNIT_TEST( Rythmos_GlobalErrorEstimator, create ) {
//  RCP<GlobalErrorEstimator<double> > gee = globalErrorEstimator<double>();
//  TEST_ASSERT( !is_null(gee) );
//}

TEUCHOS_UNIT_TEST( Rythmos_GlobalErrorEstimator, SinCos ) {
  typedef Teuchos::ScalarTraits<double> ST;
  // Forward Solve, storing data in linear interpolation buffer
  int storageLimit = 100;
  double finalTime = 0.1;
  double dt = 0.1;
  RCP<IntegratorBuilder<double> > ib = integratorBuilder<double>();
  {
    RCP<ParameterList> ibPL = Teuchos::parameterList();
    ibPL->sublist("Integrator Settings").sublist("Integrator Selection").set("Integrator Type","Default Integrator");
    ibPL->sublist("Integrator Settings").set("Final Time",finalTime);
    ibPL->sublist("Integration Control Strategy Selection").set("Integration Control Strategy Type","Simple Integration Control Strategy");
    ibPL->sublist("Integration Control Strategy Selection").sublist("Simple Integration Control Strategy").set("Take Variable Steps",false);
    ibPL->sublist("Integration Control Strategy Selection").sublist("Simple Integration Control Strategy").set("Fixed dt",dt);

    ibPL->sublist("Stepper Settings").sublist("Stepper Selection").set("Stepper Type","Backward Euler");
    //ibPL->sublist("Stepper Settings").sublist("Stepper Selection").set("Stepper Type","Implicit RK");
    //ibPL->sublist("Stepper Settings").sublist("Runge Kutta Butcher Tableau Selection").set("Runge Kutta Butcher Tableau Type","Backward Euler");
    ibPL->sublist("Interpolation Buffer Settings").sublist("Trailing Interpolation Buffer Selection").set("Interpolation Buffer Type","Interpolation Buffer");
    ibPL->sublist("Interpolation Buffer Settings").sublist("Trailing Interpolation Buffer Selection").sublist("Interpolation Buffer").set("StorageLimit",storageLimit);
    ibPL->sublist("Interpolation Buffer Settings").sublist("Interpolator Selection").set("Interpolator Type","Linear Interpolator");
    ib->setParameterList(ibPL);
  }
  RCP<SinCosModel> fwdModel = sinCosModel(true); // implicit formulation
  Thyra::ModelEvaluatorBase::InArgs<double> fwdIC = fwdModel->getNominalValues();
  RCP<Thyra::NonlinearSolverBase<double> > fwdNLSolver = timeStepNonlinearSolver<double>();
  RCP<IntegratorBase<double> > fwdIntegrator = ib->create(fwdModel,fwdIC,fwdNLSolver);
  RCP<const VectorBase<double> > x_final;
  {
    Array<double> time_vec;
    time_vec.push_back(finalTime);
    Array<RCP<const Thyra::VectorBase<double> > > x_final_array;
    fwdIntegrator->getFwdPoints(time_vec,&x_final_array,NULL,NULL);
    x_final = x_final_array[0];
  }
  // Verify x_final is correct
  {
    // Defaults from SinCos Model:
    double f = 1.0;
    double L = 1.0;
    double a = 0.0;
    double x_ic_0 = 0.0;
    double x_ic_1 = 1.0;
    double x_0 = dt/(1.0+std::pow(dt*f/L,2))*(x_ic_0/dt+x_ic_1+dt*std::pow(f/L,2)*a);
    double x_1 = dt/(1.0+std::pow(dt*f/L,2))*(-std::pow(f/L,2)*x_ic_0+x_ic_1/dt+std::pow(f/L,2)*a);
    double tol = 1.0e-10;
    Thyra::ConstDetachedVectorView<double> x_final_view( *x_final );
    TEST_FLOATING_EQUALITY( x_final_view[0], x_0, tol );
    TEST_FLOATING_EQUALITY( x_final_view[1], x_1, tol );
  }
  // Copy InterpolationBuffer data into Cubic Spline interpolation buffer for use in Adjoint Solve
  TimeRange<double> fwdTimeRange; 
  RCP<InterpolationBufferBase<double> > fwdCubicSplineInterpBuffer;
  {
    RCP<PointwiseInterpolationBufferAppender<double> > piba = pointwiseInterpolationBufferAppender<double>();
    RCP<InterpolationBuffer<double> > sinkInterpBuffer = interpolationBuffer<double>();
    sinkInterpBuffer->setStorage(storageLimit);
    RCP<CubicSplineInterpolator<double> > csi = cubicSplineInterpolator<double>();
    sinkInterpBuffer->setInterpolator(csi);
    RCP<const InterpolationBufferBase<double> > sourceInterpBuffer;
    {
      RCP<TrailingInterpolationBufferAcceptingIntegratorBase<double> > tibaib = 
        Teuchos::rcp_dynamic_cast<TrailingInterpolationBufferAcceptingIntegratorBase<double> >(fwdIntegrator,true);
      sourceInterpBuffer = tibaib->getTrailingInterpolationBuffer();
    }
    fwdTimeRange = sourceInterpBuffer->getTimeRange();
    piba->append(*sourceInterpBuffer, fwdTimeRange, Teuchos::outArg(*sinkInterpBuffer));
    fwdCubicSplineInterpBuffer = sinkInterpBuffer;

    TimeRange<double> sourceRange = sourceInterpBuffer->getTimeRange();
    TimeRange<double> sinkRange = sinkInterpBuffer->getTimeRange();
    TEST_EQUALITY( sourceRange.lower(), sinkRange.lower() );
    TEST_EQUALITY( sourceRange.upper(), sinkRange.upper() );
  }
  // Adjoint Solve, reading forward solve data from Cubic Spline interpolation buffer
  {
    RCP<ParameterList> ibPL = Teuchos::parameterList();
    ibPL->sublist("Integrator Settings").sublist("Integrator Selection").set("Integrator Type","Default Integrator");
    ibPL->sublist("Integrator Settings").set("Final Time",finalTime);
    ibPL->sublist("Integration Control Strategy Selection").set("Integration Control Strategy Type","Simple Integration Control Strategy");
    ibPL->sublist("Integration Control Strategy Selection").sublist("Simple Integration Control Strategy").set("Take Variable Steps",false);
    ibPL->sublist("Integration Control Strategy Selection").sublist("Simple Integration Control Strategy").set("Fixed dt",dt);

    ibPL->sublist("Stepper Settings").sublist("Stepper Selection").set("Stepper Type","Backward Euler");
    //ibPL->sublist("Stepper Settings").sublist("Stepper Selection").set("Stepper Type","Implicit RK");
    //ibPL->sublist("Stepper Settings").sublist("Runge Kutta Butcher Tableau Selection").set("Runge Kutta Butcher Tableau Type","Implicit 1 Stage 2nd order Gauss");
    ibPL->sublist("Interpolation Buffer Settings").sublist("Trailing Interpolation Buffer Selection").set("Interpolation Buffer Type","Interpolation Buffer");
    ibPL->sublist("Interpolation Buffer Settings").sublist("Trailing Interpolation Buffer Selection").sublist("Interpolation Buffer").set("StorageLimit",storageLimit);
    ibPL->sublist("Interpolation Buffer Settings").sublist("Interpolator Selection").set("Interpolator Type","Linear Interpolator");
    ib->setParameterList(ibPL);
  }
  RCP<Thyra::ModelEvaluator<double> > adjModel;
  {
    RCP<Rythmos::AdjointModelEvaluator<double> > model = 
      Rythmos::adjointModelEvaluator<double>(
          fwdModel, fwdTimeRange
          );
    //model->setFwdStateSolutionBuffer(fwdCubicSplineInterpBuffer);
    adjModel = model;
  }
  Thyra::ModelEvaluatorBase::InArgs<double> adjIC = adjModel->getNominalValues();
  double phi_ic_0 = 2.0;
  double phi_ic_1 = 3.0;
  {
    // Initial conditions for adjoint:
    const RCP<const Thyra::VectorSpaceBase<double> >
      f_space = fwdModel->get_f_space();
    const RCP<Thyra::VectorBase<double> > x_ic = createMember(f_space);
    {
      Thyra::DetachedVectorView<double> x_ic_view( *x_ic );
      x_ic_view[0] = phi_ic_0;
      x_ic_view[1] = phi_ic_1;
    }
    const RCP<Thyra::VectorBase<double> > xdot_ic = createMember(f_space);
    V_S( Teuchos::outArg(*xdot_ic), ST::zero() );
    adjIC.set_x(x_ic);
    adjIC.set_x_dot(xdot_ic);
  }
  RCP<Thyra::LinearNonlinearSolver<double> > adjNLSolver = Thyra::linearNonlinearSolver<double>();
  RCP<IntegratorBase<double> > adjIntegrator = ib->create(adjModel,adjIC,adjNLSolver);
  RCP<const VectorBase<double> > phi_final;
  {
    Array<double> time_vec;
    time_vec.push_back(finalTime);
    Array<RCP<const Thyra::VectorBase<double> > > phi_final_array;
    adjIntegrator->getFwdPoints(time_vec,&phi_final_array,NULL,NULL);
    phi_final = phi_final_array[0];
  }
  // Verify phi_final is correct
  {
    // Defaults from SinCos Model:
    double f = 1.0;
    double L = 1.0;
    double h = -dt;
    double phi_0 = 1.0/(1.0+std::pow(f*h/L,2.0))*(phi_ic_0+std::pow(f/L,2.0)*h*phi_ic_1);
    double phi_1 = 1.0/(1.0+std::pow(f*h/L,2.0))*(-h*phi_ic_0+phi_ic_1);
    double tol = 1.0e-10;
    Thyra::ConstDetachedVectorView<double> phi_final_view( *phi_final );
    TEST_FLOATING_EQUALITY( phi_final_view[0], phi_0, tol );
    TEST_FLOATING_EQUALITY( phi_final_view[1], phi_1, tol );
  }
  // Compute error estimate
  //TEST_ASSERT( false );
}


/*
class ABaseClass
{
  public:
  virtual ~ABaseClass() {}
  std::string print() { return "ABaseClass"; }
};
class ADerivedClass : virtual public ABaseClass
{
  public:
  ADerivedClass() {}
  ~ADerivedClass() {}
  std::string print() { return "ADerivedClass"; }
};
class VBaseClass
{
  public:
  virtual ~VBaseClass() {}
  virtual std::string print() { return "VBaseClass"; }
};
class VDerivedClass : virtual public VBaseClass
{
  public:
  VDerivedClass() {}
  ~VDerivedClass() {}
  std::string print() { return "VDerivedClass"; }
};
class ADerivedVDerivedClass : virtual public VDerivedClass
{
  public:
  ADerivedVDerivedClass() {}
  ~ADerivedVDerivedClass() {}
  std::string print() { return "ADerivedVDerivedClass"; }
};

TEUCHOS_UNIT_TEST( Rythmos_GlobalErrorEstimator, DerivedClass ) {
  {
    RCP<ABaseClass> ac;
    {
      RCP<ADerivedClass> dc = Teuchos::rcp(new ADerivedClass());
      TEST_EQUALITY( dc->print(), "ADerivedClass" );
      ac = dc;
    }
    TEST_EQUALITY( ac->print(), "ADerivedClass" );
  }

  {
    RCP<VBaseClass> av;
    {
      RCP<VDerivedClass> dv = Teuchos::rcp(new VDerivedClass());
      TEST_EQUALITY( dv->print(), "VDerivedClass" );
      av = dv;
    }
    TEST_EQUALITY( av->print(), "VDerivedClass" );
  }

  {
    RCP<VBaseClass> av;
    {
      RCP<ADerivedVDerivedClass> dv = Teuchos::rcp(new ADerivedVDerivedClass());
      TEST_EQUALITY( dv->print(), "ADerivedVDerivedClass" );
      av = dv;
    }
    TEST_EQUALITY( av->print(), "ADerivedVDerivedClass" );
  }

  {
    RCP<VDerivedClass> av;
    {
      RCP<ADerivedVDerivedClass> dv = Teuchos::rcp(new ADerivedVDerivedClass());
      TEST_EQUALITY( dv->print(), "ADerivedVDerivedClass" );
      av = dv;
    }
    TEST_EQUALITY( av->print(), "ADerivedVDerivedClass" );
  }
}
*/

} // namespace Rythmos 


