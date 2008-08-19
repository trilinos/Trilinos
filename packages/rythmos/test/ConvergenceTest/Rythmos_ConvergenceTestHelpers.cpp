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

#include "Rythmos_ConvergenceTestHelpers.hpp"
#include "../SinCos/SinCosModel.hpp"

namespace Rythmos {

// We are using the SinCosModel due to the need to get an exact solution.  They
// will fail the rcp_dynamic_cast if the underlying model is not a SinCosModel
// TODO:  Add some facility to the ModelEvaluator to provide this e.g.:
// OUT_ARG_exact_solution returns InArg object
double computeOrderByLocalErrorConvergenceStudy(const StepperFactoryBase<double>& stepperFactory)
{
//  Array<double> stepSize;
  Array<double> logStepSize;

//  Array<double> errorNorm;
  Array<double> logErrorNorm;

  double h = 0.5;
  int N = 10;
  for (int i=0 ; i<N ; ++i) {
    RCP<StepperBase<double> > stepper = stepperFactory.create();
    double stepTaken = stepper->takeStep(h,STEP_TYPE_FIXED);
    TEUCHOS_ASSERT_EQUALITY( stepTaken, h );
    // get solution 
    RCP<const VectorBase<double> > x = stepper->getStepStatus().solution;
    // get exact solution
    RCP<const VectorBase<double> > x_exact;
    {
      // This is the only place where we're using the SinCosModel specifically
      RCP<const SinCosModel> sinCosModel = Teuchos::rcp_dynamic_cast<const SinCosModel>(stepper->getModel(),true);
      x_exact = sinCosModel->getExactSolution(h).get_x();
    }
    RCP<VectorBase<double> > tmp_vec = createMember(stepper->getModel()->get_x_space());
    // take norm of difference
    Thyra::V_StVpStV(&*tmp_vec,1.0,*x_exact,-1.0,*x);
//    stepSize.push_back(h);
    logStepSize.push_back(log(h));
    double nrm = Thyra::norm_inf(*tmp_vec);
//    errorNorm.push_back(nrm);
    logErrorNorm.push_back(log(nrm));
    h = h/2.0;
  }

//  using std::cout;
//  using std::endl;
//  cout << endl;
//  cout << "stepSize = " << stepSize << endl;
//  cout << "logStepSize = " << logStepSize << endl;
//  cout << "errorNorm = " << errorNorm << endl;
//  cout << "logErrorNorm = " << logErrorNorm << endl;

  double slope = computeLinearRegressionSlope<double>(logStepSize,logErrorNorm);
//  cout << "slope = " << slope << endl;

  return slope;
}

// We are using the SinCosModel due to the need to get an exact solution.  They
// will fail the rcp_dynamic_cast if the underlying model is not a SinCosModel
// TODO:  Add some facility to the ModelEvaluator to provide this e.g.:
// OUT_ARG_exact_solution returns InArg object
double computeOrderByGlobalErrorConvergenceStudy(const StepperFactoryBase<double>& stepperFactory)
{
//  Array<double> stepSize;
  Array<double> logStepSize;

//  Array<double> errorNorm;
  Array<double> logErrorNorm;

  double t_final = 1.0;


  double h = 0.5;
  int N = 10;
  for (int i=0 ; i<N ; ++i) {
    RCP<StepperBase<double> > stepper = stepperFactory.create();
    int numSteps = int(round(t_final/h));
    for (int s=0 ; s<numSteps ; ++s) {
      double stepTaken = stepper->takeStep(h,STEP_TYPE_FIXED);
      TEUCHOS_ASSERT_EQUALITY( stepTaken, h );
    }
    // get solution 
    RCP<const VectorBase<double> > x = stepper->getStepStatus().solution;
    double t = stepper->getStepStatus().time;
    TEUCHOS_ASSERT_EQUALITY( t, t_final );
    RCP<VectorBase<double> > tmp_vec = createMember(stepper->getModel()->get_x_space());
    // get exact solution 
    RCP<const VectorBase<double> > x_exact;
    {
      // This is the only place where we're using the SinCosModel specifically
      RCP<const SinCosModel> sinCosModel = Teuchos::rcp_dynamic_cast<const SinCosModel>(stepper->getModel(),true);
      x_exact = sinCosModel->getExactSolution(t_final).get_x();
    }
    // take norm of difference
    Thyra::V_StVpStV(&*tmp_vec,1.0,*x_exact,-1.0,*x);
//    stepSize.push_back(h);
    logStepSize.push_back(log(h));
    double nrm = Thyra::norm_inf(*tmp_vec);
//    errorNorm.push_back(nrm);
    logErrorNorm.push_back(log(nrm));
    h = h/2.0;
  }

//  using std::cout;
//  using std::endl;
//  cout << endl;
//  cout << "stepSize = " << stepSize << endl;
//  cout << "logStepSize = " << logStepSize << endl;
//  cout << "errorNorm = " << errorNorm << endl;
//  cout << "logErrorNorm = " << logErrorNorm << endl;

  double slope = computeLinearRegressionSlope<double>(logStepSize,logErrorNorm);
//  cout << "slope = " << slope << endl;

  return slope;
}

} // namespace Rythmos



