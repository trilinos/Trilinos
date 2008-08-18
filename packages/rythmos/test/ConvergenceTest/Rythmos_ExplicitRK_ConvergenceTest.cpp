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

#include "Rythmos_Types.hpp"
#include "Rythmos_ConvergenceTestHelpers.hpp"

#include "Rythmos_ExplicitRKStepper.hpp"

#include "../SinCos/SinCosModel.hpp"

namespace Rythmos {

using Thyra::VectorBase;
using Thyra::VectorSpaceBase;
using Teuchos::is_null;


TEUCHOS_UNIT_TEST( Rythmos_ExplicitRKStepper, convergenceStudy ) {
  double tol = 1.0e-10;
  int order = 4;

  RCP<SinCosModel> model = sinCosModel(false);
  RCP<VectorBase<double> > tmp_vec = createMember(model->get_x_space());

  Array<double> stepSize;
  Array<double> logStepSize;

  Array<double> errorNorm;
  Array<double> logErrorNorm;

  double t_final = 1.0;

  // get exact solution 
  RCP<const VectorBase<double> > x_exact = model->getExactSolution(t_final).get_x();

  double h = 1.0;
  int N = 10;
  for (int i=0 ; i<N ; ++i) {
    RCP<ExplicitRKStepper<double> > stepper = explicitRKStepper<double>(model);
    int numSteps = int(round(t_final/h));
    for (int s=0 ; s<numSteps ; ++s) {
      double stepTaken = stepper->takeStep(h,STEP_TYPE_FIXED);
      TEST_FLOATING_EQUALITY( stepTaken, h, tol );
    }
    // get solution 
    RCP<const VectorBase<double> > x = stepper->getStepStatus().solution;
    double t = stepper->getStepStatus().time;
    TEST_FLOATING_EQUALITY( t, t_final, tol );
    // take norm of difference
    Thyra::V_StVpStV(&*tmp_vec,1.0,*x_exact,-1.0,*x);
    stepSize.push_back(h);
    logStepSize.push_back(log10(h));
    double nrm = Thyra::norm_inf(*tmp_vec);
    errorNorm.push_back(nrm);
    logErrorNorm.push_back(log10(nrm));
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

  TEST_COMPARE( slope, >=, 1.0*order ); // is slope > order?
  tol = 1.0e-3;
  TEST_FLOATING_EQUALITY( slope, 1.0*order, tol ); // is slope close to order?
}



} // namespace Rythmos

