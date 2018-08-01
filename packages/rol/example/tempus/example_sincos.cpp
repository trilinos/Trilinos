// @HEADER
// ************************************************************************
//
//               Rapid Optimization Library (ROL) Package
//                 Copyright (2014) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact lead developers:
//              Drew Kouri   (dpkouri@sandia.gov) and
//              Denis Ridzal (dridzal@sandia.gov)
//
// ************************************************************************
// @HEADER

#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_DefaultComm.hpp"

#include "example_sincos.hpp"
#include "ROL_TempusReducedObjective.hpp"
#include "ROL_OptimizationProblem.hpp"
#include "ROL_OptimizationSolver.hpp"
#include "ROL_RandomVector.hpp"

// ************************************************************
// ************************************************************
int main(int argc, char *argv[])
{
  using Teuchos::RCP;
  using Teuchos::ParameterList;
  using Teuchos::sublist;
  using Teuchos::getParametersFromXmlFile;

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  // Save ROL output to rol.out
  std::ofstream fout;
  fout.rdbuf()->pubsetbuf(0, 0); // Set's to unbuffered (must be called before open())
  fout.open("rol.out");
  ROL::Ptr<std::ostream> outStream = ROL::makePtrFromRef(fout);

  int errorFlag = 0;

  try {

    // Read Tempus and model params from .xml file
    RCP<ParameterList> pList =
      getParametersFromXmlFile("example_sincos.xml");

    // Setup the SinCosModel
    RCP<ParameterList> scm_pl = sublist(pList, "SinCosModel", true);
    RCP<SinCosModel<double> > model =
      Teuchos::rcp(new SinCosModel<double>(scm_pl));

    // Create ROL objective
    RCP<ParameterList> pl = sublist(pList, "Tempus", true);
    RCP<ParameterList> objective_params = sublist(pList, "Objective", true);
    ROL::Ptr<ROL::TempusReducedObjective<double> > objective =
      ROL::makePtr<ROL::TempusReducedObjective<double> >(
        model, pl, objective_params);

    // Create target -- do forward integration with perturbed parameter values
    RCP<ROL::Vector<double> > zs = objective->create_design_vector();
    zs->setScalar(1.5);
    RCP<ROL::Vector<double> > r = objective->create_response_vector();
    objective->run_tempus(*r, *zs);
    objective->set_target(r);

    // Check gradient
    const bool check_gradient = true;
    if (check_gradient) {
      RCP<ROL::Vector<double> > x = objective->create_design_vector();
      RCP<ROL::Vector<double> > v = objective->create_design_vector();
      srand(12345);
      ROL::RandomizeVector(*v,-1.0,1.0);
      const int numSteps = 6;
      const double largestStep = 0.1;
      const double stepReduction = 0.1;
      const int fdOrder = 1;
      std::vector<double> fdSteps(numSteps,largestStep);
      for (int i=1; i<numSteps; ++i)
        fdSteps[i] = stepReduction * fdSteps[i-1];
      objective->checkGradient(*x, *v, fdSteps, true, *outStream, fdOrder);
    }

    // Create ROL optimizer
    RCP<ParameterList> rol_pl = sublist(pList, "ROL", true);
    RCP<ROL::Vector<double> > z = objective->create_design_vector();
    ROL::OptimizationProblem<double> problem(objective, z);
    ROL::OptimizationSolver<double> solver(problem, *rol_pl);
    solver.solve(*outStream);

    RCP<ROL::ThyraVector<double> > zt =
      Teuchos::rcp_dynamic_cast<ROL::ThyraVector<double> >(z,true);
    *outStream << "Final solution: [ ";
    const int n = z->dimension();
    for (int i=0; i<n; ++i)
      *outStream << Thyra::get_ele(*(zt->getVector()),i) << " ";
    *outStream << "]" << std::endl;

    // Calculate the error from expected solution
    z->axpy(-1.0, *zs);
    double nrm = z->norm();
    *outStream << "Error in solution:  " << nrm << std::endl;
    if (nrm > 1.0e-4)
      errorFlag = 1;

  }
  catch (std::logic_error err) {
    *outStream << err.what() << "\n";
    errorFlag = -1000;
  }; // end try

  if (errorFlag != 0)
    std::cout << "End Result: TEST FAILED\n";
  else
    std::cout << "End Result: TEST PASSED\n";

  return 0;
}
