// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "example_10.hpp"

typedef double RealT;

int main(int argc, char* argv[]) {

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  ROL::Ptr<const Teuchos::Comm<int>> comm
    = ROL::toPtr(Teuchos::DefaultComm<int>::getComm());

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint = argc - 1;
  ROL::Ptr<std::ostream> outStream;
  ROL::nullstream bhs; // outputs nothing
  if (iprint > 0 && Teuchos::rank<int>(*comm)==0)
    outStream = ROL::makePtrFromRef(std::cout);
  else
    outStream = ROL::makePtrFromRef(bhs);

  int errorFlag  = 0;

  try {
    // Get ROL parameterlist
    auto parlist = ROL::getParametersFromXmlFile("input_ex10.xml");
    /**********************************************************************************************/
    /************************* CONSTRUCT VECTORS **************************************************/
    /**********************************************************************************************/
    int nx = 256;
    ROL::Ptr<ROL::Vector<RealT>> z = ROL::makePtr<ROL::StdVector<RealT>>(nx+2,0.0);
    ROL::Ptr<ROL::Vector<RealT>> u = ROL::makePtr<ROL::StdVector<RealT>>(nx,1.0);
    ROL::Ptr<ROL::Vector<RealT>> p = ROL::makePtr<ROL::StdVector<RealT>>(nx,0.0);
    /**********************************************************************************************/
    /************************* CONSTRUCT SOL COMPONENTS *******************************************/
    /**********************************************************************************************/
    // Build samplers
    int dim = 4, nSamp = parlist->sublist("Problem").get("Number of Samples",100);
    std::vector<RealT> tmp = {-1, 1};
    std::vector<std::vector<RealT>> bounds(dim,tmp);
    ROL::Ptr<ROL::BatchManager<RealT>> bman
      = ROL::makePtr<ROL::StdTeuchosBatchManager<RealT,int>>(comm);
    ROL::Ptr<ROL::SampleGenerator<RealT>> sampler
      = ROL::makePtr<ROL::MonteCarloGenerator<RealT>>(nSamp,bounds,bman);
    /**********************************************************************************************/
    /************************* CONSTRUCT OBJECTIVE FUNCTION ***************************************/
    /**********************************************************************************************/
    // Build risk-averse objective function
    RealT alpha = 1.e-3;
    ROL::Ptr<ROL::Objective_SimOpt<RealT>> objSimOpt
      = ROL::makePtr<Objective_BurgersControl<RealT>>(alpha,nx);
    ROL::Ptr<ROL::Constraint_SimOpt<RealT>> conSimOpt
      = ROL::makePtr<Constraint_BurgersControl<RealT>>(nx);
    conSimOpt->setSolveParameters(*parlist);
    ROL::Ptr<ROL::Objective<RealT>> robj
      = ROL::makePtr<ROL::Reduced_Objective_SimOpt<RealT>>(objSimOpt,conSimOpt,u,z,p);
    /**********************************************************************************************/
    /************************* SOLVE OPTIMIZATION PROBLEM *****************************************/
    /**********************************************************************************************/
    bool runBundle = parlist->sublist("Problem").get("Run Bundle",false);
    // Solve using bundle
    if (runBundle) {
      z->zero();
      ROL::Ptr<ROL::OptimizationProblem<double>> problem2
        = ROL::makePtr<ROL::OptimizationProblem<double>>(robj, z);
      problem2->setStochasticObjective(*parlist, sampler);
      parlist->sublist("Step").set("Type","Bundle");
      parlist->sublist("Step").sublist("Bundle").set("Distance Measure Coefficient",0.0);
      ROL::OptimizationSolver<double> solver2(*problem2,*parlist);
      solver2.solve(*outStream);
    }
    
    ROL::Ptr<ROL::Problem<double>> problem
      = ROL::makePtr<ROL::Problem<double>>(robj, z);
    ROL::PrimalDualRisk<double> solver(problem, sampler, *parlist);
    if (parlist->sublist("Problem").get("Run Derivative Check",false)) {
      problem->check(true,*outStream);
      solver.check(*outStream);
    }
    solver.run(*outStream);
  }
  catch (std::logic_error& err) {
    *outStream << err.what() << "\n";
    errorFlag = -1000;
  }; // end try

  if (errorFlag != 0)
    std::cout << "End Result: TEST FAILED\n";
  else
    std::cout << "End Result: TEST PASSED\n";

  return 0;
}
