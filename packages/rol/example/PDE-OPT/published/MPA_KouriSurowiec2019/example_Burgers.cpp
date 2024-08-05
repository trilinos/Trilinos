// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "example_Burgers.hpp"

int main(int argc, char* argv[]) {
  //feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
  using RealT = double;

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint = argc - 1;
  ROL::Ptr<std::ostream> outStream;
  ROL::nullstream bhs; // outputs nothing

  /*** Initialize communicator. ***/
  Teuchos::GlobalMPISession mpiSession (&argc, &argv, &bhs);
  ROL::Ptr<const Teuchos::Comm<int>> comm = Tpetra::getDefaultComm();

  const int myRank = comm->getRank();
  if (iprint > 0 && myRank == 0) outStream = ROL::makePtrFromRef(std::cout);
  else                           outStream = ROL::makePtrFromRef(bhs);
  int errorFlag  = 0;

  try {
    // Get ROL parameterlist
    auto parlist = ROL::getParametersFromXmlFile("input_Burgers.xml");
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
    std::string method = parlist->sublist("Problem").get("Method","PD-Risk");
    if (method == "Bundle") {
      ROL::Ptr<ROL::StochasticProblem<RealT>> problem
        = ROL::makePtr<ROL::StochasticProblem<RealT>>(robj, z);
      problem->makeObjectiveStochastic(*parlist, sampler);
      problem->finalize(false,true,*outStream);
      if (parlist->sublist("Problem").get("Run Derivative Check",false)) {
        problem->check(*outStream);
      }
      parlist->sublist("Step").set("Type","Bundle");
      parlist->sublist("Step").sublist("Bundle").set("Distance Measure Coefficient",0.0);
      ROL::Solver<RealT> solver(problem,*parlist);
      solver.solve(*outStream);
    }
    else if (method == "Epi-Reg") {
      ROL::ParameterList list = *parlist;
      std::string rm = list.sublist("SOL").sublist("Objective").sublist("Risk Measure").get("Name","CVaR");
      list.sublist("Step").set("Type","Trust Region");
      list.sublist("SOL").sublist("Objective").sublist("Risk Measure").sublist(rm).sublist("Distribution").set("Name","Parabolic");
      list.sublist("SOL").sublist("Objective").sublist("Risk Measure").sublist(rm).sublist("Distribution").sublist("Parabolic").set("Lower Bound",-0.5);
      list.sublist("SOL").sublist("Objective").sublist("Risk Measure").sublist(rm).sublist("Distribution").sublist("Parabolic").set("Upper Bound", 0.5);

      RealT eps0 = parlist->sublist("SOL").sublist("Primal Dual Risk").get("Initial Gradient Tolerance",1e-2);
      RealT eps1 = parlist->sublist("Status Test").get("Gradient Tolerance",1e-10);
      RealT epss = parlist->sublist("SOL").sublist("Primal Dual Risk").get("Solver Tolerance Update Scale",1e-1);
      RealT rho0 = parlist->sublist("SOL").sublist("Primal Dual Risk").get("Initial Penalty Parameter",1e1);
      RealT rho1 = parlist->sublist("SOL").sublist("Primal Dual Risk").get("Maximum Penalty Parameter",1e16);
      RealT rhos = parlist->sublist("SOL").sublist("Primal Dual Risk").get("Penalty Update Scale",1e1);
      int maxIt = static_cast<int>(std::floor(std::log(rho1/rho0)/std::log(rhos)));

      ROL::Ptr<ROL::Vector<RealT>> zp = z->clone(); zp->set(*z);
      RealT eps = eps0, rho = rho0, tau = 1.0/rho, stat(1), statp(1), err(1);
      for (int i = 0; i < maxIt; ++i) {
        // Set tolerances
        list.sublist("SOL").sublist("Objective").sublist("Risk Measure").sublist(rm).set("Smoothing Parameter",tau);
        list.sublist("Status Test").set("Gradient Tolerance",eps);
        // Set statistic
        list.sublist("SOL").sublist("Objective").set("Initial Statistic",stat);
        // Solver smoothed problem
        ROL::Ptr<ROL::StochasticProblem<RealT>> problem
          = ROL::makePtr<ROL::StochasticProblem<RealT>>(robj, z);
        problem->makeObjectiveStochastic(list, sampler);
        problem->finalize(false,true,*outStream);
        ROL::Solver<RealT> solver(problem,list);
        solver.solve(*outStream);
        // Get solution statistic
        stat = problem->getSolutionStatistic();
        // Compute iteration errors
        zp->axpy(-1.0,*z);
        err = std::sqrt(zp->dot(*zp) + std::pow(stat-statp,2));
        *outStream << std::endl << std::endl;
        *outStream << "  iter = " << i
                   << "  rho = " << rho
                   << "  eps = " << eps
                   << "  err = " << err;
        *outStream << std::endl << std::endl;
        if (eps <= eps1 && err <= eps1) break;
        zp->set(*z);
        statp = stat;
        // Update tolerances
        eps *= epss;
        rho *= rhos;
        tau  = 1.0/rho;
        if (eps < eps1) eps = eps1;
      }
    }
    else {
      ROL::Ptr<ROL::Problem<RealT>> problem
        = ROL::makePtr<ROL::Problem<RealT>>(robj, z);
      problem->finalize(false,true,*outStream);
      ROL::PrimalDualRisk<RealT> solver(problem, sampler, *parlist);
      if (parlist->sublist("Problem").get("Run Derivative Check",false)) {
        std::vector<RealT> param(4,0);
        robj->setParameter(param);
        problem->check(*outStream);
        solver.check(*outStream);
      }
      solver.run(*outStream);
    }
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
