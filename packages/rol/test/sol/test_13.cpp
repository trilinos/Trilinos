// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "ROL_Stream.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include "ROL_ParameterList.hpp"
#include "ROL_LinearRegression.hpp"
#include "ROL_Solver.hpp"
#include "ROL_MonteCarloGenerator.hpp"

typedef double RealT;

int main(int argc, char* argv[]) {

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  ROL::Ptr<std::ostream> outStream;
  ROL::nullstream bhs; // outputs nothing
  if (iprint > 0)
    outStream = ROL::makePtrFromRef(std::cout);
  else
    outStream = ROL::makePtrFromRef(bhs);

  int errorFlag  = 0;

  try {
    // Set up problem data
    int nsamp = 1000, fdim = 4;
    std::vector<ROL::Ptr<ROL::Distribution<RealT>>> dist(fdim+1);
    ROL::ParameterList list0;
    list0.sublist("Distribution").set("Name","Beta");
    list0.sublist("Distribution").sublist("Beta").set("Shape 1",5.0);
    list0.sublist("Distribution").sublist("Beta").set("Shape 2",2.0);
    dist[0] = ROL::DistributionFactory<RealT>(list0);
    ROL::ParameterList list1;
    list1.sublist("Distribution").set("Name","Exponential");
    list1.sublist("Distribution").sublist("Exponential").set("Location",0.0);
    list1.sublist("Distribution").sublist("Exponential").set("Scale",   1.0);
    dist[1] = ROL::DistributionFactory<RealT>(list1);
    ROL::ParameterList list2;
    list2.sublist("Distribution").set("Name","Gaussian");
    list2.sublist("Distribution").sublist("Gaussian").set("Mean",    1.0);
    list2.sublist("Distribution").sublist("Gaussian").set("Variance",2.0);
    dist[2] = ROL::DistributionFactory<RealT>(list2);
    ROL::ParameterList list3;
    list3.sublist("Distribution").set("Name","Uniform");
    list3.sublist("Distribution").sublist("Uniform").set("Lower Bound",0.5);
    list3.sublist("Distribution").sublist("Uniform").set("Upper Bound",0.75);
    dist[3] = ROL::DistributionFactory<RealT>(list3);
    ROL::ParameterList list4;
    list4.sublist("Distribution").set("Name","Triangle");
    list4.sublist("Distribution").sublist("Triangle").set("Lower Bound",  0.25);
    list4.sublist("Distribution").sublist("Triangle").set("Peak Location",0.5);
    list4.sublist("Distribution").sublist("Triangle").set("Upper Bound",  1.25);
    dist[4] = ROL::DistributionFactory<RealT>(list4);
    ROL::Ptr<ROL::BatchManager<RealT>> bman
      = ROL::makePtr<ROL::BatchManager<RealT>>();
    ROL::Ptr<ROL::SampleGenerator<RealT>> data
      = ROL::makePtr<ROL::MonteCarloGenerator<RealT>>(nsamp,dist,bman);
    ROL::LinearRegression<RealT> linReg(data);
    // Set up linear regression solver
    ROL::ParameterList parlist;
    ROL::Ptr<ROL::Problem<RealT>> problem;
    ROL::Ptr<ROL::Solver<RealT>>  solver;
    parlist.sublist("Status Test").set("Gradient Tolerance",1e-8);
    parlist.sublist("Status Test").set("Step Tolerance",    1e-12);
    parlist.sublist("Status Test").set("Iteration Limit",   100);
    parlist.sublist("SOL").set("Type","Error");
    std::vector<std::vector<RealT>> coeff;
    for (ROL::EErrorMeasure ed = ROL::ERRORMEASURE_MEANVARIANCEQUADRANGLE; ed != ROL::ERRORMEASURE_LAST; ed++) {
      std::string name = ROL::EErrorMeasureToString(ed);
      *outStream << name << std::endl;
      parlist.sublist("SOL").sublist("Error Measure").set("Name",name);
      if (ed == ROL::ERRORMEASURE_MEANVARIANCEQUADRANGLE) {
        parlist.sublist("Step").set("Type","Trust Region");
        parlist.sublist("Step").sublist("Trust Region").set("Subproblem Solver","Truncated CG");
      }
      else if (ed == ROL::ERRORMEASURE_TRUNCATEDMEANQUADRANGLE) {
        ROL::ParameterList &list
          = parlist.sublist("SOL").sublist("Error Measure").sublist("Huber");
        list.set("Threshold",1e-1);
        parlist.sublist("Step").set("Type","Trust Region");
        parlist.sublist("Step").sublist("Trust Region").set("Subproblem Solver","Truncated CG");
      }
      else if (ed == ROL::ERRORMEASURE_QUANTILEQUADRANGLE) {
        ROL::ParameterList &list
          = parlist.sublist("SOL").sublist("Error Measure").sublist("Koenker-Bassett");
        list.set("Confidence Level",0.75);
        list.set("Convex Combination Parameter",0.0);
        list.set("Smoothing Parameter",1e-4);
        parlist.sublist("Step").set("Type","Bundle");
        parlist.sublist("Step").sublist("Bundle").set("Distance Measure Coefficient",0.0);
      }
      else if (ed == ROL::ERRORMEASURE_MOREAUYOSIDACVAR) {
        ROL::ParameterList &list
          = parlist.sublist("SOL").sublist("Error Measure").sublist("Moreau-Yosida-Koenker-Bassett");
        list.set("Confidence Level",0.75);
        list.set("Smoothing Parameter",1e-2);
        parlist.sublist("Step").set("Type","Trust Region");
        parlist.sublist("Step").sublist("Trust Region").set("Subproblem Solver","Truncated CG");
      }
      else if (ed == ROL::ERRORMEASURE_GENMOREAUYOSIDACVAR) {
        ROL::ParameterList &list
          = parlist.sublist("SOL").sublist("Error Measure").sublist("Generalized Moreau-Yosida-Koenker-Bassett");
        list.set("Confidence Level",0.75);
        list.set("Convex Combination Parameter",0.0);
        list.set("Smoothing Parameter",1e-2);
        parlist.sublist("Step").set("Type","Trust Region");
        parlist.sublist("Step").sublist("Trust Region").set("Subproblem Solver","Truncated CG");
      }
      else if (ed == ROL::ERRORMEASURE_LOGEXPONENTIALQUADRANGLE) {
        ROL::ParameterList &list
          = parlist.sublist("SOL").sublist("Error Measure").sublist("Exponential");
        list.set("Rate",1.0);
        parlist.sublist("Step").set("Type","Trust Region");
        parlist.sublist("Step").sublist("Trust Region").set("Subproblem Solver","Truncated CG");
      }
      else if (ed == ROL::ERRORMEASURE_LOGQUANTILEQUADRANGLE) {
        ROL::ParameterList &list
          = parlist.sublist("SOL").sublist("Error Measure").sublist("Log Quantile");
        list.set("Slope for Linear Growth",0.5);
        list.set("Rate for Exponential Growth",1.0);
        list.set("Smoothing Parameter",1e-2);
        parlist.sublist("Step").set("Type","Bundle");
        parlist.sublist("Step").sublist("Bundle").set("Distance Measure Coefficient",0.0);
      }
      else if (ed == ROL::ERRORMEASURE_SMOOTHEDWORSTCASEQUADRANGLE) {
        ROL::ParameterList &list
          = parlist.sublist("SOL").sublist("Error Measure").sublist("Smoothed Worst Case");
        list.set("Smoothing Parameter",1e-4);
        parlist.sublist("Step").set("Type","Trust Region");
        parlist.sublist("Step").sublist("Trust Region").set("Subproblem Solver","Truncated CG");
      }
      linReg.setErrorMeasure(parlist,true);
      problem = linReg.getProblem();
      problem->check(true,*outStream);
      solver  = ROL::makePtr<ROL::Solver<RealT>>(problem,parlist);
      solver->solve(*outStream);
      coeff.push_back(*linReg.getCoefficients());
      linReg.print(*outStream);
      linReg.reset();
    }
    int cnt = 0;
    *outStream << std::endl << "Summary: Coefficients" << std::endl;
    *outStream << std::scientific << std::setprecision(6);
    for (ROL::EErrorMeasure ed = ROL::ERRORMEASURE_MEANVARIANCEQUADRANGLE; ed != ROL::ERRORMEASURE_LAST; ed++) {
      std::string name = ROL::EErrorMeasureToString(ed);
      *outStream << std::endl << "  " << std::setw(45) << std::left << name;
      for (int j = 0; j < fdim+1; ++j) {
        *outStream << std::setw(15) << std::left << coeff[cnt][j];
      }
      cnt++;
    }
    *outStream << std::endl << std::endl;
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
