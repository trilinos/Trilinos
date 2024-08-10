// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  example_07.cpp
    \brief Shows how to solve a steady Burgers' optimal control problem using
           full-space methods.
*/

#include "ROL_OptimizationSolver.hpp"
#include "ROL_Reduced_Objective_SimOpt.hpp"
#include "ROL_MonteCarloGenerator.hpp"
#include "ROL_ParameterList.hpp"

#include "ROL_Stream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_Comm.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_CommHelpers.hpp"

#include <iostream>
#include <fstream>
#include <algorithm>

#include "example_07.hpp"

typedef double RealT;
typedef H1VectorPrimal<RealT> PrimalStateVector;
typedef H1VectorDual<RealT> DualStateVector;
typedef L2VectorPrimal<RealT> PrimalControlVector;
typedef L2VectorDual<RealT> DualControlVector;
typedef H1VectorDual<RealT> PrimalConstraintVector;
typedef H1VectorPrimal<RealT> DualConstraintVector;

int main(int argc, char *argv[]) {

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  ROL::Ptr<const Teuchos::Comm<int>> comm
    = ROL::toPtr(Teuchos::DefaultComm<int>::getComm());

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint = argc - 1;
  bool print = (iprint>0);
  ROL::Ptr<std::ostream> outStream;
  ROL::nullstream bhs; // outputs nothing
  if (print)
    outStream = ROL::makePtrFromRef(std::cout);
  else
    outStream = ROL::makePtrFromRef(bhs);

  bool print0 = print && !(comm->getRank());
  ROL::Ptr<std::ostream> outStream0;
  if (print0)
    outStream0 = ROL::makePtrFromRef(std::cout);
  else
    outStream0 = ROL::makePtrFromRef(bhs);

  int errorFlag  = 0;

  // *** Example body.

  try {
    /*************************************************************************/
    /************* INITIALIZE BURGERS FEM CLASS ******************************/
    /*************************************************************************/
    int nx    = 512;   // Set spatial discretization.
    RealT x   = 0.0;   // Set penalty parameter.
    RealT nl  = 1.0;   // Nonlinearity parameter (1 = Burgers, 0 = linear).
    RealT cH1 = 1.0;   // Scale for derivative term in H1 norm.
    RealT cL2 = 0.0;   // Scale for mass term in H1 norm.
    ROL::Ptr<BurgersFEM<RealT>> fem
      = ROL::makePtr<BurgersFEM<RealT>>(nx,nl,cH1,cL2);
    fem->test_inverse_mass(*outStream0);
    fem->test_inverse_H1(*outStream0);
    /*************************************************************************/
    /************* INITIALIZE SIMOPT OBJECTIVE FUNCTION **********************/
    /*************************************************************************/
    ROL::Ptr<ROL::Objective_SimOpt<RealT>> pobj
      = ROL::makePtr<Objective_BurgersControl<RealT>>(fem,x);
    /*************************************************************************/
    /************* INITIALIZE SIMOPT EQUALITY CONSTRAINT *********************/
    /*************************************************************************/
    bool hess = true;
    ROL::Ptr<ROL::Constraint_SimOpt<RealT>> pcon
      = ROL::makePtr<Constraint_BurgersControl<RealT>>(fem,hess);
    /*************************************************************************/
    /************* INITIALIZE VECTOR STORAGE *********************************/
    /*************************************************************************/
    ROL::Ptr<std::vector<RealT>> z_ptr, u_ptr, c_ptr, l_ptr;
    z_ptr = ROL::makePtr<std::vector<RealT>>(nx+2, 0.0);
    u_ptr = ROL::makePtr<std::vector<RealT>>(nx, 1.0);
    c_ptr = ROL::makePtr<std::vector<RealT>>(nx, 0.0);
    l_ptr = ROL::makePtr<std::vector<RealT>>(nx, 0.0);
    ROL::Ptr<ROL::Vector<RealT>> zp, up, cp, lp;
    zp = ROL::makePtr<PrimalControlVector>(z_ptr,fem);
    up = ROL::makePtr<PrimalStateVector>(u_ptr,fem);
    cp = ROL::makePtr<PrimalConstraintVector>(c_ptr,fem);
    lp = ROL::makePtr<DualConstraintVector>(l_ptr,fem);
    /*************************************************************************/
    /************* INITIALIZE SAMPLE GENERATOR *******************************/
    /*************************************************************************/
    int dim = 4, nSamp = 1000;
    std::vector<RealT> tmp(2,0.0); tmp[0] = -1.0; tmp[1] = 1.0;
    std::vector<std::vector<RealT>> bounds(dim,tmp);
    ROL::Ptr<ROL::BatchManager<RealT>> bman
      = ROL::makePtr<L2VectorBatchManager<RealT,int>>(comm);
    ROL::Ptr<ROL::SampleGenerator<RealT>> sampler
      = ROL::makePtr<ROL::MonteCarloGenerator<RealT>>(
          nSamp,bounds,bman,false,false,100);
    /*************************************************************************/
    /************* INITIALIZE OBJECTIVE FUNCTION *****************************/
    /*************************************************************************/
    bool storage = true, fdhess = false;
    ROL::Ptr<ROL::Objective<RealT>> robj
      = ROL::makePtr<ROL::Reduced_Objective_SimOpt<RealT>>(
          pobj,pcon,up,zp,lp,storage,fdhess);
    /*************************************************************************/
    /************* INITIALIZE BOUND CONSTRAINTS ******************************/
    /*************************************************************************/
    std::vector<RealT> Zlo(nx+2,0.0), Zhi(nx+2,10.0);
    for (int i = 0; i < nx+2; i++) {
      if ( i < (int)((nx+2)/3) ) {
        Zlo[i] = -1.0;
        Zhi[i] = 1.0;
      }
      if ( i >= (int)((nx+2)/3) && i < (int)(2*(nx+2)/3) ) {
        Zlo[i] = 1.0;
        Zhi[i] = 5.0;
      }
      if ( i >= (int)(2*(nx+2)/3) ) {
        Zlo[i] = 5.0;
        Zhi[i] = 10.0;
      }
    }
    ROL::Ptr<ROL::BoundConstraint<RealT>> bnd
      = ROL::makePtr<L2BoundConstraint<RealT>>(Zlo,Zhi,fem);
    /*************************************************************************/
    /************* INITIALIZE RISK-AVERSE OPTIMIZATION PROBLEM ***************/
    /*************************************************************************/
    RealT order = 2.0, threshold = -0.85*(1.0-x);
    ROL::Ptr<ROL::ParameterList> bpoelist = ROL::makePtr<ROL::ParameterList>();
    bpoelist->sublist("SOL").set("Store Sampled Value and Gradient",true);
    bpoelist->sublist("SOL").set("Type","Probability");
    bpoelist->sublist("SOL").sublist("Probability").set("Name","bPOE");
    bpoelist->sublist("SOL").sublist("Probability").sublist("bPOE").set("Threshold",threshold);
    bpoelist->sublist("SOL").sublist("Probability").sublist("bPOE").set("Moment Order",order);
    ROL::OptimizationProblem<RealT> problem(robj,zp,bnd);
    problem.setStochasticObjective(*bpoelist,sampler);
    // CHECK OBJECTIVE DERIVATIVES
    bool derivcheck = false;
    if (derivcheck) {
      problem.check(*outStream0);
    }
    /*************************************************************************/
    /************* RUN OPTIMIZATION ******************************************/
    /*************************************************************************/
    // READ IN XML INPUT
    std::string filename = "input.xml";
    auto parlist = ROL::getParametersFromXmlFile( filename );
    // RUN OPTIMIZATION
    ROL::OptimizationSolver<RealT> solver(problem,*parlist);
    solver.solve(*outStream);
    /*************************************************************************/
    /************* PRINT CONTROL AND STATE TO SCREEN *************************/
    /*************************************************************************/
    if ( print0 ) {
      std::ofstream ofs;
      ofs.open("output_example_09.txt",std::ofstream::out);
      for ( int i = 0; i < nx+2; i++ ) {
        ofs << std::scientific << std::setprecision(10);
        ofs << std::setw(20) << std::left << (RealT)i/((RealT)nx+1.0);
        ofs << std::setw(20) << std::left << (*z_ptr)[i];
        ofs << "\n";
      }
      ofs.close();
    }
    *outStream0 << "Scalar Parameter: " << problem.getSolutionStatistic() << "\n\n";
  }
  catch (std::logic_error& err) {
    *outStream << err.what() << "\n";
    errorFlag = -1000;
  }; // end try

  comm->barrier();
  if (errorFlag != 0)
    std::cout << "End Result: TEST FAILED\n";
  else
    std::cout << "End Result: TEST PASSED\n";

  return 0;
}
