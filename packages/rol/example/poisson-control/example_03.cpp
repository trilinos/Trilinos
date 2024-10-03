// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// Burgers includes
#include "example_02.hpp"
// ROL includes
#include "ROL_StdVector.hpp"
#include "ROL_StdTeuchosBatchManager.hpp"
#include "ROL_MonteCarloGenerator.hpp"
#include "ROL_Reduced_Objective_SimOpt.hpp"
#include "ROL_ParameterList.hpp"
#include "ROL_Solver.hpp"
#include "ROL_StochasticProblem.hpp"
#include "ROL_PrimalDualRisk.hpp"

// Teuchos includes
#include "Teuchos_Time.hpp"
#include "ROL_Stream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_Comm.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_CommHelpers.hpp"

int main( int argc, char *argv[] ) {  

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  auto comm = ROL::toPtr( Teuchos::DefaultComm<int>::getComm() );

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  ROL::Ptr<std::ostream> outStream;
  ROL::nullstream bhs; // outputs nothing
  if (iprint > 0 && comm->getRank()==0)
    outStream = ROL::makePtrFromRef(std::cout);
  else
    outStream = ROL::makePtrFromRef(bhs);

  int errorFlag  = 0;

  // *** Example body.

  try {

    /***************************************************************************/
    /***************** GRAB INPUTS *********************************************/
    /***************************************************************************/
    // Get finite element parameter list
    std::string filename = "input_ex03.xml";
    auto parlist = ROL::getParametersFromXmlFile( filename );

    if ( parlist->sublist("Problem Data").get("Display Option",0) && (comm->getRank() > 0) ) {
      parlist->set("Display Option",0);
    }
  
    /***************************************************************************/
    /***************** INITIALIZE SAMPLERS *************************************/
    /***************************************************************************/
    int dim    = 2;
    int nSamp  = parlist->sublist("Problem Data").get("Number of Monte Carlo Samples",1000);
    std::vector<double> tmp(2); tmp[0] = -1.0; tmp[1] = 1.0;
    std::vector<std::vector<double> > bounds(dim,tmp);
    ROL::Ptr<ROL::BatchManager<double> > bman
      = ROL::makePtr<ROL::StdTeuchosBatchManager<double,int>>(comm);
    ROL::Ptr<ROL::SampleGenerator<double> > sampler
      = ROL::makePtr<ROL::MonteCarloGenerator<double>>(nSamp,bounds,bman,false);
  
    /***************************************************************************/
    /***************** INITIALIZE CONTROL VECTOR *******************************/
    /***************************************************************************/
    int nx = parlist->sublist("Problem Data").get("Number of Elements", 128);
    ROL::Ptr<std::vector<double> > z_ptr = ROL::makePtr<std::vector<double>>(nx+1, 0.0);
    ROL::Ptr<ROL::Vector<double> > z = ROL::makePtr<ROL::StdVector<double>>(z_ptr);
    ROL::Ptr<ROL::Vector<double> > u = ROL::makePtr<ROL::StdVector<double>>(nx-1);
    ROL::Ptr<ROL::Vector<double> > p = ROL::makePtr<ROL::StdVector<double>>(nx-1);
  
    /***************************************************************************/
    /***************** INITIALIZE OBJECTIVE FUNCTION ***************************/
    /***************************************************************************/
    double alpha = parlist->sublist("Problem Data").get("Penalty Parameter", 1.e-4);
    ROL::Ptr<FEM<double> > fem = ROL::makePtr<FEM<double>>(nx);
    ROL::Ptr<ROL::Objective_SimOpt<double> > pObj
      = ROL::makePtr<DiffusionObjective<double>>(fem, alpha);
    ROL::Ptr<ROL::Constraint_SimOpt<double> > pCon
      = ROL::makePtr<DiffusionConstraint<double>>(fem);
    ROL::Ptr<ROL::Objective<double> > robj
      = ROL::makePtr<ROL::Reduced_Objective_SimOpt<double>>(pObj,pCon,u,z,p);
    robj->setParameter({0.0,0.0});
  
    /***************************************************************************/
    /***************** INITIALIZE ROL ALGORITHM ********************************/
    /***************************************************************************/
    bool runBundle = parlist->sublist("Problem Data").get("Run Bundle",false);
    // Solve using bundle
    if (runBundle) {
      z->zero();
      ROL::Ptr<ROL::StochasticProblem<double>> problem2
        = ROL::makePtr<ROL::StochasticProblem<double>>(robj, z);
      problem2->makeObjectiveStochastic(*parlist, sampler);
      problem2->finalize(false,true,*outStream);
      parlist->sublist("Step").set("Type","Bundle");
      parlist->sublist("Step").sublist("Bundle").set("Distance Measure Coefficient",0.0);
      ROL::Solver<double> solver2(problem2,*parlist);
      solver2.solve(*outStream);
    }
    
    ROL::Ptr<ROL::Problem<double>> problem
      = ROL::makePtr<ROL::Problem<double>>(robj, z);
    ROL::PrimalDualRisk<double> solver(problem, sampler, *parlist);
    if (parlist->sublist("Problem Data").get("Run Derivative Check",false)) {
      problem->check(true,*outStream);
      solver.check(*outStream);
    }
    solver.run(*outStream);
  
    /***************************************************************************/
    /***************** PRINT RESULTS *******************************************/
    /***************************************************************************/
    int my_number_samples = sampler->numMySamples(), number_samples = 0;
    Teuchos::reduceAll<int,int>(*comm,Teuchos::REDUCE_SUM,1,&my_number_samples,&number_samples);
    int my_number_solves  = ROL::dynamicPtrCast<DiffusionConstraint<double> >(pCon)->getNumSolves(), number_solves = 0;
    Teuchos::reduceAll<int,int>(*comm,Teuchos::REDUCE_SUM,1,&my_number_solves,&number_solves);
    if (comm->getRank() == 0) {
      std::cout << "Number of Samples    = " << number_samples << "\n";
      std::cout << "Number of Solves     = " << number_solves  << "\n";
    }
  
    if ( comm->getRank() == 0 ) {
      std::ofstream file;
      file.open("control.txt");
      std::vector<double> xmesh(fem->nz(),0.0);
      fem->build_mesh(xmesh);
      for (int i = 0; i < fem->nz(); i++ ) {
        file << std::setprecision(std::numeric_limits<double>::digits10) << std::scientific << xmesh[i] << "  "  
             << std::setprecision(std::numeric_limits<double>::digits10) << std::scientific << (*z_ptr)[i] 
             << "\n";
      }
      file.close();
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




