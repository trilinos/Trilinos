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
#include "ROL_Algorithm.hpp"
#include "ROL_TrustRegionStep.hpp"
#include "ROL_LineSearchStep.hpp"
#include "ROL_StdVector.hpp"
#include "ROL_StdTeuchosBatchManager.hpp"
#include "ROL_MonteCarloGenerator.hpp"
#include "ROL_Reduced_Objective_SimOpt.hpp"
#include "ROL_RiskNeutralObjective.hpp"
#include "ROL_Vector_SimOpt.hpp"
#include "ROL_Bounds.hpp"
#include "ROL_ParameterList.hpp"

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
  if (iprint > 0)
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
    std::string filename = "example_02.xml";
    auto parlist = ROL::getParametersFromXmlFile( filename );

    if ( parlist->get("Display Option",0) && (comm->getRank() > 0) ) {
      parlist->set("Display Option",0);
    }
    // Get ROL parameter list
    filename = "input.xml";
    auto ROL_parlist = ROL::getParametersFromXmlFile( filename );
  
    /***************************************************************************/
    /***************** INITIALIZE SAMPLERS *************************************/
    /***************************************************************************/
    int dim    = 2;
    bool useSA = parlist->get("Use Stochastic Approximation",false);
    int nSamp  = 1;
    if ( !useSA ) {
      nSamp  = parlist->get("Number of Monte Carlo Samples",1000);
    }
    std::vector<double> tmp(2); tmp[0] = -1.0; tmp[1] = 1.0;
    std::vector<std::vector<double> > bounds(dim,tmp);
    ROL::Ptr<ROL::BatchManager<double> > bman
      = ROL::makePtr<ROL::StdTeuchosBatchManager<double,int>>(comm);
    ROL::Ptr<ROL::SampleGenerator<double> > sampler
      = ROL::makePtr<ROL::MonteCarloGenerator<double>>(nSamp,bounds,bman,useSA);
  
    /***************************************************************************/
    /***************** INITIALIZE CONTROL VECTOR *******************************/
    /***************************************************************************/
    int nx = parlist->get("Number of Elements", 128);
    ROL::Ptr<std::vector<double> > z_ptr = ROL::makePtr<std::vector<double>>(nx+1, 0.0);
    ROL::Ptr<ROL::Vector<double> > z = ROL::makePtr<ROL::StdVector<double>>(z_ptr);
    ROL::Ptr<std::vector<double> > u_ptr = ROL::makePtr<std::vector<double>>(nx-1, 0.0);
    ROL::Ptr<ROL::Vector<double> > u = ROL::makePtr<ROL::StdVector<double>>(u_ptr);
    ROL::Vector_SimOpt<double> x(u,z);
    ROL::Ptr<std::vector<double> > p_ptr = ROL::makePtr<std::vector<double>>(nx-1, 0.0);
    ROL::Ptr<ROL::Vector<double> > p = ROL::makePtr<ROL::StdVector<double>>(p_ptr);
    ROL::Ptr<std::vector<double> > U_ptr = ROL::makePtr<std::vector<double>>(nx+1, 35.0);
    ROL::Ptr<ROL::Vector<double> > U = ROL::makePtr<ROL::StdVector<double>>(U_ptr);
    ROL::Ptr<std::vector<double> > L_ptr = ROL::makePtr<std::vector<double>>(nx+1, -5.0);
    ROL::Ptr<ROL::Vector<double> > L = ROL::makePtr<ROL::StdVector<double>>(L_ptr);
    ROL::Bounds<double> bnd(L,U);
  
    /***************************************************************************/
    /***************** INITIALIZE OBJECTIVE FUNCTION ***************************/
    /***************************************************************************/
    double alpha = parlist->get("Penalty Parameter", 1.e-4);
    ROL::Ptr<FEM<double> > fem = ROL::makePtr<FEM<double>>(nx);
    ROL::Ptr<ROL::Objective_SimOpt<double> > pObj
      = ROL::makePtr<DiffusionObjective<double>>(fem, alpha);
    ROL::Ptr<ROL::Constraint_SimOpt<double> > pCon
      = ROL::makePtr<DiffusionConstraint<double>>(fem);
    ROL::Ptr<ROL::Objective<double> > robj
      = ROL::makePtr<ROL::Reduced_Objective_SimOpt<double>>(pObj,pCon,u,z,p);
    ROL::RiskNeutralObjective<double> obj(robj,sampler);
  
    /***************************************************************************/
    /***************** RUN DERIVATIVE CHECK ************************************/
    /***************************************************************************/
    if (parlist->get("Run Derivative Check",false)) {
      // Direction to test finite differences
      ROL::Ptr<std::vector<double> > dz_ptr = ROL::makePtr<std::vector<double>>(nx+1, 0.0);
      ROL::Ptr<ROL::Vector<double> > dz = ROL::makePtr<ROL::StdVector<double>>(dz_ptr);
      ROL::Ptr<std::vector<double> > du_ptr = ROL::makePtr<std::vector<double>>(nx-1, 0.0);
      ROL::Ptr<ROL::Vector<double> > du = ROL::makePtr<ROL::StdVector<double>>(du_ptr);
      ROL::Vector_SimOpt<double> d(du,dz);
      // Set to random vectors
      srand(12345);
      for (int i=0; i<nx+1; i++) {
        (*dz_ptr)[i] = 2.0*(double)rand()/(double)RAND_MAX - 1.0;
        (*z_ptr)[i] = 2.0*(double)rand()/(double)RAND_MAX - 1.0;
      }
      for (int i=0; i<nx-1; i++) {
        (*du_ptr)[i] = 2.0*(double)rand()/(double)RAND_MAX - 1.0;
        (*u_ptr)[i] = 2.0*(double)rand()/(double)RAND_MAX - 1.0;
      }
      // Run derivative checks
      std::vector<double> param(dim,0.0);
      robj->setParameter(param);
      if ( comm->getRank() == 0 ) {
        std::cout << "\nRUN DERIVATIVE CHECK FOR PARAMETRIZED OBJECTIVE FUNCTION SIMOPT\n";
      }
      pObj->checkGradient(x,d,(comm->getRank()==0));
      pObj->checkHessVec(x,d,(comm->getRank()==0));
      if ( comm->getRank() == 0 ) {
        std::cout << "\nRUN DERIVATIVE CHECK FOR PARAMETRIZED EQUALITY CONSTRAINT SIMOPT\n";
      }
      pCon->checkApplyJacobian(x,d,*p,(comm->getRank()==0));
      pCon->checkApplyAdjointJacobian(x,*du,*p,x,(comm->getRank()==0));
      pCon->checkApplyAdjointHessian(x,*du,d,x,(comm->getRank()==0));
      if ( comm->getRank() == 0 ) {
        std::cout << "\nRUN DERIVATIVE CHECK FOR PARAMETRIZED OBJECTIVE FUNCTION\n";
      }
      robj->checkGradient(*z,*dz,(comm->getRank()==0));
      robj->checkHessVec(*z,*dz,(comm->getRank()==0));
      // Run derivative checks
      if ( comm->getRank() == 0 ) {
        std::cout << "\nRUN DERIVATIVE CHECK FOR RISK-NEUTRAL OBJECTIVE FUNCTION\n";
      }
      obj.checkGradient(*z,*dz,(comm->getRank()==0));
      obj.checkHessVec(*z,*dz,(comm->getRank()==0));
    }
  
    /***************************************************************************/
    /***************** INITIALIZE ROL ALGORITHM ********************************/
    /***************************************************************************/
    ROL::Ptr<ROL::Algorithm<double>>  algo; 
    ROL::Ptr<ROL::Step<double>>       step;
    ROL::Ptr<ROL::StatusTest<double>> status;
    if ( useSA ) {
      ROL_parlist->sublist("General").set("Recompute Objective Function",false);
      ROL_parlist->sublist("Step").sublist("Line Search").set("Initial Step Size",0.1/alpha);
      ROL_parlist->sublist("Step").sublist("Line Search").set("User Defined Initial Step Size",true);
      ROL_parlist->sublist("Step").sublist("Line Search").sublist("Line-Search Method").set("Type","Iteration Scaling");
      ROL_parlist->sublist("Step").sublist("Line Search").sublist("Descent Method").set("Type","Steepest Descent");
      ROL_parlist->sublist("Step").sublist("Line Search").sublist("Curvature Condition").set("Type","Null Curvature Condition");
      status = ROL::makePtr<ROL::StatusTest<double>>(*ROL_parlist);
      step   = ROL::makePtr<ROL::LineSearchStep<double>>(*ROL_parlist);
      algo   = ROL::makePtr<ROL::Algorithm<double>>(step,status,false);
    } 
    else {
      status = ROL::makePtr<ROL::StatusTest<double>>(*ROL_parlist);
      step   = ROL::makePtr<ROL::TrustRegionStep<double>>(*ROL_parlist);
      algo   = ROL::makePtr<ROL::Algorithm<double>>(step,status,false);
    }
  
    /***************************************************************************/
    /***************** PERFORM OPTIMIZATION ************************************/
    /***************************************************************************/
    Teuchos::Time timer("Optimization Time",true);
    z->zero();
    algo->run(*z,obj,bnd,(comm->getRank()==0));
    double optTime = timer.stop();
  
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
      std::cout << "Optimization Time    = " << optTime        << "\n\n";
    }
  
    if ( comm->getRank() == 0 ) {
      std::ofstream file;
      if (useSA) {
        file.open("control_SA.txt");
      }
      else {
        file.open("control_SAA.txt");
      }
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




