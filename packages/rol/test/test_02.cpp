// @HEADER
// ************************************************************************
// @HEADER


/*! \file  test_02.cpp
    \brief Test Rosenbrock.
*/

#include "ROL_TestObjectives.hpp"
#include "ROL_LineSearchStep.hpp"
#include "ROL_TrustRegionStep.hpp"
#include "ROL_Algorithm.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include <iostream>

typedef double RealT;

int main(int argc, char *argv[]) {

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  Teuchos::RCP<std::ostream> outStream;
  Teuchos::oblackholestream bhs; // outputs nothing
  if (iprint > 0)
    outStream = Teuchos::rcp(&std::cout, false);
  else
    outStream = Teuchos::rcp(&bhs, false);

  int errorFlag  = 0;

  // *** Test body.

  try {

    // Initial Guess Vector 
    Teuchos::RCP<std::vector<RealT> > x0_rcp = Teuchos::rcp( new std::vector<RealT> );
    ROL::StdVector<RealT> x0(x0_rcp);

    // Exact Solution Vector
    Teuchos::RCP<std::vector<RealT> > z_rcp = Teuchos::rcp( new std::vector<RealT> );
    ROL::StdVector<RealT> z(z_rcp);

    // Get Objective Function
    ROL::ETestObjectives objFunc = ROL::TESTOBJECTIVES_ROSENBROCK;
    //ROL::ETestObjectives objFunc = ROL::TESTOBJECTIVES_SUMOFSQUARES;
    //ROL::ETestObjectives objFunc = ROL::TESTOBJECTIVES_LEASTSQUARES;
    Teuchos::RCP<ROL::Objective<RealT> > obj = Teuchos::null;
    ROL::getTestObjectives<RealT>(obj,x0,z,objFunc);

    // Get Dimension of Problem
    int dim = 
      Teuchos::rcp_const_cast<std::vector<RealT> >((Teuchos::dyn_cast<ROL::StdVector<RealT> >(x0)).getVector())->size();

    // Iteration Vector
    Teuchos::RCP<std::vector<RealT> > x_rcp = Teuchos::rcp( new std::vector<RealT> (dim, 0.0) );
    ROL::StdVector<RealT> x(x_rcp);
    x.set(x0);

    // Random Direction for Derivative Checks
    Teuchos::RCP<std::vector<RealT> > y_rcp = Teuchos::rcp( new std::vector<RealT> (dim, 0.0) );
    ROL::StdVector<RealT> y(y_rcp);
    for (int i=0; i<dim; i++) {
      (*y_rcp)[i] = ((RealT)rand())/((RealT)RAND_MAX);
    }

    // Error Vector
    Teuchos::RCP<std::vector<RealT> > e_rcp = Teuchos::rcp( new std::vector<RealT> (dim, 0.0) );
    ROL::StdVector<RealT> e(e_rcp);
    e.zero();

    // Run Finite Difference Checks
    obj->checkGradient(x,y,true);

    /* BEGIN SECANT DEFINITION */
    ROL::SecantType Stype = ROL::Secant_lDFP;
    //ROL::SecantType Stype = ROL::Secant_lSR1;
    //ROL::SecantType Stype = ROL::Secant_lBFGS;
    //ROL::SecantType Stype = ROL::Secant_BarzilaiBorwein;

    int L        = 10;
    int BBtype   = 1;
    /* END SECANT DEFINTION */

    /* BEGIN LINE SEARCH STEP DEFINTION */
    //ROL::LineSearchStepType LSStype = ROL::LineSearchStep_Newton;
    ROL::LineSearchStepType LSStype = ROL::LineSearchStep_NewtonKrylov;
    //ROL::LineSearchStepType LSStype = ROL::LineSearchStep_NewtonKrylovSecantPreconditioning;
    //ROL::LineSearchStepType LSStype = ROL::LineSearchStep_Secant;
    //ROL::LineSearchStepType LSStype = ROL::LineSearchStep_Gradient;

    ROL::LineSearchType LStype = ROL::LineSearchType_Backtracking;
    //ROL::LineSearchType LStype = ROL::LineSearchType_SimpleBacktracking;
    //ROL::LineSearchType LStype = ROL::LineSearchType_Brents;
    //ROL::LineSearchType LStype = ROL::LineSearchType_Bisection;
    //ROL::LineSearchType LStype = ROL::LineSearchType_GoldenSection;

    ROL::LineSearchCondition LScond = ROL::LineSearchCondition_Wolfe;
    //ROL::LineSearchCondition LScond = ROL::LineSearchCondition_StrongWolfe;
    //ROL::LineSearchCondition LScond = ROL::LineSearchCondition_Goldstein;
 
    int maxit    = 20;
    RealT rho    = 0.5;
    RealT c1     = 1.e-4;
    RealT c2     = 0.9;
    RealT tol    = 1.e-8;

    RealT CGtol1 = 1.e-4;
    RealT CGtol2 = 1.e-2;
    int maxitCG  = 200;

    ROL::LineSearchStep<RealT> LS_step(LStype,LScond,LSStype,maxit,c1,c2,tol,rho,
                                       Stype,L,BBtype,CGtol1,CGtol2,maxitCG);
    /* END LINE SEARCH STEP DEFINITION */

    /* BEGIN TRUST REGION STEP DEFINTION */
    //ROL::TrustRegionStepType TRStype = ROL::TrustRegionStep_Newton;
    ROL::TrustRegionStepType TRStype = ROL::TrustRegionStep_NewtonKrylov;
    //ROL::TrustRegionStepType TRStype = ROL::TrustRegionStep_NewtonKrylovSecantPreconditioning;
    //ROL::TrustRegionStepType TRStype = ROL::TrustRegionStep_Secant;
    //ROL::TrustRegionStepType TRStype = ROL::TrustRegionStep_Gradient;

    //ROL::TrustRegionType TRtype = ROL::TrustRegionType_CauchyPoint;
    ROL::TrustRegionType TRtype = ROL::TrustRegionType_TruncatedCG;
    //ROL::TrustRegionType TRtype = ROL::TrustRegionType_DoubleDogleg;  
    //ROL::TrustRegionType TRtype = ROL::TrustRegionType_Dogleg;

    maxit        = 200;
    RealT tol1   = 1.e-4;
    RealT tol2   = 1.e-2;
    RealT del    = 100.0;
    RealT delmin = 1.e-8;
    RealT delmax = 5000.0;
    RealT eta0   = 0.05;
    RealT eta1   = 0.05;
    RealT eta2   = 0.9;
    RealT gamma0 = 0.0625;
    RealT gamma1 = 0.25;
    RealT gamma2 = 2.50;
    RealT TRsafe = 1.0;

    ROL::TrustRegionStep<RealT> TR_step(TRtype,TRStype,maxit,tol1,tol2,del,delmin,delmax,
                                        eta0,eta1,eta2,gamma0,gamma1,gamma2,TRsafe,
                                        Stype,L,BBtype);
    /* END TRUST REGION STEP DEFINITION */ 

    ROL::StatusTest<RealT> status(1.e-8,1.e-16,100);    

    x.set(x0);
    ROL::DefaultAlgorithm<RealT> LS_algo(LS_step,status);
    LS_algo.run(x, *obj);
    e.set(x);
    e.axpy(-1.0,z);
    *outStream << "\nNorm of Error: " << e.norm() << "\n";

    x.set(x0);
    ROL::DefaultAlgorithm<RealT> TR_algo(TR_step,status);
    TR_algo.run(x, *obj);
    e.set(x);
    e.axpy(-1.0,z);
    *outStream << "\nNorm of Error: " << e.norm() << "\n";

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

