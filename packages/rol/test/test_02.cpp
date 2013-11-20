// @HEADER
// ************************************************************************
// @HEADER

/*! \file  test_02.cpp
    \brief Test Rosenbrock.
*/

#define USE_HESSVEC 0

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

    // Define Descent Type
    //ROL::EDescent edesc = ROL::DESCENT_STEEPEST;
    //ROL::EDescent edesc = ROL::DESCENT_SECANT;
    //ROL::EDescent edesc = ROL::DESCENT_NEWTON;
    ROL::EDescent edesc = ROL::DESCENT_NEWTONKRYLOV;
    //ROL::EDescent edesc = ROL::DESCENT_SECANTPRECOND;

    // Define Secant Type
    //ROL::ESecant esec = ROL::SECANT_LBFGS;
    ROL::ESecant esec = ROL::SECANT_LDFP;
    //ROL::ESecant esec = ROL::SECANT_LSR1;
    //ROL::ESecant esec = ROL::SECANT_BARZILAIBORWEIN;
    int L        = 10;
    int BBtype   = 1;

    /* BEGIN LINE SEARCH STEP DEFINTION */
    ROL::ELineSearch els = ROL::LINESEARCH_BACKTRACKING;
    //ROL::ELineSearch els = ROL::LINESEARCH_CUBICINTERP;
    //ROL::ELineSearch els = ROL::LINESEARCH_BISECTION;
    //ROL::ELineSearch els = ROL::LINESEARCH_GOLDENSECTION;
    //ROL::ELineSearch els = ROL::LINESEARCH_BRENTS;

    ROL::ECurvatureCondition econd = ROL::CURVATURECONDITION_WOLFE;
    //ROL::ECurvatureCondition econd = ROL::CURVATURECONDITION_STRONGWOLFE;
    //ROL::ECurvatureCondition econd = ROL::CURVATURECONDITION_GOLDSTEIN;
    /* END LINE SEARCH STEP DEFINITION */

    /* BEGIN TRUST REGION STEP DEFINTION */
    //ROL::ETrustRegion etr = ROL::TRUSTREGION_CAUCHYPOINT;
    ROL::ETrustRegion etr = ROL::TRUSTREGION_TRUNCATEDCG;
    //ROL::ETrustRegion etr = ROL::TRUSTREGION_DOGLEG;
    //ROL::ETrustRegion etr = ROL::TRUSTREGION_DOUBLEDOGLEG;
    /* END TRUST REGION STEP DEFINITION */ 

    // Define Status Test
    ROL::StatusTest<RealT> status(1.e-10,1.e-12,100);    

    for ( ROL::ETestObjectives objFunc = ROL::TESTOBJECTIVES_ROSENBROCK; objFunc < ROL::TESTOBJECTIVES_LAST; objFunc++ ) {
      *outStream << "\n\n" << ROL::ETestObjectivesToString(objFunc) << "\n\n";

      // Initial Guess Vector 
      Teuchos::RCP<std::vector<RealT> > x0_rcp = Teuchos::rcp( new std::vector<RealT> );
      ROL::StdVector<RealT> x0(x0_rcp);

      // Exact Solution Vector
      Teuchos::RCP<std::vector<RealT> > z_rcp = Teuchos::rcp( new std::vector<RealT> );
      ROL::StdVector<RealT> z(z_rcp);

      // Get Objective Function
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
      obj->checkHessVec(x,y,true);

      // RUN LINE SEARCH
      int maxit       = 20;
      RealT rho       = 0.5;
      RealT c1        = 1.e-4;
      RealT c2        = 0.9;
      RealT tol       = 1.e-8;
      RealT CGtol1    = 1.e-4;
      RealT CGtol2    = 1.e-2;
      int maxitCG     = 2*dim;
      bool useInexact = true;

      ROL::LineSearchStep<RealT> LS_step(els,econd,edesc,useInexact,maxit,c1,c2,tol,
                                         rho,esec,L,BBtype,CGtol1,CGtol2,maxitCG);
      x.set(x0);
      ROL::DefaultAlgorithm<RealT> LS_algo(LS_step,status);
      LS_algo.run(x, *obj);
      e.set(x);
      e.axpy(-1.0,z);
      *outStream << "\nNorm of Error: " << e.norm() << "\n";

      // RUN TRUST REGION
      maxit        = 2*dim;
      RealT tol1   = 1.e-4;
      RealT tol2   = 1.e-2;
      RealT del    = -1.0;
      RealT delmin = 1.e-8;
      RealT delmax = 5000.0;
      RealT eta0   = 0.05;
      RealT eta1   = 0.05;
      RealT eta2   = 0.9;
      RealT gamma0 = 0.0625;
      RealT gamma1 = 0.25;
      RealT gamma2 = 2.50;
      RealT TRsafe = 1.0;

      ROL::TrustRegionStep<RealT> TR_step(etr,edesc,maxit,tol1,tol2,del,delmin,delmax,
                                          eta0,eta1,eta2,gamma0,gamma1,gamma2,TRsafe,
                                          esec,L,BBtype);
      x.set(x0);
      ROL::DefaultAlgorithm<RealT> TR_algo(TR_step,status);
      TR_algo.run(x, *obj);
      e.set(x);
      e.axpy(-1.0,z);
      *outStream << "\nNorm of Error: " << e.norm() << "\n";
    }
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

