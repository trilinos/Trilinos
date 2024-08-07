// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "example_05.hpp"

typedef double RealT;

template<class Real>
Real random(const ROL::Ptr<const Teuchos::Comm<int> > &comm) {
  Real val = 0.0;
  if ( Teuchos::rank<int>(*comm)==0 ) {
    val = (Real)rand()/(Real)RAND_MAX;
  }
  Teuchos::broadcast<int,Real>(*comm,0,1,&val);
  return val;
}

int main(int argc, char* argv[]) {

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  ROL::Ptr<const Teuchos::Comm<int> > comm
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
    /**********************************************************************************************/
    /************************* CONSTRUCT ROL ALGORITHM ********************************************/
    /**********************************************************************************************/
    // Get ROL parameterlist
    std::string filename = "input.xml";
    auto parlist = ROL::getParametersFromXmlFile( filename );
    // Build ROL algorithm
    parlist->sublist("Status Test").set("Gradient Tolerance",1.e-7);
    parlist->sublist("Status Test").set("Step Tolerance",1.e-14);
    parlist->sublist("Status Test").set("Iteration Limit",100);
    /**********************************************************************************************/
    /************************* CONSTRUCT VECTORS **************************************************/
    /**********************************************************************************************/
    // Build control vectors
    int nx = 256;
    // Construct storage for optimal solution
    ROL::Ptr<std::vector<RealT> > z_ptr  = ROL::makePtr<std::vector<RealT>>(nx+2,0);
    ROL::Ptr<ROL::Vector<RealT> > zp  = ROL::makePtr<ROL::StdVector<RealT>>(z_ptr);
    ROL::Ptr<std::vector<RealT> > x1_ptr = ROL::makePtr<std::vector<RealT>>(nx+2,0);
    ROL::Ptr<ROL::Vector<RealT> > x1p = ROL::makePtr<ROL::StdVector<RealT>>(x1_ptr);
    ROL::Ptr<std::vector<RealT> > x2_ptr = ROL::makePtr<std::vector<RealT>>(nx+2,0);
    ROL::Ptr<ROL::Vector<RealT> > x2p = ROL::makePtr<ROL::StdVector<RealT>>(x2_ptr);
    ROL::Ptr<std::vector<RealT> > x3_ptr = ROL::makePtr<std::vector<RealT>>(nx+2,0);
    ROL::Ptr<ROL::Vector<RealT> > x3p = ROL::makePtr<ROL::StdVector<RealT>>(x3_ptr);
    std::vector<ROL::Ptr<ROL::Vector<RealT> > > xvec = {x1p, x2p, x3p};
    // Create vectors for derivative check
    ROL::Ptr<std::vector<RealT> > xr_ptr = ROL::makePtr<std::vector<RealT>>(nx+2,0);
    ROL::StdVector<RealT> xr(xr_ptr);
    ROL::Ptr<std::vector<RealT> > d_ptr  = ROL::makePtr<std::vector<RealT>>(nx+2,0);
    ROL::StdVector<RealT> d(d_ptr);
    for ( int i = 0; i < nx+2; i++ ) {
      (*xr_ptr)[i] = random<RealT>(comm);
      (*d_ptr)[i]  = random<RealT>(comm);
    }
    // Build state and adjoint vectors
    ROL::Ptr<std::vector<RealT> > u_ptr = ROL::makePtr<std::vector<RealT>>(nx,1);
    ROL::Ptr<ROL::Vector<RealT> > up = ROL::makePtr<ROL::StdVector<RealT>>(u_ptr);
    ROL::Ptr<std::vector<RealT> > p_ptr = ROL::makePtr<std::vector<RealT>>(nx,0);
    ROL::Ptr<ROL::Vector<RealT> > pp = ROL::makePtr<ROL::StdVector<RealT>>(p_ptr);
    /**********************************************************************************************/
    /************************* CONSTRUCT SOL COMPONENTS *******************************************/
    /**********************************************************************************************/
    // Build samplers
    int dim = 4, nSamp = 100;
    std::vector<RealT> tmp = {-1, 1};
    std::vector<std::vector<RealT> > bounds(dim,tmp);
    ROL::Ptr<ROL::BatchManager<RealT> > bman
      = ROL::makePtr<ROL::StdTeuchosBatchManager<RealT,int>>(comm);
    ROL::Ptr<ROL::SampleGenerator<RealT> > sampler
      = ROL::makePtr<ROL::MonteCarloGenerator<RealT>>(nSamp,bounds,bman,false,false,100);
    /**********************************************************************************************/
    /************************* CONSTRUCT OBJECTIVE FUNCTION ***************************************/
    /**********************************************************************************************/
    // Build risk-averse objective function
    RealT alpha = 1.e-3;
    ROL::Ptr<ROL::Objective_SimOpt<RealT> > pobjSimOpt
      = ROL::makePtr<Objective_BurgersControl<RealT>>(alpha,nx);
    ROL::Ptr<ROL::Constraint_SimOpt<RealT> > pconSimOpt
      = ROL::makePtr<Constraint_BurgersControl<RealT>>(nx);
    pconSimOpt->setSolveParameters(*parlist);
    ROL::Ptr<ROL::Objective<RealT> > pObj
      = ROL::makePtr<ROL::Reduced_Objective_SimOpt<RealT>>(pobjSimOpt,pconSimOpt,up,zp,pp);
    // Test parametrized objective functions
    *outStream << "Check Derivatives of Parametrized Objective Function\n";
    xvec[0]->set(xr);
    pObj->setParameter(sampler->getMyPoint(0));
    pObj->checkGradient(*xvec[0],d,true,*outStream);
    pObj->checkHessVec(*xvec[0],d,true,*outStream);
    /**********************************************************************************************/
    /************************* SMOOTHED CVAR 1.e-2, 1.e-4, 1.e-6 **********************************/
    /**********************************************************************************************/
    const RealT cl(0.9), cc(1), lb(-0.5), ub(0.5);
    const std::string ra = "Risk Averse", rm = "CVaR", dist = "Parabolic";
    const bool storage = true;
    RealT eps(1.e-2);
    std::vector<RealT> stat(3,0);
    ROL::Ptr<ROL::OptimizationProblem<RealT>> optProb;
    ROL::Ptr<ROL::OptimizationSolver<RealT>>  solver;
    for (int i = 0; i < 3; ++i) {
      *outStream << "\nSOLVE SMOOTHED CONDITIONAL VALUE AT RISK WITH TRUST REGION\n";
      // Build CVaR risk measure
      ROL::ParameterList list;
      list.sublist("SOL").set("Type",ra);
      list.sublist("SOL").set("Store Sampled Value and Gradient",storage);
      list.sublist("SOL").sublist("Risk Measure").set("Name",rm);
      list.sublist("SOL").sublist("Risk Measure").sublist(rm).set("Confidence Level",cl);
      list.sublist("SOL").sublist("Risk Measure").sublist(rm).set("Convex Combination Parameter",cc);
      list.sublist("SOL").sublist("Risk Measure").sublist(rm).set("Smoothing Parameter",eps);
      list.sublist("SOL").sublist("Risk Measure").sublist(rm).sublist("Distribution").set("Name",dist);
      list.sublist("SOL").sublist("Risk Measure").sublist(rm).sublist("Distribution").sublist(dist).set("Lower Bound",lb);
      list.sublist("SOL").sublist("Risk Measure").sublist(rm).sublist("Distribution").sublist(dist).set("Upper Bound",ub);
      // Build stochastic problem
      if ( i==0 ) { xvec[i]->zero();          }
      else        { xvec[i]->set(*xvec[i-1]); }
      optProb = ROL::makePtr<ROL::OptimizationProblem<RealT>>(pObj,xvec[i]);
      RealT init_stat(1);
      if ( i > 0 ) { init_stat = stat[i-1]; }
      list.sublist("SOL").set("Initial Statistic",init_stat);
      optProb->setStochasticObjective(list,sampler);
      optProb->check(*outStream);
      // Run ROL algorithm
      parlist->sublist("Step").set("Type","Trust Region");
      solver = ROL::makePtr<ROL::OptimizationSolver<RealT>>(*optProb,*parlist);
      clock_t start = clock();
      solver->solve(*outStream);
      *outStream << "Optimization time: " << (RealT)(clock()-start)/(RealT)CLOCKS_PER_SEC << " seconds.\n";
      // Get solution statistic
      stat[i] = optProb->getSolutionStatistic();
      // Update smoothing parameter
      eps *= static_cast<RealT>(1.e-2);
    }
    /**********************************************************************************************/
    /************************* NONSMOOTH PROBLEM **************************************************/
    /**********************************************************************************************/
    *outStream << "\nSOLVE NONSMOOTH CVAR PROBLEM WITH BUNDLE TRUST REGION\n";
    ROL::ParameterList list;
    list.sublist("SOL").set("Type",ra);
    list.sublist("SOL").set("Store Sampled Value and Gradient",storage);
    list.sublist("SOL").sublist("Risk Measure").set("Name",rm);
    list.sublist("SOL").sublist("Risk Measure").sublist(rm).set("Confidence Level",cl);
    list.sublist("SOL").sublist("Risk Measure").sublist(rm).set("Convex Combination Parameter",cc);
    list.sublist("SOL").sublist("Risk Measure").sublist(rm).set("Smoothing Parameter",0.);
    list.sublist("SOL").sublist("Risk Measure").sublist(rm).sublist("Distribution").set("Name","Dirac");
    list.sublist("SOL").sublist("Risk Measure").sublist(rm).sublist("Distribution").sublist("Dirac").set("Location",0.);
    // Build stochastic problem
    zp->set(*xvec[2]);
    optProb = ROL::makePtr<ROL::OptimizationProblem<RealT>>(pObj,zp);
    list.sublist("SOL").set("Initial Statistic",stat[2]);
    optProb->setStochasticObjective(list,sampler);
    optProb->check(*outStream);
    // Run ROL algorithm
    parlist->sublist("Status Test").set("Iteration Limit",1000);
    parlist->sublist("Step").sublist("Bundle").set("Epsilon Solution Tolerance",1.e-7);
    parlist->sublist("Step").set("Type","Bundle");
    solver = ROL::makePtr<ROL::OptimizationSolver<RealT>>(*optProb,*parlist);
    clock_t start = clock();
    solver->solve(*outStream);
    *outStream << "Optimization time: " << (RealT)(clock()-start)/(RealT)CLOCKS_PER_SEC << " seconds.\n";
    /**********************************************************************************************/
    /************************* COMPUTE ERROR ******************************************************/
    /**********************************************************************************************/
    ROL::Ptr<ROL::Vector<RealT> > cErr = zp->clone();
    RealT zstat = optProb->getSolutionStatistic();
    *outStream << "\nSUMMARY:\n";
    *outStream << "  ---------------------------------------------\n";
    *outStream << "    True Value-At-Risk    = " << zstat << "\n";
    *outStream << "  ---------------------------------------------\n";
    RealT VARerror  = std::abs(zstat-stat[0]);
    cErr->set(*xvec[0]); cErr->axpy(-1.0,*zp);
    RealT CTRLerror = cErr->norm();
    RealT TOTerror1 = std::sqrt(std::pow(VARerror,2)+std::pow(CTRLerror,2));
    *outStream << "    Value-At-Risk (1.e-2) = " <<   stat[0] << "\n";
    *outStream << "    Value-At-Risk Error   = " <<  VARerror << "\n";
    *outStream << "    Control Error         = " << CTRLerror << "\n";
    *outStream << "    Total Error           = " << TOTerror1 << "\n";
    *outStream << "  ---------------------------------------------\n";
    VARerror  = std::abs(zstat-stat[1]);
    cErr->set(*xvec[1]); cErr->axpy(-1.0,*zp);
    CTRLerror = cErr->norm();
    RealT TOTerror2 = std::sqrt(std::pow(VARerror,2)+std::pow(CTRLerror,2));
    *outStream << "    Value-At-Risk (1.e-4) = " <<   stat[1] << "\n";
    *outStream << "    Value-At-Risk Error   = " <<  VARerror << "\n";
    *outStream << "    Control Error         = " << CTRLerror << "\n";
    *outStream << "    Total Error           = " << TOTerror2 << "\n";
    *outStream << "  ---------------------------------------------\n";
    VARerror  = std::abs(zstat-stat[2]);
    cErr->set(*xvec[2]); cErr->axpy(-1.0,*zp);
    CTRLerror = cErr->norm();
    RealT TOTerror3 = std::sqrt(std::pow(VARerror,2)+std::pow(CTRLerror,2));
    *outStream << "    Value-At-Risk (1.e-6) = " <<   stat[2] << "\n";
    *outStream << "    Value-At-Risk Error   = " <<  VARerror << "\n";
    *outStream << "    Control Error         = " << CTRLerror << "\n";
    *outStream << "    Total Error           = " << TOTerror3 << "\n";
    *outStream << "  ---------------------------------------------\n\n";
    // Comparison
    errorFlag += ((TOTerror1 < 90.*TOTerror2) && (TOTerror2 < 90.*TOTerror3)) ? 1 : 0;

    // Output controls
    std::ofstream control;
    control.open("example04_control.txt");
    for (int n = 0; n < nx+2; n++) {
      control << std::scientific << std::setprecision(15)
              << std::setw(25) << static_cast<RealT>(n)/static_cast<RealT>(nx+1)
              << std::setw(25) << (*z_ptr)[n]
              << std::endl;
    }
    control.close();

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
