// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "test_01.hpp"

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
    = Teuchos::DefaultComm<int>::getComm();

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint = argc - 1;
  ROL::Ptr<std::ostream> outStream;
  ROL::nullstream bhs; // outputs nothing
  if (iprint > 0 && Teuchos::rank<int>(*comm)==0)
    ROL::makePtrFromRef(std::cout);
  else
    ROL::makePtrFromRef(bhs);

  int errorFlag  = 0;

  try {
    /**********************************************************************************************/
    /************************* CONSTRUCT ROL ALGORITHM ********************************************/
    /**********************************************************************************************/
    // Get ROL parameterlist
    std::string filename = "input.xml";
    Teuchos::RCP<Teuchos::ParameterList> parlist = Teuchos::rcp( new Teuchos::ParameterList() );
    Teuchos::updateParametersFromXmlFile( filename, parlist.ptr() );
    // Build ROL algorithm
    parlist->sublist("Status Test").set("Gradient Tolerance",1.e-8);
    parlist->sublist("Status Test").set("Step Tolerance",1.e-14);
    parlist->sublist("Status Test").set("Iteration Limit",100);
    /**********************************************************************************************/
    /************************* CONSTRUCT VECTORS **************************************************/
    /**********************************************************************************************/
    // Build control vectors
    int nx = 256;
    ROL::Ptr<std::vector<RealT> > x_ptr  = ROL::makePtr<std::vector<RealT>>(nx+2,0.0);
    ROL::StdVector<RealT> x(x_ptr);
    ROL::Ptr<std::vector<RealT> > d_ptr  = ROL::makePtr<std::vector<RealT>>(nx+2,0.0);
    ROL::StdVector<RealT> d(d_ptr);
    for ( int i = 0; i < nx+2; i++ ) {
      (*x_ptr)[i] = random<RealT>(comm);
      (*d_ptr)[i] = random<RealT>(comm);
    }
    ROL::Ptr<ROL::Vector<RealT> > xp = &x,false;
    // Build state and adjoint vectors
    ROL::Ptr<std::vector<RealT> > u_ptr  = ROL::makePtr<std::vector<RealT>>(nx,1.0);
    ROL::StdVector<RealT> u(u_ptr);
    ROL::Ptr<std::vector<RealT> > p_ptr  = ROL::makePtr<std::vector<RealT>>(nx,0.0);
    ROL::StdVector<RealT> p(p_ptr);
    ROL::Ptr<ROL::Vector<RealT> > up = &u,false;
    ROL::Ptr<ROL::Vector<RealT> > pp = &p,false;
    /**********************************************************************************************/
    /************************* CONSTRUCT SOL COMPONENTS *******************************************/
    /**********************************************************************************************/
    // Build samplers
//    int dim = 4;
//    int nSamp = 100;
//    std::vector<RealT> tmp(2,0.0); tmp[0] = -1.0; tmp[1] = 1.0;
//    std::vector<std::vector<RealT> > bounds(dim,tmp);
//    ROL::Ptr<ROL::BatchManager<RealT> > bman
//      = ROL::makePtr<ROL::StdTeuchosBatchManager<RealT,int>>(comm);
//    ROL::Ptr<ROL::SampleGenerator<RealT> > sampler
//      = ROL::makePtr<ROL::MonteCarloGenerator<RealT>>(nSamp,bounds,bman,false,false,100);
    ROL::QuadratureInfo info;
    info.dim = 4;
    info.maxLevel = 7;
    info.rule1D.clear(); info.rule1D.resize(info.dim,ROL::QUAD_CLENSHAWCURTIS);    
    info.growth1D.clear(); info.growth1D.resize(info.dim,ROL::GROWTH_DEFAULT);
    info.normalized = true;
    info.adaptive = false;
    ROL::Ptr<ROL::BatchManager<RealT> > bman
      = ROL::makePtr<ROL::StdTeuchosBatchManager<RealT,int>>(comm);
    ROL::Ptr<ROL::SampleGenerator<RealT> > sampler
      = ROL::makePtr<ROL::SparseGridGenerator<RealT>>(bman,info);
    // Print samples
    int width = 21;
    std::stringstream name;
    name << "samples_" << sampler->batchID() << ".txt";
    std::ofstream file(name.str().c_str());
    if (file.is_open()) {
      file << std::scientific << std::setprecision(12);
      for (int i = 0; i < sampler->numMySamples(); ++i) {
        std::vector<RealT> pt = sampler->getMyPoint(i);
        RealT wt = sampler->getMyWeight(i);
        for (int j = 0; j < static_cast<int>(pt.size()); ++j) {
          file << std::setw(width) << std::left << pt[j];
        }
        file << std::setw(width) << std::left << wt << std::endl;
      }
      file.close();
    }
    else {
      ROL_TEST_FOR_EXCEPTION(true, std::invalid_argument,
        ">>> (adapters/trikota/sol/test/test_01): Unable to open file!");
    }

    /**********************************************************************************************/
    /************************* CONSTRUCT OBJECTIVE FUNCTION ***************************************/
    /**********************************************************************************************/
    // Build risk-averse objective function
    RealT alpha = 1.e-3;
    ROL::Ptr<ROL::Objective_SimOpt<RealT> > pobjSimOpt
      = ROL::makePtr<Objective_BurgersControl<RealT>>(alpha,nx);
    ROL::Ptr<ROL::Constraint_SimOpt<RealT> > pconSimOpt
      = ROL::makePtr<Constraint_BurgersControl<RealT>>(nx);
    ROL::Ptr<ROL::Objective<RealT> > pObj
      = ROL::makePtr<ROL::Reduced_Objective_SimOpt<RealT>>(pobjSimOpt,pconSimOpt,up,xp,pp);
    ROL::Ptr<ROL::Objective<RealT> > obj;
    // Test parametrized objective functions
    *outStream << "Check Derivatives of Parametrized Objective Function\n";
    pObj->setParameter(sampler->getMyPoint(0));
    pObj->checkGradient(x,d,true,*outStream);
    pObj->checkHessVec(x,d,true,*outStream);
    /**********************************************************************************************/
    /************************* RISK NEUTRAL *******************************************************/
    /**********************************************************************************************/
    *outStream << "\nSOLVE RISK NEUTRAL OPTIMAL CONTROL PROBLEM WITH TRUST REGION\n";
    // Build CVaR objective function
    Teuchos::ParameterList list;
    list.sublist("SOL").set("Type","Risk Neutral");
    list.sublist("SOL").set("Store Sampled Value and Gradient",true);
    // Build stochastic problem
    ROL::OptimizationProblem<RealT> optProb(pObj,xp);
    optProb.setStochasticObjective(list,sampler);
    optProb.check(*outStream);
    // Run ROL algorithm
    parlist->sublist("Step").set("Type","Trust Region");
    ROL::OptimizationSolver<RealT> solver(optProb,*parlist);
    clock_t start = clock();
    xp->zero();
    solver.solve(*outStream);
    *outStream << "Optimization time: " << (RealT)(clock()-start)/(RealT)CLOCKS_PER_SEC << " seconds.\n";
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
