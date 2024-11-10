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
    Teuchos::RCP<Teuchos::ParameterList> parlist = Teuchos::rcp( new Teuchos::ParameterList() );
    Teuchos::updateParametersFromXmlFile( filename, parlist.ptr() );
    RealT initZ = parlist->sublist("Problem Description").get("Initial Control Guess", 0.0);
    RealT cvarLevel = parlist->sublist("Problem Description").get("CVaR Level", 0.8);
    RealT pfuncSmoothing = parlist->sublist("Problem Description").get("Plus Function Smoothing Parameter", 1e-2);
    /**********************************************************************************************/
    /************************* CONSTRUCT VECTORS **************************************************/
    /**********************************************************************************************/
    // Build control vectors
    int nx = 256;
    ROL::Ptr<std::vector<RealT> > x1_ptr  = ROL::makePtr<std::vector<RealT>>(nx+2,0.0);
    //ROL::StdVector<RealT> x1(x1_ptr);
      ROL::Ptr<ROL::StdVector<RealT> > x1 = ROL::makePtr<ROL::StdVector<RealT>>(x1_ptr);
    ROL::Ptr<std::vector<RealT> > x2_ptr  = ROL::makePtr<std::vector<RealT>>(nx+2,0.0);
    ROL::StdVector<RealT> x2(x2_ptr);
    ROL::Ptr<std::vector<RealT> > x3_ptr  = ROL::makePtr<std::vector<RealT>>(nx+2,0.0);
    ROL::StdVector<RealT> x3(x3_ptr);
    ROL::Ptr<std::vector<RealT> > z_ptr  = ROL::makePtr<std::vector<RealT>>(nx+2,0.0);
    //ROL::StdVector<RealT> z(z_ptr);
      ROL::Ptr<ROL::StdVector<RealT> > z = ROL::makePtr<ROL::StdVector<RealT>>(z_ptr);
    ROL::Ptr<std::vector<RealT> > xr_ptr = ROL::makePtr<std::vector<RealT>>(nx+2,0.0);
    ROL::StdVector<RealT> xr(xr_ptr);
    ROL::Ptr<std::vector<RealT> > d_ptr  = ROL::makePtr<std::vector<RealT>>(nx+2,0.0);
    //ROL::StdVector<RealT> d(d_ptr);
      ROL::Ptr<ROL::StdVector<RealT> > d = ROL::makePtr<ROL::StdVector<RealT>>(d_ptr);
    for ( int i = 0; i < nx+2; i++ ) {
      (*xr_ptr)[i] = random<RealT>(comm);
      (*d_ptr)[i]  = random<RealT>(comm);
      (*z_ptr)[i]  = initZ;
    }
    Teuchos::RCP<Teuchos::ParameterList> cvarlist = Teuchos::rcp( new Teuchos::ParameterList() );
    cvarlist->sublist("SOL").sublist("Risk Measure").set("Name", "CVaR");
    ROL::RiskVector<RealT> zR(cvarlist,z), x1R(cvarlist,x1), dR(cvarlist,d);
    // Build state and adjoint vectors
    ROL::Ptr<std::vector<RealT> > u_ptr  = ROL::makePtr<std::vector<RealT>>(nx,1.0);
    ROL::StdVector<RealT> u(u_ptr);
    ROL::Ptr<std::vector<RealT> > p_ptr  = ROL::makePtr<std::vector<RealT>>(nx,0.0);
    ROL::StdVector<RealT> p(p_ptr);
    ROL::Ptr<ROL::Vector<RealT> > up = ROL::makePtrFromRef(u);
    ROL::Ptr<ROL::Vector<RealT> > pp = ROL::makePtrFromRef(p);
    /**********************************************************************************************/
    /************************* CONSTRUCT SOL COMPONENTS *******************************************/
    /**********************************************************************************************/
    // Build samplers
    int dim = 4;
    int nSamp = parlist->sublist("Problem Description").get("Number of Samples", 20);
    std::vector<RealT> tmp(2,0.0); tmp[0] = -1.0; tmp[1] = 1.0;
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
    ROL::Ptr<ROL::Objective<RealT> > pObj
      = ROL::makePtr<ROL::Reduced_Objective_SimOpt<RealT>>(pobjSimOpt,pconSimOpt,up,z,pp);
    //ROL::Ptr<ROL::Objective<RealT>> obj = ROL::makePtr<ROL::RiskNeutralObjective<RealT>>(pObj, sampler, true);
    ROL::Ptr<ROL::Distribution<RealT>>      dist  = ROL::makePtr<ROL::Parabolic<RealT>>(-0.5, 0.5);
    ROL::Ptr<ROL::PlusFunction<RealT>>      pfunc = ROL::makePtr<ROL::PlusFunction<RealT>>(dist, pfuncSmoothing);
    ROL::Ptr<ROL::RandVarFunctional<RealT>> rmeas = ROL::makePtr<ROL::CVaR<RealT>>(cvarLevel, 1.0, pfunc);
    ROL::Ptr<ROL::Objective<RealT>>         obj   = ROL::makePtr<ROL::StochasticObjective<RealT>>(pObj, rmeas, sampler);
    // Test parametrized objective functions
    *outStream << "Check Derivatives of Parametrized Objective Function\n";
    //x1.set(xr);
      x1->set(xr);
    pObj->setParameter(sampler->getMyPoint(0));
    //pObj->checkGradient(x1,d,true,*outStream);
      pObj->checkGradient(*x1,*d,true,*outStream);
    //pObj->checkHessVec(x1,d,true,*outStream);
      pObj->checkHessVec(*x1,*d,true,*outStream);
    //obj->checkGradient(x1,d,true,*outStream);
      obj->checkGradient(x1R,dR,true,*outStream);
    //obj->checkHessVec(x1,d,true,*outStream);
      obj->checkHessVec(x1R,dR,true,*outStream);
    ROL::Ptr<ROL::Step<RealT>>
      steprs = ROL::makePtr<ROL::TrustRegionStep<RealT>>(*parlist);
    ROL::Ptr<ROL::StatusTest<RealT>>
      statusrs = ROL::makePtr<ROL::StatusTest<RealT>>(*parlist);
    ROL::Algorithm<RealT> algors(steprs,statusrs,false);
    //algors.run(z, *obj, true, *outStream);
    algors.run(zR, *obj, true, *outStream);
    /**********************************************************************************************/
    /****************** CONSTRUCT SIMULATED CONSTRAINT AND VECTORS ********************************/
    /**********************************************************************************************/

    // Construct SimulatedConstraint.
    int useW = parlist->sublist("Problem Description").get("Use Constraint Weights", true);
    ROL::SimulatedConstraint<RealT> simcon(sampler, pconSimOpt, useW);
    // Construct SimulatedObjective.
    ROL::SimulatedObjectiveCVaR<RealT> simobj(sampler, pobjSimOpt, pfunc, cvarLevel);
    // Simulated vectors.
    std::vector<ROL::Ptr<ROL::Vector<RealT> > > xu_ptr;
    std::vector<ROL::Ptr<ROL::Vector<RealT> > > vu_ptr;
    int nvecloc = sampler->numMySamples();
    RealT right = 1, left = 0;
    for( int k=0; k<nvecloc; ++k ) {
      ROL::Ptr<std::vector<RealT> > xuk_ptr = ROL::makePtr<std::vector<RealT>>(nx,1.0);
      ROL::Ptr<std::vector<RealT> > vuk_ptr = ROL::makePtr<std::vector<RealT>>(nx,1.0);
      ROL::Ptr<ROL::Vector<RealT> > xuk = ROL::makePtr<ROL::StdVector<RealT>>( xuk_ptr );
      ROL::Ptr<ROL::Vector<RealT> > vuk = ROL::makePtr<ROL::StdVector<RealT>>( vuk_ptr );
      for( int i=0; i<nx; ++i ) {
        (*xuk_ptr)[i] = ( (RealT)rand() / (RealT)RAND_MAX ) * (right - left) + left;
        (*vuk_ptr)[i] = ( (RealT)rand() / (RealT)RAND_MAX ) * (right - left) + left;
      }
      xu_ptr.push_back(xuk);
      vu_ptr.push_back(vuk);
    }
    ROL::Ptr<ROL::SimulatedVector<RealT> > xu = ROL::makePtr<ROL::SimulatedVector<RealT>>(xu_ptr, bman);
    ROL::Ptr<ROL::SimulatedVector<RealT> > vu = ROL::makePtr<ROL::SimulatedVector<RealT>>(vu_ptr, bman);
    // SimOpt vectors.
    ROL::Ptr<std::vector<RealT> > zvec_ptr = ROL::makePtr<std::vector<RealT>>(nx+2,0.0);
    ROL::Ptr<ROL::StdVector<RealT> > zvec = ROL::makePtr<ROL::StdVector<RealT>>(zvec_ptr);
    ROL::Ptr<std::vector<RealT> > dvec_ptr = ROL::makePtr<std::vector<RealT>>(nx+2,0.0);
    ROL::Ptr<ROL::StdVector<RealT> > dvec = ROL::makePtr<ROL::StdVector<RealT>>(dvec_ptr);
    for ( int i = 0; i < nx+2; i++ ) {
      (*zvec_ptr)[i] = random<RealT>(comm);
      (*dvec_ptr)[i] = random<RealT>(comm);
    }
    ROL::Ptr<ROL::RiskVector<RealT> > rz = ROL::makePtr<ROL::RiskVector<RealT>>(cvarlist, zvec);
    ROL::Ptr<ROL::RiskVector<RealT> > rd = ROL::makePtr<ROL::RiskVector<RealT>>(cvarlist, dvec);
    ROL::Vector_SimOpt<RealT> x(xu, rz);
    ROL::Vector_SimOpt<RealT> v(vu, rd);

    *outStream << std::endl << "TESTING SimulatedConstraint" << std::endl; 
    simcon.checkApplyJacobian(x, v, *vu, true, *outStream);
    simcon.checkAdjointConsistencyJacobian(*vu, v, x, *vu, x, true, *outStream);
    simcon.checkApplyAdjointHessian(x, *vu, v, x, true, *outStream);
    *outStream << std::endl << "TESTING SimulatedObjective" << std::endl;
    RealT tol = 1e-8;
    simobj.value(x, tol);
    simobj.checkGradient(x, v, true, *outStream);
    simobj.checkHessVec(x, v, true, *outStream);

    ROL::Ptr<ROL::Step<RealT>> step1, step2, step3, step4, step5;
    ROL::Ptr<ROL::StatusTest<RealT>> status1, status2, status3, status4, status5;
    step1 = ROL::makePtr<ROL::CompositeStep<RealT>>(*parlist);
    step2 = ROL::makePtr<ROL::CompositeStep<RealT>>(*parlist);
    step3 = ROL::makePtr<ROL::CompositeStep<RealT>>(*parlist);
    step4 = ROL::makePtr<ROL::CompositeStep<RealT>>(*parlist);
    step5 = ROL::makePtr<ROL::CompositeStep<RealT>>(*parlist);
    status1 = ROL::makePtr<ROL::ConstraintStatusTest<RealT>>(*parlist);
    status2 = ROL::makePtr<ROL::ConstraintStatusTest<RealT>>(*parlist);
    status3 = ROL::makePtr<ROL::ConstraintStatusTest<RealT>>(*parlist);
    status4 = ROL::makePtr<ROL::ConstraintStatusTest<RealT>>(*parlist);
    status5 = ROL::makePtr<ROL::ConstraintStatusTest<RealT>>(*parlist);
    ROL::Algorithm<RealT> algo(step1,status1,false);
    ROL::Algorithm<RealT> algo2(step2,status2,false);
    ROL::Algorithm<RealT> algo3(step3,status3,false);
    ROL::Algorithm<RealT> algo4(step4,status4,false);
    ROL::Algorithm<RealT> algo5(step5,status5,false);
    vu->zero();
    for ( int i = 0; i < nx+2; i++ ) {
      (*zvec_ptr)[i] = initZ;
    }
    ROL::SimulatedObjectiveCVaR<RealT> simobjExpval(sampler, pobjSimOpt, pfunc, 0.0);
    ROL::SimulatedObjectiveCVaR<RealT> simobjCVaR3(sampler, pobjSimOpt, pfunc, 0.3);
    ROL::SimulatedObjectiveCVaR<RealT> simobjCVaR6(sampler, pobjSimOpt, pfunc, 0.6);
    ROL::SimulatedObjectiveCVaR<RealT> simobjCVaR7(sampler, pobjSimOpt, pfunc, 0.7);
    algo2.run(x, *vu, simobjExpval, simcon, true, *outStream);
    algo3.run(x, *vu, simobjCVaR3, simcon, true, *outStream);
    algo4.run(x, *vu, simobjCVaR6, simcon, true, *outStream);
    algo5.run(x, *vu, simobjCVaR7, simcon, true, *outStream);
    algo.run(x, *vu, simobj, simcon, true, *outStream);

    // Output control to file.
    if (Teuchos::rank<int>(*comm)==0) {
      std::ofstream file;
      file.open("control-fs-cvar.txt");
      for ( int i = 0; i < nx+2; ++i ) {
        file << (*zvec_ptr)[i] << "\n";
      }
      file.close();
    }

    ROL::RiskVector<RealT> &rxfz = dynamic_cast<ROL::RiskVector<RealT>&>(*(x.get_2()));
    ROL::Ptr<ROL::Vector<RealT> > rfz = rxfz.getVector();
    ROL::StdVector<RealT> &rfz_std = dynamic_cast<ROL::StdVector<RealT>&>(*rfz);
    z->set(rfz_std);
    ROL::Ptr<ROL::Step<RealT>>
      steprs2 = ROL::makePtr<ROL::TrustRegionStep<RealT>>(*parlist);
    ROL::Ptr<ROL::StatusTest<RealT>>
      statusrs2 = ROL::makePtr<ROL::StatusTest<RealT>>(*parlist);
    ROL::Algorithm<RealT> algors2(steprs2,statusrs2,false);
    algors2.run(zR, *obj, true, *outStream);

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
