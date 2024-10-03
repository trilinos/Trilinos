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
    /**********************************************************************************************/
    /************************* CONSTRUCT VECTORS **************************************************/
    /**********************************************************************************************/
    // Build control vectors
    int nx = 1024;
    ROL::Ptr<std::vector<RealT> > x1_ptr  = ROL::makePtr<std::vector<RealT>>(nx+2,0.0);
    ROL::StdVector<RealT> x1(x1_ptr);
    ROL::Ptr<std::vector<RealT> > x2_ptr  = ROL::makePtr<std::vector<RealT>>(nx+2,0.0);
    ROL::StdVector<RealT> x2(x2_ptr);
    ROL::Ptr<std::vector<RealT> > x3_ptr  = ROL::makePtr<std::vector<RealT>>(nx+2,0.0);
    ROL::StdVector<RealT> x3(x3_ptr);
    ROL::Ptr<std::vector<RealT> > z_ptr  = ROL::makePtr<std::vector<RealT>>(nx+2,0.0);
    ROL::StdVector<RealT> z(z_ptr);
    ROL::Ptr<std::vector<RealT> > xr_ptr = ROL::makePtr<std::vector<RealT>>(nx+2,0.0);
    ROL::StdVector<RealT> xr(xr_ptr);
    ROL::Ptr<std::vector<RealT> > d_ptr  = ROL::makePtr<std::vector<RealT>>(nx+2,0.0);
    ROL::StdVector<RealT> d(d_ptr);
    for ( int i = 0; i < nx+2; i++ ) {
      (*xr_ptr)[i] = random<RealT>(comm);
      (*d_ptr)[i]  = random<RealT>(comm);
    }
    // Build state and adjoint vectors
    ROL::Ptr<std::vector<RealT> > u_ptr  = ROL::makePtr<std::vector<RealT>>(nx,1.0);
    ROL::StdVector<RealT> u(u_ptr);
    ROL::Ptr<std::vector<RealT> > p_ptr  = ROL::makePtr<std::vector<RealT>>(nx,0.0);
    ROL::StdVector<RealT> p(p_ptr);
    ROL::Ptr<ROL::Vector<RealT> > up = ROL::makePtrFromRef(u);
    ROL::Ptr<ROL::Vector<RealT> > zp = ROL::makePtrFromRef(z);
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
      = ROL::makePtr<ROL::Reduced_Objective_SimOpt<RealT>>(pobjSimOpt,pconSimOpt,up,zp,pp);
    ROL::Ptr<ROL::Objective<RealT> > obj = ROL::makePtr<ROL::RiskNeutralObjective<RealT>>(pObj, sampler, true);
    // Test parametrized objective functions
    *outStream << "Check Derivatives of Parametrized Objective Function\n";
    x1.set(xr);
    pObj->setParameter(sampler->getMyPoint(0));
    pObj->checkGradient(x1,d,true,*outStream);
    pObj->checkHessVec(x1,d,true,*outStream);
    obj->checkGradient(x1,d,true,*outStream);
    obj->checkHessVec(x1,d,true,*outStream);
    ROL::Ptr<ROL::Step<RealT>>
      steprs = ROL::makePtr<ROL::TrustRegionStep<RealT>>(*parlist);
    ROL::Ptr<ROL::StatusTest<RealT>>
      statusrs = ROL::makePtr<ROL::StatusTest<RealT>>(*parlist);
    ROL::Algorithm<RealT> algors(steprs,statusrs,false);
    algors.run(z, *obj, true, *outStream);
    /**********************************************************************************************/
    /****************** CONSTRUCT SIMULATED CONSTRAINT AND VECTORS ********************************/
    /**********************************************************************************************/

    // Construct SimulatedConstraint.
    int useW = parlist->sublist("Problem Description").get("Use Constraint Weights", true);
    ROL::SimulatedConstraint<RealT> simcon(sampler, pconSimOpt, useW);
    // Construct SimulatedObjective.
    ROL::SimulatedObjective<RealT> simobj(sampler, pobjSimOpt);
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
    ROL::Vector_SimOpt<RealT> x(xu, zvec);
    ROL::Vector_SimOpt<RealT> v(vu, dvec);
    *outStream << std::endl << "TESTING SimulatedConstraint" << std::endl; 
    simcon.checkApplyJacobian(x, v, *vu, true, *outStream);
    simcon.checkAdjointConsistencyJacobian(*vu, v, x, *vu, x, true, *outStream);
    simcon.checkApplyAdjointHessian(x, *vu, v, x, true, *outStream);
    *outStream << std::endl << "TESTING SimulatedObjective" << std::endl;
    simobj.checkGradient(x, v, true, *outStream);
    simobj.checkHessVec(x, v, true, *outStream);
    ROL::Ptr<ROL::Step<RealT>>
      step = ROL::makePtr<ROL::CompositeStep<RealT>>(*parlist);
    ROL::Ptr<ROL::StatusTest<RealT>>
      status = ROL::makePtr<ROL::ConstraintStatusTest<RealT>>(*parlist);
    ROL::Algorithm<RealT> algo(step,status,false);
    vu->zero();
    x.scale(1);
    algo.run(x, *vu, simobj, simcon, true, *outStream);

    // Output control to file.
    if (Teuchos::rank<int>(*comm)==0) {
      std::ofstream file;
      file.open("control-fs-expv.txt");
      for ( int i = 0; i < nx+2; ++i ) {
        file << (*zvec_ptr)[i] << "\n";
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
