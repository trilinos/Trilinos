// @HEADER
// ************************************************************************
//
//               Rapid Optimization Library (ROL) Package
//                 Copyright (2014) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact lead developers:
//              Drew Kouri   (dpkouri@sandia.gov) and
//              Denis Ridzal (dridzal@sandia.gov)
//
// ************************************************************************
// @HEADER

#include "test_01.hpp"

typedef double RealT;

template<class Real>
Real random(const Teuchos::RCP<const Teuchos::Comm<int> > &comm) {
  Real val = 0.0;
  if ( Teuchos::rank<int>(*comm)==0 ) {
    val = (Real)rand()/(Real)RAND_MAX;
  }
  Teuchos::broadcast<int,Real>(*comm,0,1,&val);
  return val;
}

int main(int argc, char* argv[]) {

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  Teuchos::RCP<const Teuchos::Comm<int> > comm
    = Teuchos::DefaultComm<int>::getComm();

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint = argc - 1;
  Teuchos::RCP<std::ostream> outStream;
  Teuchos::oblackholestream bhs; // outputs nothing
  if (iprint > 0 && Teuchos::rank<int>(*comm)==0)
    outStream = Teuchos::rcp(&std::cout, false);
  else
    outStream = Teuchos::rcp(&bhs, false);

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
    Teuchos::RCP<std::vector<RealT> > x1_rcp  = Teuchos::rcp( new std::vector<RealT>(nx+2,0.0) );
    //ROL::StdVector<RealT> x1(x1_rcp);
      Teuchos::RCP<ROL::StdVector<RealT> > x1 = Teuchos::rcp(new ROL::StdVector<RealT>(x1_rcp));
    Teuchos::RCP<std::vector<RealT> > x2_rcp  = Teuchos::rcp( new std::vector<RealT>(nx+2,0.0) );
    ROL::StdVector<RealT> x2(x2_rcp);
    Teuchos::RCP<std::vector<RealT> > x3_rcp  = Teuchos::rcp( new std::vector<RealT>(nx+2,0.0) );
    ROL::StdVector<RealT> x3(x3_rcp);
    Teuchos::RCP<std::vector<RealT> > z_rcp  = Teuchos::rcp( new std::vector<RealT>(nx+2,0.0) );
    //ROL::StdVector<RealT> z(z_rcp);
      Teuchos::RCP<ROL::StdVector<RealT> > z = Teuchos::rcp(new ROL::StdVector<RealT>(z_rcp));
    Teuchos::RCP<std::vector<RealT> > xr_rcp = Teuchos::rcp( new std::vector<RealT>(nx+2,0.0) );
    ROL::StdVector<RealT> xr(xr_rcp);
    Teuchos::RCP<std::vector<RealT> > d_rcp  = Teuchos::rcp( new std::vector<RealT>(nx+2,0.0) );
    //ROL::StdVector<RealT> d(d_rcp);
      Teuchos::RCP<ROL::StdVector<RealT> > d = Teuchos::rcp(new ROL::StdVector<RealT>(d_rcp));
    for ( int i = 0; i < nx+2; i++ ) {
      (*xr_rcp)[i] = random<RealT>(comm);
      (*d_rcp)[i]  = random<RealT>(comm);
      (*z_rcp)[i]  = initZ;
    }
    Teuchos::RCP<Teuchos::ParameterList> cvarlist = Teuchos::rcp(new Teuchos::ParameterList);
    cvarlist->sublist("SOL").sublist("Risk Measure").set("Name", "CVaR");
    ROL::RiskVector<RealT> zR(cvarlist,z), x1R(cvarlist,x1), dR(cvarlist,d);
    // Build state and adjoint vectors
    Teuchos::RCP<std::vector<RealT> > u_rcp  = Teuchos::rcp( new std::vector<RealT>(nx,1.0) );
    ROL::StdVector<RealT> u(u_rcp);
    Teuchos::RCP<std::vector<RealT> > p_rcp  = Teuchos::rcp( new std::vector<RealT>(nx,0.0) );
    ROL::StdVector<RealT> p(p_rcp);
    Teuchos::RCP<ROL::Vector<RealT> > up = Teuchos::rcp(&u,false);
    Teuchos::RCP<ROL::Vector<RealT> > pp = Teuchos::rcp(&p,false);
    /**********************************************************************************************/
    /************************* CONSTRUCT SOL COMPONENTS *******************************************/
    /**********************************************************************************************/
    // Build samplers
    int dim = 4;
    int nSamp = parlist->sublist("Problem Description").get("Number of Samples", 20);
    std::vector<RealT> tmp(2,0.0); tmp[0] = -1.0; tmp[1] = 1.0;
    std::vector<std::vector<RealT> > bounds(dim,tmp);
    Teuchos::RCP<ROL::BatchManager<RealT> > bman
      = Teuchos::rcp(new ROL::StdTeuchosBatchManager<RealT,int>(comm));
    Teuchos::RCP<ROL::SampleGenerator<RealT> > sampler
      = Teuchos::rcp(new ROL::MonteCarloGenerator<RealT>(nSamp,bounds,bman,false,false,100));
    /**********************************************************************************************/
    /************************* CONSTRUCT OBJECTIVE FUNCTION ***************************************/
    /**********************************************************************************************/
    // Build risk-averse objective function
    RealT alpha = 1.e-3;
    Teuchos::RCP<ROL::Objective_SimOpt<RealT> > pobjSimOpt
      = Teuchos::rcp(new Objective_BurgersControl<RealT>(alpha,nx));
    Teuchos::RCP<ROL::Constraint_SimOpt<RealT> > pconSimOpt
      = Teuchos::rcp(new Constraint_BurgersControl<RealT>(nx));
    Teuchos::RCP<ROL::Objective<RealT> > pObj
      = Teuchos::rcp(new ROL::Reduced_Objective_SimOpt<RealT>(pobjSimOpt,pconSimOpt,up,z,pp));
    //Teuchos::RCP<ROL::Objective<RealT> > obj = Teuchos::rcp(new ROL::RiskNeutralObjective<RealT>(pObj, sampler, true));
    Teuchos::RCP<ROL::Distribution<RealT> > dist  = Teuchos::rcp(new ROL::Parabolic<RealT>(-0.5, 0.5));
    Teuchos::RCP<ROL::PlusFunction<RealT> > pfunc = Teuchos::rcp(new ROL::PlusFunction<RealT>(dist, pfuncSmoothing));
    Teuchos::RCP<ROL::RiskMeasure<RealT> >  rmeas = Teuchos::rcp(new ROL::CVaR<RealT>(cvarLevel, 1.0, pfunc));
    Teuchos::RCP<ROL::Objective<RealT> >    obj   = Teuchos::rcp(new ROL::RiskAverseObjective<RealT>(pObj, rmeas, sampler));
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
    ROL::Algorithm<RealT> algors("Trust Region", *parlist);
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
    std::vector<Teuchos::RCP<ROL::Vector<RealT> > > xu_rcp;
    std::vector<Teuchos::RCP<ROL::Vector<RealT> > > vu_rcp;
    int nvecloc = sampler->numMySamples();
    RealT right = 1, left = 0;
    for( int k=0; k<nvecloc; ++k ) {
      Teuchos::RCP<std::vector<RealT> > xuk_rcp = Teuchos::rcp( new std::vector<RealT>(nx,1.0) );
      Teuchos::RCP<std::vector<RealT> > vuk_rcp = Teuchos::rcp( new std::vector<RealT>(nx,1.0) );
      Teuchos::RCP<ROL::Vector<RealT> > xuk = Teuchos::rcp( new ROL::StdVector<RealT>( xuk_rcp ) );
      Teuchos::RCP<ROL::Vector<RealT> > vuk = Teuchos::rcp( new ROL::StdVector<RealT>( vuk_rcp ) );
      for( int i=0; i<nx; ++i ) {
        (*xuk_rcp)[i] = ( (RealT)rand() / (RealT)RAND_MAX ) * (right - left) + left;
        (*vuk_rcp)[i] = ( (RealT)rand() / (RealT)RAND_MAX ) * (right - left) + left;
      }
      xu_rcp.push_back(xuk);
      vu_rcp.push_back(vuk);
    }
    Teuchos::RCP<ROL::SimulatedVector<RealT> > xu = Teuchos::rcp(new ROL::SimulatedVector<RealT>(xu_rcp, bman));
    Teuchos::RCP<ROL::SimulatedVector<RealT> > vu = Teuchos::rcp(new ROL::SimulatedVector<RealT>(vu_rcp, bman));
    // SimOpt vectors.
    Teuchos::RCP<std::vector<RealT> > zvec_rcp = Teuchos::rcp(new std::vector<RealT>(nx+2,0.0));
    Teuchos::RCP<ROL::StdVector<RealT> > zvec = Teuchos::rcp(new ROL::StdVector<RealT>(zvec_rcp));
    Teuchos::RCP<std::vector<RealT> > dvec_rcp = Teuchos::rcp(new std::vector<RealT>(nx+2,0.0));
    Teuchos::RCP<ROL::StdVector<RealT> > dvec = Teuchos::rcp(new ROL::StdVector<RealT>(dvec_rcp));
    for ( int i = 0; i < nx+2; i++ ) {
      (*zvec_rcp)[i] = random<RealT>(comm);
      (*dvec_rcp)[i] = random<RealT>(comm);
    }
    Teuchos::RCP<ROL::RiskVector<RealT> > rz = Teuchos::rcp(new ROL::RiskVector<RealT>(cvarlist, zvec));
    Teuchos::RCP<ROL::RiskVector<RealT> > rd = Teuchos::rcp(new ROL::RiskVector<RealT>(cvarlist, dvec));
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

    ROL::Algorithm<RealT> algo("Composite Step", *parlist);
    ROL::Algorithm<RealT> algo2("Composite Step", *parlist);
    ROL::Algorithm<RealT> algo3("Composite Step", *parlist);
    ROL::Algorithm<RealT> algo4("Composite Step", *parlist);
    ROL::Algorithm<RealT> algo5("Composite Step", *parlist);
    vu->zero();
    for ( int i = 0; i < nx+2; i++ ) {
      (*zvec_rcp)[i] = initZ;
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
        file << (*zvec_rcp)[i] << "\n";
      }
      file.close();
    }

    ROL::RiskVector<RealT> &rxfz = Teuchos::dyn_cast<ROL::RiskVector<RealT> >(*(x.get_2()));
    Teuchos::RCP<ROL::Vector<RealT> > rfz = rxfz.getVector();
    ROL::StdVector<RealT> &rfz_std = Teuchos::dyn_cast<ROL::StdVector<RealT> >(*rfz);
    z->set(rfz_std);
    ROL::Algorithm<RealT> algors2("Trust Region", *parlist);
    algors2.run(zR, *obj, true, *outStream);

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
