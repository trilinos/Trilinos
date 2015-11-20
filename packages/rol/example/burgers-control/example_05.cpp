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

#include "example_05.hpp"

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
    // Build ROL algorithm
    parlist->sublist("Status Test").set("Gradient Tolerance",1.e-8);
    parlist->sublist("Status Test").set("Step Tolerance",1.e-14);
    parlist->sublist("Status Test").set("Iteration Limit",100);
    Teuchos::RCP<ROL::Algorithm<double> > algo;
    /**********************************************************************************************/
    /************************* CONSTRUCT VECTORS **************************************************/
    /**********************************************************************************************/
    // Build control vectors
    int nx = 256;
    Teuchos::RCP<std::vector<double> > x1_rcp  = Teuchos::rcp( new std::vector<double>(nx+2,0.0) );
    ROL::StdVector<double> x1(x1_rcp);
    Teuchos::RCP<std::vector<double> > x2_rcp  = Teuchos::rcp( new std::vector<double>(nx+2,0.0) );
    ROL::StdVector<double> x2(x2_rcp);
    Teuchos::RCP<std::vector<double> > x3_rcp  = Teuchos::rcp( new std::vector<double>(nx+2,0.0) );
    ROL::StdVector<double> x3(x3_rcp);
    Teuchos::RCP<std::vector<double> > z_rcp  = Teuchos::rcp( new std::vector<double>(nx+2,0.0) );
    ROL::StdVector<double> z(z_rcp);
    Teuchos::RCP<std::vector<double> > xr_rcp = Teuchos::rcp( new std::vector<double>(nx+2,0.0) );
    ROL::StdVector<double> xr(xr_rcp);
    Teuchos::RCP<std::vector<double> > d_rcp  = Teuchos::rcp( new std::vector<double>(nx+2,0.0) );
    ROL::StdVector<double> d(d_rcp);
    for ( int i = 0; i < nx+2; i++ ) {
      (*xr_rcp)[i] = random<double>(comm);
      (*d_rcp)[i]  = random<double>(comm);
    }
    // Build state and adjoint vectors
    Teuchos::RCP<std::vector<double> > u_rcp  = Teuchos::rcp( new std::vector<double>(nx,0.0) );
    ROL::StdVector<double> u(u_rcp);
    Teuchos::RCP<std::vector<double> > p_rcp  = Teuchos::rcp( new std::vector<double>(nx,0.0) );
    ROL::StdVector<double> p(p_rcp);
    Teuchos::RCP<ROL::Vector<double> > up = Teuchos::rcp(&u,false);
    Teuchos::RCP<ROL::Vector<double> > pp = Teuchos::rcp(&p,false);
    /**********************************************************************************************/
    /************************* CONSTRUCT SOL COMPONENTS *******************************************/
    /**********************************************************************************************/
    // Build samplers
    int dim = 4;
    int nSamp = 1000;
    std::vector<double> tmp(2,0.0); tmp[0] = -1.0; tmp[1] = 1.0;
    std::vector<std::vector<double> > bounds(dim,tmp);
    Teuchos::RCP<ROL::BatchManager<double> > bman
      = Teuchos::rcp(new ROL::StdTeuchosBatchManager<double,int>(comm));
    Teuchos::RCP<ROL::SampleGenerator<double> > sampler
      = Teuchos::rcp(new ROL::MonteCarloGenerator<double>(nSamp,bounds,bman,false,false,100));
    /**********************************************************************************************/
    /************************* CONSTRUCT OBJECTIVE FUNCTION ***************************************/
    /**********************************************************************************************/
    // Build risk-averse objective function
    double alpha = 1.e-3;
    Teuchos::RCP<ROL::ParametrizedObjective_SimOpt<double> > pobjSimOpt
      = Teuchos::rcp(new Objective_BurgersControl<double>(alpha,nx));
    Teuchos::RCP<ROL::ParametrizedEqualityConstraint_SimOpt<double> > pconSimOpt
      = Teuchos::rcp(new EqualityConstraint_BurgersControl<double>(nx));
    Teuchos::RCP<ROL::ParametrizedObjective<double> > pObj
      = Teuchos::rcp(new ROL::Reduced_ParametrizedObjective_SimOpt<double>(pobjSimOpt,pconSimOpt,up,pp));
    Teuchos::RCP<ROL::Objective<double> > obj;
    // Test parametrized objective functions
    *outStream << "Check Derivatives of Parametrized Objective Function\n";
    x1.set(xr);
    pObj->setParameter(sampler->getMyPoint(0));
    pObj->checkGradient(x1,d,true,*outStream);
    pObj->checkHessVec(x1,d,true,*outStream);
    /**********************************************************************************************/
    /************************* SMOOTHED CVAR 1.e-2 ************************************************/
    /**********************************************************************************************/
    *outStream << "\nSOLVE SMOOTHED CONDITIONAL VALUE AT RISK WITH TRUST REGION\n";
    // Build CVaR objective function
    Teuchos::ParameterList list;
    list.sublist("SOL").set("Store Sampled Value and Gradient",true);
    list.sublist("SOL").sublist("Risk Measure").set("Name","CVaR");
    list.sublist("SOL").sublist("Risk Measure").sublist("CVaR").set("Confidence Level",0.99);
    list.sublist("SOL").sublist("Risk Measure").sublist("CVaR").set("Convex Combination Parameter",1.0);
    list.sublist("SOL").sublist("Risk Measure").sublist("CVaR").set("Smoothing Parameter",1.e-2);
    list.sublist("SOL").sublist("Risk Measure").sublist("CVaR").sublist("Distribution").set("Name","Parabolic");
    list.sublist("SOL").sublist("Risk Measure").sublist("CVaR").sublist("Distribution").sublist("Parabolic").set("Lower Bound",-0.5);
    list.sublist("SOL").sublist("Risk Measure").sublist("CVaR").sublist("Distribution").sublist("Parabolic").set("Upper Bound", 0.5);
    obj = Teuchos::rcp( new ROL::RiskAverseObjective<double>(pObj,list,sampler) );
    // Build CVaR vectors
    double x1v = 10.0*random<double>(comm)-5.0;
    Teuchos::RCP<ROL::Vector<double> > x1p = Teuchos::rcp(&x1,false);
    ROL::RiskVector<double> x1c(x1p,true,x1v);
    // Run ROL algorithm
    algo = Teuchos::rcp(new ROL::Algorithm<double>("Trust Region",*parlist,false));
    x1c.zero();
    clock_t start = clock();
    algo->run(x1c,*obj,true,*outStream);
    *outStream << "Optimization time: " << (double)(clock()-start)/(double)CLOCKS_PER_SEC << " seconds.\n";
    /**********************************************************************************************/
    /************************* SMOOTHED CVAR 1.e-4 ************************************************/
    /**********************************************************************************************/
    *outStream << "\nSOLVE SMOOTHED CONDITIONAL VALUE AT RISK WITH TRUST REGION\n";
    // Build CVaR objective function
    list.sublist("SOL").sublist("Risk Measure").sublist("CVaR").set("Smoothing Parameter",1.e-4);
    obj = Teuchos::rcp( new ROL::RiskAverseObjective<double>(pObj,list,sampler) );
    // Build CVaR vectors
    double x2v = 10.0*random<double>(comm)-5.0;
    Teuchos::RCP<ROL::Vector<double> > x2p = Teuchos::rcp(&x2,false);
    ROL::RiskVector<double> x2c(x2p,true,x2v);
    // Run ROL algorithm
    algo = Teuchos::rcp(new ROL::Algorithm<double>("Trust Region",*parlist,false));
    x2c.set(x1c);
    start = clock();
    algo->run(x2c,*obj,true,*outStream);
    *outStream << "Optimization time: " << (double)(clock()-start)/(double)CLOCKS_PER_SEC << " seconds.\n";
    /**********************************************************************************************/
    /************************* SMOOTHED CVAR 1.e-6 ************************************************/
    /**********************************************************************************************/
    *outStream << "\nSOLVE SMOOTHED CONDITIONAL VALUE AT RISK WITH TRUST REGION\n";
    // Build CVaR objective function
    list.sublist("SOL").sublist("Risk Measure").sublist("CVaR").set("Smoothing Parameter",1.e-6);
    obj = Teuchos::rcp( new ROL::RiskAverseObjective<double>(pObj,list,sampler) );
    // Build CVaR vectors
    double x3v = 10.0*random<double>(comm)-5.0;
    Teuchos::RCP<ROL::Vector<double> > x3p = Teuchos::rcp(&x3,false);
    ROL::RiskVector<double> x3c(x3p,true,x3v);
    // Run ROL algorithm
    algo = Teuchos::rcp(new ROL::Algorithm<double>("Trust Region",*parlist,false));
    x3c.set(x2c);
    start = clock();
    algo->run(x3c,*obj,true,*outStream);
    *outStream << "Optimization time: " << (double)(clock()-start)/(double)CLOCKS_PER_SEC << " seconds.\n";
    /**********************************************************************************************/
    /************************* NONSMOOTH PROBLEM **************************************************/
    /**********************************************************************************************/
    *outStream << "\nSOLVE NONSMOOTH CVAR PROBLEM WITH BUNDLE TRUST REGION\n";
    // Build CVaR objective function
    list.sublist("SOL").sublist("Risk Measure").sublist("CVaR").set("Smoothing Parameter",0.);
    list.sublist("SOL").sublist("Risk Measure").sublist("CVaR").sublist("Distribution").set("Name","Dirac");
    list.sublist("SOL").sublist("Risk Measure").sublist("CVaR").sublist("Distribution").sublist("Dirac").set("Location",0.);
    obj = Teuchos::rcp( new ROL::RiskAverseObjective<double>(pObj,list,sampler) );
    // Build CVaR vector
    double zv = 10.0*random<double>(comm)-5.0;
    Teuchos::RCP<ROL::Vector<double> > zp = Teuchos::rcp(&z,false);
    ROL::RiskVector<double> zc(zp,true,zv);
    // Run ROL algorithm
    parlist->sublist("Status Test").set("Iteration Limit",1000);
    parlist->sublist("Step").sublist("Bundle").set("Epsilon Solution Tolerance",1.e-8);
    algo = Teuchos::rcp(new ROL::Algorithm<double>("Bundle",*parlist,false));
    zc.set(x3c);
    start = clock();
    algo->run(zc,*obj,true,*outStream);
    *outStream << "Optimization time: " << (double)(clock()-start)/(double)CLOCKS_PER_SEC << " seconds.\n";
    /**********************************************************************************************/
    /************************* COMPUTE ERROR ******************************************************/
    /**********************************************************************************************/
    *outStream << "\nSUMMARY:\n";
    *outStream << "  ---------------------------------------------\n";
    *outStream << "    True Value-At-Risk    = " << zc.getStatistic() << "\n";
    *outStream << "  ---------------------------------------------\n";
    double VARerror  = std::abs(zc.getStatistic()-x1c.getStatistic());
    Teuchos::RCP<ROL::Vector<double> > cErr = x1.clone();
    cErr->set(x1); cErr->axpy(-1.0,z);
    double CTRLerror = cErr->norm();
    cErr = x1c.clone();
    cErr->set(x1c); cErr->axpy(-1.0,zc);
    double TOTerror1 = cErr->norm();
    *outStream << "    Value-At-Risk (1.e-2) = " << x1c.getStatistic() << "\n";
    *outStream << "    Value-At-Risk Error   = " <<  VARerror << "\n";
    *outStream << "    Control Error         = " << CTRLerror << "\n";
    *outStream << "    Total Error           = " << TOTerror1 << "\n";
    *outStream << "  ---------------------------------------------\n";
    VARerror  = std::abs(zc.getStatistic()-x2c.getStatistic());
    cErr = x2.clone();
    cErr->set(x2); cErr->axpy(-1.0,z);
    CTRLerror = cErr->norm();
    cErr = x2c.clone();
    cErr->set(x2c); cErr->axpy(-1.0,zc);
    double TOTerror2 = cErr->norm();
    *outStream << "    Value-At-Risk (1.e-4) = " << x2c.getStatistic() << "\n";
    *outStream << "    Value-At-Risk Error   = " <<  VARerror << "\n";
    *outStream << "    Control Error         = " << CTRLerror << "\n";
    *outStream << "    Total Error           = " << TOTerror2 << "\n";
    *outStream << "  ---------------------------------------------\n";
    VARerror  = std::abs(zc.getStatistic()-x3c.getStatistic());
    cErr = x3.clone();
    cErr->set(x3); cErr->axpy(-1.0,z);
    CTRLerror = cErr->norm();
    cErr = x3c.clone();
    cErr->set(x3c); cErr->axpy(-1.0,zc);
    double TOTerror3 = cErr->norm();
    *outStream << "    Value-At-Risk (1.e-6) = " << x3c.getStatistic() << "\n";
    *outStream << "    Value-At-Risk Error   = " <<  VARerror << "\n";
    *outStream << "    Control Error         = " << CTRLerror << "\n";
    *outStream << "    Total Error           = " << TOTerror3 << "\n";
    *outStream << "  ---------------------------------------------\n\n";
    // Comparison
    errorFlag += ((TOTerror1 < 90.*TOTerror2) && (TOTerror2 < 90.*TOTerror3)) ? 1 : 0;
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
