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
    Teuchos::RCP<std::vector<RealT> > z_rcp  = Teuchos::rcp(new std::vector<RealT>(nx+2,0));
    Teuchos::RCP<ROL::Vector<RealT> > zp  = Teuchos::rcp(new ROL::StdVector<RealT>(z_rcp));
    Teuchos::RCP<std::vector<RealT> > x1_rcp = Teuchos::rcp(new std::vector<RealT>(nx+2,0));
    Teuchos::RCP<ROL::Vector<RealT> > x1p = Teuchos::rcp(new ROL::StdVector<RealT>(x1_rcp));
    Teuchos::RCP<std::vector<RealT> > x2_rcp = Teuchos::rcp(new std::vector<RealT>(nx+2,0));
    Teuchos::RCP<ROL::Vector<RealT> > x2p = Teuchos::rcp(new ROL::StdVector<RealT>(x2_rcp));
    Teuchos::RCP<std::vector<RealT> > x3_rcp = Teuchos::rcp(new std::vector<RealT>(nx+2,0));
    Teuchos::RCP<ROL::Vector<RealT> > x3p = Teuchos::rcp(new ROL::StdVector<RealT>(x3_rcp));
    std::vector<Teuchos::RCP<ROL::Vector<RealT> > > xvec = {x1p, x2p, x3p};
    // Create vectors for derivative check
    Teuchos::RCP<std::vector<RealT> > xr_rcp = Teuchos::rcp(new std::vector<RealT>(nx+2,0));
    ROL::StdVector<RealT> xr(xr_rcp);
    Teuchos::RCP<std::vector<RealT> > d_rcp  = Teuchos::rcp(new std::vector<RealT>(nx+2,0));
    ROL::StdVector<RealT> d(d_rcp);
    for ( int i = 0; i < nx+2; i++ ) {
      (*xr_rcp)[i] = random<RealT>(comm);
      (*d_rcp)[i]  = random<RealT>(comm);
    }
    // Build state and adjoint vectors
    Teuchos::RCP<std::vector<RealT> > u_rcp = Teuchos::rcp(new std::vector<RealT>(nx,1));
    Teuchos::RCP<ROL::Vector<RealT> > up = Teuchos::rcp(new ROL::StdVector<RealT>(u_rcp));
    Teuchos::RCP<std::vector<RealT> > p_rcp = Teuchos::rcp(new std::vector<RealT>(nx,0));
    Teuchos::RCP<ROL::Vector<RealT> > pp = Teuchos::rcp(new ROL::StdVector<RealT>(p_rcp));
    /**********************************************************************************************/
    /************************* CONSTRUCT SOL COMPONENTS *******************************************/
    /**********************************************************************************************/
    // Build samplers
    int dim = 4, nSamp = 100;
    std::vector<RealT> tmp = {-1, 1};
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
    pconSimOpt->setSolveParameters(*parlist);
    Teuchos::RCP<ROL::Objective<RealT> > pObj
      = Teuchos::rcp(new ROL::Reduced_Objective_SimOpt<RealT>(pobjSimOpt,pconSimOpt,up,zp,pp));
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
    Teuchos::RCP<ROL::Algorithm<RealT> > algo;
    Teuchos::RCP<ROL::OptimizationProblem<RealT> > optProb;
    for (int i = 0; i < 3; ++i) {
      *outStream << "\nSOLVE SMOOTHED CONDITIONAL VALUE AT RISK WITH TRUST REGION\n";
      // Build CVaR risk measure
      Teuchos::ParameterList list;
      list.sublist("SOL").set("Stochastic Component Type",ra);
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
      optProb = Teuchos::rcp(new ROL::OptimizationProblem<RealT>(pObj,xvec[i]));
      RealT init_stat(1);
      if ( i > 0 ) { init_stat = stat[i-1]; }
      list.sublist("SOL").set("Initial Statistic",init_stat);
      optProb->setStochasticObjective(list,sampler);
      optProb->check(*outStream);
      // Run ROL algorithm
      algo = Teuchos::rcp(new ROL::Algorithm<RealT>("Trust Region",*parlist,false));
      clock_t start = clock();
      algo->run(*optProb,true,*outStream);
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
    Teuchos::ParameterList list;
    list.sublist("SOL").set("Stochastic Component Type",ra);
    list.sublist("SOL").set("Store Sampled Value and Gradient",storage);
    list.sublist("SOL").sublist("Risk Measure").set("Name",rm);
    list.sublist("SOL").sublist("Risk Measure").sublist(rm).set("Confidence Level",cl);
    list.sublist("SOL").sublist("Risk Measure").sublist(rm).set("Convex Combination Parameter",cc);
    list.sublist("SOL").sublist("Risk Measure").sublist(rm).set("Smoothing Parameter",0.);
    list.sublist("SOL").sublist("Risk Measure").sublist(rm).sublist("Distribution").set("Name","Dirac");
    list.sublist("SOL").sublist("Risk Measure").sublist(rm).sublist("Distribution").sublist("Dirac").set("Location",0.);
    // Build stochastic problem
    zp->set(*xvec[2]);
    optProb = Teuchos::rcp(new ROL::OptimizationProblem<RealT>(pObj,zp));
    list.sublist("SOL").set("Initial Statistic",stat[2]);
    optProb->setStochasticObjective(list,sampler);
    optProb->check(*outStream);
    // Run ROL algorithm
    parlist->sublist("Status Test").set("Iteration Limit",1000);
    parlist->sublist("Step").sublist("Bundle").set("Epsilon Solution Tolerance",1.e-7);
    algo = Teuchos::rcp(new ROL::Algorithm<RealT>("Bundle",*parlist,false));
    clock_t start = clock();
    algo->run(*optProb,true,*outStream);
    *outStream << "Optimization time: " << (RealT)(clock()-start)/(RealT)CLOCKS_PER_SEC << " seconds.\n";
    /**********************************************************************************************/
    /************************* COMPUTE ERROR ******************************************************/
    /**********************************************************************************************/
    Teuchos::RCP<ROL::Vector<RealT> > cErr = zp->clone();
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
              << std::setw(25) << (*z_rcp)[n]
              << std::endl;
    }
    control.close();

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
