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
    Teuchos::RCP<ROL::Algorithm<RealT> > algo;
    /**********************************************************************************************/
    /************************* CONSTRUCT VECTORS **************************************************/
    /**********************************************************************************************/
    // Build control vectors
    int nx = 256;
    Teuchos::RCP<std::vector<RealT> > x1_rcp  = Teuchos::rcp( new std::vector<RealT>(nx+2,0.0) );
    ROL::StdVector<RealT> x1(x1_rcp);
    Teuchos::RCP<std::vector<RealT> > x2_rcp  = Teuchos::rcp( new std::vector<RealT>(nx+2,0.0) );
    ROL::StdVector<RealT> x2(x2_rcp);
    Teuchos::RCP<std::vector<RealT> > x3_rcp  = Teuchos::rcp( new std::vector<RealT>(nx+2,0.0) );
    ROL::StdVector<RealT> x3(x3_rcp);
    Teuchos::RCP<std::vector<RealT> > z_rcp  = Teuchos::rcp( new std::vector<RealT>(nx+2,0.0) );
    ROL::StdVector<RealT> z(z_rcp);
    Teuchos::RCP<std::vector<RealT> > xr_rcp = Teuchos::rcp( new std::vector<RealT>(nx+2,0.0) );
    ROL::StdVector<RealT> xr(xr_rcp);
    Teuchos::RCP<std::vector<RealT> > d_rcp  = Teuchos::rcp( new std::vector<RealT>(nx+2,0.0) );
    ROL::StdVector<RealT> d(d_rcp);
    for ( int i = 0; i < nx+2; i++ ) {
      (*xr_rcp)[i] = random<RealT>(comm);
      (*d_rcp)[i]  = random<RealT>(comm);
    }
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
    int nSamp = 100;
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
    Teuchos::RCP<ROL::ParametrizedObjective_SimOpt<RealT> > pobjSimOpt
      = Teuchos::rcp(new Objective_BurgersControl<RealT>(alpha,nx));
    Teuchos::RCP<ROL::ParametrizedEqualityConstraint_SimOpt<RealT> > pconSimOpt
      = Teuchos::rcp(new EqualityConstraint_BurgersControl<RealT>(nx));
    Teuchos::RCP<ROL::ParametrizedObjective<RealT> > pObj
      = Teuchos::rcp(new ROL::Reduced_ParametrizedObjective_SimOpt<RealT>(pobjSimOpt,pconSimOpt,up,pp));
    Teuchos::RCP<ROL::Objective<RealT> > obj;
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
    Teuchos::ParameterList list1;
    list1.sublist("SOL").set("Stochastic Optimization Type","Risk Averse");
    list1.sublist("SOL").set("Store Sampled Value and Gradient",true);
    list1.sublist("SOL").sublist("Risk Measure").set("Name","CVaR");
    list1.sublist("SOL").sublist("Risk Measure").sublist("CVaR").set("Confidence Level",0.99);
    list1.sublist("SOL").sublist("Risk Measure").sublist("CVaR").set("Convex Combination Parameter",1.0);
    list1.sublist("SOL").sublist("Risk Measure").sublist("CVaR").set("Smoothing Parameter",1.e-2);
    list1.sublist("SOL").sublist("Risk Measure").sublist("CVaR").sublist("Distribution").set("Name","Parabolic");
    list1.sublist("SOL").sublist("Risk Measure").sublist("CVaR").sublist("Distribution").sublist("Parabolic").set("Lower Bound",-0.5);
    list1.sublist("SOL").sublist("Risk Measure").sublist("CVaR").sublist("Distribution").sublist("Parabolic").set("Upper Bound", 0.5);
    // Build stochastic problem
    Teuchos::RCP<ROL::Vector<RealT> > x1p = Teuchos::rcp(&x1,false);
    x1p->zero();
    ROL::StochasticProblem<RealT> optProb1(list1,pObj,sampler,x1p);
    optProb1.checkObjectiveGradient(d,true,*outStream);
    optProb1.checkObjectiveHessVec(d,true,*outStream);
    // Run ROL algorithm
    algo = Teuchos::rcp(new ROL::Algorithm<RealT>("Trust Region",*parlist,false));
    clock_t start = clock();
    algo->run(optProb1,true,*outStream);
    *outStream << "Optimization time: " << (RealT)(clock()-start)/(RealT)CLOCKS_PER_SEC << " seconds.\n";
    /**********************************************************************************************/
    /************************* SMOOTHED CVAR 1.e-4 ************************************************/
    /**********************************************************************************************/
    *outStream << "\nSOLVE SMOOTHED CONDITIONAL VALUE AT RISK WITH TRUST REGION\n";
    Teuchos::ParameterList list2;
    list2.sublist("SOL").set("Stochastic Optimization Type","Risk Averse");
    list2.sublist("SOL").set("Store Sampled Value and Gradient",true);
    list2.sublist("SOL").sublist("Risk Measure").set("Name","CVaR");
    list2.sublist("SOL").sublist("Risk Measure").sublist("CVaR").set("Confidence Level",0.99);
    list2.sublist("SOL").sublist("Risk Measure").sublist("CVaR").set("Convex Combination Parameter",1.0);
    list2.sublist("SOL").sublist("Risk Measure").sublist("CVaR").set("Smoothing Parameter",1.e-4);
    list2.sublist("SOL").sublist("Risk Measure").sublist("CVaR").sublist("Distribution").set("Name","Parabolic");
    list2.sublist("SOL").sublist("Risk Measure").sublist("CVaR").sublist("Distribution").sublist("Parabolic").set("Lower Bound",-0.5);
    list2.sublist("SOL").sublist("Risk Measure").sublist("CVaR").sublist("Distribution").sublist("Parabolic").set("Upper Bound", 0.5);
    // Build stochastic problem
    Teuchos::RCP<ROL::Vector<RealT> > x2p = Teuchos::rcp(&x2,false);
    x2p->set(*x1p);
    ROL::StochasticProblem<RealT> optProb2(list2,pObj,sampler,x2p);
    optProb2.setSolutionStatistic(optProb1.getSolutionStatistic());
    optProb2.checkObjectiveGradient(d,true,*outStream);
    optProb2.checkObjectiveHessVec(d,true,*outStream);
    // Run ROL algorithm
    algo = Teuchos::rcp(new ROL::Algorithm<RealT>("Trust Region",*parlist,false));
    start = clock();
    algo->run(optProb2,true,*outStream);
    *outStream << "Optimization time: " << (RealT)(clock()-start)/(RealT)CLOCKS_PER_SEC << " seconds.\n";
    /**********************************************************************************************/
    /************************* SMOOTHED CVAR 1.e-6 ************************************************/
    /**********************************************************************************************/
    *outStream << "\nSOLVE SMOOTHED CONDITIONAL VALUE AT RISK WITH TRUST REGION\n";
    Teuchos::ParameterList list3;
    list3.sublist("SOL").set("Stochastic Optimization Type","Risk Averse");
    list3.sublist("SOL").set("Store Sampled Value and Gradient",true);
    list3.sublist("SOL").sublist("Risk Measure").set("Name","CVaR");
    list3.sublist("SOL").sublist("Risk Measure").sublist("CVaR").set("Confidence Level",0.99);
    list3.sublist("SOL").sublist("Risk Measure").sublist("CVaR").set("Convex Combination Parameter",1.0);
    list3.sublist("SOL").sublist("Risk Measure").sublist("CVaR").set("Smoothing Parameter",1.e-6);
    list3.sublist("SOL").sublist("Risk Measure").sublist("CVaR").sublist("Distribution").set("Name","Parabolic");
    list3.sublist("SOL").sublist("Risk Measure").sublist("CVaR").sublist("Distribution").sublist("Parabolic").set("Lower Bound",-0.5);
    list3.sublist("SOL").sublist("Risk Measure").sublist("CVaR").sublist("Distribution").sublist("Parabolic").set("Upper Bound", 0.5);
    // Build stochastic problem
    Teuchos::RCP<ROL::Vector<RealT> > x3p = Teuchos::rcp(&x3,false);
    x3p->set(*x2p);
    ROL::StochasticProblem<RealT> optProb3(list3,pObj,sampler,x3p);
    optProb3.setSolutionStatistic(optProb2.getSolutionStatistic());
    optProb3.checkObjectiveGradient(d,true,*outStream);
    optProb3.checkObjectiveHessVec(d,true,*outStream);
    // Run ROL algorithm
    algo = Teuchos::rcp(new ROL::Algorithm<RealT>("Trust Region",*parlist,false));
    start = clock();
    algo->run(optProb3,true,*outStream);
    *outStream << "Optimization time: " << (RealT)(clock()-start)/(RealT)CLOCKS_PER_SEC << " seconds.\n";
    /**********************************************************************************************/
    /************************* NONSMOOTH PROBLEM **************************************************/
    /**********************************************************************************************/
    *outStream << "\nSOLVE NONSMOOTH CVAR PROBLEM WITH BUNDLE TRUST REGION\n";
    Teuchos::ParameterList list;
    list.sublist("SOL").set("Stochastic Optimization Type","Risk Averse");
    list.sublist("SOL").set("Store Sampled Value and Gradient",true);
    list.sublist("SOL").sublist("Risk Measure").set("Name","CVaR");
    list.sublist("SOL").sublist("Risk Measure").sublist("CVaR").set("Confidence Level",0.99);
    list.sublist("SOL").sublist("Risk Measure").sublist("CVaR").set("Convex Combination Parameter",1.0);
    list.sublist("SOL").sublist("Risk Measure").sublist("CVaR").set("Smoothing Parameter",0.);
    list.sublist("SOL").sublist("Risk Measure").sublist("CVaR").sublist("Distribution").set("Name","Dirac");
    list.sublist("SOL").sublist("Risk Measure").sublist("CVaR").sublist("Distribution").sublist("Dirac").set("Location",0.);
    // Build stochastic problem
    Teuchos::RCP<ROL::Vector<RealT> > zp = Teuchos::rcp(&z,false);
    zp->set(*x3p);
    ROL::StochasticProblem<RealT> optProb(list,pObj,sampler,zp);
    optProb.setSolutionStatistic(optProb3.getSolutionStatistic());
    optProb.checkObjectiveGradient(d,true,*outStream);
    optProb.checkObjectiveHessVec(d,true,*outStream);
    // Run ROL algorithm
    parlist->sublist("Status Test").set("Iteration Limit",1000);
    parlist->sublist("Step").sublist("Bundle").set("Epsilon Solution Tolerance",1.e-7);
    algo = Teuchos::rcp(new ROL::Algorithm<RealT>("Bundle",*parlist,false));
    start = clock();
    algo->run(optProb,true,*outStream);
    *outStream << "Optimization time: " << (RealT)(clock()-start)/(RealT)CLOCKS_PER_SEC << " seconds.\n";
    /**********************************************************************************************/
    /************************* COMPUTE ERROR ******************************************************/
    /**********************************************************************************************/
    *outStream << "\nSUMMARY:\n";
    *outStream << "  ---------------------------------------------\n";
    *outStream << "    True Value-At-Risk    = " << optProb.getSolutionStatistic() << "\n";
    *outStream << "  ---------------------------------------------\n";
    RealT VARerror  = std::abs(optProb.getSolutionStatistic()-optProb1.getSolutionStatistic());
    Teuchos::RCP<ROL::Vector<RealT> > cErr = x1.clone();
    cErr->set(x1); cErr->axpy(-1.0,z);
    RealT CTRLerror = cErr->norm();
    RealT TOTerror1 = std::sqrt(std::pow(VARerror,2)+std::pow(CTRLerror,2));
    *outStream << "    Value-At-Risk (1.e-2) = " << optProb1.getSolutionStatistic() << "\n";
    *outStream << "    Value-At-Risk Error   = " <<  VARerror << "\n";
    *outStream << "    Control Error         = " << CTRLerror << "\n";
    *outStream << "    Total Error           = " << TOTerror1 << "\n";
    *outStream << "  ---------------------------------------------\n";
    VARerror  = std::abs(optProb.getSolutionStatistic()-optProb2.getSolutionStatistic());
    cErr = x2.clone();
    cErr->set(x2); cErr->axpy(-1.0,z);
    CTRLerror = cErr->norm();
    RealT TOTerror2 = std::sqrt(std::pow(VARerror,2)+std::pow(CTRLerror,2));
    *outStream << "    Value-At-Risk (1.e-4) = " << optProb2.getSolutionStatistic() << "\n";
    *outStream << "    Value-At-Risk Error   = " <<  VARerror << "\n";
    *outStream << "    Control Error         = " << CTRLerror << "\n";
    *outStream << "    Total Error           = " << TOTerror2 << "\n";
    *outStream << "  ---------------------------------------------\n";
    VARerror  = std::abs(optProb.getSolutionStatistic()-optProb3.getSolutionStatistic());
    cErr = x3.clone();
    cErr->set(x3); cErr->axpy(-1.0,z);
    CTRLerror = cErr->norm();
    RealT TOTerror3 = std::sqrt(std::pow(VARerror,2)+std::pow(CTRLerror,2));
    *outStream << "    Value-At-Risk (1.e-6) = " << optProb3.getSolutionStatistic() << "\n";
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
