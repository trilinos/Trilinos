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
    // Build ROL algorithm
    parlist->sublist("Status Test").set("Gradient Tolerance",1.e-8);
    parlist->sublist("Status Test").set("Step Tolerance",1.e-14);
    parlist->sublist("Status Test").set("Iteration Limit",100);
    /**********************************************************************************************/
    /************************* CONSTRUCT VECTORS **************************************************/
    /**********************************************************************************************/
    // Build control vectors
    int nx = 256;
    Teuchos::RCP<std::vector<RealT> > x_rcp  = Teuchos::rcp( new std::vector<RealT>(nx+2,0.0) );
    ROL::StdVector<RealT> x(x_rcp);
    Teuchos::RCP<std::vector<RealT> > d_rcp  = Teuchos::rcp( new std::vector<RealT>(nx+2,0.0) );
    ROL::StdVector<RealT> d(d_rcp);
    for ( int i = 0; i < nx+2; i++ ) {
      (*x_rcp)[i] = random<RealT>(comm);
      (*d_rcp)[i] = random<RealT>(comm);
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
//    int dim = 4;
//    int nSamp = 100;
//    std::vector<RealT> tmp(2,0.0); tmp[0] = -1.0; tmp[1] = 1.0;
//    std::vector<std::vector<RealT> > bounds(dim,tmp);
//    Teuchos::RCP<ROL::BatchManager<RealT> > bman
//      = Teuchos::rcp(new ROL::StdTeuchosBatchManager<RealT,int>(comm));
//    Teuchos::RCP<ROL::SampleGenerator<RealT> > sampler
//      = Teuchos::rcp(new ROL::MonteCarloGenerator<RealT>(nSamp,bounds,bman,false,false,100));
    ROL::QuadratureInfo info;
    info.dim = 4;
    info.maxLevel = 7;
    info.rule1D.clear(); info.rule1D.resize(info.dim,ROL::QUAD_CLENSHAWCURTIS);    
    info.growth1D.clear(); info.growth1D.resize(info.dim,ROL::GROWTH_DEFAULT);
    info.normalized = true;
    info.adaptive = false;
    Teuchos::RCP<ROL::BatchManager<RealT> > bman
      = Teuchos::rcp(new ROL::StdTeuchosBatchManager<RealT,int>(comm));
    Teuchos::RCP<ROL::SampleGenerator<RealT> > sampler
      = Teuchos::rcp(new ROL::SparseGridGenerator<RealT>(bman,info));
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
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument,
        ">>> (adapters/trikota/sol/test/test_01): Unable to open file!");
    }

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
    pObj->setParameter(sampler->getMyPoint(0));
    pObj->checkGradient(x,d,true,*outStream);
    pObj->checkHessVec(x,d,true,*outStream);
    /**********************************************************************************************/
    /************************* RISK NEUTRAL *******************************************************/
    /**********************************************************************************************/
    *outStream << "\nSOLVE RISK NEUTRAL OPTIMAL CONTROL PROBLEM WITH TRUST REGION\n";
    // Build CVaR objective function
    Teuchos::ParameterList list;
    list.sublist("SOL").set("Stochastic Optimization Type","Risk Neutral");
    list.sublist("SOL").set("Store Sampled Value and Gradient",true);
    // Build stochastic problem
    Teuchos::RCP<ROL::Vector<RealT> > xp = Teuchos::rcp(&x,false);
    ROL::StochasticProblem<RealT> optProb(list,pObj,sampler,xp);
    optProb.checkObjectiveGradient(d,true,*outStream);
    optProb.checkObjectiveHessVec(d,true,*outStream);
    // Run ROL algorithm
    ROL::Algorithm<RealT> algo("Trust Region",*parlist,false);
    clock_t start = clock();
    xp->zero();
    algo.run(optProb,true,*outStream);
    *outStream << "Optimization time: " << (RealT)(clock()-start)/(RealT)CLOCKS_PER_SEC << " seconds.\n";
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
