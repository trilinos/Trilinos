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

/*! \file  example_03.cpp
    \brief Shows how to solve a steady Burgers' optimal control problem using
           full-space methods.
*/

#include "ROL_Algorithm.hpp"

#include "ROL_Reduced_Objective_SimOpt.hpp"
#include "ROL_HMCRObjective.hpp"
#include "ROL_RiskVector.hpp"

#include "ROL_MonteCarloGenerator.hpp"

#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_Comm.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_CommHelpers.hpp"

#include <iostream>
#include <algorithm>

#include "example_06.hpp"

typedef double RealT;
typedef H1VectorPrimal<RealT> PrimalStateVector;
typedef H1VectorDual<RealT> DualStateVector;
typedef L2VectorPrimal<RealT> PrimalControlVector;
typedef L2VectorDual<RealT> DualControlVector;
typedef H1VectorDual<RealT> PrimalConstraintVector;
typedef H1VectorPrimal<RealT> DualConstraintVector;

int main(int argc, char *argv[]) {

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  Teuchos::RCP<const Teuchos::Comm<int> > comm
    = Teuchos::DefaultComm<int>::getComm();

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint = argc - 1;
  bool print = (iprint>0);
  Teuchos::RCP<std::ostream> outStream;
  Teuchos::oblackholestream bhs; // outputs nothing
  if (print)
    outStream = Teuchos::rcp(&std::cout, false);
  else
    outStream = Teuchos::rcp(&bhs, false);
    
  bool print0 = print && !comm->getRank();
  Teuchos::RCP<std::ostream> outStream0;
  if (print0)
    outStream0 = Teuchos::rcp(&std::cout, false); 
  else
    outStream0 = Teuchos::rcp(&bhs, false);

  int errorFlag  = 0;

  // *** Example body.

  try {
    /*************************************************************************/
    /************* INITIALIZE BURGERS FEM CLASS ******************************/
    /*************************************************************************/
    int nx      = 256;   // Set spatial discretization.
    RealT alpha = 1.e-3; // Set penalty parameter.
    RealT nl    = 1.0;   // Nonlinearity parameter (1 = Burgers, 0 = linear).
    RealT cH1   = 1.0;   // Scale for derivative term in H1 norm.
    RealT cL2   = 0.0;   // Scale for mass term in H1 norm.
    Teuchos::RCP<BurgersFEM<RealT> > fem
      = Teuchos::rcp(new BurgersFEM<RealT>(nx,nl,cH1,cL2));
    fem->test_inverse_mass(*outStream0);
    fem->test_inverse_H1(*outStream0);
    /*************************************************************************/
    /************* INITIALIZE SIMOPT OBJECTIVE FUNCTION **********************/
    /*************************************************************************/
    Teuchos::RCP<std::vector<RealT> > ud_rcp
      = Teuchos::rcp( new std::vector<RealT> (nx, 1.0) );
    Teuchos::RCP<ROL::Vector<RealT> > ud
      = Teuchos::rcp(new L2VectorPrimal<RealT>(ud_rcp,fem));
    Teuchos::RCP<ROL::Objective_SimOpt<RealT> > pobj
      = Teuchos::rcp(new Objective_BurgersControl<RealT>(fem,ud,alpha));
    /*************************************************************************/
    /************* INITIALIZE SIMOPT EQUALITY CONSTRAINT *********************/
    /*************************************************************************/
    bool hess = true;
    Teuchos::RCP<ROL::Constraint_SimOpt<RealT> > pcon
      = Teuchos::rcp(new Constraint_BurgersControl<RealT>(fem,hess));
    /*************************************************************************/
    /************* INITIALIZE VECTOR STORAGE *********************************/
    /*************************************************************************/
    // INITIALIZE CONTROL VECTORS
    Teuchos::RCP<std::vector<RealT> > z_rcp
      = Teuchos::rcp( new std::vector<RealT> (nx+2, 1.0) );
    Teuchos::RCP<std::vector<RealT> > gz_rcp
      = Teuchos::rcp( new std::vector<RealT> (nx+2, 1.0) );
    Teuchos::RCP<std::vector<RealT> > yz_rcp
      = Teuchos::rcp( new std::vector<RealT> (nx+2, 1.0) );
    for (int i=0; i<nx+2; i++) {
      (*yz_rcp)[i] = 2.0*random<RealT>(comm)-1.0;
    }
    Teuchos::RCP<ROL::Vector<RealT> > zp
      = Teuchos::rcp(new PrimalControlVector(z_rcp,fem));
    Teuchos::RCP<ROL::Vector<RealT> > gzp
      = Teuchos::rcp(new DualControlVector(gz_rcp,fem));
    Teuchos::RCP<ROL::Vector<RealT> > yzp
      = Teuchos::rcp(new PrimalControlVector(yz_rcp,fem));
    RealT zvar = 0.0*random<RealT>(comm);
    RealT gvar = random<RealT>(comm);
    RealT yvar = random<RealT>(comm);
    Teuchos::RCP<Teuchos::ParameterList> hmcrlist = Teuchos::rcp(new Teuchos::ParameterList);
    hmcrlist->sublist("SOL").sublist("Risk Measure").set("Name","HMCR");
    ROL::RiskVector<RealT> z(hmcrlist,zp,zvar), g(hmcrlist,gzp,gvar), y(hmcrlist,yzp,yvar);
    // INITIALIZE STATE VECTORS
    Teuchos::RCP<std::vector<RealT> > u_rcp
      = Teuchos::rcp( new std::vector<RealT> (nx, 1.0) );
    Teuchos::RCP<std::vector<RealT> > gu_rcp
      = Teuchos::rcp( new std::vector<RealT> (nx, 1.0) );
    Teuchos::RCP<ROL::Vector<RealT> > up
      = Teuchos::rcp(new PrimalStateVector(u_rcp,fem));
    Teuchos::RCP<ROL::Vector<RealT> > gup
      = Teuchos::rcp(new DualStateVector(gu_rcp,fem));
    // INITIALIZE CONSTRAINT VECTORS
    Teuchos::RCP<std::vector<RealT> > c_rcp
      = Teuchos::rcp( new std::vector<RealT> (nx, 1.0) );
    Teuchos::RCP<std::vector<RealT> > l_rcp
      = Teuchos::rcp( new std::vector<RealT> (nx, 1.0) );
    for (int i=0; i<nx; i++) {
      (*l_rcp)[i] = random<RealT>(comm);
    }
    Teuchos::RCP<ROL::Vector<RealT> > cp
      = Teuchos::rcp(new PrimalConstraintVector(c_rcp,fem));
    Teuchos::RCP<ROL::Vector<RealT> > lp
      = Teuchos::rcp(new DualConstraintVector(l_rcp,fem));
    /*************************************************************************/
    /************* INITIALIZE SAMPLE GENERATOR *******************************/
    /*************************************************************************/
    int dim = 4, nSamp = 1000;
    std::vector<RealT> tmp(2,0.0); tmp[0] = -1.0; tmp[1] = 1.0;
    std::vector<std::vector<RealT> > bounds(dim,tmp);
    Teuchos::RCP<ROL::BatchManager<RealT> > bman
      = Teuchos::rcp(new L2VectorBatchManager<RealT,int>(comm));
    Teuchos::RCP<ROL::SampleGenerator<RealT> > sampler
      = Teuchos::rcp(new ROL::MonteCarloGenerator<RealT>(
          nSamp,bounds,bman,false,false,100));
    /*************************************************************************/
    /************* INITIALIZE RISK-AVERSE OBJECTIVE FUNCTION *****************/
    /*************************************************************************/
    bool storage = true, fdhess = false;
    Teuchos::RCP<ROL::Objective<RealT> > robj
      = Teuchos::rcp(new ROL::Reduced_Objective_SimOpt<RealT>(
          pobj,pcon,up,zp,lp,gup,gzp,cp,storage,fdhess));
    RealT order = 2.0, prob = 0.95;
    Teuchos::RCP<ROL::Objective<RealT> > obj
      = Teuchos::rcp(new ROL::HMCRObjective<RealT>(
          robj,order,prob,sampler,storage));
    /*************************************************************************/
    /************* CHECK DERIVATIVES AND CONSISTENCY *************************/
    /*************************************************************************/
    // CHECK OBJECTIVE DERIVATIVES
    bool derivcheck = false;
    if (derivcheck) {
      int nranks = sampler->numBatches();
      for (int pid = 0; pid < nranks; pid++) {
        if ( pid == sampler->batchID() ) {
          for (int i = sampler->start(); i < sampler->numMySamples(); i++) {
            *outStream << "Sample " << i << "  Rank " << sampler->batchID() << "\n";
            *outStream << "(" << sampler->getMyPoint(i)[0] << ", "
                              << sampler->getMyPoint(i)[1] << ", "
                              << sampler->getMyPoint(i)[2] << ", "
                              << sampler->getMyPoint(i)[3] << ")\n";
            pcon->setParameter(sampler->getMyPoint(i));
            pcon->checkSolve(*up,*zp,*cp,print,*outStream);
            robj->setParameter(sampler->getMyPoint(i));
            *outStream << "\n";
            robj->checkGradient(*zp,*gzp,*yzp,print,*outStream);
            robj->checkHessVec(*zp,*gzp,*yzp,print,*outStream);
            *outStream << "\n\n";
          }
        }
        comm->barrier();
      }
    }
    obj->checkGradient(z,g,y,print0,*outStream0);
    obj->checkHessVec(z,g,y,print0,*outStream0);
    /*************************************************************************/
    /************* RUN OPTIMIZATION ******************************************/
    /*************************************************************************/
    // READ IN XML INPUT
    std::string filename = "input.xml";
    Teuchos::RCP<Teuchos::ParameterList> parlist
      = Teuchos::rcp( new Teuchos::ParameterList() );
    Teuchos::updateParametersFromXmlFile( filename, parlist.ptr() );
    // DEFINE ALGORITHM
    ROL::Algorithm<RealT> algo("Trust Region",*parlist,false);
    // RUN OPTIMIZATION
    z.zero();
    algo.run(z, g, *obj, print0, *outStream0);
    /*************************************************************************/
    /************* PRINT CONTROL AND STATE TO SCREEN *************************/
    /*************************************************************************/
    *outStream0 << "\n";
    for ( int i = 0; i < nx+2; i++ ) {
      *outStream0 << std::scientific << std::setprecision(10);
      *outStream0 << std::setw(20) << std::left << (RealT)i/((RealT)nx+1.0);
      *outStream0 << std::setw(20) << std::left << (*z_rcp)[i];
      *outStream0 << "\n";
    }
    *outStream0 << "\n";
    *outStream0 << "Scalar Parameter: " << z.getStatistic(0) << "\n";
  }
  catch (std::logic_error err) {
    *outStream << err.what() << "\n";
    errorFlag = -1000;
  }; // end try

  comm->barrier();
  if (errorFlag != 0)
    std::cout << "End Result: TEST FAILED\n";
  else
    std::cout << "End Result: TEST PASSED\n";

  return 0;
}
