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
#include "ROL_MonteCarloGenerator.hpp"
#include "ROL_OptimizationProblem.hpp"

#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_Comm.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_CommHelpers.hpp"

#include <iostream>
#include <fstream>
#include <algorithm>

#include "example_08.hpp"

typedef double RealT;
typedef H1VectorPrimal<RealT> PrimalStateVector;
typedef H1VectorDual<RealT> DualStateVector;
typedef L2VectorPrimal<RealT> PrimalControlVector;
typedef L2VectorDual<RealT> DualControlVector;
typedef H1VectorDual<RealT> PrimalConstraintVector;
typedef H1VectorPrimal<RealT> DualConstraintVector;

int main(int argc, char *argv[]) {

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  ROL::Ptr<const Teuchos::Comm<int> > comm
    = Teuchos::DefaultComm<int>::getComm();

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint = argc - 1;
  bool print = (iprint>0); // && !(comm->getRank());
  ROL::Ptr<std::ostream> outStream;
  Teuchos::oblackholestream bhs; // outputs nothing
  if (print)
    outStream = ROL::makePtrFromRef(std::cout);
  else
    outStream = ROL::makePtrFromRef(bhs);

  bool print0 = print && !comm->getRank();
  ROL::Ptr<std::ostream> outStream0;
  if (print0)
    outStream0 = ROL::makePtrFromRef(std::cout);
  else
    outStream0 = ROL::makePtrFromRef(bhs);

  int errorFlag  = 0;

  // *** Example body.

  try {
    /*************************************************************************/
    /************* INITIALIZE BURGERS FEM CLASS ******************************/
    /*************************************************************************/
    int nx      = 512;   // Set spatial discretization.
    RealT alpha = 1.e-3; // Set penalty parameter.
    RealT nl    = 1.0;   // Nonlinearity parameter (1 = Burgers, 0 = linear).
    RealT cH1   = 1.0;   // Scale for derivative term in H1 norm.
    RealT cL2   = 0.0;   // Scale for mass term in H1 norm.
    ROL::Ptr<BurgersFEM<RealT> > fem
      = ROL::makePtr<BurgersFEM<RealT>>(nx,nl,cH1,cL2);
    fem->test_inverse_mass(*outStream0);
    fem->test_inverse_H1(*outStream0);
    /*************************************************************************/
    /************* INITIALIZE SIMOPT OBJECTIVE FUNCTION **********************/
    /*************************************************************************/
    ROL::Ptr<std::vector<RealT> > ud_ptr
      = ROL::makePtr<std::vector<RealT>>(nx, 1.0);
    ROL::Ptr<ROL::Vector<RealT> > ud
      = ROL::makePtr<L2VectorPrimal<RealT>>(ud_ptr,fem);
    ROL::Ptr<ROL::Objective_SimOpt<RealT> > pobj
      = ROL::makePtr<Objective_BurgersControl<RealT>>(fem,ud,alpha);
    /*************************************************************************/
    /************* INITIALIZE SIMOPT EQUALITY CONSTRAINT *********************/
    /*************************************************************************/
    bool hess = true;
    ROL::Ptr<ROL::Constraint_SimOpt<RealT> > pcon
      = ROL::makePtr<Constraint_BurgersControl<RealT>>(fem,hess);
    /*************************************************************************/
    /************* INITIALIZE VECTOR STORAGE *********************************/
    /*************************************************************************/
    // INITIALIZE CONTROL VECTORS
    ROL::Ptr<std::vector<RealT> > z_ptr
      = ROL::makePtr<std::vector<RealT>>(nx+2, 1.0);
    ROL::Ptr<std::vector<RealT> > gz_ptr
      = ROL::makePtr<std::vector<RealT>>(nx+2, 1.0);
    ROL::Ptr<std::vector<RealT> > yz_ptr
      = ROL::makePtr<std::vector<RealT>>(nx+2, 1.0);
    for (int i=0; i<nx+2; i++) {
      (*yz_ptr)[i] = 2.0*random<RealT>(comm)-1.0;
    }
    ROL::Ptr<ROL::Vector<RealT> > zp
      = ROL::makePtr<PrimalControlVector>(z_ptr,fem);
    ROL::Ptr<ROL::Vector<RealT> > gzp
      = ROL::makePtr<DualControlVector>(gz_ptr,fem);
    ROL::Ptr<ROL::Vector<RealT> > yzp
      = ROL::makePtr<PrimalControlVector>(yz_ptr,fem);
    // INITIALIZE STATE VECTORS
    ROL::Ptr<std::vector<RealT> > u_ptr
      = ROL::makePtr<std::vector<RealT>>(nx, 1.0);
    ROL::Ptr<std::vector<RealT> > gu_ptr
      = ROL::makePtr<std::vector<RealT>>(nx, 1.0);
    ROL::Ptr<ROL::Vector<RealT> > up
      = ROL::makePtr<PrimalStateVector>(u_ptr,fem);
    ROL::Ptr<ROL::Vector<RealT> > gup
      = ROL::makePtr<DualStateVector>(gu_ptr,fem);
    // INITIALIZE CONSTRAINT VECTORS
    ROL::Ptr<std::vector<RealT> > c_ptr
      = ROL::makePtr<std::vector<RealT>>(nx, 1.0);
    ROL::Ptr<std::vector<RealT> > l_ptr
      = ROL::makePtr<std::vector<RealT>>(nx, 1.0);
    for (int i=0; i<nx; i++) {
      (*l_ptr)[i] = random<RealT>(comm);
    }
    ROL::Ptr<ROL::Vector<RealT> > cp
      = ROL::makePtr<PrimalConstraintVector>(c_ptr,fem);
    ROL::Ptr<ROL::Vector<RealT> > lp
      = ROL::makePtr<DualConstraintVector>(l_ptr,fem);
    /*************************************************************************/
    /************* INITIALIZE SAMPLE GENERATOR *******************************/
    /*************************************************************************/
    int dim = 4, nSamp = 1000;
    std::vector<RealT> tmp(2,0.0); tmp[0] = -1.0; tmp[1] = 1.0;
    std::vector<std::vector<RealT> > bounds(dim,tmp);
    ROL::Ptr<ROL::BatchManager<RealT> > bman
      = ROL::makePtr<L2VectorBatchManager<RealT,int>>(comm);
    ROL::Ptr<ROL::SampleGenerator<RealT> > sampler
      = ROL::makePtr<ROL::MonteCarloGenerator<RealT>>(
          nSamp,bounds,bman,false,false,100);
    /*************************************************************************/
    /************* INITIALIZE REDUCED OBJECTIVE FUNCTION *********************/
    /*************************************************************************/
    bool storage = true, fdhess = false;
    ROL::Ptr<ROL::Objective<RealT> > robj
      = ROL::makePtr<ROL::Reduced_Objective_SimOpt<RealT>>(
          pobj,pcon,up,zp,lp,gup,gzp,cp,storage,fdhess);
    /*************************************************************************/
    /************* INITIALIZE BOUND CONSTRAINTS ******************************/
    /*************************************************************************/
    std::vector<RealT> Zlo(nx+2,0.0), Zhi(nx+2,10.0);
    for (int i = 0; i < nx+2; i++) {
      if ( i < (int)((nx+2)/3) ) {
        Zlo[i] = -1.0;
        Zhi[i] = 1.0;
      }
      if ( i >= (int)((nx+2)/3) && i < (int)(2*(nx+2)/3) ) {
        Zlo[i] = 1.0;
        Zhi[i] = 5.0;
      }
      if ( i >= (int)(2*(nx+2)/3) ) {
        Zlo[i] = 5.0;
        Zhi[i] = 10.0;
      }
    }
    ROL::Ptr<ROL::BoundConstraint<RealT> > Zbnd
      = ROL::makePtr<L2BoundConstraint<RealT>>(Zlo,Zhi,fem);
    /*************************************************************************/
    /************* INITIALIZE OPTIMIZATION PROBLEM ***************************/
    /*************************************************************************/
    Teuchos::ParameterList SOLlist;
    SOLlist.sublist("SOL").set("Stochastic Component Type","Risk Averse");
    SOLlist.sublist("SOL").set("Store Sampled Value and Gradient",storage);
    SOLlist.sublist("SOL").sublist("Risk Measure").set("Name","KL Divergence");
    SOLlist.sublist("SOL").sublist("Risk Measure").sublist("KL Divergence").set("Threshold",1.e-2);
    ROL::OptimizationProblem<RealT> optProb(robj,zp,Zbnd);
    optProb.setStochasticObjective(SOLlist,sampler);
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
    optProb.check(*outStream0);
    /*************************************************************************/
    /************* RUN OPTIMIZATION ******************************************/
    /*************************************************************************/
    // READ IN XML INPUT
    std::string filename = "input.xml";
    Teuchos::RCP<Teuchos::ParameterList> parlist
      = Teuchos::rcp( new Teuchos::ParameterList() );
    Teuchos::updateParametersFromXmlFile( filename, parlist.ptr() );
    // RUN OPTIMIZATION
    ROL::Algorithm<RealT> algo("Trust Region",*parlist,false);
    zp->zero();
    algo.run(optProb,print0,*outStream0);
    /*************************************************************************/
    /************* PRINT CONTROL AND STATE TO SCREEN *************************/
    /*************************************************************************/
    if ( print0 ) {
      std::ofstream ofs;
      ofs.open("output_example_08.txt",std::ofstream::out);
      for ( int i = 0; i < nx+2; i++ ) {
        ofs << std::scientific << std::setprecision(10);
        ofs << std::setw(20) << std::left << (RealT)i/((RealT)nx+1.0);
        ofs << std::setw(20) << std::left << (*z_ptr)[i];
        ofs << "\n";
      }
      ofs.close();
    }
    *outStream0 << "Scalar Parameter: " << optProb.getSolutionStatistic() << "\n\n";
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
