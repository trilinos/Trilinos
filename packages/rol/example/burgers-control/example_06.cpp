// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  example_03.cpp
    \brief Shows how to solve a steady Burgers' optimal control problem using
           full-space methods.
*/

#include "ROL_Algorithm.hpp"
#include "ROL_TrustRegionStep.hpp"
#include "ROL_StatusTest.hpp"

#include "ROL_Reduced_Objective_SimOpt.hpp"
//#include "ROL_HMCRObjective.hpp"
#include "ROL_RiskVector.hpp"
#include "ROL_StochasticObjective.hpp"
#include "ROL_ParameterList.hpp"

#include "ROL_MonteCarloGenerator.hpp"

#include "ROL_Stream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_Comm.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_CommHelpers.hpp"

#include <iostream>
#include <algorithm>

#include "example_06.hpp"

typedef double RealT;
typedef H1VectorPrimal<RealT> PrimalStateVector;
typedef H1VectorDual<RealT>   DualStateVector;
typedef L2VectorPrimal<RealT> PrimalControlVector;
typedef L2VectorDual<RealT>   DualControlVector;
typedef H1VectorDual<RealT>   PrimalConstraintVector;
typedef H1VectorPrimal<RealT> DualConstraintVector;

int main(int argc, char *argv[]) {

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  auto comm = ROL::toPtr(Teuchos::DefaultComm<int>::getComm());

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint = argc - 1;
  bool print = (iprint>0);
  ROL::Ptr<std::ostream> outStream;
  ROL::nullstream bhs; // outputs nothing
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
    int nx      = 256;   // Set spatial discretization.
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
    RealT zvar = 0.0*random<RealT>(comm);
    RealT gvar = random<RealT>(comm);
    RealT yvar = random<RealT>(comm);
    ROL::Ptr<ROL::ParameterList> hmcrlist = ROL::makePtr<ROL::ParameterList>();
    hmcrlist->sublist("SOL").sublist("Risk Measure").set("Name","HMCR");
    ROL::RiskVector<RealT> z(hmcrlist,zp,zvar), g(hmcrlist,gzp,gvar), y(hmcrlist,yzp,yvar);
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
    /************* INITIALIZE RISK-AVERSE OBJECTIVE FUNCTION *****************/
    /*************************************************************************/
    bool storage = true, fdhess = false;
    ROL::Ptr<ROL::Objective<RealT> > robj
      = ROL::makePtr<ROL::Reduced_Objective_SimOpt<RealT>>(
          pobj,pcon,up,zp,lp,gup,gzp,cp,storage,fdhess);
    //RealT order = 2.0, prob = 0.95;
    //ROL::Ptr<ROL::Objective<RealT> > obj
    //  = ROL::makePtr<ROL::HMCRObjective<RealT>>(
    //      robj,order,prob,sampler,storage);
    hmcrlist->sublist("SOL").sublist("Risk Measure").sublist("HMCR").set("Order",2);
    hmcrlist->sublist("SOL").sublist("Risk Measure").sublist("HMCR").set("Confidence Level",0.95);
    hmcrlist->sublist("SOL").sublist("Risk Measure").sublist("HMCR").set("Convex Combination Parameter",0.0);
    ROL::Ptr<ROL::Objective<RealT> > obj
      = ROL::makePtr<ROL::StochasticObjective<RealT> >(robj,*hmcrlist,sampler);
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
    auto parlist = ROL::getParametersFromXmlFile( filename );
    // DEFINE ALGORITHM
    ROL::Ptr<ROL::Step<RealT>>
      step = ROL::makePtr<ROL::TrustRegionStep<RealT>>(*parlist);
    ROL::Ptr<ROL::StatusTest<RealT>>
      status = ROL::makePtr<ROL::StatusTest<RealT>>(*parlist);
    ROL::Algorithm<RealT> algo(step,status,false);
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
      *outStream0 << std::setw(20) << std::left << (*z_ptr)[i];
      *outStream0 << "\n";
    }
    *outStream0 << "\n";
    *outStream0 << "Scalar Parameter: " << z.getStatistic(0) << "\n";
  }
  catch (std::logic_error& err) {
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
