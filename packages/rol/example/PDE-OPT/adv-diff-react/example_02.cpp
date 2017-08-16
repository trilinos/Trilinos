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

/*! \file  example_02.cpp
    \brief Shows how to solve the stochastic advection-diffusion problem.
*/

#include "Teuchos_Comm.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "Tpetra_DefaultPlatform.hpp"
#include "Tpetra_Version.hpp"

#include <iostream>
#include <algorithm>
//#include <fenv.h>

#include "ROL_Algorithm.hpp"
#include "ROL_Bounds.hpp"
#include "ROL_Reduced_Objective_SimOpt.hpp"
#include "ROL_MonteCarloGenerator.hpp"
#include "ROL_OptimizationProblem.hpp"
#include "ROL_TpetraTeuchosBatchManager.hpp"

#include "../TOOLS/meshmanager.hpp"
#include "../TOOLS/pdeconstraint.hpp"
#include "../TOOLS/pdeobjective.hpp"
#include "../TOOLS/pdevector.hpp"
#include "../TOOLS/batchmanager.hpp"
#include "pde_stoch_adv_diff.hpp"
#include "obj_stoch_adv_diff.hpp"
#include "mesh_stoch_adv_diff.hpp"

typedef double RealT;

template<class Real>
Real random(const Teuchos::Comm<int> &comm,
            const Real a = -1, const Real b = 1) {
  Real val(0), u(0);
  if ( Teuchos::rank<int>(comm)==0 ) {
    u   = static_cast<Real>(rand())/static_cast<Real>(RAND_MAX);
    val = (b-a)*u + a;
  }
  Teuchos::broadcast<int,Real>(comm,0,1,&val);
  return val;
}

int main(int argc, char *argv[]) {
//  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  Teuchos::RCP<std::ostream> outStream;
  Teuchos::oblackholestream bhs; // outputs nothing

  /*** Initialize communicator. ***/
  Teuchos::GlobalMPISession mpiSession (&argc, &argv, &bhs);
  Teuchos::RCP<const Teuchos::Comm<int> > comm
    = Tpetra::DefaultPlatform::getDefaultPlatform().getComm();
  Teuchos::RCP<const Teuchos::Comm<int> > serial_comm
    = Teuchos::rcp(new Teuchos::SerialComm<int>());
  const int myRank = comm->getRank();
  if ((iprint > 0) && (myRank == 0)) {
    outStream = Teuchos::rcp(&std::cout, false);
  }
  else {
    outStream = Teuchos::rcp(&bhs, false);
  }
  int errorFlag  = 0;

  // *** Example body.
  try {

    /*** Read in XML input ***/
    std::string filename = "input.xml";
    Teuchos::RCP<Teuchos::ParameterList> parlist = Teuchos::rcp( new Teuchos::ParameterList() );
    Teuchos::updateParametersFromXmlFile( filename, parlist.ptr() );

    // Problem dimensions
    const int controlDim = 9, stochDim = 37;
    const RealT one(1); 

    /*************************************************************************/
    /***************** BUILD GOVERNING PDE ***********************************/
    /*************************************************************************/
    /*** Initialize main data structure. ***/
    Teuchos::RCP<MeshManager<RealT> > meshMgr
      = Teuchos::rcp(new MeshManager_stoch_adv_diff<RealT>(*parlist));
    // Initialize PDE describing advection-diffusion equation
    Teuchos::RCP<PDE_stoch_adv_diff<RealT> > pde
      = Teuchos::rcp(new PDE_stoch_adv_diff<RealT>(*parlist));
    Teuchos::RCP<ROL::Constraint_SimOpt<RealT> > con
      = Teuchos::rcp(new PDE_Constraint<RealT>(pde,meshMgr,serial_comm,*parlist,*outStream));
    Teuchos::RCP<PDE_Constraint<RealT> > pdeCon
      = Teuchos::rcp_dynamic_cast<PDE_Constraint<RealT> >(con);
    pdeCon->getAssembler()->printMeshData(*outStream);
    con->setSolveParameters(*parlist);

    /*************************************************************************/
    /***************** BUILD VECTORS *****************************************/
    /*************************************************************************/
    Teuchos::RCP<Tpetra::MultiVector<> >  u_rcp = pdeCon->getAssembler()->createStateVector();
    Teuchos::RCP<Tpetra::MultiVector<> >  p_rcp = pdeCon->getAssembler()->createStateVector();
    Teuchos::RCP<Tpetra::MultiVector<> > du_rcp = pdeCon->getAssembler()->createStateVector();
    u_rcp->randomize();  //u_rcp->putScalar(static_cast<RealT>(1));
    p_rcp->randomize();  //p_rcp->putScalar(static_cast<RealT>(1));
    du_rcp->randomize(); //du_rcp->putScalar(static_cast<RealT>(0));
    Teuchos::RCP<ROL::Vector<RealT> > up
      = Teuchos::rcp(new PDE_PrimalSimVector<RealT>(u_rcp,pde,pdeCon->getAssembler()));
    Teuchos::RCP<ROL::Vector<RealT> > pp
      = Teuchos::rcp(new PDE_PrimalSimVector<RealT>(p_rcp,pde,pdeCon->getAssembler()));
    Teuchos::RCP<ROL::Vector<RealT> > dup
      = Teuchos::rcp(new PDE_PrimalSimVector<RealT>(du_rcp,pde,pdeCon->getAssembler()));
    // Create residual vectors
    Teuchos::RCP<Tpetra::MultiVector<> > r_rcp = pdeCon->getAssembler()->createResidualVector();
    r_rcp->randomize(); //r_rcp->putScalar(static_cast<RealT>(1));
    Teuchos::RCP<ROL::Vector<RealT> > rp
      = Teuchos::rcp(new PDE_DualSimVector<RealT>(r_rcp,pde,pdeCon->getAssembler()));
    // Create control vector and set to ones
    Teuchos::RCP<std::vector<RealT> >  z_rcp = Teuchos::rcp(new std::vector<RealT>(controlDim));
    Teuchos::RCP<std::vector<RealT> > dz_rcp = Teuchos::rcp(new std::vector<RealT>(controlDim));
    Teuchos::RCP<std::vector<RealT> > yz_rcp = Teuchos::rcp(new std::vector<RealT>(controlDim));
    // Create control direction vector and set to random
    for (int i = 0; i < controlDim; ++i) {
      (*z_rcp)[i]  = random<RealT>(*comm);
      (*dz_rcp)[i] = random<RealT>(*comm);
      (*yz_rcp)[i] = random<RealT>(*comm);
    }
    Teuchos::RCP<ROL::Vector<RealT> > zp
      = Teuchos::rcp(new PDE_OptVector<RealT>(Teuchos::rcp(new ROL::StdVector<RealT>(z_rcp))));
    Teuchos::RCP<ROL::Vector<RealT> > dzp
      = Teuchos::rcp(new PDE_OptVector<RealT>(Teuchos::rcp(new ROL::StdVector<RealT>(dz_rcp))));
    Teuchos::RCP<ROL::Vector<RealT> > yzp
      = Teuchos::rcp(new PDE_OptVector<RealT>(Teuchos::rcp(new ROL::StdVector<RealT>(yz_rcp))));
    // Create ROL SimOpt vectors
    ROL::Vector_SimOpt<RealT> x(up,zp);
    ROL::Vector_SimOpt<RealT> d(dup,dzp);

    Teuchos::RCP<Tpetra::MultiVector<> > dualu_rcp = pdeCon->getAssembler()->createStateVector();
    Teuchos::RCP<ROL::Vector<RealT> > dualup
      = Teuchos::rcp(new PDE_DualSimVector<RealT>(dualu_rcp,pde,pdeCon->getAssembler()));

    /*************************************************************************/
    /***************** BUILD COST FUNCTIONAL *********************************/
    /*************************************************************************/
    std::vector<Teuchos::RCP<QoI<RealT> > > qoi_vec(2,Teuchos::null);
    qoi_vec[0] = Teuchos::rcp(new QoI_State_Cost_stoch_adv_diff<RealT>(pde->getFE()));
    qoi_vec[1] = Teuchos::rcp(new QoI_Control_Cost_stoch_adv_diff<RealT>());
    RealT stateCost   = parlist->sublist("Problem").get("State Cost",1.e5);
    RealT controlCost = parlist->sublist("Problem").get("Control Cost",1.e0);
    std::vector<RealT> wts = {stateCost, controlCost};
    Teuchos::RCP<ROL::Objective_SimOpt<RealT> > obj
      = Teuchos::rcp(new PDE_Objective<RealT>(qoi_vec,wts,pdeCon->getAssembler()));
    bool storage = parlist->sublist("Problem").get("Use State and Adjoint Storage",true);
    Teuchos::RCP<ROL::Reduced_Objective_SimOpt<RealT> > objReduced
      = Teuchos::rcp(new ROL::Reduced_Objective_SimOpt<RealT>(obj, con, up, zp, pp, storage, false));

    /*************************************************************************/
    /***************** BUILD BOUND CONSTRAINT ********************************/
    /*************************************************************************/
    Teuchos::RCP<std::vector<RealT> > zlo_rcp = Teuchos::rcp(new std::vector<RealT>(controlDim,0));
    Teuchos::RCP<std::vector<RealT> > zhi_rcp = Teuchos::rcp(new std::vector<RealT>(controlDim,1));
    Teuchos::RCP<ROL::Vector<RealT> > zlop
      = Teuchos::rcp(new PDE_OptVector<RealT>(Teuchos::rcp(new ROL::StdVector<RealT>(zlo_rcp))));
    Teuchos::RCP<ROL::Vector<RealT> > zhip
      = Teuchos::rcp(new PDE_OptVector<RealT>(Teuchos::rcp(new ROL::StdVector<RealT>(zhi_rcp))));
    Teuchos::RCP<ROL::BoundConstraint<RealT> > bnd
      = Teuchos::rcp(new ROL::Bounds<RealT>(zlop,zhip));

    /*************************************************************************/
    /***************** BUILD SAMPLER *****************************************/
    /*************************************************************************/
    int nsamp = parlist->sublist("Problem").get("Number of Samples",100);
    std::vector<RealT> tmp = {-one,one};
    std::vector<std::vector<RealT> > bounds(stochDim,tmp);
    Teuchos::RCP<ROL::BatchManager<RealT> > bman
      = Teuchos::rcp(new PDE_OptVector_BatchManager<RealT>(comm));
    Teuchos::RCP<ROL::SampleGenerator<RealT> > sampler
      = Teuchos::rcp(new ROL::MonteCarloGenerator<RealT>(nsamp,bounds,bman));

    /*************************************************************************/
    /***************** BUILD STOCHASTIC PROBLEM ******************************/
    /*************************************************************************/
    ROL::OptimizationProblem<RealT> opt(objReduced,zp,bnd);
    parlist->sublist("SOL").set("Initial Statistic",one);
    opt.setStochasticObjective(*parlist,sampler);

    /*************************************************************************/
    /***************** RUN VECTOR AND DERIVATIVE CHECKS **********************/
    /*************************************************************************/
    bool checkDeriv = parlist->sublist("Problem").get("Check Derivatives",false);
    if ( checkDeriv ) {
      up->checkVector(*pp,*dup,true,*outStream);
      zp->checkVector(*yzp,*dzp,true,*outStream);
      std::vector<RealT> param(stochDim,0);
      objReduced->setParameter(param);
      *outStream << "\n\nCheck Gradient of Full Objective Function\n";
      obj->checkGradient(x,d,true,*outStream);
      *outStream << "\n\nCheck Hessian of Full Objective Function\n";
      obj->checkHessVec(x,d,true,*outStream);
      *outStream << "\n\nCheck Full Jacobian of PDE Constraint\n";
      con->checkApplyJacobian(x,d,*rp,true,*outStream);
      *outStream << "\n\nCheck Jacobian_1 of PDE Constraint\n";
      con->checkApplyJacobian_1(*up,*zp,*dup,*rp,true,*outStream);
      *outStream << "\n\nCheck Jacobian_2 of PDE Constraint\n";
      con->checkApplyJacobian_2(*up,*zp,*dzp,*rp,true,*outStream);
      *outStream << "\n\nCheck Full Hessian of PDE Constraint\n";
      con->checkApplyAdjointHessian(x,*pp,d,x,true,*outStream);
      *outStream << "\n\nCheck Hessian_11 of PDE Constraint\n";
      con->checkApplyAdjointHessian_11(*up,*zp,*pp,*dup,*dualup,true,*outStream);
      *outStream << "\n\nCheck Hessian_21 of PDE Constraint\n";
      con->checkApplyAdjointHessian_21(*up,*zp,*pp,*dzp,*dualup,true,*outStream);
      *outStream << "\n\nCheck Hessian_12 of PDE Constraint\n";
      con->checkApplyAdjointHessian_12(*up,*zp,*pp,*dup,*dzp,true,*outStream);
      *outStream << "\n\nCheck Hessian_22 of PDE Constraint\n";
      con->checkApplyAdjointHessian_22(*up,*zp,*pp,*dzp,*dzp,true,*outStream);
      *outStream << "\n";
      con->checkAdjointConsistencyJacobian(*dup,d,x,true,*outStream);
      *outStream << "\n";
      con->checkInverseJacobian_1(*up,*up,*up,*zp,true,*outStream);
      *outStream << "\n";
      *outStream << "\n\nCheck Gradient of Reduced Objective Function\n";
      objReduced->checkGradient(*zp,*dzp,true,*outStream);
      *outStream << "\n\nCheck Hessian of Reduced Objective Function\n";
      objReduced->checkHessVec(*zp,*dzp,true,*outStream);

      opt.check(*outStream);
    }

    /*************************************************************************/
    /***************** SOLVE OPTIMIZATION PROBLEM ****************************/
    /*************************************************************************/
    ROL::Algorithm<RealT> algo("Trust Region",*parlist,false);
    zp->zero();
    std::clock_t timer = std::clock();
    algo.run(opt,true,*outStream);
    *outStream << "Optimization time: "
               << static_cast<RealT>(std::clock()-timer)/static_cast<RealT>(CLOCKS_PER_SEC)
               << " seconds." << std::endl << std::endl;

    /*************************************************************************/
    /***************** OUTPUT RESULTS ****************************************/
    /*************************************************************************/
    std::clock_t timer_print = std::clock();
    // Output control to file
    if ( myRank == 0 ) {
      std::ofstream zfile;
      zfile.open("control.txt");
      for (int i = 0; i < controlDim; i++) {
        zfile << (*z_rcp)[i] << "\n";
      }
      zfile.close();
    }
    // Output expected state and samples to file
    up->zero(); pp->zero(); dup->zero();
    RealT tol(1.e-8);
    Teuchos::RCP<ROL::BatchManager<RealT> > bman_Eu
      = Teuchos::rcp(new ROL::TpetraTeuchosBatchManager<RealT>(comm));
    std::vector<RealT> sample(stochDim);
    std::stringstream name_samp;
    name_samp << "samples_" << bman->batchID() << ".txt";
    std::ofstream file_samp;
    file_samp.open(name_samp.str());
    for (int i = 0; i < sampler->numMySamples(); ++i) {
      sample = sampler->getMyPoint(i);
      con->setParameter(sample);
      con->solve(*rp,*dup,*zp,tol);
      up->axpy(sampler->getMyWeight(i),*dup);
      for (int j = 0; j < stochDim; ++j) {
        file_samp << sample[j] << "  ";
      }
      file_samp << "\n";
    }
    file_samp.close();
    bman_Eu->sumAll(*up,*pp);
    pdeCon->getAssembler()->outputTpetraVector(p_rcp,"mean_state.txt");
    // Build objective function distribution
    RealT val(0);
    int nsamp_dist = parlist->sublist("Problem").get("Number of Output Samples",100);
    Teuchos::RCP<ROL::SampleGenerator<RealT> > sampler_dist
      = Teuchos::rcp(new ROL::MonteCarloGenerator<RealT>(nsamp_dist,bounds,bman));
    std::stringstream name;
    name << "obj_samples_" << bman->batchID() << ".txt";
    std::ofstream file;
    file.open(name.str());
    for (int i = 0; i < sampler_dist->numMySamples(); ++i) {
      sample = sampler_dist->getMyPoint(i);
      objReduced->setParameter(sample);
      val = objReduced->value(*zp,tol);
      for (int j = 0; j < stochDim; ++j) {
        file << sample[j] << "  ";
      }
      file << val << "\n";
    }
    file.close();
    *outStream << "Output time: "
               << static_cast<RealT>(std::clock()-timer_print)/static_cast<RealT>(CLOCKS_PER_SEC)
               << " seconds." << std::endl << std::endl;

    Teuchos::Array<RealT> res(1,0);
    pdeCon->solve(*rp,*up,*zp,tol);
    r_rcp->norm2(res.view(0,1));

    /*************************************************************************/
    /***************** CHECK RESIDUAL NORM ***********************************/
    /*************************************************************************/
    *outStream << "Residual Norm: " << res[0] << std::endl << std::endl;
    errorFlag += (res[0] > 1.e-6 ? 1 : 0);
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
