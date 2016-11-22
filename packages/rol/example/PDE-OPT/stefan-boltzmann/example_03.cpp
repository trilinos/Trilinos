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

/*! \file  example_01.cpp
    \brief Shows how to solve the stochastic Stefan-Boltzmann problem.
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
#include "ROL_BoundConstraint.hpp"
#include "ROL_Reduced_Objective_SimOpt.hpp"
#include "ROL_MonteCarloGenerator.hpp"
#include "ROL_StochasticProblem.hpp"
#include "ROL_TpetraTeuchosBatchManager.hpp"

#include "../TOOLS/meshmanager.hpp"
#include "../TOOLS/pdeconstraint.hpp"
#include "../TOOLS/pdeobjective.hpp"
#include "../TOOLS/pdevector.hpp"
#include "../TOOLS/batchmanager.hpp"
#include "pde_stoch_stefan_boltzmann.hpp"
#include "obj_stoch_stefan_boltzmann.hpp"
#include "mesh_stoch_stefan_boltzmann.hpp"

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
    const int stochDim = 32, controlDim = 1;
    const RealT one(1); 

    /*************************************************************************/
    /***************** BUILD GOVERNING PDE ***********************************/
    /*************************************************************************/
    /*** Initialize main data structure. ***/
    Teuchos::RCP<MeshManager<RealT> > meshMgr
      = Teuchos::rcp(new MeshManager_BackwardFacingStepChannel<RealT>(*parlist));
    //  = Teuchos::rcp(new StochasticStefanBoltzmannMeshManager<RealT>(*parlist));
    // Initialize PDE describing advection-diffusion equation
    Teuchos::RCP<StochasticStefanBoltzmannPDE<RealT> > pde
      = Teuchos::rcp(new StochasticStefanBoltzmannPDE<RealT>(*parlist));
    Teuchos::RCP<ROL::EqualityConstraint_SimOpt<RealT> > con
      = Teuchos::rcp(new PDE_Constraint<RealT>(pde,meshMgr,serial_comm,*parlist,*outStream));
    // Cast the constraint and get the assembler.
    Teuchos::RCP<PDE_Constraint<RealT> > pdecon
      = Teuchos::rcp_dynamic_cast<PDE_Constraint<RealT> >(con);
    Teuchos::RCP<Assembler<RealT> > assembler = pdecon->getAssembler();
    con->setSolveParameters(*parlist);

    /*************************************************************************/
    /***************** BUILD VECTORS *****************************************/
    /*************************************************************************/
    Teuchos::RCP<Tpetra::MultiVector<> >  u_rcp = assembler->createStateVector();
    Teuchos::RCP<Tpetra::MultiVector<> >  p_rcp = assembler->createStateVector();
    Teuchos::RCP<Tpetra::MultiVector<> > du_rcp = assembler->createStateVector();
    u_rcp->randomize();  //u_rcp->putScalar(static_cast<RealT>(1));
    p_rcp->randomize();  //p_rcp->putScalar(static_cast<RealT>(1));
    du_rcp->randomize(); //du_rcp->putScalar(static_cast<RealT>(0));
    Teuchos::RCP<ROL::Vector<RealT> > up
      = Teuchos::rcp(new PDE_PrimalSimVector<RealT>(u_rcp,pde,assembler));
    Teuchos::RCP<ROL::Vector<RealT> > pp
      = Teuchos::rcp(new PDE_PrimalSimVector<RealT>(p_rcp,pde,assembler));
    Teuchos::RCP<ROL::Vector<RealT> > dup
      = Teuchos::rcp(new PDE_PrimalSimVector<RealT>(du_rcp,pde,assembler));
    // Create residual vectors
    Teuchos::RCP<Tpetra::MultiVector<> > r_rcp = assembler->createResidualVector();
    r_rcp->randomize(); //r_rcp->putScalar(static_cast<RealT>(1));
    Teuchos::RCP<ROL::Vector<RealT> > rp
      = Teuchos::rcp(new PDE_DualSimVector<RealT>(r_rcp,pde,assembler));
    // Create control vector --- FIELD
    Teuchos::RCP<Tpetra::MultiVector<> >  zbc_rcp
      = assembler->createControlVector();
    Teuchos::RCP<Tpetra::MultiVector<> > dzbc_rcp
      = assembler->createControlVector();
    Teuchos::RCP<Tpetra::MultiVector<> > yzbc_rcp
      = assembler->createControlVector();
    zbc_rcp->randomize();  zbc_rcp->putScalar(static_cast<RealT>(280));
    dzbc_rcp->randomize(); //dzbc_rcp->putScalar(static_cast<RealT>(0));
    yzbc_rcp->randomize(); //yzbc_rcp->putScalar(static_cast<RealT>(0));
    Teuchos::RCP<ROL::TpetraMultiVector<RealT> > zbc
      = Teuchos::rcp(new PDE_PrimalOptVector<RealT>(zbc_rcp,pde,assembler));
    Teuchos::RCP<ROL::TpetraMultiVector<RealT> > dzbc
      = Teuchos::rcp(new PDE_PrimalOptVector<RealT>(dzbc_rcp,pde,assembler));
    Teuchos::RCP<ROL::TpetraMultiVector<RealT> > yzbc
      = Teuchos::rcp(new PDE_PrimalOptVector<RealT>(yzbc_rcp,pde,assembler));
    // Create control vector --- PARAMETER
    Teuchos::RCP<std::vector<RealT> >  zp_rcp
      = Teuchos::rcp(new std::vector<RealT>(controlDim));
    Teuchos::RCP<std::vector<RealT> > dzp_rcp
      = Teuchos::rcp(new std::vector<RealT>(controlDim));
    Teuchos::RCP<std::vector<RealT> > yzp_rcp
      = Teuchos::rcp(new std::vector<RealT>(controlDim));
    for (int i = 0; i < controlDim; ++i) {
      (*zp_rcp)[i]  = random<RealT>(*comm);
      (*dzp_rcp)[i] = random<RealT>(*comm);
      (*yzp_rcp)[i] = random<RealT>(*comm);
    }
    Teuchos::RCP<ROL::StdVector<RealT> > zparam
      = Teuchos::rcp(new ROL::StdVector<RealT>(zp_rcp));
    Teuchos::RCP<ROL::StdVector<RealT> > dzparam
      = Teuchos::rcp(new ROL::StdVector<RealT>(dzp_rcp));
    Teuchos::RCP<ROL::StdVector<RealT> > yzparam
      = Teuchos::rcp(new ROL::StdVector<RealT>(yzp_rcp));
    // Create OptVectors
    Teuchos::RCP<ROL::Vector<RealT> > zp
      = Teuchos::rcp(new PDE_OptVector<RealT>(zbc, zparam));
    Teuchos::RCP<ROL::Vector<RealT> > dzp
      = Teuchos::rcp(new PDE_OptVector<RealT>(dzbc, dzparam));
    Teuchos::RCP<ROL::Vector<RealT> > yzp
      = Teuchos::rcp(new PDE_OptVector<RealT>(yzbc, yzparam));
    // Create ROL SimOpt vectors
    ROL::Vector_SimOpt<RealT> x(up,zp);
    ROL::Vector_SimOpt<RealT> d(dup,dzp);

    /*************************************************************************/
    /***************** BUILD COST FUNCTIONAL *********************************/
    /*************************************************************************/
    std::vector<Teuchos::RCP<QoI<RealT> > > qoi_vec(3,Teuchos::null);
    qoi_vec[0] = Teuchos::rcp(new QoI_StateCost<RealT>(pde->getVolFE(),*parlist));
    qoi_vec[1] = Teuchos::rcp(new QoI_ControlCost<RealT>(
      pde->getVolFE(),pde->getBdryFE(0),pde->getBdryCellLocIds(0),*parlist));
    qoi_vec[2] = Teuchos::rcp(new QoI_AdvectionCost<RealT>());
    Teuchos::RCP<StochasticStefanBoltzmannStdObjective3<RealT> > std_obj
      = Teuchos::rcp(new StochasticStefanBoltzmannStdObjective3<RealT>(*parlist));
    Teuchos::RCP<ROL::Objective_SimOpt<RealT> > obj
      = Teuchos::rcp(new PDE_Objective<RealT>(qoi_vec,std_obj,assembler));
    Teuchos::RCP<ROL::Reduced_Objective_SimOpt<RealT> > objReduced
      = Teuchos::rcp(new ROL::Reduced_Objective_SimOpt<RealT>(obj, con, up, pp, true, false));

    /*************************************************************************/
    /***************** BUILD BOUND CONSTRAINT ********************************/
    /*************************************************************************/
    // Bounds for boundary control
    RealT lower_bc = parlist->sublist("Problem").get("Lower Control Bound", 280.0);
    RealT upper_bc = parlist->sublist("Problem").get("Upper Control Bound", 370.0);
    Teuchos::RCP<Tpetra::MultiVector<> >  zlo_bc_rcp
      = assembler->createControlVector();
    Teuchos::RCP<Tpetra::MultiVector<> >  zhi_bc_rcp
      = assembler->createControlVector();
    zlo_bc_rcp->putScalar(static_cast<RealT>(lower_bc));
    zhi_bc_rcp->putScalar(static_cast<RealT>(upper_bc));
    Teuchos::RCP<ROL::TpetraMultiVector<RealT> > zlo_bc
      = Teuchos::rcp(new PDE_PrimalOptVector<RealT>(zlo_bc_rcp,pde,assembler));
    Teuchos::RCP<ROL::TpetraMultiVector<RealT> > zhi_bc
      = Teuchos::rcp(new PDE_PrimalOptVector<RealT>(zhi_bc_rcp,pde,assembler));
    // Bounds for advection control
    RealT lower = parlist->sublist("Problem").get("Lower Advection Bound",-100.0);
    RealT upper = parlist->sublist("Problem").get("Upper Advection Bound", 100.0);
    Teuchos::RCP<std::vector<RealT> > zlo_param_rcp
      = Teuchos::rcp(new std::vector<RealT>(controlDim,lower));
    Teuchos::RCP<std::vector<RealT> > zhi_param_rcp
      = Teuchos::rcp(new std::vector<RealT>(controlDim,upper));
    Teuchos::RCP<ROL::StdVector<RealT> > zlo_adv
      = Teuchos::rcp(new ROL::StdVector<RealT>(zlo_param_rcp));
    Teuchos::RCP<ROL::StdVector<RealT> > zhi_adv
      = Teuchos::rcp(new ROL::StdVector<RealT>(zhi_param_rcp));
    // Combined bounds
    Teuchos::RCP<ROL::Vector<RealT> > zlop
      = Teuchos::rcp(new PDE_OptVector<RealT>(zlo_bc, zlo_adv));
    Teuchos::RCP<ROL::Vector<RealT> > zhip
      = Teuchos::rcp(new PDE_OptVector<RealT>(zhi_bc, zhi_adv));
    Teuchos::RCP<ROL::BoundConstraint<RealT> > bnd
      = Teuchos::rcp(new ROL::BoundConstraint<RealT>(zlop,zhip));

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
    ROL::StochasticProblem<RealT> opt(*parlist,objReduced,sampler,zp,bnd);
    opt.setSolutionStatistic(one);

    /*************************************************************************/
    /***************** RUN VECTOR AND DERIVATIVE CHECKS **********************/
    /*************************************************************************/
    bool checkDeriv = parlist->sublist("Problem").get("Check Derivatives",false);
    if ( checkDeriv ) {
      up->checkVector(*pp,*dup,true,*outStream);
      zp->checkVector(*yzp,*dzp,true,*outStream);
      obj->checkGradient(x,d,true,*outStream);
      obj->checkHessVec(x,d,true,*outStream);
      con->checkApplyJacobian(x,d,*up,true,*outStream);
      con->checkApplyAdjointHessian(x,*dup,d,x,true,*outStream);
      con->checkAdjointConsistencyJacobian(*dup,d,x,true,*outStream);
      con->checkInverseJacobian_1(*up,*up,*up,*zp,true,*outStream);
      con->checkInverseAdjointJacobian_1(*up,*up,*up,*zp,true,*outStream);
      objReduced->checkGradient(*zp,*dzp,true,*outStream);
      objReduced->checkHessVec(*zp,*dzp,true,*outStream);
      opt.checkObjectiveGradient(*dzp,true,*outStream);
      opt.checkObjectiveHessVec(*dzp,true,*outStream);
    }

    /*************************************************************************/
    /***************** SOLVE OPTIMIZATION PROBLEM ****************************/
    /*************************************************************************/
    ROL::Algorithm<RealT> algo("Trust Region",*parlist,false);
    (*zp_rcp)[0] = parlist->sublist("Problem").get("Advection Magnitude",0.0);
    u_rcp->putScalar(450.0);
    std::clock_t timer = std::clock();
    algo.run(opt,true,*outStream);
    *outStream << "Optimization time: "
               << static_cast<RealT>(std::clock()-timer)/static_cast<RealT>(CLOCKS_PER_SEC)
               << " seconds." << std::endl << std::endl;

    /*************************************************************************/
    /***************** OUTPUT RESULTS ****************************************/
    /*************************************************************************/
    std::clock_t timer_print = std::clock();
    assembler->printMeshData(*outStream);
    // Output control to file
    pdecon->outputTpetraVector(zbc_rcp,"control.txt");
    *outStream << std::endl << "Advection value: " << (*zp_rcp)[0] << std::endl;
    // Output expected state and samples to file
    *outStream << std::endl << "Print Expected Value of State" << std::endl;
    up->zero(); pp->zero(); dup->zero();
    RealT tol(1.e-8);
    Teuchos::RCP<ROL::BatchManager<RealT> > bman_Eu
      = Teuchos::rcp(new ROL::TpetraTeuchosBatchManager<RealT>(comm));
    std::vector<RealT> sample(stochDim);
    std::stringstream name_samp;
    name_samp << "samples_" << bman->batchID() << ".txt";
    std::ofstream file_samp;
    file_samp.open(name_samp.str());
    file_samp << std::scientific << std::setprecision(15);
    for (int i = 0; i < sampler->numMySamples(); ++i) {
      *outStream << "Sample i = " << i << std::endl;
      sample = sampler->getMyPoint(i);
      con->setParameter(sample);
      con->solve(*rp,*dup,*zp,tol);
      up->axpy(sampler->getMyWeight(i),*dup);
      for (int j = 0; j < stochDim; ++j) {
        file_samp << std::setw(25) << std::left << sample[j];
      }
      file_samp << std::endl;
    }
    file_samp.close();
    bman_Eu->sumAll(*up,*pp);
    pdecon->outputTpetraVector(p_rcp,"mean_state.txt");
    // Build objective function distribution
    *outStream << std::endl << "Print Objective CDF" << std::endl;
    RealT val1(0), val2(0);
    int nsamp_dist = parlist->sublist("Problem").get("Number of Output Samples",100);
    Teuchos::RCP<ROL::Objective_SimOpt<RealT> > stateCost
      = Teuchos::rcp(new IntegralObjective<RealT>(qoi_vec[0],assembler));
    Teuchos::RCP<ROL::Reduced_Objective_SimOpt<RealT> > redStateCost
      = Teuchos::rcp(new ROL::Reduced_Objective_SimOpt<RealT>(stateCost, con, up, pp, true, false));
    Teuchos::RCP<ROL::Objective_SimOpt<RealT> > ctrlCost
      = Teuchos::rcp(new IntegralObjective<RealT>(qoi_vec[1],assembler));
    Teuchos::RCP<ROL::Reduced_Objective_SimOpt<RealT> > redCtrlCost
      = Teuchos::rcp(new ROL::Reduced_Objective_SimOpt<RealT>(ctrlCost, con, up, pp, true, false));
    Teuchos::RCP<ROL::SampleGenerator<RealT> > sampler_dist
      = Teuchos::rcp(new ROL::MonteCarloGenerator<RealT>(nsamp_dist,bounds,bman));
    std::stringstream name;
    name << "obj_samples_" << bman->batchID() << ".txt";
    std::ofstream file;
    file.open(name.str());
    file << std::scientific << std::setprecision(15);
    for (int i = 0; i < sampler_dist->numMySamples(); ++i) {
      sample = sampler_dist->getMyPoint(i);
      redStateCost->setParameter(sample);
      val1 = redStateCost->value(*zp,tol);
      redCtrlCost->setParameter(sample);
      val2 = redCtrlCost->value(*zp,tol);
      for (int j = 0; j < stochDim; ++j) {
        file << std::setw(25) << std::left << sample[j];
      }
      file << std::setw(25) << std::left << val1;
      file << std::setw(25) << std::left << val2 << std::endl;
    }
    file.close();
    *outStream << "Output time: "
               << static_cast<RealT>(std::clock()-timer_print)/static_cast<RealT>(CLOCKS_PER_SEC)
               << " seconds." << std::endl << std::endl;

    Teuchos::Array<RealT> res(1,0);
    con->value(*rp,*up,*zp,tol);
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
