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
#include "ROL_Stream.hpp"
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
  ROL::Ptr<std::ostream> outStream;
  ROL::nullstream bhs; // outputs nothing

  /*** Initialize communicator. ***/
  Teuchos::GlobalMPISession mpiSession (&argc, &argv, &bhs);
  ROL::Ptr<const Teuchos::Comm<int> > comm
    = Tpetra::DefaultPlatform::getDefaultPlatform().getComm();
  ROL::Ptr<const Teuchos::Comm<int> > serial_comm
    = ROL::makePtr<Teuchos::SerialComm<int>>();
  const int myRank = comm->getRank();
  if ((iprint > 0) && (myRank == 0)) {
    outStream = ROL::makePtrFromRef(std::cout);
  }
  else {
    outStream = ROL::makePtrFromRef(bhs);
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
    ROL::Ptr<MeshManager<RealT> > meshMgr
      = ROL::makePtr<MeshManager_BackwardFacingStepChannel<RealT>>(*parlist);
    //  = ROL::makePtr<StochasticStefanBoltzmannMeshManager<RealT>>(*parlist);
    // Initialize PDE describing advection-diffusion equation
    ROL::Ptr<StochasticStefanBoltzmannPDE<RealT> > pde
      = ROL::makePtr<StochasticStefanBoltzmannPDE<RealT>>(*parlist);
    ROL::Ptr<ROL::Constraint_SimOpt<RealT> > con
      = ROL::makePtr<PDE_Constraint<RealT>>(pde,meshMgr,serial_comm,*parlist,*outStream);
    // Cast the constraint and get the assembler.
    ROL::Ptr<PDE_Constraint<RealT> > pdecon
      = ROL::dynamicPtrCast<PDE_Constraint<RealT> >(con);
    ROL::Ptr<Assembler<RealT> > assembler = pdecon->getAssembler();
    con->setSolveParameters(*parlist);

    /*************************************************************************/
    /***************** BUILD VECTORS *****************************************/
    /*************************************************************************/
    ROL::Ptr<Tpetra::MultiVector<> >  u_ptr = assembler->createStateVector();
    ROL::Ptr<Tpetra::MultiVector<> >  p_ptr = assembler->createStateVector();
    ROL::Ptr<Tpetra::MultiVector<> > du_ptr = assembler->createStateVector();
    u_ptr->randomize();  //u_ptr->putScalar(static_cast<RealT>(1));
    p_ptr->randomize();  //p_ptr->putScalar(static_cast<RealT>(1));
    du_ptr->randomize(); //du_ptr->putScalar(static_cast<RealT>(0));
    ROL::Ptr<ROL::Vector<RealT> > up
      = ROL::makePtr<PDE_PrimalSimVector<RealT>>(u_ptr,pde,assembler);
    ROL::Ptr<ROL::Vector<RealT> > pp
      = ROL::makePtr<PDE_PrimalSimVector<RealT>>(p_ptr,pde,assembler);
    ROL::Ptr<ROL::Vector<RealT> > dup
      = ROL::makePtr<PDE_PrimalSimVector<RealT>>(du_ptr,pde,assembler);
    // Create residual vectors
    ROL::Ptr<Tpetra::MultiVector<> > r_ptr = assembler->createResidualVector();
    r_ptr->randomize(); //r_ptr->putScalar(static_cast<RealT>(1));
    ROL::Ptr<ROL::Vector<RealT> > rp
      = ROL::makePtr<PDE_DualSimVector<RealT>>(r_ptr,pde,assembler);
    // Create control vector --- FIELD
    ROL::Ptr<Tpetra::MultiVector<> >  zbc_ptr
      = assembler->createControlVector();
    ROL::Ptr<Tpetra::MultiVector<> > dzbc_ptr
      = assembler->createControlVector();
    ROL::Ptr<Tpetra::MultiVector<> > yzbc_ptr
      = assembler->createControlVector();
    zbc_ptr->randomize();  zbc_ptr->putScalar(static_cast<RealT>(280));
    dzbc_ptr->randomize(); //dzbc_ptr->putScalar(static_cast<RealT>(0));
    yzbc_ptr->randomize(); //yzbc_ptr->putScalar(static_cast<RealT>(0));
    ROL::Ptr<ROL::TpetraMultiVector<RealT> > zbc
      = ROL::makePtr<PDE_PrimalOptVector<RealT>>(zbc_ptr,pde,assembler);
    ROL::Ptr<ROL::TpetraMultiVector<RealT> > dzbc
      = ROL::makePtr<PDE_PrimalOptVector<RealT>>(dzbc_ptr,pde,assembler);
    ROL::Ptr<ROL::TpetraMultiVector<RealT> > yzbc
      = ROL::makePtr<PDE_PrimalOptVector<RealT>>(yzbc_ptr,pde,assembler);
    // Create control vector --- PARAMETER
    ROL::Ptr<std::vector<RealT> >  zp_ptr
      = ROL::makePtr<std::vector<RealT>>(controlDim);
    ROL::Ptr<std::vector<RealT> > dzp_ptr
      = ROL::makePtr<std::vector<RealT>>(controlDim);
    ROL::Ptr<std::vector<RealT> > yzp_ptr
      = ROL::makePtr<std::vector<RealT>>(controlDim);
    for (int i = 0; i < controlDim; ++i) {
      (*zp_ptr)[i]  = random<RealT>(*comm);
      (*dzp_ptr)[i] = random<RealT>(*comm);
      (*yzp_ptr)[i] = random<RealT>(*comm);
    }
    ROL::Ptr<ROL::StdVector<RealT> > zparam
      = ROL::makePtr<ROL::StdVector<RealT>>(zp_ptr);
    ROL::Ptr<ROL::StdVector<RealT> > dzparam
      = ROL::makePtr<ROL::StdVector<RealT>>(dzp_ptr);
    ROL::Ptr<ROL::StdVector<RealT> > yzparam
      = ROL::makePtr<ROL::StdVector<RealT>>(yzp_ptr);
    // Create OptVectors
    ROL::Ptr<ROL::Vector<RealT> > zp
      = ROL::makePtr<PDE_OptVector<RealT>>(zbc, zparam);
    ROL::Ptr<ROL::Vector<RealT> > dzp
      = ROL::makePtr<PDE_OptVector<RealT>>(dzbc, dzparam);
    ROL::Ptr<ROL::Vector<RealT> > yzp
      = ROL::makePtr<PDE_OptVector<RealT>>(yzbc, yzparam);
    // Create ROL SimOpt vectors
    ROL::Vector_SimOpt<RealT> x(up,zp);
    ROL::Vector_SimOpt<RealT> d(dup,dzp);

    /*************************************************************************/
    /***************** BUILD COST FUNCTIONAL *********************************/
    /*************************************************************************/
    std::vector<ROL::Ptr<QoI<RealT> > > qoi_vec(3,ROL::nullPtr);
    qoi_vec[0] = ROL::makePtr<QoI_StateCost<RealT>>(pde->getVolFE(),*parlist);
    qoi_vec[1] = ROL::makePtr<QoI_ControlCost<RealT>>(
      pde->getVolFE(),pde->getBdryFE(0),pde->getBdryCellLocIds(0),*parlist);
    qoi_vec[2] = ROL::makePtr<QoI_AdvectionCost<RealT>>();
    ROL::Ptr<StochasticStefanBoltzmannStdObjective3<RealT> > std_obj
      = ROL::makePtr<StochasticStefanBoltzmannStdObjective3<RealT>>(*parlist);
    ROL::Ptr<ROL::Objective_SimOpt<RealT> > obj
      = ROL::makePtr<PDE_Objective<RealT>>(qoi_vec,std_obj,assembler);
    ROL::Ptr<ROL::Reduced_Objective_SimOpt<RealT> > objReduced
      = ROL::makePtr<ROL::Reduced_Objective_SimOpt<RealT>>(obj, con, up, zp, pp, true, false);

    /*************************************************************************/
    /***************** BUILD BOUND CONSTRAINT ********************************/
    /*************************************************************************/
    // Bounds for boundary control
    RealT lower_bc = parlist->sublist("Problem").get("Lower Control Bound", 280.0);
    RealT upper_bc = parlist->sublist("Problem").get("Upper Control Bound", 370.0);
    ROL::Ptr<Tpetra::MultiVector<> >  zlo_bc_ptr
      = assembler->createControlVector();
    ROL::Ptr<Tpetra::MultiVector<> >  zhi_bc_ptr
      = assembler->createControlVector();
    zlo_bc_ptr->putScalar(static_cast<RealT>(lower_bc));
    zhi_bc_ptr->putScalar(static_cast<RealT>(upper_bc));
    ROL::Ptr<ROL::TpetraMultiVector<RealT> > zlo_bc
      = ROL::makePtr<PDE_PrimalOptVector<RealT>>(zlo_bc_ptr,pde,assembler);
    ROL::Ptr<ROL::TpetraMultiVector<RealT> > zhi_bc
      = ROL::makePtr<PDE_PrimalOptVector<RealT>>(zhi_bc_ptr,pde,assembler);
    // Bounds for advection control
    RealT lower = parlist->sublist("Problem").get("Lower Advection Bound",-100.0);
    RealT upper = parlist->sublist("Problem").get("Upper Advection Bound", 100.0);
    ROL::Ptr<std::vector<RealT> > zlo_param_ptr
      = ROL::makePtr<std::vector<RealT>>(controlDim,lower);
    ROL::Ptr<std::vector<RealT> > zhi_param_ptr
      = ROL::makePtr<std::vector<RealT>>(controlDim,upper);
    ROL::Ptr<ROL::StdVector<RealT> > zlo_adv
      = ROL::makePtr<ROL::StdVector<RealT>>(zlo_param_ptr);
    ROL::Ptr<ROL::StdVector<RealT> > zhi_adv
      = ROL::makePtr<ROL::StdVector<RealT>>(zhi_param_ptr);
    // Combined bounds
    ROL::Ptr<ROL::Vector<RealT> > zlop
      = ROL::makePtr<PDE_OptVector<RealT>>(zlo_bc, zlo_adv);
    ROL::Ptr<ROL::Vector<RealT> > zhip
      = ROL::makePtr<PDE_OptVector<RealT>>(zhi_bc, zhi_adv);
    ROL::Ptr<ROL::BoundConstraint<RealT> > bnd
      = ROL::makePtr<ROL::Bounds<RealT>>(zlop,zhip);

    /*************************************************************************/
    /***************** BUILD SAMPLER *****************************************/
    /*************************************************************************/
    int nsamp = parlist->sublist("Problem").get("Number of Samples",100);
    std::vector<RealT> tmp = {-one,one};
    std::vector<std::vector<RealT> > bounds(stochDim,tmp);
    ROL::Ptr<ROL::BatchManager<RealT> > bman
      = ROL::makePtr<PDE_OptVector_BatchManager<RealT>>(comm);
    ROL::Ptr<ROL::SampleGenerator<RealT> > sampler
      = ROL::makePtr<ROL::MonteCarloGenerator<RealT>>(nsamp,bounds,bman);

    /*************************************************************************/
    /***************** BUILD STOCHASTIC PROBLEM ******************************/
    /*************************************************************************/
    ROL::OptimizationProblem<RealT> opt(objReduced,zp,bnd);
    parlist->sublist("SOL").set("Initial Statistic", one);
    opt.setStochasticObjective(*parlist,sampler);

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
      opt.check(*outStream);
    }

    /*************************************************************************/
    /***************** SOLVE OPTIMIZATION PROBLEM ****************************/
    /*************************************************************************/
    ROL::Algorithm<RealT> algo("Trust Region",*parlist,false);
    (*zp_ptr)[0] = parlist->sublist("Problem").get("Advection Magnitude",0.0);
    u_ptr->putScalar(450.0);
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
    pdecon->outputTpetraVector(zbc_ptr,"control.txt");
    *outStream << std::endl << "Advection value: " << (*zp_ptr)[0] << std::endl;
    // Output expected state and samples to file
    *outStream << std::endl << "Print Expected Value of State" << std::endl;
    up->zero(); pp->zero(); dup->zero();
    RealT tol(1.e-8);
    ROL::Ptr<ROL::BatchManager<RealT> > bman_Eu
      = ROL::makePtr<ROL::TpetraTeuchosBatchManager<RealT>>(comm);
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
    pdecon->outputTpetraVector(p_ptr,"mean_state.txt");
    // Build objective function distribution
    *outStream << std::endl << "Print Objective CDF" << std::endl;
    RealT val1(0), val2(0);
    int nsamp_dist = parlist->sublist("Problem").get("Number of Output Samples",100);
    ROL::Ptr<ROL::Objective_SimOpt<RealT> > stateCost
      = ROL::makePtr<IntegralObjective<RealT>>(qoi_vec[0],assembler);
    ROL::Ptr<ROL::Reduced_Objective_SimOpt<RealT> > redStateCost
      = ROL::makePtr<ROL::Reduced_Objective_SimOpt<RealT>>(stateCost, con, up, zp, pp, true, false);
    ROL::Ptr<ROL::Objective_SimOpt<RealT> > ctrlCost
      = ROL::makePtr<IntegralObjective<RealT>>(qoi_vec[1],assembler);
    ROL::Ptr<ROL::Reduced_Objective_SimOpt<RealT> > redCtrlCost
      = ROL::makePtr<ROL::Reduced_Objective_SimOpt<RealT>>(ctrlCost, con, up, zp, pp, true, false);
    ROL::Ptr<ROL::SampleGenerator<RealT> > sampler_dist
      = ROL::makePtr<ROL::MonteCarloGenerator<RealT>>(nsamp_dist,bounds,bman);
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
    con->solve(*rp,*up,*zp,tol);
    r_ptr->norm2(res.view(0,1));

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
