// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  example_01.cpp
    \brief Shows how to solve the stochastic Stefan-Boltzmann problem.
*/

#include "Teuchos_Comm.hpp"
#include "ROL_Stream.hpp"
#include "ROL_GlobalMPISession.hpp"

#include "Tpetra_Core.hpp"
#include "Tpetra_Version.hpp"

#include <iostream>
#include <algorithm>
//#include <fenv.h>

#include "ROL_Bounds.hpp"
#include "ROL_Reduced_Objective_SimOpt.hpp"
#include "ROL_MonteCarloGenerator.hpp"
#include "ROL_StochasticProblem.hpp"
#include "ROL_Solver.hpp"
#include "ROL_TpetraTeuchosBatchManager.hpp"

#include "../TOOLS/meshmanagerK.hpp"
#include "../TOOLS/pdeconstraintK.hpp"
#include "../TOOLS/pdeobjectiveK.hpp"
#include "../TOOLS/pdevectorK.hpp"
#include "../TOOLS/batchmanagerK.hpp"
#include "pde_navier-stokesK.hpp"
#include "obj_navier-stokesK.hpp"

using RealT = double;
using DeviceT = Kokkos::HostSpace;

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
  //feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  ROL::Ptr<std::ostream> outStream;
  ROL::nullstream bhs; // outputs nothing

  /*** Initialize communicator. ***/
  ROL::GlobalMPISession mpiSession (&argc, &argv, &bhs);
  Kokkos::ScopeGuard kokkosScope (argc, argv);
  ROL::Ptr<const Teuchos::Comm<int> > comm
    = Tpetra::getDefaultComm();
  ROL::Ptr<const Teuchos::Comm<int> > serial_comm
    = ROL::makePtr<Teuchos::SerialComm<int>>();
  const int myRank = comm->getRank();
  if ((iprint > 0) && (myRank == 0))
    outStream = ROL::makePtrFromRef(std::cout);
  else
    outStream = ROL::makePtrFromRef(bhs);
  int errorFlag  = 0;

  // *** Example body.
  try {

    /*** Read in XML input ***/
    std::string filename = "input.xml";
    auto parlist = ROL::getParametersFromXmlFile(filename);

    // Problem dimensions
    const int stochDim = 2;
    const RealT one(1); 

    /*************************************************************************/
    /***************** BUILD GOVERNING PDE ***********************************/
    /*************************************************************************/
    /*** Initialize main data structure. ***/
    auto meshMgr = ROL::makePtr<MeshManager_BackwardFacingStepChannel<RealT,DeviceT>>(*parlist);
    // Initialize PDE describing advection-diffusion equation
    auto pde = ROL::makePtr<PDE_NavierStokes<RealT,DeviceT>>(*parlist);
    auto con = ROL::makePtr<PDE_Constraint<RealT,DeviceT>>(pde,meshMgr,serial_comm,*parlist,*outStream);
    // Cast the constraint and get the assembler.
    auto assembler = con->getAssembler();
    con->setSolveParameters(*parlist);

    /*************************************************************************/
    /***************** BUILD VECTORS *****************************************/
    /*************************************************************************/
    auto  u_ptr = assembler->createStateVector();    u_ptr->randomize();
    auto du_ptr = assembler->createStateVector();    du_ptr->randomize();
    auto  p_ptr = assembler->createStateVector();    p_ptr->randomize();
    auto  r_ptr = assembler->createResidualVector(); r_ptr->randomize();
    auto  z_ptr = assembler->createControlVector();  z_ptr->randomize();
    auto up  = ROL::makePtr<PDE_PrimalSimVector<RealT,DeviceT>>(u_ptr,pde,assembler,*parlist);
    auto dup = ROL::makePtr<PDE_PrimalSimVector<RealT,DeviceT>>(u_ptr,pde,assembler,*parlist);
    auto pp  = ROL::makePtr<PDE_PrimalSimVector<RealT,DeviceT>>(p_ptr,pde,assembler,*parlist);
    auto rp  = ROL::makePtr<PDE_DualSimVector<RealT,DeviceT>>(r_ptr,pde,assembler,*parlist);
    ROL::Ptr<ROL::TpetraMultiVector<RealT>> zpde
      = ROL::makePtr<PDE_PrimalOptVector<RealT,DeviceT>>(z_ptr,pde,assembler,*parlist);
    auto zp = ROL::makePtr<PDE_OptVector<RealT>>(zpde);
    // Create ROL SimOpt vectors
    ROL::Vector_SimOpt<RealT> x(up,zp);

    /*************************************************************************/
    /***************** BUILD COST FUNCTIONAL *********************************/
    /*************************************************************************/
    std::vector<ROL::Ptr<QoI<RealT,DeviceT>>> qoi_vec(2,ROL::nullPtr);
    qoi_vec[0] = ROL::makePtr<QoI_State_NavierStokes<RealT,DeviceT>>(*parlist,
                                                             pde->getVelocityFE(),
                                                             pde->getPressureFE(),
                                                             pde->getFieldHelper());
    qoi_vec[1] = ROL::makePtr<QoI_L2Penalty_NavierStokes<RealT,DeviceT>>(pde->getVelocityFE(),
                                                                 pde->getPressureFE(),
                                                                 pde->getVelocityBdryFE(),
                                                                 pde->getBdryCellLocIds(),
                                                                 pde->getFieldHelper());
    auto std_obj = ROL::makePtr<StdObjective_NavierStokes<RealT>>(*parlist);
    auto obj = ROL::makePtr<PDE_Objective<RealT,DeviceT>>(qoi_vec,std_obj,assembler);
    auto objReduced = ROL::makePtr<ROL::Reduced_Objective_SimOpt<RealT>>(obj, con, up, zp, pp, true, false);

    /*************************************************************************/
    /***************** BUILD BOUND CONSTRAINT ********************************/
    /*************************************************************************/
    auto zlo_ptr = assembler->createControlVector();
    auto zhi_ptr = assembler->createControlVector();
    zlo_ptr->putScalar(static_cast<RealT>(0));
    zhi_ptr->putScalar(ROL::ROL_INF<RealT>());
    ROL::Ptr<ROL::TpetraMultiVector<RealT>> zlopde, zhipde;
    zlopde = ROL::makePtr<PDE_PrimalOptVector<RealT,DeviceT>>(zlo_ptr,pde,assembler,*parlist);
    zhipde = ROL::makePtr<PDE_PrimalOptVector<RealT,DeviceT>>(zhi_ptr,pde,assembler,*parlist);
    auto zlop = ROL::makePtr<PDE_OptVector<RealT>>(zlopde);
    auto zhip = ROL::makePtr<PDE_OptVector<RealT>>(zhipde);
    auto bnd = ROL::makePtr<ROL::Bounds<RealT>>(zlop,zhip);
    bool useBounds = parlist->sublist("Problem").get("Use bounds", false);
    if (!useBounds) bnd->deactivate();

    /*************************************************************************/
    /***************** BUILD SAMPLER *****************************************/
    /*************************************************************************/
    int nsamp = parlist->sublist("Problem").get("Number of samples",100);
    std::vector<RealT> tmp = {-one,one};
    std::vector<std::vector<RealT>> bounds(stochDim,tmp);
    auto bman = ROL::makePtr<PDE_OptVector_BatchManager<RealT,DeviceT>>(comm);
    auto sampler = ROL::makePtr<ROL::MonteCarloGenerator<RealT>>(nsamp,bounds,bman);

    /*************************************************************************/
    /***************** BUILD STOCHASTIC PROBLEM ******************************/
    /*************************************************************************/
    auto opt = ROL::makePtr<ROL::StochasticProblem<RealT>>(objReduced,zp);
    opt->addBoundConstraint(bnd);
    parlist->sublist("SOL").set("Initial Statistic",one);
    opt->makeObjectiveStochastic(*parlist,sampler);
    opt->finalize(false,true,*outStream);

    /*************************************************************************/
    /***************** RUN VECTOR AND DERIVATIVE CHECKS **********************/
    /*************************************************************************/
    bool checkDeriv = parlist->sublist("Problem").get("Check derivatives",false);
    if ( checkDeriv ) {
      auto dzp = zp->clone(); dzp->randomize(-one,one);
      auto yzp = zp->clone(); yzp->randomize(-one,one);
      ROL::Vector_SimOpt<RealT> d(dup,dzp);
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
      opt->check(true,*outStream);
    }

    /*************************************************************************/
    /***************** SOLVE OPTIMIZATION PROBLEM ****************************/
    /*************************************************************************/
    parlist->sublist("Step").set("Type","Trust Region");
    ROL::Solver<RealT> solver(opt,*parlist);
    up->zero(); zp->zero();
    std::clock_t timer = std::clock();
    solver.solve(*outStream);
    *outStream << "Optimization time: "
               << static_cast<RealT>(std::clock()-timer)/static_cast<RealT>(CLOCKS_PER_SEC)
               << " seconds." << std::endl << std::endl;

    /*************************************************************************/
    /***************** OUTPUT RESULTS ****************************************/
    /*************************************************************************/
    std::clock_t timer_print = std::clock();
    assembler->printMeshData(*outStream);
    // Output control to file
    con->outputTpetraVector(z_ptr,"control.txt");
    // Output expected state and samples to file
    *outStream << std::endl << "Print Expected Value of State" << std::endl;
    up->zero(); pp->zero(); dup->zero();
    RealT tol(1.e-8);
    auto bman_Eu = ROL::makePtr<ROL::TpetraTeuchosBatchManager<RealT>>(comm);
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
      for (int j = 0; j < stochDim; ++j)
        file_samp << std::setw(25) << std::left << sample[j];
      file_samp << std::endl;
    }
    file_samp.close();
    bman_Eu->sumAll(*up,*pp);
    con->outputTpetraVector(p_ptr,"mean_state.txt");
    // Output extreme cases
    std::vector<RealT> param(stochDim,0);
    param[0] = -one; param[1] = -one;
    con->setParameter(param);
    con->solve(*rp,*dup,*zp,tol);
    con->outputTpetraVector(du_ptr,"state_m1_m1.txt");
    param[0] =  one; param[1] = -one;
    con->setParameter(param);
    con->solve(*rp,*dup,*zp,tol);
    con->outputTpetraVector(du_ptr,"state_p1_m1.txt");
    param[0] = -one; param[1] =  one;
    con->setParameter(param);
    con->solve(*rp,*dup,*zp,tol);
    con->outputTpetraVector(du_ptr,"state_m1_p1.txt");
    param[0] =  one; param[1] =  one;
    con->setParameter(param);
    con->solve(*rp,*dup,*zp,tol);
    con->outputTpetraVector(du_ptr,"state_p1_p1.txt");
    // Build objective function distribution
    *outStream << std::endl << "Print Objective CDF" << std::endl;
    RealT val1(0), val2(0);
    int nsamp_dist = parlist->sublist("Problem").get("Number of output samples",100);
    auto stateCost = ROL::makePtr<IntegralObjective<RealT,DeviceT>>(qoi_vec[0],assembler);
    auto redStateCost = ROL::makePtr<ROL::Reduced_Objective_SimOpt<RealT>>(stateCost, con, up, zp, pp, true, false);
    auto ctrlCost = ROL::makePtr<IntegralObjective<RealT,DeviceT>>(qoi_vec[1],assembler);
    auto redCtrlCost = ROL::makePtr<ROL::Reduced_Objective_SimOpt<RealT>>(ctrlCost, con, up, zp, pp, true, false);
    auto sampler_dist = ROL::makePtr<ROL::MonteCarloGenerator<RealT>>(nsamp_dist,bounds,bman);
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
      for (int j = 0; j < stochDim; ++j)
        file << std::setw(25) << std::left << sample[j];
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
  catch (std::logic_error& err) {
    *outStream << err.what() << "\n";
    errorFlag = -1000;
  }; // end try

  if (errorFlag != 0)
    std::cout << "End Result: TEST FAILED\n";
  else
    std::cout << "End Result: TEST PASSED\n";

  return 0;
}
