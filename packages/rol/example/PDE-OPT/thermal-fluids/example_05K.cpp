// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  example_01.cpp
    \brief Shows how to solve the Navier-Stokes control problem.
*/

#include "Teuchos_Comm.hpp"
#include "ROL_Stream.hpp"
#include "ROL_GlobalMPISession.hpp"

#include "Tpetra_Core.hpp"
#include "Tpetra_Version.hpp"

#include <iostream>
#include <algorithm>
//#include <fenv.h>

#include "ROL_TpetraMultiVector.hpp"
#include "ROL_Solver.hpp"
#include "ROL_Reduced_Objective_SimOpt.hpp"
#include "ROL_MonteCarloGenerator.hpp"
#include "ROL_SimulatedConstraint.hpp"
#include "ROL_SimulatedObjectiveCVaR.hpp"
#include "ROL_SimulatedObjective.hpp"
#include "ROL_TpetraTeuchosBatchManager.hpp"

#include "../TOOLS/pdeconstraintK.hpp"
#include "../TOOLS/pdeobjectiveK.hpp"
#include "../TOOLS/pdevectorK.hpp"
#include "pde_thermal-fluids_ex03K.hpp"
#include "obj_thermal-fluids_ex03K.hpp"
#include "mesh_thermal-fluidsK.hpp"

using RealT = double;
using DeviceT = Kokkos::HostSpace;

template<class Real>
void print(ROL::Objective<Real> &obj,
           const ROL::Vector<Real> &z,
           ROL::SampleGenerator<Real> &sampler,
           int ngsamp,
           const ROL::Ptr<const Teuchos::Comm<int>> &comm,
           std::string filename) {
  Real tol(1e-8);
  // Build objective function distribution
  int nsamp = sampler.numMySamples();
  std::vector<Real> myvalues(nsamp), myzerovec(nsamp, 0);
  std::vector<double> gvalues(ngsamp), gzerovec(ngsamp, 0);
  std::vector<Real> sample = sampler.getMyPoint(0);
  int sdim = sample.size();
  std::vector<std::vector<Real> > mysamples(sdim, myzerovec);
  std::vector<std::vector<double> > gsamples(sdim, gzerovec);
  for (int i = 0; i < nsamp; ++i) {
    sample = sampler.getMyPoint(i);
    obj.setParameter(sample);
    myvalues[i] = static_cast<double>(obj.value(z,tol));
    for (int j = 0; j < sdim; ++j)
      mysamples[j][i] = static_cast<double>(sample[j]);
  }

  // Send data to root processor
#ifdef HAVE_MPI
  auto mpicomm = ROL::dynamicPtrCast<const Teuchos::MpiComm<int>>(comm);
  int nproc = Teuchos::size<int>(*mpicomm);
  std::vector<int> sampleCounts(nproc, 0), sampleDispls(nproc, 0);
  MPI_Gather(&nsamp,1,MPI_INT,&sampleCounts[0],1,MPI_INT,0,*(mpicomm->getRawMpiComm())());
  for (int i = 1; i < nproc; ++i)
    sampleDispls[i] = sampleDispls[i-1] + sampleCounts[i-1];
  MPI_Gatherv(&myvalues[0],nsamp,MPI_DOUBLE,&gvalues[0],&sampleCounts[0],&sampleDispls[0],MPI_DOUBLE,0,*(mpicomm->getRawMpiComm())());
  for (int j = 0; j < sdim; ++j)
    MPI_Gatherv(&mysamples[j][0],nsamp,MPI_DOUBLE,&gsamples[j][0],&sampleCounts[0],&sampleDispls[0],MPI_DOUBLE,0,*(mpicomm->getRawMpiComm())());
#else
  gvalues.assign(myvalues.begin(),myvalues.end());
  for (int j = 0; j < sdim; ++j)
    gsamples[j].assign(mysamples[j].begin(),mysamples[j].end());
#endif

  // Print
  int rank  = Teuchos::rank<int>(*comm);
  if ( rank==0 ) {
    std::ofstream file;
    file.open(filename);
    file << std::scientific << std::setprecision(15);
    for (int i = 0; i < ngsamp; ++i) {
      for (int j = 0; j < sdim; ++j)
        file << std::setw(25) << std::left << gsamples[j][i];
      file << std::setw(25) << std::left << gvalues[i] << std::endl;
    }
    file.close();
  }
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
  auto comm = Tpetra::getDefaultComm();
  auto serial_comm = ROL::makePtr<Teuchos::SerialComm<int>>();
  const int myRank = comm->getRank();
  if ((iprint > 0) && (myRank == 0))
    outStream = ROL::makePtrFromRef(std::cout);
  else
    outStream = ROL::makePtrFromRef(bhs);
  int errorFlag  = 0;

  // *** Example body.
  try {

    /*** Read in XML input ***/
    std::string filename = "input_ex05.xml";
    auto parlist = ROL::getParametersFromXmlFile(filename);
    parlist->sublist("SimOpt").sublist("Solve").set("Output Iteration History",myRank==0);

    /*** Initialize main data structure. ***/
    auto meshMgr = ROL::makePtr<MeshManager_ThermalFluids<RealT,DeviceT>>(*parlist);
    // Initialize PDE describing Navier-Stokes equations.
    auto pde = ROL::makePtr<PDE_ThermalFluids_ex03<RealT,DeviceT>>(*parlist);
    auto con = ROL::makePtr<PDE_Constraint<RealT,DeviceT>>(pde,meshMgr,serial_comm,*parlist,*outStream);
    // Cast the constraint and get the assembler.
    auto assembler = con->getAssembler();
    con->setSolveParameters(*parlist);
    con->outputTpetraData();

    // Create state vector and set to zeroes
    auto u_ptr  = assembler->createStateVector();     u_ptr->randomize();
    auto p_ptr  = assembler->createStateVector();     p_ptr->randomize();
    auto y_ptr  = assembler->createStateVector();     y_ptr->randomize();
    auto r_ptr  = assembler->createResidualVector();  r_ptr->randomize();
    auto z_ptr  = assembler->createControlVector();   z_ptr->putScalar(1.234); //z_ptr->randomize();
    auto s_ptr  = assembler->createControlVector();   s_ptr->putScalar(2.345); //s_ptr->randomize();
    auto t_ptr  = assembler->createControlVector();   t_ptr->putScalar(3.456); //t_ptr->randomize();
    auto up  = ROL::makePtr<PDE_PrimalSimVector<RealT,DeviceT>>(u_ptr,pde,assembler);
    auto pp  = ROL::makePtr<PDE_PrimalSimVector<RealT,DeviceT>>(p_ptr,pde,assembler);
    auto yp  = ROL::makePtr<PDE_PrimalSimVector<RealT,DeviceT>>(y_ptr,pde,assembler);
    auto rp  = ROL::makePtr<PDE_DualSimVector<RealT,DeviceT>>(r_ptr,pde,assembler);
    auto zp  = ROL::makePtr<PDE_PrimalOptVector<RealT,DeviceT>>(z_ptr,pde,assembler);
    auto sp  = ROL::makePtr<PDE_PrimalOptVector<RealT,DeviceT>>(s_ptr,pde,assembler);
    auto tp  = ROL::makePtr<PDE_PrimalOptVector<RealT,DeviceT>>(t_ptr,pde,assembler);

    // Initialize objective function.
    std::vector<ROL::Ptr<QoI<RealT,DeviceT> > > qoi_vec(2,ROL::nullPtr);
    qoi_vec[0] = ROL::makePtr<QoI_State_ThermalFluids<RealT,DeviceT>>(*parlist,
                                                                 pde->getVelocityFE(),
                                                                 pde->getPressureFE(),
                                                                 pde->getThermalFE(),
                                                                 pde->getFieldHelper());
    qoi_vec[1] = ROL::makePtr<QoI_L2Penalty_ThermalFluids<RealT,DeviceT>>(pde->getVelocityFE(),
                                                                     pde->getPressureFE(),
                                                                     pde->getThermalFE(),
                                                                     pde->getThermalBdryFE(),
                                                                     pde->getBdryCellLocIds(),
                                                                     pde->getFieldHelper());
    auto std_obj = ROL::makePtr<StdObjective_ThermalFluids<RealT>>(*parlist);
    auto obj = ROL::makePtr<PDE_Objective<RealT,DeviceT>>(qoi_vec,std_obj,assembler);
    auto stateStore = ROL::makePtr<ROL::VectorController<RealT>>();
    auto robj = ROL::makePtr<ROL::Reduced_Objective_SimOpt<RealT>>(obj, con, stateStore, up, zp, pp, true, false);

    /*************************************************************************/
    /***************** BUILD SAMPLER *****************************************/
    /*************************************************************************/
    int Nbottom = parlist->sublist("Problem").get("Bottom KL Truncation Order",5);
    int Nleft   = parlist->sublist("Problem").get("Left KL Truncation Order",5);
    int Nright  = parlist->sublist("Problem").get("Right KL Truncation Order",5);
    int stochDim = Nbottom + Nleft + Nright + 3;
    int nsamp = parlist->sublist("Problem").get("Number of samples",100);
    // Build vector of distributions
    std::vector<ROL::Ptr<ROL::Distribution<RealT> > > distVec(stochDim);
    Teuchos::ParameterList UList;
    UList.sublist("Distribution").set("Name","Uniform");
    UList.sublist("Distribution").sublist("Uniform").set("Lower Bound",-1.0);
    UList.sublist("Distribution").sublist("Uniform").set("Upper Bound", 1.0);
    for (int i = 0; i < stochDim; ++i)
      distVec[i] = ROL::DistributionFactory<RealT>(UList);
    // Sampler
    auto bman = ROL::makePtr<ROL::TpetraTeuchosBatchManager<RealT>>(comm);
    auto sampler = ROL::makePtr<ROL::MonteCarloGenerator<RealT>>(nsamp,distVec,bman);

    /*************************************************************************/
    /***************** BUILD STOCHASTIC PROBLEM ******************************/
    /*************************************************************************/
    bool useW    = parlist->sublist("SOL").sublist("Simulated").get("Use Constraint Weights", true);
    bool useCVaR = parlist->sublist("SOL").sublist("Simulated").get("Use CVaR", false);
    auto simcon  = ROL::makePtr<ROL::SimulatedConstraint<RealT>>(sampler, con, useW);
    ROL::Ptr<ROL::Objective<RealT>> simobj;
    if (useCVaR) {
      ROL::ParameterList list = parlist->sublist("SOL").sublist("Simulated");
      auto pf = ROL::makePtr<ROL::PlusFunction<RealT>>(list);
      RealT alpha = parlist->sublist("SOL").sublist("Simulated").get("CVaR Confidence Level", 0.9);
      simobj = ROL::makePtr<ROL::SimulatedObjectiveCVaR<RealT>>(sampler, obj, pf, alpha);
    }
    else {
      simobj = ROL::makePtr<ROL::SimulatedObjective<RealT>>(sampler, obj);
    }
    std::vector<ROL::Ptr<ROL::Vector<RealT>>> vuvec, vpvec, vyvec;
    for (int i = 0; i < sampler->numMySamples(); ++i) {
      auto vu_ptr  = assembler->createStateVector(); vu_ptr->putScalar(4.567); //vu_ptr->randomize();
      auto vp_ptr  = assembler->createStateVector(); vp_ptr->putScalar(5.678); //vp_ptr->randomize();
      auto vy_ptr  = assembler->createStateVector(); vy_ptr->putScalar(6.789); //vy_ptr->randomize();
      auto vup  = ROL::makePtr<PDE_PrimalSimVector<RealT,DeviceT>>(vu_ptr,pde,assembler);
      auto vpp  = ROL::makePtr<PDE_PrimalSimVector<RealT,DeviceT>>(vp_ptr,pde,assembler);
      auto vyp  = ROL::makePtr<PDE_PrimalSimVector<RealT,DeviceT>>(vy_ptr,pde,assembler);
      vuvec.push_back(vup);
      vpvec.push_back(vpp);
      vyvec.push_back(vyp);
    }
    auto vu = ROL::makePtr<ROL::SimulatedVector<RealT>>(vuvec,bman);
    auto vp = ROL::makePtr<ROL::SimulatedVector<RealT>>(vpvec,bman);
    auto vy = ROL::makePtr<ROL::SimulatedVector<RealT>>(vyvec,bman);
    ROL::Ptr<ROL::Vector<RealT>> rz, rs, rt;
    if (useCVaR) {
      Teuchos::RCP<Teuchos::ParameterList> cvarlist = Teuchos::rcp( new Teuchos::ParameterList() );
      cvarlist->sublist("SOL").sublist("Risk Measure").set("Name","CVaR");
      rz = ROL::makePtr<ROL::RiskVector<RealT>>(cvarlist, zp);
      rs = ROL::makePtr<ROL::RiskVector<RealT>>(cvarlist, sp);
      rt = ROL::makePtr<ROL::RiskVector<RealT>>(cvarlist, tp);
    }
    else {
      rz = zp;
      rs = sp;
      rt = tp;
    }
    auto x = ROL::makePtr<ROL::Vector_SimOpt<RealT>>(vu,rz);
    auto p = ROL::makePtr<ROL::Vector_SimOpt<RealT>>(vp,rs);
    auto y = ROL::makePtr<ROL::Vector_SimOpt<RealT>>(vy,rt);
    x->checkVector(*p,*y,true,*outStream);

    bool derivCheck = parlist->sublist("Problem").get("Check derivatives",false);
    if (derivCheck) {
      *outStream << std::endl << "TESTING SimulatedConstraint" << std::endl;
      simcon->checkApplyJacobian(*x, *p, *vu, true, *outStream);
      simcon->checkAdjointConsistencyJacobian(*vu, *p, *x, true, *outStream);
      simcon->checkApplyAdjointHessian(*x, *vu, *p, *x, true, *outStream);
      *outStream << std::endl << "TESTING SimulatedObjective" << std::endl;
      RealT tol = 1e-8;
      simobj->value(*x, tol);
      simobj->checkGradient(*x, *p, true, *outStream);
      simobj->checkHessVec(*x, *p, true, *outStream);
    }

    zp->zero();
    auto vusim = ROL::dynamicPtrCast<ROL::SimulatedVector<RealT>>(vu);
    for (int i = 0; i < sampler->numMySamples(); ++i) {
      RealT tol = 1e-8;
      std::vector<RealT> param = sampler->getMyPoint(i);
      con->setParameter(param);
      vusim->get(i)->zero();
      con->update(*(vusim->get(i)),*zp);
      con->solve(*rp,*(vusim->get(i)),*zp,tol);
    }

    bool zeroInit = parlist->sublist("Problem").get("Zero initial guess",true);
    if (zeroInit) {
      x->zero();
      vp->zero();
    }

    /*************************************************************************/
    /***************** SOLVE PROBLEM *****************************************/
    /*************************************************************************/
    auto opt = ROL::makePtr<ROL::Problem<RealT>>(simobj, x);
    opt->addConstraint("PDE", simcon, vp);
    opt->finalize(false,true,*outStream);
    ROL::Solver<RealT> solver(opt,*parlist);
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
    up->zero(); pp->zero();
    for (int i = 0; i < sampler->numMySamples(); ++i)
      up->axpy(sampler->getMyWeight(i),*(vusim->get(i)));
    bman->sumAll(*up,*pp);
    con->outputTpetraVector(p_ptr,"mean_state.txt");
    // Build full objective function distribution
    *outStream << std::endl << "Print Objective CDF" << std::endl;
    int nsamp_dist = parlist->sublist("Problem").get("Number of output samples",100);
    auto sampler_dist = ROL::makePtr<ROL::MonteCarloGenerator<RealT>>(nsamp_dist,distVec,bman);
    print<RealT>(*robj,*zp,*sampler_dist,nsamp_dist,comm,"obj_samples.txt");
    // Build vorticity objective function distribution
    auto obj0 = ROL::makePtr<IntegralObjective<RealT,DeviceT>>(qoi_vec[0],assembler);
    auto stateStore0 = ROL::makePtr<ROL::VectorController<RealT>>();
    auto robj0 = ROL::makePtr<ROL::Reduced_Objective_SimOpt<RealT>>(obj0, con, stateStore0, up, zp, pp, true, false);
    print<RealT>(*robj0,*zp,*sampler_dist,nsamp_dist,comm,"vort_samples.txt");

    *outStream << "Output time: "
               << static_cast<RealT>(std::clock()-timer_print)/static_cast<RealT>(CLOCKS_PER_SEC)
               << " seconds." << std::endl << std::endl;
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
