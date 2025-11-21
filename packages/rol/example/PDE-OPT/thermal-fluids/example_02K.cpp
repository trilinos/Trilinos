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
#include "ROL_Reduced_Objective_SimOpt.hpp"
#include "ROL_MonteCarloGenerator.hpp"
#include "ROL_StochasticProblem.hpp"
#include "ROL_Solver.hpp"
#include "ROL_DistributionFactory.hpp"
#include "ROL_TpetraTeuchosBatchManager.hpp"

#include "../TOOLS/meshmanagerK.hpp"
#include "../TOOLS/pdeconstraintK.hpp"
#include "../TOOLS/pdeobjectiveK.hpp"
#include "../TOOLS/pdevectorK.hpp"
#include "pde_thermal-fluidsK.hpp"
#include "obj_thermal-fluidsK.hpp"

using RealT = double;
using DeviceT = Kokkos::HostSpace;

template<class Real>
void print(ROL::Objective<Real> &obj,
           const ROL::Vector<Real> &z,
           ROL::SampleGenerator<Real> &sampler,
           const int ngsamp,
           const ROL::Ptr<const Teuchos::Comm<int> > &comm,
           const std::string &filename) {
  Real tol(1e-8);
  // Build objective function distribution
  int nsamp = sampler.numMySamples();
  std::vector<Real> myvalues(nsamp), myzerovec(nsamp, 0);
  std::vector<double> gvalues(ngsamp), gzerovec(ngsamp, 0);
  std::vector<Real> sample = sampler.getMyPoint(0);
  int sdim = sample.size();
  std::vector<std::vector<Real>> mysamples(sdim, myzerovec);
  std::vector<std::vector<double>> gsamples(sdim, gzerovec);
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
    std::string filename = "input.xml";
    auto parlist = ROL::getParametersFromXmlFile(filename);

    /*** Initialize main data structure. ***/
    auto meshMgr = ROL::makePtr<MeshManager_BackwardFacingStepChannel<RealT,DeviceT>>(*parlist);
    // Initialize PDE describing Navier-Stokes equations.
    auto pde = ROL::makePtr<PDE_ThermalFluids<RealT,DeviceT>>(*parlist);
    auto con = ROL::makePtr<PDE_Constraint<RealT,DeviceT>>(pde,meshMgr,serial_comm,*parlist,*outStream);
    // Cast the constraint and get the assembler.
    auto assembler = con->getAssembler();
    con->setSolveParameters(*parlist);
    con->outputTpetraData();

    // Create state vector and set to zeroes
    auto u_ptr  = assembler->createStateVector();     u_ptr->randomize();
    auto p_ptr  = assembler->createStateVector();     p_ptr->randomize();
    auto du_ptr = assembler->createStateVector();     du_ptr->randomize();
    auto yu_ptr = assembler->createStateVector();     yu_ptr->randomize();
    auto yp_ptr = assembler->createStateVector();     yp_ptr->randomize();
    auto r_ptr  = assembler->createResidualVector();  r_ptr->randomize();
    auto z_ptr  = assembler->createControlVector();   z_ptr->randomize();
    auto dz_ptr = assembler->createControlVector();   dz_ptr->randomize();
    auto up  = ROL::makePtr<PDE_PrimalSimVector<RealT,DeviceT>>(u_ptr,pde,assembler);
    auto pp  = ROL::makePtr<PDE_PrimalSimVector<RealT,DeviceT>>(p_ptr,pde,assembler);
    auto dup = ROL::makePtr<PDE_PrimalSimVector<RealT,DeviceT>>(du_ptr,pde,assembler);
    auto yup = ROL::makePtr<PDE_PrimalSimVector<RealT,DeviceT>>(yu_ptr,pde,assembler);
    auto ypp = ROL::makePtr<PDE_PrimalSimVector<RealT,DeviceT>>(yp_ptr,pde,assembler);
    auto rp  = ROL::makePtr<PDE_DualSimVector<RealT,DeviceT>>(r_ptr,pde,assembler);
    auto zp  = ROL::makePtr<PDE_PrimalOptVector<RealT,DeviceT>>(z_ptr,pde,assembler);
    auto dzp = ROL::makePtr<PDE_PrimalOptVector<RealT,DeviceT>>(dz_ptr,pde,assembler);
    // Create ROL SimOpt vectors
    ROL::Vector_SimOpt<RealT> x(up,zp);
    ROL::Vector_SimOpt<RealT> d(dup,dzp);

    // Initialize objective function.
    std::vector<ROL::Ptr<QoI<RealT,DeviceT>>> qoi_vec(2,ROL::nullPtr);
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
    int stochDim = 6;
    int nsamp = parlist->sublist("Problem").get("Number of samples",100);
    std::vector<ROL::Ptr<ROL::Distribution<RealT>>> distVec(stochDim);
    // Reynolds Number
    Teuchos::ParameterList ReList;
    ReList.sublist("Distribution").set("Name","Uniform");
    ReList.sublist("Distribution").sublist("Uniform").set("Lower Bound",0.0);
    ReList.sublist("Distribution").sublist("Uniform").set("Upper Bound",1.0);
    distVec[0] = ROL::DistributionFactory<RealT>(ReList);
    // Prandtl Number
    Teuchos::ParameterList PrList;
    PrList.sublist("Distribution").set("Name","Uniform");
    PrList.sublist("Distribution").sublist("Uniform").set("Lower Bound",0.0);
    PrList.sublist("Distribution").sublist("Uniform").set("Upper Bound",1.0);
    distVec[1] = ROL::DistributionFactory<RealT>(PrList);
    // Grashof Number
    Teuchos::ParameterList GrList;
    GrList.sublist("Distribution").set("Name","Beta");
    GrList.sublist("Distribution").sublist("Beta").set("Shape 1",5.0);
    GrList.sublist("Distribution").sublist("Beta").set("Shape 2",2.0);
    distVec[2] = ROL::DistributionFactory<RealT>(GrList);
    // Thickness Number
    Teuchos::ParameterList ThList;
    ThList.sublist("Distribution").set("Name","Beta");
    ThList.sublist("Distribution").sublist("Beta").set("Shape 1",5.0);
    ThList.sublist("Distribution").sublist("Beta").set("Shape 2",3.0);
    distVec[3] = ROL::DistributionFactory<RealT>(ThList);
    // Inflow/Outflow
    Teuchos::ParameterList BdList;
    BdList.sublist("Distribution").set("Name","Uniform");
    BdList.sublist("Distribution").sublist("Uniform").set("Lower Bound",-1.0);
    BdList.sublist("Distribution").sublist("Uniform").set("Upper Bound", 1.0);
    distVec[4] = ROL::DistributionFactory<RealT>(BdList);
    distVec[5] = ROL::DistributionFactory<RealT>(BdList);
    // Sampler
    auto bman = ROL::makePtr<ROL::TpetraTeuchosBatchManager<RealT>>(comm);
    auto sampler = ROL::makePtr<ROL::MonteCarloGenerator<RealT>>(nsamp,distVec,bman);

    /*************************************************************************/
    /***************** BUILD STOCHASTIC PROBLEM ******************************/
    /*************************************************************************/
    auto problem = ROL::makePtr<ROL::StochasticProblem<RealT>>(robj, zp);
    parlist->sublist("SOL").sublist("Objective").sublist("Risk Measure").set("Name","CVaR");
    parlist->sublist("SOL").set("Initial Statistic", 1.0);
    problem->makeObjectiveStochastic(*parlist, sampler);

    // Run derivative checks
    bool checkDeriv = parlist->sublist("Problem").get("Check derivatives",false);
    if ( checkDeriv ) {
      *outStream << "Check Gradient of Full Objective Function" << std::endl;
      obj->checkGradient(x,d,true,*outStream);
      *outStream << std::endl << "Check Hessian of Full Objective Function" << std::endl;
      obj->checkHessVec(x,d,true,*outStream);
      *outStream << std::endl << "Check Jacobian of Constraint" << std::endl;
      con->checkApplyJacobian(x,d,*up,true,*outStream);
      *outStream << std::endl << "Check Jacobian_1 of Constraint" << std::endl;
      con->checkApplyJacobian_1(*up,*zp,*dup,*rp,true,*outStream);
      *outStream << std::endl << "Check Jacobian_2 of Constraint" << std::endl;
      con->checkApplyJacobian_2(*up,*zp,*dzp,*rp,true,*outStream);
      *outStream << std::endl << "Check Hessian of Constraint" << std::endl;
      con->checkApplyAdjointHessian(x,*dup,d,x,true,*outStream);
      *outStream << std::endl << "Check Hessian_11 of Constraint" << std::endl;
      con->checkApplyAdjointHessian_11(*up,*zp,*pp,*dup,*rp,true,*outStream);
      *outStream << std::endl << "Check Hessian_12 of Constraint" << std::endl;
      con->checkApplyAdjointHessian_12(*up,*zp,*pp,*dup,*dzp,true,*outStream);
      *outStream << std::endl << "Check Hessian_21 of Constraint" << std::endl;
      con->checkApplyAdjointHessian_21(*up,*zp,*pp,*dzp,*rp,true,*outStream);
      *outStream << std::endl << "Check Hessian_22 of Constraint" << std::endl;
      con->checkApplyAdjointHessian_22(*up,*zp,*pp,*dzp,*dzp,true,*outStream);

      *outStream << std::endl << "Check Adjoint Jacobian of Constraint" << std::endl;
      con->checkAdjointConsistencyJacobian(*dup,d,x,true,*outStream);
      *outStream << std::endl << "Check Adjoint Jacobian_1 of Constraint" << std::endl;
      con->checkAdjointConsistencyJacobian_1(*pp,*dup,*up,*zp,true,*outStream);
      *outStream << std::endl << "Check Adjoint Jacobian_2 of Constraint" << std::endl;
      con->checkAdjointConsistencyJacobian_2(*pp,*dzp,*up,*zp,true,*outStream);

      *outStream << std::endl << "Check Constraint Solve" << std::endl;
      con->checkSolve(*up,*zp,*rp,true,*outStream);
      *outStream << std::endl << "Check Inverse Jacobian_1 of Constraint" << std::endl;
      con->checkInverseJacobian_1(*rp,*dup,*up,*zp,true,*outStream);
      *outStream << std::endl << "Check Inverse Adjoint Jacobian_1 of Constraint" << std::endl;
      con->checkInverseAdjointJacobian_1(*rp,*pp,*up,*zp,true,*outStream);

      *outStream << std::endl << "Check Gradient of Reduced Objective Function" << std::endl;
      robj->checkGradient(*zp,*dzp,true,*outStream);
      *outStream << std::endl << "Check Hessian of Reduced Objective Function" << std::endl;
      robj->checkHessVec(*zp,*dzp,true,*outStream);
    }

    up->zero(); zp->zero();
    parlist->sublist("Step").set("Type","Trust Region");
    ROL::Solver<RealT> solver(problem,*parlist);
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
    up->zero(); pp->zero(); dup->zero(); dzp->zero(); yup->zero();
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
      con->solve(*rp,*dup,*dzp,tol);
      yup->axpy(sampler->getMyWeight(i),*dup);
      for (int j = 0; j < stochDim; ++j)
        file_samp << std::setw(25) << std::left << sample[j];
      file_samp << std::endl;
    }
    file_samp.close();
    bman_Eu->sumAll(*up,*pp);
    bman_Eu->sumAll(*yup,*ypp);
    con->outputTpetraVector(p_ptr,"mean_state.txt");
    con->outputTpetraVector(yp_ptr,"mean_uncontrolled_state.txt");
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
