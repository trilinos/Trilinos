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
#include "ROL_TpetraTeuchosBatchManager.hpp"
#include "ROL_StochasticProblem.hpp"
#include "ROL_Solver.hpp"

#include "../TOOLS/pdeconstraintK.hpp"
#include "../TOOLS/pdeobjectiveK.hpp"
#include "../TOOLS/pdevectorK.hpp"
#include "../TOOLS/meshreaderK.hpp"
#include "pde_thermal-fluids_ex07K.hpp"
#include "obj_thermal-fluids_ex03K.hpp"

using RealT = double;
using DeviceT = Kokkos::HostSpace;

template<class Real>
void print(ROL::Objective<Real> &obj,
           const ROL::Vector<Real> &z,
           ROL::SampleGenerator<Real> &sampler,
           int ngsamp,
           const ROL::Ptr<const Teuchos::Comm<int> > &comm,
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
  auto mpicomm = ROL::dynamicPtrCast<const Teuchos::MpiComm<int> >(comm);
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
  int rank = Teuchos::rank<int>(*comm);
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
    std::string filename = "input_ex07.xml";
    auto parlist = ROL::getParametersFromXmlFile(filename);

    /*** Initialize main data structure. ***/
    auto meshMgr = ROL::makePtr<MeshReader<RealT,DeviceT>>(*parlist);
    // Initialize PDE describing Navier-Stokes equations.
    auto pde = ROL::makePtr<PDE_ThermalFluids_ex07<RealT,DeviceT>>(*parlist);
    auto con = ROL::makePtr<PDE_Constraint<RealT,DeviceT>>(pde,meshMgr,serial_comm,*parlist,*outStream);
    // Cast the constraint and get the assembler.
    auto assembler = con->getAssembler();
    con->setSolveParameters(*parlist);
    con->outputTpetraData();

    // Create state vector and set to zeroes
    auto u_ptr  = assembler->createStateVector();
    auto p_ptr  = assembler->createStateVector();
    auto r_ptr  = assembler->createResidualVector();
    auto z_ptr  = assembler->createControlVector();
    auto up  = ROL::makePtr<PDE_PrimalSimVector<RealT,DeviceT>>(u_ptr,pde,assembler);
    auto pp  = ROL::makePtr<PDE_PrimalSimVector<RealT,DeviceT>>(p_ptr,pde,assembler);
    auto rp  = ROL::makePtr<PDE_DualSimVector<RealT,DeviceT>>(r_ptr,pde,assembler);
    auto zp  = ROL::makePtr<PDE_PrimalOptVector<RealT,DeviceT>>(z_ptr,pde,assembler);

    // Create ROL SimOpt vectors
    //ROL::Vector_SimOpt<RealT> x(up,zp);

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
    int Nbottom = parlist->sublist("Problem").get("Bottom KL Truncation Order",5);
    int N0      = parlist->sublist("Problem").get("Side 0 KL Truncation Order",5);
    int N1      = parlist->sublist("Problem").get("Side 1 KL Truncation Order",5);
    int N2      = parlist->sublist("Problem").get("Side 2 KL Truncation Order",5);
    int N3      = parlist->sublist("Problem").get("Side 3 KL Truncation Order",5);
    int stochDim = 2*(Nbottom + N0 + N1 + N2 + N3) + 3;
    int nsamp = parlist->sublist("Problem").get("Number of samples",100);
    int nsamp_dist = parlist->sublist("Problem").get("Number of output samples",100);
    // Build vector of distributions
    std::vector<ROL::Ptr<ROL::Distribution<RealT>>> distVec(stochDim);
    Teuchos::ParameterList UList;
    UList.sublist("Distribution").set("Name","Uniform");
    UList.sublist("Distribution").sublist("Uniform").set("Lower Bound",-1.0);
    UList.sublist("Distribution").sublist("Uniform").set("Upper Bound", 1.0);
    for (int i = 0; i < stochDim; ++i)
      distVec[i] = ROL::DistributionFactory<RealT>(UList);
    // Sampler
    auto bman = ROL::makePtr<ROL::TpetraTeuchosBatchManager<RealT>>(comm);
    auto sampler = ROL::makePtr<ROL::MonteCarloGenerator<RealT>>(nsamp,distVec,bman);
    auto sampler_dist = ROL::makePtr<ROL::MonteCarloGenerator<RealT>>(nsamp_dist,distVec,bman);

    // Set up optimization problem, check derivatives and solve
    zp->zero();
    auto optProb = ROL::makePtr<ROL::StochasticProblem<RealT>>(robj,zp);
    bool isStoch = parlist->sublist("Problem").get("Is stochastic?",false);
    if (isStoch) {
      optProb->makeObjectiveStochastic(*parlist,sampler);
    }
    if ( parlist->sublist("Problem").get("Check derivatives",false) ) {
      optProb->check(true,*outStream);
      auto xp = ROL::makePtr<ROL::Vector_SimOpt<RealT>>(up,zp);
      ROL::Problem<RealT> op(obj,xp);
      op.addConstraint("PDE",con,pp);
      op.check(true,*outStream);
    }
    optProb->finalize(false,true,*outStream);
    ROL::Solver<RealT> optSolver(optProb,*parlist);
    optSolver.solve(*outStream);

    // Output.
    RealT tol = std::sqrt(ROL::ROL_EPSILON<RealT>());
    if (isStoch) {
      std::vector<RealT> pt;
      RealT wt;
      up->zero(); pp->zero();
      for (int i = 0; i < sampler->numMySamples(); ++i) {
        pt = sampler->getMyPoint(i);
        wt = sampler->getMyWeight(i);
        con->setParameter(pt);
        con->solve(*rp,*up,*zp,tol);
        pp->axpy(wt,*up);
      }
      up->zero();
      sampler->sumAll(*pp,*up);
    } else {
      con->solve(*rp,*up,*zp,tol);
    }
    con->outputTpetraVector(u_ptr,"state.txt");

    con->outputTpetraVector(z_ptr,"control.txt");
    assembler->printMeshData(*outStream);
    if (isStoch) {
      print<RealT>(*robj,*zp,*sampler_dist,nsamp_dist,comm,"obj_samples.txt");
    }

    // Get a summary from the time monitor.
    Teuchos::TimeMonitor::summarize();
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
