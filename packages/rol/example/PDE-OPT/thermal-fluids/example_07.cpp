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
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "Tpetra_Core.hpp"
#include "Tpetra_Version.hpp"

#include <iostream>
#include <algorithm>
//#include <fenv.h>

#include "ROL_TpetraMultiVector.hpp"
#include "ROL_Reduced_Objective_SimOpt.hpp"
#include "ROL_OptimizationSolver.hpp"
#include "ROL_MonteCarloGenerator.hpp"
#include "ROL_TpetraTeuchosBatchManager.hpp"

#include "../TOOLS/pdeconstraint.hpp"
#include "../TOOLS/pdeobjective.hpp"
#include "../TOOLS/pdevector.hpp"
#include "../TOOLS/meshreader.hpp"
#include "pde_thermal-fluids_ex07.hpp"
#include "obj_thermal-fluids_ex03.hpp"

typedef double RealT;

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
  std::vector<std::vector<Real> > mysamples(sdim, myzerovec);
  std::vector<std::vector<double> > gsamples(sdim, gzerovec);
  for (int i = 0; i < nsamp; ++i) {
    sample = sampler.getMyPoint(i);
    obj.setParameter(sample);
    myvalues[i] = static_cast<double>(obj.value(z,tol));
    for (int j = 0; j < sdim; ++j) {
      mysamples[j][i] = static_cast<double>(sample[j]);
    }
  }

  // Send data to root processor
#ifdef HAVE_MPI
  ROL::Ptr<const Teuchos::MpiComm<int> > mpicomm
    = ROL::dynamicPtrCast<const Teuchos::MpiComm<int> >(comm);
  int nproc = Teuchos::size<int>(*mpicomm);
  std::vector<int> sampleCounts(nproc, 0), sampleDispls(nproc, 0);
  MPI_Gather(&nsamp,1,MPI_INT,&sampleCounts[0],1,MPI_INT,0,*(mpicomm->getRawMpiComm())());
  for (int i = 1; i < nproc; ++i) {
    sampleDispls[i] = sampleDispls[i-1] + sampleCounts[i-1];
  }
  MPI_Gatherv(&myvalues[0],nsamp,MPI_DOUBLE,&gvalues[0],&sampleCounts[0],&sampleDispls[0],MPI_DOUBLE,0,*(mpicomm->getRawMpiComm())());
  for (int j = 0; j < sdim; ++j) {
    MPI_Gatherv(&mysamples[j][0],nsamp,MPI_DOUBLE,&gsamples[j][0],&sampleCounts[0],&sampleDispls[0],MPI_DOUBLE,0,*(mpicomm->getRawMpiComm())());
  }
#else
  gvalues.assign(myvalues.begin(),myvalues.end());
  for (int j = 0; j < sdim; ++j) {
    gsamples[j].assign(mysamples[j].begin(),mysamples[j].end());
  }
#endif

  // Print
  int rank  = Teuchos::rank<int>(*comm);
  if ( rank==0 ) {
    std::ofstream file;
    file.open(filename);
    file << std::scientific << std::setprecision(15);
    for (int i = 0; i < ngsamp; ++i) {
      for (int j = 0; j < sdim; ++j) {
        file << std::setw(25) << std::left << gsamples[j][i];
      }
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
  Teuchos::GlobalMPISession mpiSession (&argc, &argv, &bhs);
  ROL::Ptr<const Teuchos::Comm<int> > comm
    = Tpetra::getDefaultComm();
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
    std::string filename = "input_ex07.xml";
    Teuchos::RCP<Teuchos::ParameterList> parlist = Teuchos::rcp( new Teuchos::ParameterList() );
    Teuchos::updateParametersFromXmlFile( filename, parlist.ptr() );

    /*** Initialize main data structure. ***/
    ROL::Ptr<MeshManager<RealT> > meshMgr
      = ROL::makePtr<MeshReader<RealT>>(*parlist);
    // Initialize PDE describing Navier-Stokes equations.
    ROL::Ptr<PDE_ThermalFluids_ex07<RealT> > pde
      = ROL::makePtr<PDE_ThermalFluids_ex07<RealT>>(*parlist);
    ROL::Ptr<ROL::Constraint_SimOpt<RealT> > con
      = ROL::makePtr<PDE_Constraint<RealT>>(pde,meshMgr,serial_comm,*parlist,*outStream);
    // Cast the constraint and get the assembler.
    ROL::Ptr<PDE_Constraint<RealT> > pdecon
      = ROL::dynamicPtrCast<PDE_Constraint<RealT> >(con);
    ROL::Ptr<Assembler<RealT> > assembler = pdecon->getAssembler();
    con->setSolveParameters(*parlist);
    pdecon->outputTpetraData();

    // Create state vector and set to zeroes
    ROL::Ptr<Tpetra::MultiVector<> > u_ptr, p_ptr, r_ptr, z_ptr;
    u_ptr  = assembler->createStateVector();
    p_ptr  = assembler->createStateVector();
    r_ptr  = assembler->createResidualVector();
    z_ptr  = assembler->createControlVector();
    ROL::Ptr<ROL::Vector<RealT> > up, pp, rp, zp;
    up  = ROL::makePtr<PDE_PrimalSimVector<RealT>>(u_ptr,pde,assembler);
    pp  = ROL::makePtr<PDE_PrimalSimVector<RealT>>(p_ptr,pde,assembler);
    rp  = ROL::makePtr<PDE_DualSimVector<RealT>>(r_ptr,pde,assembler);
    zp  = ROL::makePtr<PDE_PrimalOptVector<RealT>>(z_ptr,pde,assembler);

    // Create ROL SimOpt vectors
    //ROL::Vector_SimOpt<RealT> x(up,zp);

    // Initialize objective function.
    std::vector<ROL::Ptr<QoI<RealT> > > qoi_vec(2,ROL::nullPtr);
    qoi_vec[0] = ROL::makePtr<QoI_State_ThermalFluids<RealT>>(*parlist,
                                                                 pde->getVelocityFE(),
                                                                 pde->getPressureFE(),
                                                                 pde->getThermalFE(),
                                                                 pde->getFieldHelper());
    qoi_vec[1] = ROL::makePtr<QoI_L2Penalty_ThermalFluids<RealT>>(pde->getVelocityFE(),
                                                                     pde->getPressureFE(),
                                                                     pde->getThermalFE(),
                                                                     pde->getThermalBdryFE(),
                                                                     pde->getBdryCellLocIds(),
                                                                     pde->getFieldHelper());
    ROL::Ptr<StdObjective_ThermalFluids<RealT> > std_obj
      = ROL::makePtr<StdObjective_ThermalFluids<RealT>>(*parlist);
    ROL::Ptr<ROL::Objective_SimOpt<RealT> > obj
      = ROL::makePtr<PDE_Objective<RealT>>(qoi_vec,std_obj,assembler);
    ROL::Ptr<ROL::VectorController<RealT> > stateStore
      = ROL::makePtr<ROL::VectorController<RealT>>();
    ROL::Ptr<ROL::Reduced_Objective_SimOpt<RealT> > robj
      = ROL::makePtr<ROL::Reduced_Objective_SimOpt<RealT>>(obj, con, stateStore, up, zp, pp, true, false);

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
    std::vector<ROL::Ptr<ROL::Distribution<RealT> > > distVec(stochDim);
    Teuchos::ParameterList UList;
    UList.sublist("Distribution").set("Name","Uniform");
    UList.sublist("Distribution").sublist("Uniform").set("Lower Bound",-1.0);
    UList.sublist("Distribution").sublist("Uniform").set("Upper Bound", 1.0);
    for (int i = 0; i < stochDim; ++i) {
      distVec[i] = ROL::DistributionFactory<RealT>(UList);
    }
    // Sampler
    ROL::Ptr<ROL::BatchManager<RealT> > bman
      = ROL::makePtr<ROL::TpetraTeuchosBatchManager<RealT>>(comm);
    ROL::Ptr<ROL::SampleGenerator<RealT> > sampler
      = ROL::makePtr<ROL::MonteCarloGenerator<RealT>>(nsamp,distVec,bman);
    ROL::Ptr<ROL::SampleGenerator<RealT> > sampler_dist
      = ROL::makePtr<ROL::MonteCarloGenerator<RealT>>(nsamp_dist,distVec,bman);

    // Set up optimization problem, check derivatives and solve
    zp->zero();
    ROL::OptimizationProblem<RealT> optProb(robj,zp);
    bool isStoch = parlist->sublist("Problem").get("Is stochastic?",false);
    if (isStoch) {
      optProb.setStochasticObjective(*parlist,sampler);
    }
    if ( parlist->sublist("Problem").get("Check derivatives",false) ) {
      optProb.check(*outStream);
      ROL::Ptr<ROL::Vector<RealT> > xp
        = ROL::makePtr<ROL::Vector_SimOpt<RealT>>(up,zp);
      ROL::OptimizationProblem<RealT> op(obj,xp,con,pp);
      op.check(*outStream);
    }
    ROL::OptimizationSolver<RealT> optSolver(optProb,*parlist);
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
    pdecon->outputTpetraVector(u_ptr,"state.txt");

    pdecon->outputTpetraVector(z_ptr,"control.txt");
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
