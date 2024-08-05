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
#include "ROL_Algorithm.hpp"
#include "ROL_Reduced_Objective_SimOpt.hpp"
#include "ROL_MonteCarloGenerator.hpp"
#include "ROL_OptimizationSolver.hpp"
#include "ROL_TpetraTeuchosBatchManager.hpp"

#include "../TOOLS/pdeconstraint.hpp"
#include "../TOOLS/pdeobjective.hpp"
#include "../TOOLS/pdevector.hpp"
#include "../TOOLS/batchmanager.hpp"
#include "pde_thermal-fluids_ex03.hpp"
#include "obj_thermal-fluids_ex03.hpp"
#include "mesh_thermal-fluids.hpp"

typedef double RealT;

template<class Real>
void setUpAndSolve(ROL::OptimizationProblem<Real> &opt,
                   Teuchos::ParameterList &parlist,
                   std::ostream &outStream) {
  parlist.sublist("Step").set("Type","Trust Region");
  ROL::OptimizationSolver<RealT> solver(opt,parlist);
  Teuchos::Time timer("Optimization Time", true);
  solver.solve(outStream);
  timer.stop();
  outStream << "Total optimization time = " << timer.totalElapsedTime() << " seconds." << std::endl;
}

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
    std::string filename = "input_ex06.xml";
    Teuchos::RCP<Teuchos::ParameterList> parlist = Teuchos::rcp( new Teuchos::ParameterList() );
    Teuchos::updateParametersFromXmlFile( filename, parlist.ptr() );
    parlist->sublist("SimOpt").sublist("Solve").set("Output Iteration History",myRank==0);

    /*** Initialize main data structure. ***/
    ROL::Ptr<MeshManager<RealT> > meshMgr
      = ROL::makePtr<MeshManager_ThermalFluids<RealT>>(*parlist);
    // Initialize PDE describing Navier-Stokes equations.
    ROL::Ptr<PDE_ThermalFluids_ex03<RealT> > pde
      = ROL::makePtr<PDE_ThermalFluids_ex03<RealT>>(*parlist);
    ROL::Ptr<ROL::Constraint_SimOpt<RealT> > con
      = ROL::makePtr<PDE_Constraint<RealT>>(pde,meshMgr,serial_comm,*parlist,*outStream);
    // Cast the constraint and get the assembler.
    ROL::Ptr<PDE_Constraint<RealT> > pdecon
      = ROL::dynamicPtrCast<PDE_Constraint<RealT> >(con);
    ROL::Ptr<Assembler<RealT> > assembler = pdecon->getAssembler();
    con->setSolveParameters(*parlist);
    pdecon->outputTpetraData();

    // Create state vector and set to zeroes
    ROL::Ptr<Tpetra::MultiVector<> > u_ptr, p_ptr, du_ptr, dy_ptr, yu_ptr, yp_ptr, r_ptr, z_ptr, dz_ptr;
    u_ptr  = assembler->createStateVector();     u_ptr->randomize();
    p_ptr  = assembler->createStateVector();     p_ptr->randomize();
    du_ptr = assembler->createStateVector();     du_ptr->randomize();
    dy_ptr = assembler->createStateVector();     dy_ptr->randomize();
    yu_ptr = assembler->createStateVector();     yu_ptr->randomize();
    yp_ptr = assembler->createStateVector();     yp_ptr->randomize();
    r_ptr  = assembler->createResidualVector();  r_ptr->randomize();
    z_ptr  = assembler->createControlVector();   z_ptr->randomize();
    dz_ptr = assembler->createControlVector();   dz_ptr->randomize();
    ROL::Ptr<ROL::Vector<RealT> > up, pp, dup, dyp, yup, ypp, rp, zp, dzp;
    up  = ROL::makePtr<PDE_PrimalSimVector<RealT>>(u_ptr,pde,assembler);
    pp  = ROL::makePtr<PDE_PrimalSimVector<RealT>>(p_ptr,pde,assembler);
    dup = ROL::makePtr<PDE_PrimalSimVector<RealT>>(du_ptr,pde,assembler);
    dyp = ROL::makePtr<PDE_PrimalSimVector<RealT>>(dy_ptr,pde,assembler);
    yup = ROL::makePtr<PDE_PrimalSimVector<RealT>>(yu_ptr,pde,assembler);
    ypp = ROL::makePtr<PDE_PrimalSimVector<RealT>>(yp_ptr,pde,assembler);
    rp  = ROL::makePtr<PDE_DualSimVector<RealT>>(r_ptr,pde,assembler);
    zp  = ROL::makePtr<PDE_PrimalOptVector<RealT>>(z_ptr,pde,assembler);
    dzp = ROL::makePtr<PDE_PrimalOptVector<RealT>>(dz_ptr,pde,assembler);
    // Create ROL SimOpt vectors
    ROL::Vector_SimOpt<RealT> x(up,zp);
    ROL::Vector_SimOpt<RealT> d(dup,dzp);

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
    ROL::Ptr<ROL::Objective_SimOpt<RealT> > objState
      = ROL::makePtr<PDE_Objective<RealT>>(qoi_vec[0],assembler);
    ROL::Ptr<ROL::Objective_SimOpt<RealT> > objCtrl
      = ROL::makePtr<PDE_Objective<RealT>>(qoi_vec[1],assembler);
    ROL::Ptr<ROL::VectorController<RealT> > stateStore
      = ROL::makePtr<ROL::VectorController<RealT>>();
    ROL::Ptr<ROL::Reduced_Objective_SimOpt<RealT> > robj
      = ROL::makePtr<ROL::Reduced_Objective_SimOpt<RealT>>(obj, con, stateStore, up, zp, pp, true, false);

    /*************************************************************************/
    /***************** BUILD SAMPLER *****************************************/
    /*************************************************************************/
    int Nbottom = parlist->sublist("Problem").get("Bottom KL Truncation Order",5);
    int Nleft   = parlist->sublist("Problem").get("Left KL Truncation Order",5);
    int Nright  = parlist->sublist("Problem").get("Right KL Truncation Order",5);
    int stochDim = Nbottom + Nleft + Nright + 3;
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
    //  = ROL::makePtr<PDE_OptVector_BatchManager<RealT>>(comm);
    ROL::Ptr<ROL::SampleGenerator<RealT> > sampler
      = ROL::makePtr<ROL::MonteCarloGenerator<RealT>>(nsamp,distVec,bman);
    ROL::Ptr<ROL::SampleGenerator<RealT> > sampler_dist
      = ROL::makePtr<ROL::MonteCarloGenerator<RealT>>(nsamp_dist,distVec,bman);

    /*************************************************************************/
    /***************** BUILD STOCHASTIC PROBLEM ******************************/
    /*************************************************************************/
    ROL::Ptr<ROL::OptimizationProblem<RealT> > opt;
    std::vector<RealT> ctrl, stat;

    Teuchos::Array<RealT> lambdaArray
      = Teuchos::getArrayFromStringParameter<RealT>(parlist->sublist("Problem"),"Convex Combination Parameters");
    std::vector<RealT> lambda = lambdaArray.toVector();
    std::sort(lambda.begin(),lambda.end());
    int N = lambda.size();
    RealT alpha = parlist->sublist("Problem").get("CVaR Confidence Level",0.9);

    /*************************************************************************/
    /***************** SOLVE RISK NEUTRAL ************************************/
    /*************************************************************************/
    zp->zero(); up->zero();
    RealT tol(1e-8), zero(0), one(1);
    bool lambdaZero = (lambda[0] == zero);
    if ( lambdaZero ) {
      lambda.erase(lambda.begin()); --N;
      // Solve.
      parlist->sublist("SOL").set("Type","Risk Neutral");
      opt = ROL::makePtr<ROL::OptimizationProblem<RealT>>(robj,zp);
      parlist->sublist("SOL").set("Initial Statistic",one);
      opt->setStochasticObjective(*parlist,sampler);
      setUpAndSolve<RealT>(*opt,*parlist,*outStream);
      // Output.
      ctrl.push_back(objCtrl->value(*up,*zp,tol));
      stat.push_back(opt->getSolutionStatistic());
      std::string CtrlRN = "control_RN.txt";
      pdecon->outputTpetraVector(z_ptr,CtrlRN);
      std::string ObjRN = "obj_samples_RN.txt";
      print<RealT>(*robj,*zp,*sampler_dist,nsamp_dist,comm,ObjRN);
    }

    /*************************************************************************/
    /***************** SOLVE MEAN PLUS CVAR **********************************/
    /*************************************************************************/
    parlist->sublist("SOL").set("Type","Risk Averse");
    parlist->sublist("SOL").sublist("Risk Measure").set("Name","CVaR");
    parlist->sublist("SOL").sublist("Risk Measure").sublist("CVaR").set("Confidence Level",alpha);
    parlist->sublist("SOL").sublist("Risk Measure").sublist("CVaR").set("Smoothing Parameter",1e-4);
    parlist->sublist("SOL").sublist("Risk Measure").sublist("CVaR").sublist("Distribution").set("Name","Parabolic");
    parlist->sublist("SOL").sublist("Risk Measure").sublist("CVaR").sublist("Distribution").sublist("Parabolic").set("Lower Bound",0.0);
    parlist->sublist("SOL").sublist("Risk Measure").sublist("CVaR").sublist("Distribution").sublist("Parabolic").set("Upper Bound",1.0);
    for (int i = 0; i < N; ++i) {
      parlist->sublist("SOL").sublist("Risk Measure").sublist("CVaR").set("Convex Combination Parameter",one-lambda[i]);
      // Solve.
      opt = ROL::makePtr<ROL::OptimizationProblem<RealT>>(robj,zp);
      parlist->sublist("SOL").set("Initial Statistic",stat[i]);
      opt->setStochasticObjective(*parlist,sampler);
      setUpAndSolve<RealT>(*opt,*parlist,*outStream);
      // Output.
      ctrl.push_back(objCtrl->value(*up,*zp,tol));
      stat.push_back(opt->getSolutionStatistic());
      std::stringstream nameCtrl;
      nameCtrl << "control_MCVaR_" << i << ".txt";
      pdecon->outputTpetraVector(z_ptr,nameCtrl.str().c_str());
      std::stringstream nameObj;
      nameObj << "obj_samples_MCVaR_" << i << ".txt";
      print<RealT>(*robj,*zp,*sampler_dist,nsamp_dist,comm,nameObj.str());
    }
    
    /*************************************************************************/
    /***************** OUTPUT RESULTS ****************************************/
    /*************************************************************************/
    const int rank = Teuchos::rank<int>(*comm);
    if ( rank==0 ) {
      std::stringstream nameCTRL, nameVAR;
      nameCTRL << "ctrl.txt";
      nameVAR << "stat.txt";
      std::ofstream fileCTRL, fileVAR;
      fileCTRL.open(nameCTRL.str());
      fileVAR.open(nameVAR.str());
      fileCTRL << std::scientific << std::setprecision(15);
      fileVAR << std::scientific << std::setprecision(15);
      int size = stat.size();
      for (int i = 0; i < size; ++i) {
        fileCTRL << std::setw(25) << std::left << ctrl[i] << std::endl;
        fileVAR << std::setw(25) << std::left << stat[i] << std::endl;
      }
      fileCTRL.close();
      fileVAR.close();
    }
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
