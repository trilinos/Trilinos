// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  example_03.cpp
    \brief Shows how to solve the stochastic Navier-Stokes problem.
*/

#include "Teuchos_Comm.hpp"
#ifdef HAVE_MPI
#include "Teuchos_DefaultMpiComm.hpp"
#endif
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

template<class Real>
void setUpAndSolve(const ROL::Ptr<ROL::Problem<Real>> &opt,
                   Teuchos::ParameterList &parlist,
                   std::ostream &outStream) {
  parlist.sublist("Step").set("Type","Trust Region");
  ROL::Solver<Real> solver(opt,parlist);
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
           const ROL::Ptr<const Teuchos::Comm<int>> &comm,
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
    for (int j = 0; j < sdim; ++j) {
      mysamples[j][i] = static_cast<double>(sample[j]);
    }
  }

  // Send data to root processor
#ifdef HAVE_MPI
  ROL::Ptr<const Teuchos::MpiComm<int>> mpicomm
    = ROL::dynamicPtrCast<const Teuchos::MpiComm<int>>(comm);
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
  ROL::GlobalMPISession mpiSession (&argc, &argv, &bhs);
  Kokkos::ScopeGuard kokkosScope (argc, argv);
  ROL::Ptr<const Teuchos::Comm<int>> comm
    = Tpetra::getDefaultComm();
  ROL::Ptr<const Teuchos::Comm<int>> serial_comm
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
    std::string filename = "input_ex03.xml";
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
    con = ROL::makePtr<PDE_Constraint<RealT,DeviceT>>(pde,meshMgr,serial_comm,*parlist,*outStream);
    // Cast the constraint and get the assembler.
    auto assembler = pdecon->getAssembler();
    assembler->printMeshData(*outStream);
    con->setSolveParameters(*parlist);

    /*************************************************************************/
    /***************** BUILD VECTORS *****************************************/
    /*************************************************************************/
    auto  u_ptr = assembler->createStateVector();
    auto  p_ptr = assembler->createStateVector();
    auto du_ptr = assembler->createStateVector();
    u_ptr->randomize();  //u_ptr->putScalar(static_cast<RealT>(1));
    p_ptr->randomize();  //p_ptr->putScalar(static_cast<RealT>(1));
    du_ptr->randomize(); //du_ptr->putScalar(static_cast<RealT>(0));
    auto up  = ROL::makePtr<PDE_PrimalSimVector<RealT,DeviceT>>(u_ptr,pde,assembler,*parlist);
    auto pp  = ROL::makePtr<PDE_PrimalSimVector<RealT,DeviceT>>(p_ptr,pde,assembler,*parlist);
    auto dup = ROL::makePtr<PDE_PrimalSimVector<RealT,DeviceT>>(du_ptr,pde,assembler,*parlist);
    // Create residual vectors
    auto r_ptr = assembler->createResidualVector();
    r_ptr->randomize(); //r_ptr->putScalar(static_cast<RealT>(1));
    auto rp = ROL::makePtr<PDE_DualSimVector<RealT,DeviceT>>(r_ptr,pde,assembler,*parlist);
    // Create control vector and set to ones
    auto  z_ptr = assembler->createControlVector();
    auto dz_ptr = assembler->createControlVector();
    auto yz_ptr = assembler->createControlVector();
    z_ptr->randomize();  z_ptr->putScalar(static_cast<RealT>(0));
    dz_ptr->randomize(); //dz_ptr->putScalar(static_cast<RealT>(0));
    yz_ptr->randomize(); //yz_ptr->putScalar(static_cast<RealT>(0));
    ROL::Ptr<ROL::TpetraMultiVector<RealT>> zpde
      = ROL::makePtr<PDE_PrimalOptVector<RealT,DeviceT>>(z_ptr,pde,assembler,*parlist);
    ROL::Ptr<ROL::TpetraMultiVector<RealT>> dzpde
      = ROL::makePtr<PDE_PrimalOptVector<RealT,DeviceT>>(dz_ptr,pde,assembler,*parlist);
    ROL::Ptr<ROL::TpetraMultiVector<RealT>> yzpde
      = ROL::makePtr<PDE_PrimalOptVector<RealT,DeviceT>>(yz_ptr,pde,assembler,*parlist);
    auto zp  = ROL::makePtr<PDE_OptVector<RealT>>(zpde);
    auto dzp = ROL::makePtr<PDE_OptVector<RealT>>(dzpde);
    auto yzp = ROL::makePtr<PDE_OptVector<RealT>>(yzpde);
    // Create ROL SimOpt vectors
    ROL::Vector_SimOpt<RealT> x(up,zp);
    ROL::Vector_SimOpt<RealT> d(dup,dzp);

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
    auto objState = ROL::makePtr<PDE_Objective<RealT,DeviceT>>(qoi_vec[0],assembler);
    auto objCtrl = ROL::makePtr<PDE_Objective<RealT,DeviceT>>(qoi_vec[1],assembler);
    auto objRed = ROL::makePtr<ROL::Reduced_Objective_SimOpt<RealT>>(obj, con, up, zp, pp, true, false);

    /*************************************************************************/
    /***************** BUILD BOUND CONSTRAINT ********************************/
    /*************************************************************************/
    auto zlo_ptr = assembler->createControlVector();
    auto zhi_ptr = assembler->createControlVector();
    zlo_ptr->putScalar(static_cast<RealT>(0));
    zhi_ptr->putScalar(ROL::ROL_INF<RealT>());
    ROL::Ptr<ROL::TpetraMultiVector<RealT>> zlopde
      = ROL::makePtr<PDE_PrimalOptVector<RealT,DeviceT>>(zlo_ptr,pde,assembler,*parlist);
    ROL::Ptr<ROL::TpetraMultiVector<RealT>> zhipde
      = ROL::makePtr<PDE_PrimalOptVector<RealT,DeviceT>>(zhi_ptr,pde,assembler,*parlist);
    auto zlop = ROL::makePtr<PDE_OptVector<RealT>>(zlopde);
    auto zhip = ROL::makePtr<PDE_OptVector<RealT>>(zhipde);
    auto bnd = ROL::makePtr<ROL::Bounds<RealT>>(zlop,zhip);
    bool useBounds = parlist->sublist("Problem").get("Use bounds", false);
    if (!useBounds) bnd->deactivate();

    /*************************************************************************/
    /***************** BUILD SAMPLER *****************************************/
    /*************************************************************************/
    int nsamp = parlist->sublist("Problem").get("Number of samples",100);
    int nsamp_dist = parlist->sublist("Problem").get("Number of output samples",100);
    std::vector<RealT> tmp = {-one,one};
    std::vector<std::vector<RealT>> bounds(stochDim,tmp);
    auto bman = ROL::makePtr<PDE_OptVector_BatchManager<RealT>>(comm);
    auto sampler = ROL::makePtr<ROL::MonteCarloGenerator<RealT>>(nsamp,bounds,bman);
    auto sampler_dist = ROL::makePtr<ROL::MonteCarloGenerator<RealT>>(nsamp_dist,bounds,bman);

    /*************************************************************************/
    /***************** BUILD STOCHASTIC PROBLEM ******************************/
    /*************************************************************************/
    ROL::Ptr<ROL::StochasticProblem<RealT>> opt;
    std::vector<RealT> ctrl;
    std::vector<RealT> var;

    std::vector<RealT> alphaArray
      = ROL::getArrayFromStringParameter<RealT>(parlist->sublist("Problem"),"Confidence Levels");
    std::vector<RealT> alpha = alphaArray.toVector();
    std::sort(alpha.begin(),alpha.end());
    int N = alpha.size();

    /*************************************************************************/
    /***************** SOLVE RISK NEUTRAL ************************************/
    /*************************************************************************/
    RealT tol(1e-8);
    bool alphaZero = (alpha[0] == static_cast<RealT>(0));
    if ( alphaZero ) {
      alpha.erase(alpha.begin()); --N;
      // Solve.
      parlist->sublist("SOL").set("Type","Risk Neutral");
      opt = ROL::makePtr<ROL::StochasticProblem<RealT>>(objRed,zp);
      opt->addBoundConstraint(bnd);
      parlist->sublist("SOL").set("Initial Statistic",one);
      opt->makeObjectiveStochastic(*parlist,sampler);
      opt->finalize(false,true,*outStream);
      setUpAndSolve<RealT>(opt,*parlist,*outStream);
      // Output.
      ctrl.push_back(objCtrl->value(*up,*zp,tol));
      try                       { var.push_back(opt->getSolutionStatistic()); }
      catch (std::exception &e) { var.push_back(0.0); }
      std::string CtrlRN = "control_RN.txt";
      con->outputTpetraVector(z_ptr,CtrlRN);
      std::string ObjRN = "obj_samples_RN.txt";
      print<RealT>(*objRed,*zp,*sampler_dist,nsamp_dist,comm,ObjRN);
    }

    /*************************************************************************/
    /***************** SOLVE MEAN PLUS CVAR **********************************/
    /*************************************************************************/
    parlist->sublist("SOL").set("Type","Risk Averse");
    parlist->sublist("SOL").sublist("Risk Measure").set("Name","CVaR");
    parlist->sublist("SOL").sublist("Risk Measure").sublist("CVaR").set("Convex Combination Parameter",1.0);
    parlist->sublist("SOL").sublist("Risk Measure").sublist("CVaR").set("Smoothing Parameter",1e-4);
    parlist->sublist("SOL").sublist("Risk Measure").sublist("CVaR").sublist("Distribution").set("Name","Parabolic");
    parlist->sublist("SOL").sublist("Risk Measure").sublist("CVaR").sublist("Distribution").sublist("Parabolic").set("Lower Bound",0.0);
    parlist->sublist("SOL").sublist("Risk Measure").sublist("CVaR").sublist("Distribution").sublist("Parabolic").set("Upper Bound",1.0);
    for (int i = 0; i < N; ++i) {
      // Solve.
      parlist->sublist("SOL").sublist("Risk Measure").sublist("CVaR").set("Confidence Level",alpha[i]);
      opt = ROL::makePtr<ROL::StochasticProblem<RealT>>(objRed,zp);
      opt->addBoundConstraint(bnd);
      parlist->sublist("SOL").set("Initial Statistic",one);
      opt->makeObjectiveStochastic(*parlist,sampler);
      opt->finalize(false,true,*outStream);
      setUpAndSolve<RealT>(opt,*parlist,*outStream);
      // Output.
      ctrl.push_back(objCtrl->value(*up,*zp,tol));
      try                       { var.push_back(opt->getSolutionStatistic()); }
      catch (std::exception &e) { var.push_back(0.0); }
      std::stringstream nameCtrl;
      nameCtrl << "control_CVaR_" << i+1 << ".txt";
      con->outputTpetraVector(z_ptr,nameCtrl.str().c_str());
      std::stringstream nameObj;
      nameObj << "obj_samples_CVaR_" << i+1 << ".txt";
      print<RealT>(*objRed,*zp,*sampler_dist,nsamp_dist,comm,nameObj.str());
    }

    /*************************************************************************/
    /***************** PRINT CONTROL OBJ AND VAR *****************************/
    /*************************************************************************/
    const int rank = Teuchos::rank<int>(*comm);
    if ( rank==0 ) {
      std::stringstream nameCTRL, nameVAR;
      nameCTRL << "ctrl.txt";
      nameVAR << "var.txt";
      std::ofstream fileCTRL, fileVAR;
      fileCTRL.open(nameCTRL.str());
      fileVAR.open(nameVAR.str());
      fileCTRL << std::scientific << std::setprecision(15);
      fileVAR << std::scientific << std::setprecision(15);
      int size = var.size();
      for (int i = 0; i < size; ++i) {
        fileCTRL << std::setw(25) << std::left << ctrl[i] << std::endl;
        fileVAR << std::setw(25) << std::left << var[i] << std::endl;
      }
      fileCTRL.close();
      fileVAR.close();
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
