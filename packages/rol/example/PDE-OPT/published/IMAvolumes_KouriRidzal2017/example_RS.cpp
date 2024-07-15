// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  example_RS.cpp
    \brief Reduced-space solution of a thermal-fluids problem with random coefficients,
           using a risk-neutral formulation and adaptive sparse-grid sampling.
*/

#include "Teuchos_Comm.hpp"
#include "ROL_Stream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "Tpetra_Core.hpp"
#include "Tpetra_Version.hpp"

#include <iostream>
#include <algorithm>

#include "ROL_TpetraMultiVector.hpp"
#include "ROL_Reduced_Objective_SimOpt.hpp"
#include "ROL_MonteCarloGenerator.hpp"
#include "ROL_SparseGridGenerator.hpp"
#include "ROL_OptimizationSolver.hpp"
#include "ROL_TpetraTeuchosBatchManager.hpp"

#include "../../TOOLS/pdeconstraint.hpp"
#include "../../TOOLS/pdeobjective.hpp"
#include "../../TOOLS/pdevector.hpp"
#include "../../TOOLS/batchmanager.hpp"
#include "pde_thermal-fluids.hpp"
#include "obj_thermal-fluids.hpp"
#include "mesh_thermal-fluids.hpp"

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
    std::string filename = "input_RS.xml";
    Teuchos::RCP<Teuchos::ParameterList> parlist = Teuchos::rcp( new Teuchos::ParameterList() );
    Teuchos::updateParametersFromXmlFile( filename, parlist.ptr() );
    bool printSimOpt = parlist->sublist("SimOpt").sublist("Solve").get("Output Iteration History",false);
    parlist->sublist("SimOpt").sublist("Solve").set("Output Iteration History",printSimOpt&&(myRank==0));
    int verbosity = parlist->sublist("General").get("Print Verbosity",0);
    verbosity = (myRank != 0) ? 0 : verbosity;
    parlist->sublist("General").set("Print Verbosity",verbosity);

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
    int maxLevel   = parlist->sublist("Problem").get("Maximum Sparse Grid Level",7);
    bool adaptiveV = parlist->sublist("Problem").get("Use Value Adaptive Sparse Grids",false);
    bool adaptiveG = parlist->sublist("Problem").get("Use Gradient Adaptive Sparse Grids",true);
    bool printSG   = parlist->sublist("Problem").get("Print Sparse Grid Size",false);
    ROL::QuadratureInfo info, vinfo, ginfo;
    info.dim        = stochDim;
    info.maxLevel   = maxLevel;
    info.normalized = true;
    info.adaptive   = false;
    info.print      = (printSG&&(myRank==0));
    info.name       = "Full";
    info.rule1D.clear();   info.rule1D.resize(info.dim,ROL::QUAD_CLENSHAWCURTIS);
    info.growth1D.clear(); info.growth1D.resize(info.dim,ROL::GROWTH_DEFAULT);
    vinfo           = info;
    vinfo.adaptive  = adaptiveV;
    vinfo.name      = "Value";
    ginfo           = info;
    ginfo.adaptive  = adaptiveG;
    ginfo.name      = "Gradient";
    ROL::Ptr<ROL::BatchManager<RealT> > bman
      = ROL::makePtr<ROL::TpetraTeuchosBatchManager<RealT>>(comm);
    //  = ROL::makePtr<PDE_OptVector_BatchManager<RealT>>(comm);
    ROL::Ptr<ROL::SampleGenerator<RealT> > sampler
      = ROL::makePtr<ROL::SparseGridGenerator<RealT>>(bman,info,false);
    ROL::Ptr<ROL::SampleGenerator<RealT> > vsampler
      = ROL::makePtr<ROL::SparseGridGenerator<RealT>>(bman,vinfo,adaptiveV);
    ROL::Ptr<ROL::SampleGenerator<RealT> > gsampler
      = ROL::makePtr<ROL::SparseGridGenerator<RealT>>(bman,ginfo,adaptiveG);
    /*************************************************************************/
    /***************** BUILD STOCHASTIC PROBLEM ******************************/
    /*************************************************************************/
    ROL::OptimizationProblem<RealT> opt(robj,zp);
    parlist->sublist("SOL").set("Initial Statistic", static_cast<RealT>(1));
    opt.setStochasticObjective(*parlist,vsampler,gsampler);

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

    up->zero();
    zp->zero();
    parlist->sublist("Step").set("Type","Trust Region");
    ROL::OptimizationSolver<RealT> solver(opt,*parlist);
    std::clock_t timer = std::clock();
    solver.solve(*outStream);
    *outStream << "Optimization time: "
               << static_cast<RealT>(std::clock()-timer)/static_cast<RealT>(CLOCKS_PER_SEC)
               << " seconds." << std::endl << std::endl;

    parlist->sublist("General").set("Inexact Objective Function", false);
    parlist->sublist("General").set("Inexact Gradient", false);
    ROL::OptimizationProblem<RealT> optFull(robj,zp);
    parlist->sublist("SOL").set("Initial Statistic", static_cast<RealT>(1));
    opt.setStochasticObjective(*parlist,sampler);
    ROL::Algorithm<RealT> algoFull("Trust Region",*parlist,false);
    std::clock_t timerFull = std::clock();
    algoFull.run(optFull,true,*outStream);
    *outStream << "Optimization time: "
               << static_cast<RealT>(std::clock()-timerFull)/static_cast<RealT>(CLOCKS_PER_SEC)
               << " seconds." << std::endl << std::endl;
    
    /*************************************************************************/
    /***************** OUTPUT RESULTS ****************************************/
    /*************************************************************************/
    vsampler->print("Value_samples");
    gsampler->print("Gradient_samples");
    ROL::dynamicPtrCast<ROL::SparseGridGenerator<RealT> >(vsampler)->printIndexSet();
    ROL::dynamicPtrCast<ROL::SparseGridGenerator<RealT> >(gsampler)->printIndexSet();
    std::clock_t timer_print = std::clock();
    assembler->printMeshData(*outStream);
    // Output control to file
    pdecon->outputTpetraVector(z_ptr,"control.txt");
    // Output expected state and samples to file
    *outStream << std::endl << "Print Expected Value of State" << std::endl;
    up->zero(); pp->zero(); dup->zero(); dzp->zero(); yup->zero(); dyp->zero();
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
      con->solve(*rp,*dyp,*dzp,tol);
      yup->axpy(sampler->getMyWeight(i),*dyp);
      for (int j = 0; j < stochDim; ++j) {
        file_samp << std::setw(25) << std::left << sample[j];
      }
      file_samp << std::endl;
    }
    file_samp.close();
    bman_Eu->sumAll(*up,*pp);
    bman_Eu->sumAll(*yup,*ypp);
    pdecon->outputTpetraVector(p_ptr,"mean_state.txt");
    pdecon->outputTpetraVector(yp_ptr,"mean_uncontrolled_state.txt");
    // Build full objective function distribution
    *outStream << std::endl << "Print Objective CDF" << std::endl;
    int nsamp_dist = parlist->sublist("Problem").get("Number of output samples",100);
    ROL::Ptr<ROL::SampleGenerator<RealT> > sampler_dist
      = ROL::makePtr<ROL::MonteCarloGenerator<RealT>>(nsamp_dist,distVec,bman);
    print<RealT>(*robj,*zp,*sampler_dist,nsamp_dist,comm,"obj_samples.txt");
    // Build vorticity objective function distribution
    ROL::Ptr<ROL::Objective_SimOpt<RealT> > obj0
      = ROL::makePtr<IntegralObjective<RealT>>(qoi_vec[0],assembler);
    ROL::Ptr<ROL::VectorController<RealT> > stateStore0
      = ROL::makePtr<ROL::VectorController<RealT>>();
    ROL::Ptr<ROL::Reduced_Objective_SimOpt<RealT> > robj0
      = ROL::makePtr<ROL::Reduced_Objective_SimOpt<RealT>>(obj0, con, stateStore0, up, zp, pp, true, false);
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
