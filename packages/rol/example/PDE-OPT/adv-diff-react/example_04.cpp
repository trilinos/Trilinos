// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  example_02.cpp
    \brief Shows how to solve the stochastic advection-diffusion problem.
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

#include "ROL_Bounds.hpp"
#include "ROL_Reduced_Objective_SimOpt.hpp"
#include "ROL_MonteCarloGenerator.hpp"
#include "ROL_StochasticProblem.hpp"
#include "ROL_Solver.hpp"
#include "ROL_TpetraTeuchosBatchManager.hpp"

#include "../TOOLS/meshmanager.hpp"
#include "../TOOLS/pdeconstraint.hpp"
#include "../TOOLS/pdeobjective.hpp"
#include "../TOOLS/pdevector.hpp"
#include "../TOOLS/batchmanager.hpp"
#include "pde_stoch_adv_diff.hpp"
#include "obj_stoch_adv_diff.hpp"
#include "mesh_stoch_adv_diff.hpp"

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
  Teuchos::GlobalMPISession mpiSession (&argc, &argv, &bhs);
  ROL::Ptr<const Teuchos::Comm<int>> comm
    = Tpetra::getDefaultComm();
  ROL::Ptr<const Teuchos::Comm<int>> serial_comm
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
    std::string filename = "input_ex04.xml";
    Teuchos::RCP<Teuchos::ParameterList> parlist = Teuchos::rcp( new Teuchos::ParameterList() );
    Teuchos::updateParametersFromXmlFile( filename, parlist.ptr() );

    // Problem dimensions
    const int controlDim = 9, stochDim = 37;
    const RealT one(1); 

    /*************************************************************************/
    /***************** BUILD GOVERNING PDE ***********************************/
    /*************************************************************************/
    /*** Initialize main data structure. ***/
    ROL::Ptr<MeshManager<RealT>> meshMgr
      = ROL::makePtr<MeshManager_stoch_adv_diff<RealT>>(*parlist);
    // Initialize PDE describing advection-diffusion equation
    ROL::Ptr<PDE_stoch_adv_diff<RealT>> pde
      = ROL::makePtr<PDE_stoch_adv_diff<RealT>>(*parlist);
    ROL::Ptr<ROL::Constraint_SimOpt<RealT>> con
      = ROL::makePtr<PDE_Constraint<RealT>>(pde,meshMgr,serial_comm,*parlist,*outStream);
    ROL::Ptr<PDE_Constraint<RealT>> pdeCon
      = ROL::dynamicPtrCast<PDE_Constraint<RealT>>(con);
    const ROL::Ptr<Assembler<RealT>> assembler = pdeCon->getAssembler();

    /*************************************************************************/
    /***************** BUILD VECTORS *****************************************/
    /*************************************************************************/
    ROL::Ptr<Tpetra::MultiVector<>> u_ptr, p_ptr, du_ptr, r_ptr;
    ROL::Ptr<std::vector<RealT>> z_ptr, dz_ptr, ez_ptr;
    ROL::Ptr<ROL::Vector<RealT>> up, pp, dup, rp, zp, dzp, ezp;
    // Create state vectors
    u_ptr  = assembler->createStateVector();
    p_ptr  = assembler->createStateVector();
    du_ptr = assembler->createStateVector();
    u_ptr->randomize();
    p_ptr->randomize();
    du_ptr->randomize();
    up  = ROL::makePtr<PDE_PrimalSimVector<RealT>>(u_ptr,pde,assembler,*parlist);
    pp  = ROL::makePtr<PDE_PrimalSimVector<RealT>>(p_ptr,pde,assembler,*parlist);
    dup = ROL::makePtr<PDE_PrimalSimVector<RealT>>(du_ptr,pde,assembler,*parlist);
    // Create residual vector
    r_ptr  = assembler->createResidualVector();
    rp  = ROL::makePtr<PDE_DualSimVector<RealT>>(r_ptr,pde,assembler,*parlist);
    // Create control vectors
    z_ptr  = ROL::makePtr<std::vector<RealT>>(controlDim);
    dz_ptr = ROL::makePtr<std::vector<RealT>>(controlDim);
    ez_ptr = ROL::makePtr<std::vector<RealT>>(controlDim);
    for (int i = 0; i < controlDim; ++i) {
      (*z_ptr)[i]  = random<RealT>(*comm);
      (*dz_ptr)[i] = random<RealT>(*comm);
      (*ez_ptr)[i] = random<RealT>(*comm);
    }
    zp  = ROL::makePtr<PDE_OptVector<RealT>>(ROL::makePtr<ROL::StdVector<RealT>>(z_ptr));
    dzp = ROL::makePtr<PDE_OptVector<RealT>>(ROL::makePtr<ROL::StdVector<RealT>>(dz_ptr));
    ezp = ROL::makePtr<PDE_OptVector<RealT>>(ROL::makePtr<ROL::StdVector<RealT>>(ez_ptr));

    /*************************************************************************/
    /***************** BUILD COST FUNCTIONAL *********************************/
    /*************************************************************************/
    std::vector<ROL::Ptr<QoI<RealT>>> qoi_vec(2,ROL::nullPtr);
    qoi_vec[0] = ROL::makePtr<QoI_State_Cost_stoch_adv_diff<RealT>>(pde->getFE());
    qoi_vec[1] = ROL::makePtr<QoI_Control_Cost_stoch_adv_diff<RealT>>();
    ROL::Ptr<StdObjective_stoch_adv_diff<RealT>> std_obj
      = ROL::makePtr<StdObjective_stoch_adv_diff<RealT>>(*parlist);
    ROL::Ptr<ROL::Objective_SimOpt<RealT>> obj
      = ROL::makePtr<PDE_Objective<RealT>>(qoi_vec,std_obj,assembler);
    ROL::Ptr<ROL::Reduced_Objective_SimOpt<RealT>> objReduced
      = ROL::makePtr<ROL::Reduced_Objective_SimOpt<RealT>>(obj, con, up, zp, pp, true, false);

    /*************************************************************************/
    /***************** BUILD BOUND CONSTRAINT ********************************/
    /*************************************************************************/
    ROL::Ptr<std::vector<RealT>> zlo_ptr = ROL::makePtr<std::vector<RealT>>(controlDim,0);
    ROL::Ptr<std::vector<RealT>> zhi_ptr = ROL::makePtr<std::vector<RealT>>(controlDim,1);
    ROL::Ptr<ROL::Vector<RealT>> zlop, zhip;
    zlop = ROL::makePtr<PDE_OptVector<RealT>>(ROL::makePtr<ROL::StdVector<RealT>>(zlo_ptr));
    zhip = ROL::makePtr<PDE_OptVector<RealT>>(ROL::makePtr<ROL::StdVector<RealT>>(zhi_ptr));
    ROL::Ptr<ROL::BoundConstraint<RealT>> bnd
      = ROL::makePtr<ROL::Bounds<RealT>>(zlop,zhip);

    /*************************************************************************/
    /***************** BUILD SAMPLER *****************************************/
    /*************************************************************************/
    int nsamp = parlist->sublist("Problem").get("Number of Samples",100);
    std::vector<RealT> tmp = {-one,one};
    std::vector<std::vector<RealT>> bounds(stochDim,tmp);
    ROL::Ptr<ROL::BatchManager<RealT>> bman
      = ROL::makePtr<PDE_OptVector_BatchManager<RealT>>(comm);
    ROL::Ptr<ROL::SampleGenerator<RealT>> sampler
      = ROL::makePtr<ROL::MonteCarloGenerator<RealT>>(nsamp,bounds,bman);

    /*************************************************************************/
    /***************** SOLVE OPTIMIZATION PROBLEM ****************************/
    /*************************************************************************/
    ROL::Ptr<ROL::StochasticProblem<RealT>> opt;
    ROL::Ptr<ROL::Solver<RealT>>  solver;

    const int nruns = 7;
    bool checkDeriv = parlist->sublist("Problem").get("Check Derivatives",false);
    bool outputSol  = parlist->sublist("Problem").get("Output Solution",false);

    // Run names
    std::vector<std::string> nvec = {"MV","RN","CVaR","MCVaR","ER","BP","KL"};

    // Parameter lists
    std::vector<Teuchos::ParameterList> plvec(7,*parlist);
    // Mean value
    plvec[0].sublist("SOL").sublist("Objective").set("Type", "Mean Value");
    // Risk neutral
    plvec[1].sublist("SOL").sublist("Objective").set("Type", "Risk Neutral");
    // CVaR
    plvec[2].sublist("SOL").sublist("Objective").set("Type", "Risk Averse");
    plvec[2].sublist("SOL").sublist("Objective").sublist("Risk Measure").set("Name","CVaR");
    plvec[2].sublist("SOL").sublist("Objective").sublist("Risk Measure").sublist("CVaR").set("Confidence Level", 0.95);
    plvec[2].sublist("SOL").sublist("Objective").sublist("Risk Measure").sublist("CVaR").set("Convex Combination Parameter", 1.0);
    plvec[2].sublist("SOL").sublist("Objective").sublist("Risk Measure").sublist("CVaR").set("Smoothing Parameter", 1e-4);
    plvec[2].sublist("SOL").sublist("Objective").sublist("Risk Measure").sublist("CVaR").sublist("Distribution").set("Name", "Parabolic");
    plvec[2].sublist("SOL").sublist("Objective").sublist("Risk Measure").sublist("CVaR").sublist("Distribution").sublist("Parabolic").set("Lower Bound", 0.0);
    plvec[2].sublist("SOL").sublist("Objective").sublist("Risk Measure").sublist("CVaR").sublist("Distribution").sublist("Parabolic").set("Upper Bound", 1.0);
    // Mixture of expectation and CVaR
    plvec[3].sublist("SOL").sublist("Objective").set("Type", "Risk Averse");
    plvec[3].sublist("SOL").sublist("Objective").sublist("Risk Measure").set("Name","CVaR");
    plvec[3].sublist("SOL").sublist("Objective").sublist("Risk Measure").sublist("CVaR").set("Confidence Level", 0.95);
    plvec[3].sublist("SOL").sublist("Objective").sublist("Risk Measure").sublist("CVaR").set("Convex Combination Parameter", 0.5);
    plvec[3].sublist("SOL").sublist("Objective").sublist("Risk Measure").sublist("CVaR").set("Smoothing Parameter", 1e-4);
    plvec[3].sublist("SOL").sublist("Objective").sublist("Risk Measure").sublist("CVaR").sublist("Distribution").set("Name", "Parabolic");
    plvec[3].sublist("SOL").sublist("Objective").sublist("Risk Measure").sublist("CVaR").sublist("Distribution").sublist("Parabolic").set("Lower Bound", 0.0);
    plvec[3].sublist("SOL").sublist("Objective").sublist("Risk Measure").sublist("CVaR").sublist("Distribution").sublist("Parabolic").set("Upper Bound", 1.0);
    // Entropic risk
    plvec[4].sublist("SOL").sublist("Objective").set("Type", "Risk Averse");
    plvec[4].sublist("SOL").sublist("Objective").sublist("Risk Measure").set("Name","Entropic Risk");
    plvec[4].sublist("SOL").sublist("Objective").sublist("Risk Measure").sublist("Entropic Risk").set("Rate", 1.0);
    // bPOE
    plvec[5].sublist("SOL").sublist("Objective").set("Type", "Probability");
    plvec[5].sublist("SOL").sublist("Objective").sublist("Probability").set("Name","bPOE");
    plvec[5].sublist("SOL").sublist("Objective").sublist("Probability").sublist("bPOE").set("Moment Order", 2.0);
    plvec[5].sublist("SOL").sublist("Objective").sublist("Probability").sublist("bPOE").set("Threshold", 6.0);
    // KL-divergence distributionally robust optimization
    plvec[6].sublist("SOL").sublist("Objective").set("Type", "Risk Averse");
    plvec[6].sublist("SOL").sublist("Objective").sublist("Risk Measure").set("Name","KL Divergence");
    plvec[6].sublist("SOL").sublist("Objective").sublist("Risk Measure").sublist("KL Divergence").set("Threshold", 0.1);
    
    RealT stat(1);
    for (int i = 0; i < nruns; ++i) {
      //zp->zero();

      // Build stochastic optimization problem
      opt = ROL::makePtr<ROL::StochasticProblem<RealT>>(objReduced,zp);
      opt->addBoundConstraint(bnd);
      plvec[i].sublist("SOL").sublist("Objective").set("Initial Statistic", stat);
      opt->makeObjectiveStochastic(plvec[i],sampler);
      opt->finalize(false,true,*outStream);
      if (checkDeriv) opt->check(true,*outStream);

      // Solve optimization problem
      plvec[i].sublist("Step").set("Type","Trust Region");
      solver = ROL::makePtr<ROL::Solver<RealT>>(opt,plvec[i]);
      std::clock_t timer = std::clock();
      solver->solve(*outStream);
      try                       { stat = opt->getSolutionStatistic(); }
      catch (std::exception &e) { stat = RealT(1);                    }
      *outStream << "Optimization time: "
                 << static_cast<RealT>(std::clock()-timer)/static_cast<RealT>(CLOCKS_PER_SEC)
                 << " seconds." << std::endl << std::endl;

      // Output control to file
      if ( myRank == 0 && outputSol ) {
        std::stringstream zname;
        zname << nvec[i] << "_control.txt";
        std::ofstream zfile;
        zfile.open(zname.str());
        for (int i = 0; i < controlDim; i++) {
          zfile << std::scientific << std::setprecision(15)
                << std::setw(25) << (*z_ptr)[i]
                << std::endl;
        }
        zfile.close();
      }

      // Build objective function distribution
      if (outputSol) {
        int nsamp_dist = plvec[i].sublist("Problem").get("Number of Output Samples",100);
        ROL::Ptr<ROL::SampleGenerator<RealT>> sampler_dist
          = ROL::makePtr<ROL::MonteCarloGenerator<RealT>>(nsamp_dist,bounds,bman);
        std::stringstream name;
        name << "obj_samples_" << nvec[i] << ".txt";
        print<RealT>(*objReduced,*zp,*sampler_dist,nsamp_dist,comm,name.str());
      }
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
