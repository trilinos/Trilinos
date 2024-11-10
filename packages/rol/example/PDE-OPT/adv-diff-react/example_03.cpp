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

#include "ROL_Solver.hpp"
#include "ROL_Bounds.hpp"
#include "ROL_Reduced_Objective_SimOpt.hpp"
#include "ROL_MonteCarloGenerator.hpp"
#include "ROL_StochasticProblem.hpp"
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
    std::string filename = "input_ex03.xml";
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
    // Create state vectors
    ROL::Ptr<Tpetra::MultiVector<>>  u_ptr = assembler->createStateVector();
    ROL::Ptr<Tpetra::MultiVector<>>  p_ptr = assembler->createStateVector();
    ROL::Ptr<Tpetra::MultiVector<>> du_ptr = assembler->createStateVector();
    ROL::Ptr<Tpetra::MultiVector<>>  r_ptr = assembler->createResidualVector();
    u_ptr->randomize();  //u_ptr->putScalar(static_cast<RealT>(1));
    p_ptr->randomize();  //p_ptr->putScalar(static_cast<RealT>(1));
    du_ptr->randomize(); //du_ptr->putScalar(static_cast<RealT>(0));
    ROL::Ptr<ROL::Vector<RealT>> up, pp, dup, rp;
    up  = ROL::makePtr<PDE_PrimalSimVector<RealT>>(u_ptr,pde,assembler,*parlist);
    pp  = ROL::makePtr<PDE_PrimalSimVector<RealT>>(p_ptr,pde,assembler,*parlist);
    dup = ROL::makePtr<PDE_PrimalSimVector<RealT>>(du_ptr,pde,assembler,*parlist);
    rp  = ROL::makePtr<PDE_DualSimVector<RealT>>(r_ptr,pde,assembler,*parlist);
    // Create control vectors
    ROL::Ptr<std::vector<RealT>>  z_ptr = ROL::makePtr<std::vector<RealT>>(controlDim);
    ROL::Ptr<std::vector<RealT>> dz_ptr = ROL::makePtr<std::vector<RealT>>(controlDim);
    ROL::Ptr<std::vector<RealT>> ez_ptr = ROL::makePtr<std::vector<RealT>>(controlDim);
    for (int i = 0; i < controlDim; ++i) {
      (*z_ptr)[i]  = random<RealT>(*comm);
      (*dz_ptr)[i] = random<RealT>(*comm);
      (*ez_ptr)[i] = random<RealT>(*comm);
    }
    ROL::Ptr<ROL::Vector<RealT>> zp, dzp, ezp;
    zp  = ROL::makePtr<PDE_OptVector<RealT>(ROL::makePtr<ROL::StdVector<RealT>>>(z_ptr));
    dzp = ROL::makePtr<PDE_OptVector<RealT>(ROL::makePtr<ROL::StdVector<RealT>>>(dz_ptr));
    ezp = ROL::makePtr<PDE_OptVector<RealT>(ROL::makePtr<ROL::StdVector<RealT>>>(ez_ptr));

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
      = ROL::makePtr<ROL::Reduced_Objective_SimOpt<RealT>>(obj, con, up, pp, true, false);

    /*************************************************************************/
    /***************** BUILD BOUND CONSTRAINT ********************************/
    /*************************************************************************/
    ROL::Ptr<std::vector<RealT>> zlo_ptr = ROL::makePtr<std::vector<RealT>>(controlDim,0);
    ROL::Ptr<std::vector<RealT>> zhi_ptr = ROL::makePtr<std::vector<RealT>>(controlDim,1);
    ROL::Ptr<ROL::Vector<RealT>> zlop, zhip;
    zlop = ROL::makePtr<PDE_OptVector<RealT>(ROL::makePtr<ROL::StdVector<RealT>>>(zlo_ptr));
    zhip = ROL::makePtr<PDE_OptVector<RealT>(ROL::makePtr<ROL::StdVector<RealT>>>(zhi_ptr));
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
    ROL::Ptr<ROL::Solver<RealT>> solver;
    zp->zero();

    int nQuad = 11, nSmooth = 1, N(2);
    RealT eps(1), eps0(1.e-3), stat(0);
    std::vector<RealT> statVec(nQuad);
    std::vector<RealT>  errVec(nQuad);
    std::vector<RealT> normVec(nQuad);
    bool checkDeriv = parlist->sublist("Problem").get("Check Derivatives",false);

    std::string rm = "Risk Measure", sr = "Spectral Risk";
    std::string dist = "Distribution", pf = "Plus Function";
    parlist->sublist("SOL").set("Type","Risk Averse");
    parlist->sublist("SOL").sublist("Objective").sublist(rm).set("Name",sr);
    parlist->sublist("SOL").sublist("Objective").sublist(rm).sublist(sr).set("Print Quadrature to Screen",!myRank);
    parlist->sublist("SOL").sublist("Objective").sublist(rm).sublist(sr).sublist(dist).set("Name","Beta");
    parlist->sublist("SOL").sublist("Objective").sublist(rm).sublist(sr).sublist(dist).sublist("Beta").set("Shape 1",5.0);
    parlist->sublist("SOL").sublist("Objective").sublist(rm).sublist(sr).sublist(dist).sublist("Beta").set("Shape 2",2.0);
    parlist->sublist("SOL").sublist("Objective").sublist(rm).sublist(sr).sublist(pf).sublist(dist).set("Name","Parabolic");
    parlist->sublist("SOL").sublist("Objective").sublist(rm).sublist(sr).sublist(pf).sublist(dist).sublist("Parabolic").set("Lower Bound",-0.5);
    parlist->sublist("SOL").sublist("Objective").sublist(rm).sublist(sr).sublist(pf).sublist(dist).sublist("Parabolic").set("Upper Bound", 0.5);
    
    for (int i = 0; i < nQuad; ++i) {
      eps = eps0;
      parlist->sublist("SOL").sublist("Objective").sublist(rm).sublist(sr).set("Number of Quadrature Points",N);
      for (int j = 0; j < nSmooth; ++j) {
        parlist->sublist("SOL").sublist("Objective").sublist(rm).sublist(sr).sublist(pf).set("Smoothing Parameter",eps);
        // Build stochastic optimization problem
        opt = ROL::makePtr<ROL::StochasticProblem<RealT>>(objReduced,zp);
        opt->addBoundConstraint(bnd);
        parlist->sublist("SOL").sublist("Objective").set("Initial Statisitic", stat);
        opt->makeObjectiveStochastic(*parlist,sampler);
        opt->finalize(false,true,*outStream);
        if (checkDeriv) opt->check(true,*outStream);
        // Solve optimization problem
        solver = ROL::makePtr<ROL::Solver<RealT>>(opt,*parlist);
        std::clock_t timer = std::clock();
        solver->solve(*outStream);
        *outStream << "Optimization time: "
                   << static_cast<RealT>(std::clock()-timer)/static_cast<RealT>(CLOCKS_PER_SEC)
                   << " seconds." << std::endl << std::endl;
        stat = opt->getSolutionStatistic();
        // Print control and statistic values to screen
        *outStream << std::endl << std::endl;
        *outStream << std::scientific << std::setprecision(15);
        *outStream << "N = " << N << ", eps = " << eps << std::endl;
        *outStream << "Control:" << std::endl;
        for (int k = 0; k < controlDim; ++k) {
          *outStream << std::scientific << std::setprecision(15)
                     << std::setw(25) << (*z_ptr)[k]
                     << std::endl; 
        }
        *outStream << "Statistic: " << std::endl;
        *outStream << std::scientific << std::setprecision(15)
                   << std::setw(25) << stat
                   << std::endl << std::endl;
        // Update smoothing parameters
        eps *= static_cast<RealT>(0.1);
      }
      // Update number of quadrature points
      N *= 2;
      // Store control errors, control norms and statistic values
      ezp->set(*zp); ezp->axpy(-one,*dzp);
      normVec[i] = zp->norm();
      errVec[i]  = ezp->norm();
      statVec[i] = stat;
      // Store previous control
      dzp->set(*zp);
    }

    *outStream << std::endl;
    *outStream << std::setw(25) << std::left << "Control Error"
               << std::setw(25) << std::left << "Control Norm"
               << std::setw(25) << std::left << "Statistic"
               << std::endl;
    for (int i = 0; i < nQuad; ++i) {
      *outStream << std::scientific << std::setprecision(15)
                 << std::setw(25) << std::left <<  errVec[i]
                 << std::setw(25) << std::left << normVec[i]
                 << std::setw(25) << std::left << statVec[i]
                 << std::endl;
    }

    /*************************************************************************/
    /***************** OUTPUT RESULTS ****************************************/
    /*************************************************************************/
    std::clock_t timer_print = std::clock();
    assembler->printMeshData(*outStream);
    // Output control to file
    if ( myRank == 0 ) {
      std::ofstream zfile;
      zfile.open("control.txt");
      for (int i = 0; i < controlDim; i++) {
        zfile << (*z_ptr)[i] << "\n";
      }
      zfile.close();
    }
    // Output expected state and samples to file
    up->zero(); pp->zero(); dup->zero();
    RealT tol(1.e-8);
    ROL::Ptr<ROL::BatchManager<RealT>> bman_Eu
      = ROL::makePtr<ROL::TpetraTeuchosBatchManager<RealT>>(comm);
    std::vector<RealT> sample(stochDim);
    std::stringstream name_samp;
    name_samp << "samples_" << bman->batchID() << ".txt";
    std::ofstream file_samp;
    file_samp.open(name_samp.str());
    for (int i = 0; i < sampler->numMySamples(); ++i) {
      sample = sampler->getMyPoint(i);
      con->setParameter(sample);
      con->solve(*rp,*dup,*zp,tol);
      up->axpy(sampler->getMyWeight(i),*dup);
      for (int j = 0; j < stochDim; ++j) {
        file_samp << sample[j] << "  ";
      }
      file_samp << "\n";
    }
    file_samp.close();
    bman_Eu->sumAll(*up,*pp);
    assembler->outputTpetraVector(p_ptr,"mean_state.txt");
    // Build objective function distribution
    RealT val(0);
    int nsamp_dist = parlist->sublist("Problem").get("Number of Output Samples",100);
    ROL::Ptr<ROL::SampleGenerator<RealT>> sampler_dist
      = ROL::makePtr<ROL::MonteCarloGenerator<RealT>>(nsamp_dist,bounds,bman);
    std::stringstream name;
    name << "obj_samples_" << bman->batchID() << ".txt";
    std::ofstream file;
    file.open(name.str());
    for (int i = 0; i < sampler_dist->numMySamples(); ++i) {
      sample = sampler_dist->getMyPoint(i);
      objReduced->setParameter(sample);
      val = objReduced->value(*zp,tol);
      for (int j = 0; j < stochDim; ++j) {
        file << sample[j] << "  ";
      }
      file << val << "\n";
    }
    file.close();
    *outStream << "Output time: "
               << static_cast<RealT>(std::clock()-timer_print)/static_cast<RealT>(CLOCKS_PER_SEC)
               << " seconds." << std::endl << std::endl;

    Teuchos::Array<RealT> res(1,0);
    pdeCon->value(*rp,*up,*zp,tol);
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
