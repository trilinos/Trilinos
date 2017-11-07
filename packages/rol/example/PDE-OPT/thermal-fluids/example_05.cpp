// @HEADER
// ************************************************************************
//
//               Rapid Optimization Library (ROL) Package
//                 Copyright (2014) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact lead developers:
//              Drew Kouri   (dpkouri@sandia.gov) and
//              Denis Ridzal (dridzal@sandia.gov)
//
// ************************************************************************
// @HEADER

/*! \file  example_01.cpp
    \brief Shows how to solve the Navier-Stokes control problem.
*/

#include "Teuchos_Comm.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "Tpetra_DefaultPlatform.hpp"
#include "Tpetra_Version.hpp"

#include <iostream>
#include <algorithm>
//#include <fenv.h>

#include "ROL_TpetraMultiVector.hpp"
#include "ROL_Algorithm.hpp"
#include "ROL_Reduced_Objective_SimOpt.hpp"
#include "ROL_MonteCarloGenerator.hpp"
#include "ROL_SimulatedConstraint.hpp"
#include "ROL_SimulatedObjectiveCVaR.hpp"
#include "ROL_SimulatedObjective.hpp"
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
void print(ROL::Objective<Real> &obj,
           const ROL::Vector<Real> &z,
           ROL::SampleGenerator<Real> &sampler,
           const int ngsamp,
           const Teuchos::RCP<const Teuchos::Comm<int> > &comm,
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
  Teuchos::RCP<const Teuchos::MpiComm<int> > mpicomm
    = Teuchos::rcp_dynamic_cast<const Teuchos::MpiComm<int> >(comm);
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
//  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  Teuchos::RCP<std::ostream> outStream;
  Teuchos::oblackholestream bhs; // outputs nothing

  /*** Initialize communicator. ***/
  Teuchos::GlobalMPISession mpiSession (&argc, &argv, &bhs);
  Teuchos::RCP<const Teuchos::Comm<int> > comm
    = Tpetra::DefaultPlatform::getDefaultPlatform().getComm();
  Teuchos::RCP<const Teuchos::Comm<int> > serial_comm
    = Teuchos::rcp(new Teuchos::SerialComm<int>());
  const int myRank = comm->getRank();
  if ((iprint > 0) && (myRank == 0)) {
    outStream = Teuchos::rcp(&std::cout, false);
  }
  else {
    outStream = Teuchos::rcp(&bhs, false);
  }
  int errorFlag  = 0;

  // *** Example body.
  try {

    /*** Read in XML input ***/
    std::string filename = "input_ex05.xml";
    Teuchos::RCP<Teuchos::ParameterList> parlist = Teuchos::rcp( new Teuchos::ParameterList() );
    Teuchos::updateParametersFromXmlFile( filename, parlist.ptr() );
    parlist->sublist("SimOpt").sublist("Solve").set("Output Iteration History",myRank==0);

    /*** Initialize main data structure. ***/
    Teuchos::RCP<MeshManager<RealT> > meshMgr
      = Teuchos::rcp(new MeshManager_ThermalFluids<RealT>(*parlist));
    // Initialize PDE describing Navier-Stokes equations.
    Teuchos::RCP<PDE_ThermalFluids_ex03<RealT> > pde
      = Teuchos::rcp(new PDE_ThermalFluids_ex03<RealT>(*parlist));
    Teuchos::RCP<ROL::Constraint_SimOpt<RealT> > con
      = Teuchos::rcp(new PDE_Constraint<RealT>(pde,meshMgr,serial_comm,*parlist,*outStream));
    // Cast the constraint and get the assembler.
    Teuchos::RCP<PDE_Constraint<RealT> > pdecon
      = Teuchos::rcp_dynamic_cast<PDE_Constraint<RealT> >(con);
    Teuchos::RCP<Assembler<RealT> > assembler = pdecon->getAssembler();
    con->setSolveParameters(*parlist);
    pdecon->outputTpetraData();

    // Create state vector and set to zeroes
    Teuchos::RCP<Tpetra::MultiVector<> > u_rcp, p_rcp, y_rcp, r_rcp, z_rcp, s_rcp, t_rcp;
    u_rcp  = assembler->createStateVector();     u_rcp->randomize();
    p_rcp  = assembler->createStateVector();     p_rcp->randomize();
    y_rcp  = assembler->createStateVector();     y_rcp->randomize();
    r_rcp  = assembler->createResidualVector();  r_rcp->randomize();
    z_rcp  = assembler->createControlVector();   z_rcp->putScalar(1.234); //z_rcp->randomize();
    s_rcp  = assembler->createControlVector();   s_rcp->putScalar(2.345); //s_rcp->randomize();
    t_rcp  = assembler->createControlVector();   t_rcp->putScalar(3.456); //t_rcp->randomize();
    Teuchos::RCP<ROL::Vector<RealT> > up, pp, yp, rp, zp, sp, tp;
    up  = Teuchos::rcp(new PDE_PrimalSimVector<RealT>(u_rcp,pde,assembler));
    pp  = Teuchos::rcp(new PDE_PrimalSimVector<RealT>(p_rcp,pde,assembler));
    yp  = Teuchos::rcp(new PDE_PrimalSimVector<RealT>(y_rcp,pde,assembler));
    rp  = Teuchos::rcp(new PDE_DualSimVector<RealT>(r_rcp,pde,assembler));
    zp  = Teuchos::rcp(new PDE_PrimalOptVector<RealT>(z_rcp,pde,assembler));
    sp  = Teuchos::rcp(new PDE_PrimalOptVector<RealT>(s_rcp,pde,assembler));
    tp  = Teuchos::rcp(new PDE_PrimalOptVector<RealT>(t_rcp,pde,assembler));

    // Initialize objective function.
    std::vector<Teuchos::RCP<QoI<RealT> > > qoi_vec(2,Teuchos::null);
    qoi_vec[0] = Teuchos::rcp(new QoI_State_ThermalFluids<RealT>(*parlist,
                                                                 pde->getVelocityFE(),
                                                                 pde->getPressureFE(),
                                                                 pde->getThermalFE(),
                                                                 pde->getFieldHelper()));
    qoi_vec[1] = Teuchos::rcp(new QoI_L2Penalty_ThermalFluids<RealT>(pde->getVelocityFE(),
                                                                     pde->getPressureFE(),
                                                                     pde->getThermalFE(),
                                                                     pde->getThermalBdryFE(),
                                                                     pde->getBdryCellLocIds(),
                                                                     pde->getFieldHelper()));
    Teuchos::RCP<StdObjective_ThermalFluids<RealT> > std_obj
      = Teuchos::rcp(new StdObjective_ThermalFluids<RealT>(*parlist));
    Teuchos::RCP<ROL::Objective_SimOpt<RealT> > obj
      = Teuchos::rcp(new PDE_Objective<RealT>(qoi_vec,std_obj,assembler));
    Teuchos::RCP<ROL::SimController<RealT> > stateStore
      = Teuchos::rcp(new ROL::SimController<RealT>());
    Teuchos::RCP<ROL::Reduced_Objective_SimOpt<RealT> > robj
      = Teuchos::rcp(new ROL::Reduced_Objective_SimOpt<RealT>(obj, con, stateStore, up, zp, pp, true, false));

    /*************************************************************************/
    /***************** BUILD SAMPLER *****************************************/
    /*************************************************************************/
    int Nbottom = parlist->sublist("Problem").get("Bottom KL Truncation Order",5);
    int Nleft   = parlist->sublist("Problem").get("Left KL Truncation Order",5);
    int Nright  = parlist->sublist("Problem").get("Right KL Truncation Order",5);
    int stochDim = Nbottom + Nleft + Nright + 3;
    int nsamp = parlist->sublist("Problem").get("Number of samples",100);
    // Build vector of distributions
    std::vector<Teuchos::RCP<ROL::Distribution<RealT> > > distVec(stochDim);
    Teuchos::ParameterList UList;
    UList.sublist("Distribution").set("Name","Uniform");
    UList.sublist("Distribution").sublist("Uniform").set("Lower Bound",-1.0);
    UList.sublist("Distribution").sublist("Uniform").set("Upper Bound", 1.0);
    for (int i = 0; i < stochDim; ++i) {
      distVec[i] = ROL::DistributionFactory<RealT>(UList);
    }
    // Sampler
    Teuchos::RCP<ROL::BatchManager<RealT> > bman
      = Teuchos::rcp(new ROL::TpetraTeuchosBatchManager<RealT>(comm));
    //  = Teuchos::rcp(new PDE_OptVector_BatchManager<RealT>(comm));
    Teuchos::RCP<ROL::SampleGenerator<RealT> > sampler
      = Teuchos::rcp(new ROL::MonteCarloGenerator<RealT>(nsamp,distVec,bman));

    /*************************************************************************/
    /***************** BUILD STOCHASTIC PROBLEM ******************************/
    /*************************************************************************/
    bool useW    = parlist->sublist("SOL").sublist("Simulated").get("Use Constraint Weights", true);
    bool useCVaR = parlist->sublist("SOL").sublist("Simulated").get("Use CVaR", false);
    Teuchos::RCP<ROL::Constraint<RealT> > simcon
      = Teuchos::rcp(new ROL::SimulatedConstraint<RealT>(sampler, con, useW));
    Teuchos::RCP<ROL::Objective<RealT> > simobj;
    if (useCVaR) {
      Teuchos::ParameterList list = parlist->sublist("SOL").sublist("Simulated");
      Teuchos::RCP<ROL::PlusFunction<RealT> > pf
        = Teuchos::rcp(new ROL::PlusFunction<RealT>(list));
      RealT alpha = parlist->sublist("SOL").sublist("Simulated").get("CVaR Confidence Level", 0.9);
      simobj = Teuchos::rcp(new ROL::SimulatedObjectiveCVaR<RealT>(sampler, obj, pf, alpha));
    }
    else {
      simobj = Teuchos::rcp(new ROL::SimulatedObjective<RealT>(sampler, obj));
    }
    std::vector<Teuchos::RCP<ROL::Vector<RealT> > > vuvec, vpvec, vyvec;
    for (int i = 0; i < sampler->numMySamples(); ++i) {
      Teuchos::RCP<Tpetra::MultiVector<> > vu_rcp, vp_rcp, vy_rcp;
      vu_rcp  = assembler->createStateVector(); vu_rcp->putScalar(4.567); //vu_rcp->randomize();
      vp_rcp  = assembler->createStateVector(); vp_rcp->putScalar(5.678); //vp_rcp->randomize();
      vy_rcp  = assembler->createStateVector(); vy_rcp->putScalar(6.789); //vy_rcp->randomize();
      Teuchos::RCP<ROL::Vector<RealT> > vup, vpp, vyp;
      vup  = Teuchos::rcp(new PDE_PrimalSimVector<RealT>(vu_rcp,pde,assembler));
      vpp  = Teuchos::rcp(new PDE_PrimalSimVector<RealT>(vp_rcp,pde,assembler));
      vyp  = Teuchos::rcp(new PDE_PrimalSimVector<RealT>(vy_rcp,pde,assembler));
      vuvec.push_back(vup);
      vpvec.push_back(vpp);
      vyvec.push_back(vyp);
    }
    Teuchos::RCP<ROL::Vector<RealT> > vu, vp, vy;
    vu = Teuchos::rcp(new ROL::SimulatedVector<RealT>(vuvec,bman));
    vp = Teuchos::rcp(new ROL::SimulatedVector<RealT>(vpvec,bman));
    vy = Teuchos::rcp(new ROL::SimulatedVector<RealT>(vyvec,bman));
    Teuchos::RCP<ROL::Vector<RealT> > rz, rs, rt;
    if (useCVaR) {
      Teuchos::RCP<Teuchos::ParameterList> cvarlist = Teuchos::rcp(new Teuchos::ParameterList);
      cvarlist->sublist("SOL").sublist("Risk Measure").set("Name","CVaR");
      rz = Teuchos::rcp(new ROL::RiskVector<RealT>(cvarlist, zp));
      rs = Teuchos::rcp(new ROL::RiskVector<RealT>(cvarlist, sp));
      rt = Teuchos::rcp(new ROL::RiskVector<RealT>(cvarlist, tp));
    }
    else {
      rz = zp;
      rs = sp;
      rt = tp;
    }
    ROL::Vector_SimOpt<RealT> x(vu,rz);
    ROL::Vector_SimOpt<RealT> p(vp,rs);
    ROL::Vector_SimOpt<RealT> y(vy,rt);
    x.checkVector(p,y,true,*outStream);

    bool derivCheck = parlist->sublist("Problem").get("Check derivatives",false);
    if (derivCheck) {
      *outStream << std::endl << "TESTING SimulatedConstraint" << std::endl;
      simcon->checkApplyJacobian(x, p, *vu, true, *outStream);
      simcon->checkAdjointConsistencyJacobian(*vu, p, x, *vu, x, true, *outStream);
      simcon->checkApplyAdjointHessian(x, *vu, p, x, true, *outStream);
      *outStream << std::endl << "TESTING SimulatedObjective" << std::endl;
      RealT tol = 1e-8;
      simobj->value(x, tol);
      simobj->checkGradient(x, p, true, *outStream);
      simobj->checkHessVec(x, p, true, *outStream);
    }

    zp->zero();
    Teuchos::RCP<ROL::SimulatedVector<RealT> > vusim
      = Teuchos::rcp_dynamic_cast<ROL::SimulatedVector<RealT> >(vu);
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
      x.zero();
      vp->zero();
    }

    /*************************************************************************/
    /***************** SOLVE PROBLEM *****************************************/
    /*************************************************************************/
    ROL::Algorithm<RealT> algo("Composite Step",*parlist,false);
    std::clock_t timer = std::clock();
    algo.run(x, *vp, *simobj, *simcon, true, *outStream);
    *outStream << "Optimization time: "
               << static_cast<RealT>(std::clock()-timer)/static_cast<RealT>(CLOCKS_PER_SEC)
               << " seconds." << std::endl << std::endl;
    
    /*************************************************************************/
    /***************** OUTPUT RESULTS ****************************************/
    /*************************************************************************/
    std::clock_t timer_print = std::clock();
    assembler->printMeshData(*outStream);
    // Output control to file
    pdecon->outputTpetraVector(z_rcp,"control.txt");
    // Output expected state and samples to file
    *outStream << std::endl << "Print Expected Value of State" << std::endl;
    up->zero(); pp->zero();
    for (int i = 0; i < sampler->numMySamples(); ++i) {
      up->axpy(sampler->getMyWeight(i),*(vusim->get(i)));
    }
    bman->sumAll(*up,*pp);
    pdecon->outputTpetraVector(p_rcp,"mean_state.txt");
    // Build full objective function distribution
    *outStream << std::endl << "Print Objective CDF" << std::endl;
    int nsamp_dist = parlist->sublist("Problem").get("Number of output samples",100);
    Teuchos::RCP<ROL::SampleGenerator<RealT> > sampler_dist
      = Teuchos::rcp(new ROL::MonteCarloGenerator<RealT>(nsamp_dist,distVec,bman));
    print<RealT>(*robj,*zp,*sampler_dist,nsamp_dist,comm,"obj_samples.txt");
    // Build vorticity objective function distribution
    Teuchos::RCP<ROL::Objective_SimOpt<RealT> > obj0
      = Teuchos::rcp(new IntegralObjective<RealT>(qoi_vec[0],assembler));
    Teuchos::RCP<ROL::SimController<RealT> > stateStore0
      = Teuchos::rcp(new ROL::SimController<RealT>());
    Teuchos::RCP<ROL::Reduced_Objective_SimOpt<RealT> > robj0
      = Teuchos::rcp(new ROL::Reduced_Objective_SimOpt<RealT>(obj0, con, stateStore0, up, zp, pp, true, false));
    print<RealT>(*robj0,*zp,*sampler_dist,nsamp_dist,comm,"vort_samples.txt");

    *outStream << "Output time: "
               << static_cast<RealT>(std::clock()-timer_print)/static_cast<RealT>(CLOCKS_PER_SEC)
               << " seconds." << std::endl << std::endl;
  }
  catch (std::logic_error err) {
    *outStream << err.what() << "\n";
    errorFlag = -1000;
  }; // end try

  if (errorFlag != 0)
    std::cout << "End Result: TEST FAILED\n";
  else
    std::cout << "End Result: TEST PASSED\n";

  return 0;
}
