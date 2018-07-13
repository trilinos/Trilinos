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
#include "ROL_Stream.hpp"
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
//  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  ROL::Ptr<std::ostream> outStream;
  ROL::nullstream bhs; // outputs nothing

  /*** Initialize communicator. ***/
  Teuchos::GlobalMPISession mpiSession (&argc, &argv, &bhs);
  ROL::Ptr<const Teuchos::Comm<int> > comm
    = Tpetra::DefaultPlatform::getDefaultPlatform().getComm();
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
    std::string filename = "input_ex05.xml";
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
    ROL::Ptr<Tpetra::MultiVector<> > u_ptr, p_ptr, y_ptr, r_ptr, z_ptr, s_ptr, t_ptr;
    u_ptr  = assembler->createStateVector();     u_ptr->randomize();
    p_ptr  = assembler->createStateVector();     p_ptr->randomize();
    y_ptr  = assembler->createStateVector();     y_ptr->randomize();
    r_ptr  = assembler->createResidualVector();  r_ptr->randomize();
    z_ptr  = assembler->createControlVector();   z_ptr->putScalar(1.234); //z_ptr->randomize();
    s_ptr  = assembler->createControlVector();   s_ptr->putScalar(2.345); //s_ptr->randomize();
    t_ptr  = assembler->createControlVector();   t_ptr->putScalar(3.456); //t_ptr->randomize();
    ROL::Ptr<ROL::Vector<RealT> > up, pp, yp, rp, zp, sp, tp;
    up  = ROL::makePtr<PDE_PrimalSimVector<RealT>>(u_ptr,pde,assembler);
    pp  = ROL::makePtr<PDE_PrimalSimVector<RealT>>(p_ptr,pde,assembler);
    yp  = ROL::makePtr<PDE_PrimalSimVector<RealT>>(y_ptr,pde,assembler);
    rp  = ROL::makePtr<PDE_DualSimVector<RealT>>(r_ptr,pde,assembler);
    zp  = ROL::makePtr<PDE_PrimalOptVector<RealT>>(z_ptr,pde,assembler);
    sp  = ROL::makePtr<PDE_PrimalOptVector<RealT>>(s_ptr,pde,assembler);
    tp  = ROL::makePtr<PDE_PrimalOptVector<RealT>>(t_ptr,pde,assembler);

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
    ROL::Ptr<ROL::SimController<RealT> > stateStore
      = ROL::makePtr<ROL::SimController<RealT>>();
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

    /*************************************************************************/
    /***************** BUILD STOCHASTIC PROBLEM ******************************/
    /*************************************************************************/
    bool useW    = parlist->sublist("SOL").sublist("Simulated").get("Use Constraint Weights", true);
    bool useCVaR = parlist->sublist("SOL").sublist("Simulated").get("Use CVaR", false);
    ROL::Ptr<ROL::Constraint<RealT> > simcon
      = ROL::makePtr<ROL::SimulatedConstraint<RealT>>(sampler, con, useW);
    ROL::Ptr<ROL::Objective<RealT> > simobj;
    if (useCVaR) {
      Teuchos::ParameterList list = parlist->sublist("SOL").sublist("Simulated");
      ROL::Ptr<ROL::PlusFunction<RealT> > pf
        = ROL::makePtr<ROL::PlusFunction<RealT>>(list);
      RealT alpha = parlist->sublist("SOL").sublist("Simulated").get("CVaR Confidence Level", 0.9);
      simobj = ROL::makePtr<ROL::SimulatedObjectiveCVaR<RealT>>(sampler, obj, pf, alpha);
    }
    else {
      simobj = ROL::makePtr<ROL::SimulatedObjective<RealT>>(sampler, obj);
    }
    std::vector<ROL::Ptr<ROL::Vector<RealT> > > vuvec, vpvec, vyvec;
    for (int i = 0; i < sampler->numMySamples(); ++i) {
      ROL::Ptr<Tpetra::MultiVector<> > vu_ptr, vp_ptr, vy_ptr;
      vu_ptr  = assembler->createStateVector(); vu_ptr->putScalar(4.567); //vu_ptr->randomize();
      vp_ptr  = assembler->createStateVector(); vp_ptr->putScalar(5.678); //vp_ptr->randomize();
      vy_ptr  = assembler->createStateVector(); vy_ptr->putScalar(6.789); //vy_ptr->randomize();
      ROL::Ptr<ROL::Vector<RealT> > vup, vpp, vyp;
      vup  = ROL::makePtr<PDE_PrimalSimVector<RealT>>(vu_ptr,pde,assembler);
      vpp  = ROL::makePtr<PDE_PrimalSimVector<RealT>>(vp_ptr,pde,assembler);
      vyp  = ROL::makePtr<PDE_PrimalSimVector<RealT>>(vy_ptr,pde,assembler);
      vuvec.push_back(vup);
      vpvec.push_back(vpp);
      vyvec.push_back(vyp);
    }
    ROL::Ptr<ROL::Vector<RealT> > vu, vp, vy;
    vu = ROL::makePtr<ROL::SimulatedVector<RealT>>(vuvec,bman);
    vp = ROL::makePtr<ROL::SimulatedVector<RealT>>(vpvec,bman);
    vy = ROL::makePtr<ROL::SimulatedVector<RealT>>(vyvec,bman);
    ROL::Ptr<ROL::Vector<RealT> > rz, rs, rt;
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
    ROL::Ptr<ROL::SimulatedVector<RealT> > vusim
      = ROL::dynamicPtrCast<ROL::SimulatedVector<RealT> >(vu);
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
    pdecon->outputTpetraVector(z_ptr,"control.txt");
    // Output expected state and samples to file
    *outStream << std::endl << "Print Expected Value of State" << std::endl;
    up->zero(); pp->zero();
    for (int i = 0; i < sampler->numMySamples(); ++i) {
      up->axpy(sampler->getMyWeight(i),*(vusim->get(i)));
    }
    bman->sumAll(*up,*pp);
    pdecon->outputTpetraVector(p_ptr,"mean_state.txt");
    // Build full objective function distribution
    *outStream << std::endl << "Print Objective CDF" << std::endl;
    int nsamp_dist = parlist->sublist("Problem").get("Number of output samples",100);
    ROL::Ptr<ROL::SampleGenerator<RealT> > sampler_dist
      = ROL::makePtr<ROL::MonteCarloGenerator<RealT>>(nsamp_dist,distVec,bman);
    print<RealT>(*robj,*zp,*sampler_dist,nsamp_dist,comm,"obj_samples.txt");
    // Build vorticity objective function distribution
    ROL::Ptr<ROL::Objective_SimOpt<RealT> > obj0
      = ROL::makePtr<IntegralObjective<RealT>>(qoi_vec[0],assembler);
    ROL::Ptr<ROL::SimController<RealT> > stateStore0
      = ROL::makePtr<ROL::SimController<RealT>>();
    ROL::Ptr<ROL::Reduced_Objective_SimOpt<RealT> > robj0
      = ROL::makePtr<ROL::Reduced_Objective_SimOpt<RealT>>(obj0, con, stateStore0, up, zp, pp, true, false);
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
