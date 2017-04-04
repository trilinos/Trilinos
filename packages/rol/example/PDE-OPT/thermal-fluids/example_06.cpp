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
#include "ROL_StochasticProblem.hpp"
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
void setUpAndSolve(ROL::StochasticProblem<Real> &opt,
                   Teuchos::ParameterList &parlist,
                   std::ostream &outStream) {
  ROL::Algorithm<RealT> algo("Trust Region",parlist,false);
  Teuchos::Time timer("Optimization Time", true);
  algo.run(opt,true,outStream);
  timer.stop();
  outStream << "Total optimization time = " << timer.totalElapsedTime() << " seconds." << std::endl;
}

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
    std::string filename = "input_ex06.xml";
    Teuchos::RCP<Teuchos::ParameterList> parlist = Teuchos::rcp( new Teuchos::ParameterList() );
    Teuchos::updateParametersFromXmlFile( filename, parlist.ptr() );
    parlist->sublist("SimOpt").sublist("Solve").set("Output Iteration History",myRank==0);

    /*** Initialize main data structure. ***/
    Teuchos::RCP<MeshManager<RealT> > meshMgr
      = Teuchos::rcp(new MeshManager_ThermalFluids<RealT>(*parlist));
    // Initialize PDE describing Navier-Stokes equations.
    Teuchos::RCP<PDE_ThermalFluids_ex03<RealT> > pde
      = Teuchos::rcp(new PDE_ThermalFluids_ex03<RealT>(*parlist));
    Teuchos::RCP<ROL::EqualityConstraint_SimOpt<RealT> > con
      = Teuchos::rcp(new PDE_Constraint<RealT>(pde,meshMgr,serial_comm,*parlist,*outStream));
    // Cast the constraint and get the assembler.
    Teuchos::RCP<PDE_Constraint<RealT> > pdecon
      = Teuchos::rcp_dynamic_cast<PDE_Constraint<RealT> >(con);
    Teuchos::RCP<Assembler<RealT> > assembler = pdecon->getAssembler();
    con->setSolveParameters(*parlist);
    pdecon->outputTpetraData();

    // Create state vector and set to zeroes
    Teuchos::RCP<Tpetra::MultiVector<> > u_rcp, p_rcp, du_rcp, dy_rcp, yu_rcp, yp_rcp, r_rcp, z_rcp, dz_rcp;
    u_rcp  = assembler->createStateVector();     u_rcp->randomize();
    p_rcp  = assembler->createStateVector();     p_rcp->randomize();
    du_rcp = assembler->createStateVector();     du_rcp->randomize();
    dy_rcp = assembler->createStateVector();     dy_rcp->randomize();
    yu_rcp = assembler->createStateVector();     yu_rcp->randomize();
    yp_rcp = assembler->createStateVector();     yp_rcp->randomize();
    r_rcp  = assembler->createResidualVector();  r_rcp->randomize();
    z_rcp  = assembler->createControlVector();   z_rcp->randomize();
    dz_rcp = assembler->createControlVector();   dz_rcp->randomize();
    Teuchos::RCP<ROL::Vector<RealT> > up, pp, dup, dyp, yup, ypp, rp, zp, dzp;
    up  = Teuchos::rcp(new PDE_PrimalSimVector<RealT>(u_rcp,pde,assembler));
    pp  = Teuchos::rcp(new PDE_PrimalSimVector<RealT>(p_rcp,pde,assembler));
    dup = Teuchos::rcp(new PDE_PrimalSimVector<RealT>(du_rcp,pde,assembler));
    dyp = Teuchos::rcp(new PDE_PrimalSimVector<RealT>(dy_rcp,pde,assembler));
    yup = Teuchos::rcp(new PDE_PrimalSimVector<RealT>(yu_rcp,pde,assembler));
    ypp = Teuchos::rcp(new PDE_PrimalSimVector<RealT>(yp_rcp,pde,assembler));
    rp  = Teuchos::rcp(new PDE_DualSimVector<RealT>(r_rcp,pde,assembler));
    zp  = Teuchos::rcp(new PDE_PrimalOptVector<RealT>(z_rcp,pde,assembler));
    dzp = Teuchos::rcp(new PDE_PrimalOptVector<RealT>(dz_rcp,pde,assembler));
    // Create ROL SimOpt vectors
    ROL::Vector_SimOpt<RealT> x(up,zp);
    ROL::Vector_SimOpt<RealT> d(dup,dzp);

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
    Teuchos::RCP<ROL::Objective_SimOpt<RealT> > objState
      = Teuchos::rcp(new PDE_Objective<RealT>(qoi_vec[0],assembler));
    Teuchos::RCP<ROL::Objective_SimOpt<RealT> > objCtrl
      = Teuchos::rcp(new PDE_Objective<RealT>(qoi_vec[1],assembler));
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
    int nsamp_dist = parlist->sublist("Problem").get("Number of output samples",100);
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
    Teuchos::RCP<ROL::SampleGenerator<RealT> > sampler_dist
      = Teuchos::rcp(new ROL::MonteCarloGenerator<RealT>(nsamp_dist,distVec,bman));

    /*************************************************************************/
    /***************** BUILD STOCHASTIC PROBLEM ******************************/
    /*************************************************************************/
    Teuchos::RCP<ROL::StochasticProblem<RealT> > opt;
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
      parlist->sublist("SOL").set("Stochastic Optimization Type","Risk Neutral");
      opt = Teuchos::rcp(new ROL::StochasticProblem<RealT>(*parlist,robj,sampler,zp));
      opt->setSolutionStatistic(one);
      setUpAndSolve<RealT>(*opt,*parlist,*outStream);
      // Output.
      ctrl.push_back(objCtrl->value(*up,*zp,tol));
      stat.push_back(opt->getSolutionStatistic());
      std::string CtrlRN = "control_RN.txt";
      pdecon->outputTpetraVector(z_rcp,CtrlRN);
      std::string ObjRN = "obj_samples_RN.txt";
      print<RealT>(*robj,*zp,*sampler_dist,nsamp_dist,comm,ObjRN);
    }

    /*************************************************************************/
    /***************** SOLVE MEAN PLUS CVAR **********************************/
    /*************************************************************************/
    parlist->sublist("SOL").set("Stochastic Optimization Type","Risk Averse");
    parlist->sublist("SOL").sublist("Risk Measure").set("Name","Quantile-Based Quadrangle");
    parlist->sublist("SOL").sublist("Risk Measure").sublist("Quantile-Based Quadrangle").set("Confidence Level",alpha);
    parlist->sublist("SOL").sublist("Risk Measure").sublist("Quantile-Based Quadrangle").set("Smoothing Parameter",1e-4);
    parlist->sublist("SOL").sublist("Risk Measure").sublist("Quantile-Based Quadrangle").sublist("Distribution").set("Name","Parabolic");
    parlist->sublist("SOL").sublist("Risk Measure").sublist("Quantile-Based Quadrangle").sublist("Distribution").sublist("Parabolic").set("Lower Bound",0.0);
    parlist->sublist("SOL").sublist("Risk Measure").sublist("Quantile-Based Quadrangle").sublist("Distribution").sublist("Parabolic").set("Upper Bound",1.0);
    for (int i = 0; i < N; ++i) {
      parlist->sublist("SOL").sublist("Risk Measure").sublist("Quantile-Based Quadrangle").set("Convex Combination Parameter",one-lambda[i]);
      // Solve.
      opt = Teuchos::rcp(new ROL::StochasticProblem<RealT>(*parlist,robj,sampler,zp));
      opt->setSolutionStatistic(stat[i]);
      setUpAndSolve<RealT>(*opt,*parlist,*outStream);
      // Output.
      ctrl.push_back(objCtrl->value(*up,*zp,tol));
      stat.push_back(opt->getSolutionStatistic());
      std::stringstream nameCtrl;
      nameCtrl << "control_MCVaR_" << i << ".txt";
      pdecon->outputTpetraVector(z_rcp,nameCtrl.str().c_str());
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
