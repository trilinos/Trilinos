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

/*! \file  example_02.cpp
    \brief Shows how to solve the stochastic advection-diffusion problem.
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

#include "ROL_Algorithm.hpp"
#include "ROL_Bounds.hpp"
#include "ROL_Reduced_Objective_SimOpt.hpp"
#include "ROL_MonteCarloGenerator.hpp"
#include "ROL_OptimizationProblem.hpp"
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
           const ROL::SharedPointer<const Teuchos::Comm<int> > &comm,
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
  ROL::SharedPointer<const Teuchos::MpiComm<int> > mpicomm
    = ROL::dynamicPointerCast<const Teuchos::MpiComm<int> >(comm);
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
  ROL::SharedPointer<std::ostream> outStream;
  Teuchos::oblackholestream bhs; // outputs nothing

  /*** Initialize communicator. ***/
  Teuchos::GlobalMPISession mpiSession (&argc, &argv, &bhs);
  ROL::SharedPointer<const Teuchos::Comm<int> > comm
    = Tpetra::DefaultPlatform::getDefaultPlatform().getComm();
  ROL::SharedPointer<const Teuchos::Comm<int> > serial_comm
    = ROL::makeShared<Teuchos::SerialComm<int>>();
  const int myRank = comm->getRank();
  if ((iprint > 0) && (myRank == 0)) {
    outStream = ROL::makeSharedFromRef(std::cout);
  }
  else {
    outStream = ROL::makeSharedFromRef(bhs);
  }
  int errorFlag  = 0;

  // *** Example body.
  try {

    /*** Read in XML input ***/
    std::string filename = "input_ex04.xml";
    ROL::SharedPointer<Teuchos::ParameterList> parlist = ROL::makeShared<Teuchos::ParameterList>();
    Teuchos::updateParametersFromXmlFile( filename, parlist.ptr() );

    // Problem dimensions
    const int controlDim = 9, stochDim = 37;
    const RealT one(1); 

    /*************************************************************************/
    /***************** BUILD GOVERNING PDE ***********************************/
    /*************************************************************************/
    /*** Initialize main data structure. ***/
    ROL::SharedPointer<MeshManager<RealT> > meshMgr
      = ROL::makeShared<MeshManager_stoch_adv_diff<RealT>>(*parlist);
    // Initialize PDE describing advection-diffusion equation
    ROL::SharedPointer<PDE_stoch_adv_diff<RealT> > pde
      = ROL::makeShared<PDE_stoch_adv_diff<RealT>>(*parlist);
    ROL::SharedPointer<ROL::Constraint_SimOpt<RealT> > con
      = ROL::makeShared<PDE_Constraint<RealT>>(pde,meshMgr,serial_comm,*parlist,*outStream);
    ROL::SharedPointer<PDE_Constraint<RealT> > pdeCon
      = ROL::dynamicPointerCast<PDE_Constraint<RealT> >(con);
    const ROL::SharedPointer<Assembler<RealT> > assembler = pdeCon->getAssembler();

    /*************************************************************************/
    /***************** BUILD VECTORS *****************************************/
    /*************************************************************************/
    ROL::SharedPointer<Tpetra::MultiVector<> > u_ptr, p_ptr, du_ptr, r_ptr;
    ROL::SharedPointer<std::vector<RealT> > z_ptr, dz_ptr, ez_ptr;
    ROL::SharedPointer<ROL::Vector<RealT> > up, pp, dup, rp, zp, dzp, ezp;
    // Create state vectors
    u_ptr  = assembler->createStateVector();
    p_ptr  = assembler->createStateVector();
    du_ptr = assembler->createStateVector();
    u_ptr->randomize();
    p_ptr->randomize();
    du_ptr->randomize();
    up  = ROL::makeShared<PDE_PrimalSimVector<RealT>>(u_ptr,pde,assembler,*parlist);
    pp  = ROL::makeShared<PDE_PrimalSimVector<RealT>>(p_ptr,pde,assembler,*parlist);
    dup = ROL::makeShared<PDE_PrimalSimVector<RealT>>(du_ptr,pde,assembler,*parlist);
    // Create residual vector
    r_ptr  = assembler->createResidualVector();
    rp  = ROL::makeShared<PDE_DualSimVector<RealT>>(r_ptr,pde,assembler,*parlist);
    // Create control vectors
    z_ptr  = ROL::makeShared<std::vector<RealT>>(controlDim);
    dz_ptr = ROL::makeShared<std::vector<RealT>>(controlDim);
    ez_ptr = ROL::makeShared<std::vector<RealT>>(controlDim);
    for (int i = 0; i < controlDim; ++i) {
      (*z_ptr)[i]  = random<RealT>(*comm);
      (*dz_ptr)[i] = random<RealT>(*comm);
      (*ez_ptr)[i] = random<RealT>(*comm);
    }
    zp  = ROL::makeShared<PDE_OptVector<RealT>>(ROL::makeShared<ROL::StdVector<RealT>>(z_ptr));
    dzp = ROL::makeShared<PDE_OptVector<RealT>>(ROL::makeShared<ROL::StdVector<RealT>>(dz_ptr));
    ezp = ROL::makeShared<PDE_OptVector<RealT>>(ROL::makeShared<ROL::StdVector<RealT>>(ez_ptr));

    /*************************************************************************/
    /***************** BUILD COST FUNCTIONAL *********************************/
    /*************************************************************************/
    std::vector<ROL::SharedPointer<QoI<RealT> > > qoi_vec(2,ROL::nullPointer);
    qoi_vec[0] = ROL::makeShared<QoI_State_Cost_stoch_adv_diff<RealT>>(pde->getFE());
    qoi_vec[1] = ROL::makeShared<QoI_Control_Cost_stoch_adv_diff<RealT>>();
    ROL::SharedPointer<StdObjective_stoch_adv_diff<RealT> > std_obj
      = ROL::makeShared<StdObjective_stoch_adv_diff<RealT>>(*parlist);
    ROL::SharedPointer<ROL::Objective_SimOpt<RealT> > obj
      = ROL::makeShared<PDE_Objective<RealT>>(qoi_vec,std_obj,assembler);
    ROL::SharedPointer<ROL::Reduced_Objective_SimOpt<RealT> > objReduced
      = ROL::makeShared<ROL::Reduced_Objective_SimOpt<RealT>>(obj, con, up, zp, pp, true, false);

    /*************************************************************************/
    /***************** BUILD BOUND CONSTRAINT ********************************/
    /*************************************************************************/
    ROL::SharedPointer<std::vector<RealT> > zlo_ptr = ROL::makeShared<std::vector<RealT>>(controlDim,0);
    ROL::SharedPointer<std::vector<RealT> > zhi_ptr = ROL::makeShared<std::vector<RealT>>(controlDim,1);
    ROL::SharedPointer<ROL::Vector<RealT> > zlop, zhip;
    zlop = ROL::makeShared<PDE_OptVector<RealT>>(ROL::makeShared<ROL::StdVector<RealT>>(zlo_ptr));
    zhip = ROL::makeShared<PDE_OptVector<RealT>>(ROL::makeShared<ROL::StdVector<RealT>>(zhi_ptr));
    ROL::SharedPointer<ROL::BoundConstraint<RealT> > bnd
      = ROL::makeShared<ROL::Bounds<RealT>>(zlop,zhip);

    /*************************************************************************/
    /***************** BUILD SAMPLER *****************************************/
    /*************************************************************************/
    int nsamp = parlist->sublist("Problem").get("Number of Samples",100);
    std::vector<RealT> tmp = {-one,one};
    std::vector<std::vector<RealT> > bounds(stochDim,tmp);
    ROL::SharedPointer<ROL::BatchManager<RealT> > bman
      = ROL::makeShared<PDE_OptVector_BatchManager<RealT>>(comm);
    ROL::SharedPointer<ROL::SampleGenerator<RealT> > sampler
      = ROL::makeShared<ROL::MonteCarloGenerator<RealT>>(nsamp,bounds,bman);

    /*************************************************************************/
    /***************** SOLVE OPTIMIZATION PROBLEM ****************************/
    /*************************************************************************/
    ROL::SharedPointer<ROL::OptimizationProblem<RealT> > opt;
    ROL::SharedPointer<ROL::Algorithm<RealT> > algo;

    const int nruns = 7;
    bool checkDeriv = parlist->sublist("Problem").get("Check Derivatives",false);

    // Run names
    std::vector<std::string> nvec = {"MV","RN","CVaR","MCVaR","ER","BP","KL"};

    // Parameter lists
    std::vector<Teuchos::ParameterList> plvec(7,*parlist);
    // Mean value
    plvec[0].sublist("SOL").set("Stochastic Component Type", "Mean Value");
    // Risk neutral
    plvec[1].sublist("SOL").set("Stochastic Component Type", "Risk Neutral");
    // CVaR
    plvec[2].sublist("SOL").set("Stochastic Component Type", "Risk Averse");
    plvec[2].sublist("SOL").sublist("Risk Measure").set("Name","Quantile-Based Quadrangle");
    plvec[2].sublist("SOL").sublist("Risk Measure").sublist("Quantile-Based Quadrangle").set("Confidence Level", 0.95);
    plvec[2].sublist("SOL").sublist("Risk Measure").sublist("Quantile-Based Quadrangle").set("Convex Combination Parameter", 0.0);
    plvec[2].sublist("SOL").sublist("Risk Measure").sublist("Quantile-Based Quadrangle").set("Smoothing Parameter", 1e-4);
    plvec[2].sublist("SOL").sublist("Risk Measure").sublist("Quantile-Based Quadrangle").sublist("Distribution").set("Name", "Parabolic");
    plvec[2].sublist("SOL").sublist("Risk Measure").sublist("Quantile-Based Quadrangle").sublist("Distribution").sublist("Parabolic").set("Lower Bound", 0.0);
    plvec[2].sublist("SOL").sublist("Risk Measure").sublist("Quantile-Based Quadrangle").sublist("Distribution").sublist("Parabolic").set("Upper Bound", 1.0);
    // Mixture of expectation and CVaR
    plvec[3].sublist("SOL").set("Stochastic Component Type", "Risk Averse");
    plvec[3].sublist("SOL").sublist("Risk Measure").set("Name","Quantile-Based Quadrangle");
    plvec[3].sublist("SOL").sublist("Risk Measure").sublist("Quantile-Based Quadrangle").set("Confidence Level", 0.95);
    plvec[3].sublist("SOL").sublist("Risk Measure").sublist("Quantile-Based Quadrangle").set("Convex Combination Parameter", 0.5);
    plvec[3].sublist("SOL").sublist("Risk Measure").sublist("Quantile-Based Quadrangle").set("Smoothing Parameter", 1e-4);
    plvec[3].sublist("SOL").sublist("Risk Measure").sublist("Quantile-Based Quadrangle").sublist("Distribution").set("Name", "Parabolic");
    plvec[3].sublist("SOL").sublist("Risk Measure").sublist("Quantile-Based Quadrangle").sublist("Distribution").sublist("Parabolic").set("Lower Bound", 0.0);
    plvec[3].sublist("SOL").sublist("Risk Measure").sublist("Quantile-Based Quadrangle").sublist("Distribution").sublist("Parabolic").set("Upper Bound", 1.0);
    // Entropic risk
    plvec[4].sublist("SOL").set("Stochastic Component Type", "Risk Averse");
    plvec[4].sublist("SOL").sublist("Risk Measure").set("Name","Exponential Utility");
    plvec[4].sublist("SOL").sublist("Risk Measure").sublist("Exponential Utility").set("Rate", 1.0);
    // bPOE
    plvec[5].sublist("SOL").set("Stochastic Component Type", "Risk Averse");
    plvec[5].sublist("SOL").sublist("Risk Measure").set("Name","bPOE");
    plvec[5].sublist("SOL").sublist("Risk Measure").sublist("bPOE").set("Moment Order", 2.0);
    plvec[5].sublist("SOL").sublist("Risk Measure").sublist("bPOE").set("Threshold", 6.0);
    // KL-divergence distributionally robust optimization
    plvec[6].sublist("SOL").set("Stochastic Component Type", "Risk Averse");
    plvec[6].sublist("SOL").sublist("Risk Measure").set("Name","KL Divergence");
    plvec[6].sublist("SOL").sublist("Risk Measure").sublist("KL Divergence").set("Threshold", 0.1);
    
    RealT stat(1);
    for (int i = 0; i < nruns; ++i) {
      //zp->zero();

      // Build stochastic optimization problem
      opt = ROL::makeShared<ROL::OptimizationProblem<RealT>>(objReduced,zp,bnd);
      plvec[i].sublist("SOL").set("Initial Statistic", stat);
      opt->setStochasticObjective(plvec[i],sampler);
      if (checkDeriv) {
        opt->check(*outStream);
      }

      // Solve optimization problem
      algo = ROL::makeShared<ROL::Algorithm<RealT>>("Trust Region",plvec[i],false);
      std::clock_t timer = std::clock();
      algo->run(*opt,true,*outStream);
      stat = opt->getSolutionStatistic();
      *outStream << "Optimization time: "
                 << static_cast<RealT>(std::clock()-timer)/static_cast<RealT>(CLOCKS_PER_SEC)
                 << " seconds." << std::endl << std::endl;

      // Output control to file
      if ( myRank == 0 ) {
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
      int nsamp_dist = plvec[i].sublist("Problem").get("Number of Output Samples",100);
      ROL::SharedPointer<ROL::SampleGenerator<RealT> > sampler_dist
        = ROL::makeShared<ROL::MonteCarloGenerator<RealT>>(nsamp_dist,bounds,bman);
      std::stringstream name;
      name << "obj_samples_" << nvec[i] << ".txt";
      print<RealT>(*objReduced,*zp,*sampler_dist,nsamp_dist,comm,name.str());
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
