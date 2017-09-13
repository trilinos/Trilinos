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
    Teuchos::RCP<MeshManager<RealT> > meshMgr
      = Teuchos::rcp(new MeshManager_stoch_adv_diff<RealT>(*parlist));
    // Initialize PDE describing advection-diffusion equation
    Teuchos::RCP<PDE_stoch_adv_diff<RealT> > pde
      = Teuchos::rcp(new PDE_stoch_adv_diff<RealT>(*parlist));
    Teuchos::RCP<ROL::Constraint_SimOpt<RealT> > con
      = Teuchos::rcp(new PDE_Constraint<RealT>(pde,meshMgr,serial_comm,*parlist,*outStream));
    Teuchos::RCP<PDE_Constraint<RealT> > pdeCon
      = Teuchos::rcp_dynamic_cast<PDE_Constraint<RealT> >(con);
    const Teuchos::RCP<Assembler<RealT> > assembler = pdeCon->getAssembler();

    /*************************************************************************/
    /***************** BUILD VECTORS *****************************************/
    /*************************************************************************/
    // Create state vectors
    Teuchos::RCP<Tpetra::MultiVector<> >  u_rcp = assembler->createStateVector();
    Teuchos::RCP<Tpetra::MultiVector<> >  p_rcp = assembler->createStateVector();
    Teuchos::RCP<Tpetra::MultiVector<> > du_rcp = assembler->createStateVector();
    u_rcp->randomize();  //u_rcp->putScalar(static_cast<RealT>(1));
    p_rcp->randomize();  //p_rcp->putScalar(static_cast<RealT>(1));
    du_rcp->randomize(); //du_rcp->putScalar(static_cast<RealT>(0));
    Teuchos::RCP<ROL::Vector<RealT> > up
      = Teuchos::rcp(new PDE_PrimalSimVector<RealT>(u_rcp,pde,assembler,*parlist));
    Teuchos::RCP<ROL::Vector<RealT> > pp
      = Teuchos::rcp(new PDE_PrimalSimVector<RealT>(p_rcp,pde,assembler,*parlist));
    Teuchos::RCP<ROL::Vector<RealT> > dup
      = Teuchos::rcp(new PDE_PrimalSimVector<RealT>(du_rcp,pde,assembler,*parlist));
    // Create residual vector
    Teuchos::RCP<Tpetra::MultiVector<> >  r_rcp = assembler->createResidualVector();
    Teuchos::RCP<ROL::Vector<RealT> > rp
      = Teuchos::rcp(new PDE_DualSimVector<RealT>(r_rcp,pde,assembler,*parlist));
    // Create control vectors
    Teuchos::RCP<std::vector<RealT> >  z_rcp = Teuchos::rcp(new std::vector<RealT>(controlDim));
    Teuchos::RCP<std::vector<RealT> > dz_rcp = Teuchos::rcp(new std::vector<RealT>(controlDim));
    Teuchos::RCP<std::vector<RealT> > ez_rcp = Teuchos::rcp(new std::vector<RealT>(controlDim));
    for (int i = 0; i < controlDim; ++i) {
      (*z_rcp)[i]  = random<RealT>(*comm);
      (*dz_rcp)[i] = random<RealT>(*comm);
      (*ez_rcp)[i] = random<RealT>(*comm);
    }
    Teuchos::RCP<ROL::Vector<RealT> > zp
      = Teuchos::rcp(new PDE_OptVector<RealT>(Teuchos::rcp(new ROL::StdVector<RealT>(z_rcp))));
    Teuchos::RCP<ROL::Vector<RealT> > dzp
      = Teuchos::rcp(new PDE_OptVector<RealT>(Teuchos::rcp(new ROL::StdVector<RealT>(dz_rcp))));
    Teuchos::RCP<ROL::Vector<RealT> > ezp
      = Teuchos::rcp(new PDE_OptVector<RealT>(Teuchos::rcp(new ROL::StdVector<RealT>(ez_rcp))));

    /*************************************************************************/
    /***************** BUILD COST FUNCTIONAL *********************************/
    /*************************************************************************/
    std::vector<Teuchos::RCP<QoI<RealT> > > qoi_vec(2,Teuchos::null);
    qoi_vec[0] = Teuchos::rcp(new QoI_State_Cost_stoch_adv_diff<RealT>(pde->getFE()));
    qoi_vec[1] = Teuchos::rcp(new QoI_Control_Cost_stoch_adv_diff<RealT>());
    Teuchos::RCP<StdObjective_stoch_adv_diff<RealT> > std_obj
      = Teuchos::rcp(new StdObjective_stoch_adv_diff<RealT>(*parlist));
    Teuchos::RCP<ROL::Objective_SimOpt<RealT> > obj
      = Teuchos::rcp(new PDE_Objective<RealT>(qoi_vec,std_obj,assembler));
    Teuchos::RCP<ROL::Reduced_Objective_SimOpt<RealT> > objReduced
      = Teuchos::rcp(new ROL::Reduced_Objective_SimOpt<RealT>(obj, con, up, pp, true, false));

    /*************************************************************************/
    /***************** BUILD BOUND CONSTRAINT ********************************/
    /*************************************************************************/
    Teuchos::RCP<std::vector<RealT> > zlo_rcp = Teuchos::rcp(new std::vector<RealT>(controlDim,0));
    Teuchos::RCP<std::vector<RealT> > zhi_rcp = Teuchos::rcp(new std::vector<RealT>(controlDim,1));
    Teuchos::RCP<ROL::Vector<RealT> > zlop
      = Teuchos::rcp(new PDE_OptVector<RealT>(Teuchos::rcp(new ROL::StdVector<RealT>(zlo_rcp))));
    Teuchos::RCP<ROL::Vector<RealT> > zhip
      = Teuchos::rcp(new PDE_OptVector<RealT>(Teuchos::rcp(new ROL::StdVector<RealT>(zhi_rcp))));
    Teuchos::RCP<ROL::BoundConstraint<RealT> > bnd
      = Teuchos::rcp(new ROL::Bounds<RealT>(zlop,zhip));

    /*************************************************************************/
    /***************** BUILD SAMPLER *****************************************/
    /*************************************************************************/
    int nsamp = parlist->sublist("Problem").get("Number of Samples",100);
    std::vector<RealT> tmp = {-one,one};
    std::vector<std::vector<RealT> > bounds(stochDim,tmp);
    Teuchos::RCP<ROL::BatchManager<RealT> > bman
      = Teuchos::rcp(new PDE_OptVector_BatchManager<RealT>(comm));
    Teuchos::RCP<ROL::SampleGenerator<RealT> > sampler
      = Teuchos::rcp(new ROL::MonteCarloGenerator<RealT>(nsamp,bounds,bman));

    /*************************************************************************/
    /***************** SOLVE OPTIMIZATION PROBLEM ****************************/
    /*************************************************************************/
    Teuchos::RCP<ROL::OptimizationProblem<RealT> > opt;
    Teuchos::RCP<ROL::Algorithm<RealT> > algo;
    zp->zero();

    int nQuad = 11, nSmooth = 1, N(2);
    RealT eps(1), eps0(1.e-3), stat(0);
    std::vector<RealT> statVec(nQuad);
    std::vector<RealT>  errVec(nQuad);
    std::vector<RealT> normVec(nQuad);
    bool checkDeriv = parlist->sublist("Problem").get("Check Derivatives",false);

    std::string rm = "Risk Measure", sr = "Spectral Risk";
    std::string dist = "Distribution", pf = "Plus Function";
    parlist->sublist("SOL").set("Stochastic Component Type","Risk Averse");
    parlist->sublist("SOL").sublist(rm).set("Name",sr);
    parlist->sublist("SOL").sublist(rm).sublist(sr).set("Print Quadrature to Screen",!myRank);
    parlist->sublist("SOL").sublist(rm).sublist(sr).sublist(dist).set("Name","Beta");
    parlist->sublist("SOL").sublist(rm).sublist(sr).sublist(dist).sublist("Beta").set("Shape 1",5.0);
    parlist->sublist("SOL").sublist(rm).sublist(sr).sublist(dist).sublist("Beta").set("Shape 2",2.0);
    parlist->sublist("SOL").sublist(rm).sublist(sr).sublist(pf).sublist(dist).set("Name","Parabolic");
    parlist->sublist("SOL").sublist(rm).sublist(sr).sublist(pf).sublist(dist).sublist("Parabolic").set("Lower Bound",-0.5);
    parlist->sublist("SOL").sublist(rm).sublist(sr).sublist(pf).sublist(dist).sublist("Parabolic").set("Upper Bound", 0.5);
    
    for (int i = 0; i < nQuad; ++i) {
      eps = eps0;
      parlist->sublist("SOL").sublist(rm).sublist(sr).set("Number of Quadrature Points",N);
      for (int j = 0; j < nSmooth; ++j) {
        parlist->sublist("SOL").sublist(rm).sublist(sr).sublist(pf).set("Smoothing Parameter",eps);
        // Build stochastic optimization problem
        opt = Teuchos::rcp(new ROL::OptimizationProblem<RealT>(objReduced,zp,bnd));
        parlist->sublist("SOL").set("Initial Statisitic", stat);
        opt->setStochasticObjective(*parlist,sampler);
        if (checkDeriv) {
          opt->check(*outStream);
        }
        // Solve optimization problem
        algo = Teuchos::rcp(new ROL::Algorithm<RealT>("Trust Region",*parlist,false));
        std::clock_t timer = std::clock();
        algo->run(*opt,true,*outStream);
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
                     << std::setw(25) << (*z_rcp)[k]
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
        zfile << (*z_rcp)[i] << "\n";
      }
      zfile.close();
    }
    // Output expected state and samples to file
    up->zero(); pp->zero(); dup->zero();
    RealT tol(1.e-8);
    Teuchos::RCP<ROL::BatchManager<RealT> > bman_Eu
      = Teuchos::rcp(new ROL::TpetraTeuchosBatchManager<RealT>(comm));
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
    assembler->outputTpetraVector(p_rcp,"mean_state.txt");
    // Build objective function distribution
    RealT val(0);
    int nsamp_dist = parlist->sublist("Problem").get("Number of Output Samples",100);
    Teuchos::RCP<ROL::SampleGenerator<RealT> > sampler_dist
      = Teuchos::rcp(new ROL::MonteCarloGenerator<RealT>(nsamp_dist,bounds,bman));
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
    r_rcp->norm2(res.view(0,1));

    /*************************************************************************/
    /***************** CHECK RESIDUAL NORM ***********************************/
    /*************************************************************************/
    *outStream << "Residual Norm: " << res[0] << std::endl << std::endl;
    errorFlag += (res[0] > 1.e-6 ? 1 : 0);
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
