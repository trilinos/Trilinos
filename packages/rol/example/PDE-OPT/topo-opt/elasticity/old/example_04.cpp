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

/*! \file  example_04.cpp
    \brief Shows how to solve the stuctural topology optimization problem
           under uncertainty.
*/

#include "Teuchos_Comm.hpp"
#ifdef HAVE_MPI
#include "Teuchos_DefaultMpiComm.hpp"
#endif
#include "Teuchos_Time.hpp"
#include "ROL_Stream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "Tpetra_DefaultPlatform.hpp"
#include "Tpetra_Version.hpp"

#include <iostream>
#include <algorithm>

#include "ROL_Algorithm.hpp"
#include "ROL_AugmentedLagrangian.hpp"
#include "ROL_ScaledStdVector.hpp"
#include "ROL_Reduced_Objective_SimOpt.hpp"
#include "ROL_Reduced_Constraint_SimOpt.hpp"
#include "ROL_Bounds.hpp"
#include "ROL_CompositeConstraint_SimOpt.hpp"
#include "ROL_MonteCarloGenerator.hpp"
#include "ROL_OptimizationProblem.hpp"
#include "ROL_TpetraTeuchosBatchManager.hpp"

#include "../../TOOLS/pdeconstraint.hpp"
#include "../../TOOLS/linearpdeconstraint.hpp"
#include "../../TOOLS/pdeobjective.hpp"
#include "../../TOOLS/pdevector.hpp"
#include "../../TOOLS/integralconstraint.hpp"
#include "obj_compliance.hpp"
#include "obj_volume.hpp"
#include "mesh_topo-opt.hpp"
#include "pde_elasticity.hpp"
#include "pde_filter.hpp"

typedef double RealT;

template<class Real>
void setUpAndSolve(ROL::OptimizationProblem<Real> &opt,
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
    RealT tol(1e-8), one(1);

    /*** Read in XML input ***/
    std::string filename = "input_ex04.xml";
    Teuchos::RCP<Teuchos::ParameterList> parlist = Teuchos::rcp( new Teuchos::ParameterList() );
    Teuchos::updateParametersFromXmlFile( filename, parlist.ptr() );

    // Retrieve parameters.
    const RealT volFraction  = parlist->sublist("Problem").get("Volume Fraction", 0.4);
    const RealT objFactor    = parlist->sublist("Problem").get("Objective Scaling", 1e-4);

    /*** Initialize main data structure. ***/
    ROL::Ptr<MeshManager<RealT> > meshMgr
      = ROL::makePtr<MeshManager_TopoOpt<RealT>>(*parlist);
    // Initialize PDE describing elasticity equations.
    ROL::Ptr<PDE_Elasticity<RealT> > pde
      = ROL::makePtr<PDE_Elasticity<RealT>>(*parlist);
    ROL::Ptr<ROL::Constraint_SimOpt<RealT> > con
      = ROL::makePtr<PDE_Constraint<RealT>>(pde,meshMgr,serial_comm,*parlist,*outStream);
    // Initialize the filter PDE.
    ROL::Ptr<PDE_Filter<RealT> > pdeFilter
      = ROL::makePtr<PDE_Filter<RealT>>(*parlist);
    ROL::Ptr<ROL::Constraint_SimOpt<RealT> > conFilter
      = ROL::makePtr<Linear_PDE_Constraint<RealT>>(pdeFilter,meshMgr,serial_comm,*parlist,*outStream);
    // Cast the constraint and get the assembler.
    ROL::Ptr<PDE_Constraint<RealT> > pdecon
      = ROL::dynamicPtrCast<PDE_Constraint<RealT> >(con);
    ROL::Ptr<Assembler<RealT> > assembler = pdecon->getAssembler();
    pdecon->printMeshData(*outStream);
    con->setSolveParameters(*parlist);

    // Create state vector.
    ROL::Ptr<Tpetra::MultiVector<> > u_ptr = assembler->createStateVector();
    u_ptr->randomize();
    ROL::Ptr<ROL::Vector<RealT> > up
      = ROL::makePtr<PDE_PrimalSimVector<RealT>>(u_ptr,pde,assembler,*parlist);
    ROL::Ptr<Tpetra::MultiVector<> > p_ptr = assembler->createStateVector();
    p_ptr->randomize();
    ROL::Ptr<ROL::Vector<RealT> > pp
      = ROL::makePtr<PDE_PrimalSimVector<RealT>>(p_ptr,pde,assembler,*parlist);
    // Create control vector.
    ROL::Ptr<Tpetra::MultiVector<> > z_ptr = assembler->createControlVector();
    //z_ptr->randomize();
    z_ptr->putScalar(volFraction);
    //z_ptr->putScalar(0);
    ROL::Ptr<ROL::Vector<RealT> > zp
      = ROL::makePtr<PDE_PrimalOptVector<RealT>>(z_ptr,pde,assembler,*parlist);
    // Create residual vector.
    ROL::Ptr<Tpetra::MultiVector<> > r_ptr = assembler->createResidualVector();
    r_ptr->putScalar(0.0);
    ROL::Ptr<ROL::Vector<RealT> > rp
      = ROL::makePtr<PDE_DualSimVector<RealT>>(r_ptr,pde,assembler,*parlist);
    // Create state direction vector.
    ROL::Ptr<Tpetra::MultiVector<> > du_ptr = assembler->createStateVector();
    du_ptr->randomize();
    //du_ptr->putScalar(0);
    ROL::Ptr<ROL::Vector<RealT> > dup
      = ROL::makePtr<PDE_PrimalSimVector<RealT>>(du_ptr,pde,assembler,*parlist);
    // Create control direction vector.
    ROL::Ptr<Tpetra::MultiVector<> > dz_ptr = assembler->createControlVector();
    dz_ptr->randomize();
    dz_ptr->scale(0.01);
    ROL::Ptr<ROL::Vector<RealT> > dzp
      = ROL::makePtr<PDE_PrimalOptVector<RealT>>(dz_ptr,pde,assembler,*parlist);
    // Create control test vector.
    ROL::Ptr<Tpetra::MultiVector<> > rz_ptr = assembler->createControlVector();
    rz_ptr->randomize();
    ROL::Ptr<ROL::Vector<RealT> > rzp
      = ROL::makePtr<PDE_PrimalOptVector<RealT>>(rz_ptr,pde,assembler,*parlist);

    ROL::Ptr<Tpetra::MultiVector<> > dualu_ptr = assembler->createStateVector();
    ROL::Ptr<ROL::Vector<RealT> > dualup
      = ROL::makePtr<PDE_DualSimVector<RealT>>(dualu_ptr,pde,assembler,*parlist);
    ROL::Ptr<Tpetra::MultiVector<> > dualz_ptr = assembler->createControlVector();
    ROL::Ptr<ROL::Vector<RealT> > dualzp
      = ROL::makePtr<PDE_DualOptVector<RealT>>(dualz_ptr,pde,assembler,*parlist);

    // Create ROL SimOpt vectors.
    ROL::Vector_SimOpt<RealT> x(up,zp);
    ROL::Vector_SimOpt<RealT> d(dup,dzp);

    // Initialize "filtered" or "unfiltered" constraint.
    ROL::Ptr<ROL::Constraint_SimOpt<RealT> > pdeWithFilter;
    bool useFilter  = parlist->sublist("Problem").get("Use Filter", true);
    if (useFilter) {
      pdeWithFilter
        = ROL::makePtr<ROL::CompositeConstraint_SimOpt<RealT>>(con, conFilter, *rp, *rp, *up, *zp, *zp);
    }
    else {
      pdeWithFilter = con;
    }
    pdeWithFilter->setSolveParameters(*parlist);

    // Initialize compliance objective function.
    std::vector<ROL::Ptr<QoI<RealT> > > qoi_vec(2,ROL::nullPtr);
    qoi_vec[0] = ROL::makePtr<QoI_Compliance_TopoOpt<RealT>>(pde->getFE(),
                                                             pde->getLoad(),
                                                             pde->getFieldHelper(),
                                                             objFactor);
    qoi_vec[1] = ROL::makePtr<QoI_Volume_TopoOpt<RealT>>(pde->getFE(),
                                                         pde->getFieldHelper(),
                                                         volFraction);
    RealT lambda = parlist->sublist("Problem").get("Volume Cost Parameter",1.0);
    ROL::Ptr<StdObjective_TopoOpt<RealT> > std_obj
      = ROL::makePtr<StdObjective_TopoOpt<RealT>>(lambda);
    ROL::Ptr<ROL::Objective_SimOpt<RealT> > obj
      = ROL::makePtr<PDE_Objective<RealT>>(qoi_vec,std_obj,assembler);
    // Initialize volume objective
    ROL::Ptr<IntegralObjective<RealT> > volObj
      = ROL::makePtr<IntegralObjective<RealT>>(qoi_vec[1],assembler);

    // Initialize reduced compliance function.
    bool storage = parlist->sublist("Problem").get("Use state storage",true);
    ROL::Ptr<ROL::SimController<RealT> > stateStore
      = ROL::makePtr<ROL::SimController<RealT>>();
    ROL::Ptr<ROL::Reduced_Objective_SimOpt<RealT> > objRed
      = ROL::makePtr<ROL::Reduced_Objective_SimOpt<RealT>>(obj,pdeWithFilter,
                                                              stateStore,up,zp,pp,
                                                              storage);

    // Initialize bound constraints.
    ROL::Ptr<Tpetra::MultiVector<> > lo_ptr = assembler->createControlVector();
    ROL::Ptr<Tpetra::MultiVector<> > hi_ptr = assembler->createControlVector();
    lo_ptr->putScalar(0.0); hi_ptr->putScalar(1.0);
    ROL::Ptr<ROL::Vector<RealT> > lop
      = ROL::makePtr<PDE_PrimalOptVector<RealT>>(lo_ptr,pde,assembler);
    ROL::Ptr<ROL::Vector<RealT> > hip
      = ROL::makePtr<PDE_PrimalOptVector<RealT>>(hi_ptr,pde,assembler);
    ROL::Ptr<ROL::BoundConstraint<RealT> > bnd
      = ROL::makePtr<ROL::Bounds<RealT>>(lop,hip);

    /*************************************************************************/
    /***************** BUILD SAMPLER *****************************************/
    /*************************************************************************/
    int nsamp      = parlist->sublist("Problem").get("Number of samples", 4);
    int nsamp_dist = parlist->sublist("Problem").get("Number of Output Samples",100);
    Teuchos::Array<RealT> loadMag
      = Teuchos::getArrayFromStringParameter<double>(parlist->sublist("Problem").sublist("Load"), "Magnitude");
    int nLoads = loadMag.size();
    int dim    = 2;
    std::vector<ROL::Ptr<ROL::Distribution<RealT> > > distVec(dim*nLoads);
    for (int i = 0; i < nLoads; ++i) {
      std::stringstream sli;
      sli << "Stochastic Load " << i;
      Teuchos::ParameterList magList;
      magList.sublist("Distribution") = parlist->sublist("Problem").sublist(sli.str()).sublist("Magnitude");
      //magList.print(*outStream);
      distVec[i*dim + 0] = ROL::DistributionFactory<RealT>(magList);
      Teuchos::ParameterList angList;
      angList.sublist("Distribution") = parlist->sublist("Problem").sublist(sli.str()).sublist("Polar Angle");
      //angList.print(*outStream);
      distVec[i*dim + 1] = ROL::DistributionFactory<RealT>(angList);
    }
    ROL::Ptr<ROL::BatchManager<RealT> > bman
      = ROL::makePtr<ROL::TpetraTeuchosBatchManager<RealT>>(comm);
    ROL::Ptr<ROL::SampleGenerator<RealT> > sampler
      = ROL::makePtr<ROL::MonteCarloGenerator<RealT>>(nsamp,distVec,bman);
    ROL::Ptr<ROL::SampleGenerator<RealT> > sampler_dist
      = ROL::makePtr<ROL::MonteCarloGenerator<RealT>>(nsamp_dist,distVec,bman);

    /*************************************************************************/
    /***************** SOLVE OPTIMIZATION PROBLEMS ***************************/
    /*************************************************************************/
    ROL::Ptr<ROL::OptimizationProblem<RealT> > opt;
    std::vector<RealT> vol;
    std::vector<RealT> var;

    /*************************************************************************/
    /***************** SOLVE RISK NEUTRAL ************************************/
    /*************************************************************************/
    // Solve.
    parlist->sublist("SOL").set("Stochastic Component Type","Risk Neutral");
    opt = ROL::makePtr<ROL::OptimizationProblem<RealT>>(objRed,zp,bnd);
    parlist->sublist("SOL").set("Initial Statistic",one);
    opt->setStochasticObjective(*parlist,sampler);
    setUpAndSolve<RealT>(*opt,*parlist,*outStream);
    // Output.
    vol.push_back(volObj->value(*up,*zp,tol));
    var.push_back(opt->getSolutionStatistic());
    std::string DensRN = "density_RN.txt";
    pdecon->outputTpetraVector(z_ptr,DensRN);
    std::string ObjRN = "obj_samples_RN.txt";
    print<RealT>(*objRed,*zp,*sampler_dist,nsamp_dist,comm,ObjRN);

    /*************************************************************************/
    /***************** SOLVE MEAN PLUS CVAR **********************************/
    /*************************************************************************/
    const int N = parlist->sublist("Problem").get("Denominator for Convex Combination",8);
    RealT mu(0);
    parlist->sublist("SOL").set("Stochastic Component Type","Risk Averse");
    parlist->sublist("SOL").sublist("Risk Measure").set("Name","Quantile-Based Quadrangle");
    parlist->sublist("SOL").sublist("Risk Measure").sublist("Quantile-Based Quadrangle").set("Confidence Level",0.9);
    parlist->sublist("SOL").sublist("Risk Measure").sublist("Quantile-Based Quadrangle").set("Smoothing Parameter",1e-4);
    parlist->sublist("SOL").sublist("Risk Measure").sublist("Quantile-Based Quadrangle").sublist("Distribution").set("Name","Parabolic");
    parlist->sublist("SOL").sublist("Risk Measure").sublist("Quantile-Based Quadrangle").sublist("Distribution").sublist("Parabolic").set("Lower Bound",0.0);
    parlist->sublist("SOL").sublist("Risk Measure").sublist("Quantile-Based Quadrangle").sublist("Distribution").sublist("Parabolic").set("Upper Bound",1.0);
    for (int i = 0; i < N; ++i) {
      mu = static_cast<RealT>(N-i-1)/static_cast<RealT>(N);
      // Solve.
      parlist->sublist("SOL").sublist("Risk Measure").sublist("Quantile-Based Quadrangle").set("Convex Combination Parameter",mu);
      opt = ROL::makePtr<ROL::OptimizationProblem<RealT>>(objRed,zp,bnd);
      parlist->sublist("SOL").set("Initial Statistic",var[i]);
      opt->setStochasticObjective(*parlist,sampler);
      setUpAndSolve<RealT>(*opt,*parlist,*outStream);
      // Output.
      vol.push_back(volObj->value(*up,*zp,tol));
      var.push_back(opt->getSolutionStatistic());
      std::stringstream nameDens;
      nameDens << "density_CVaR_" << N-i-1 << "_" << N << ".txt";
      pdecon->outputTpetraVector(z_ptr,nameDens.str().c_str());
      std::stringstream nameObj;
      nameObj << "obj_samples_CVaR_" << N-i-1 << "_" << N << ".txt";
      print<RealT>(*objRed,*zp,*sampler_dist,nsamp_dist,comm,nameObj.str());
    }

    /*************************************************************************/
    /***************** PRINT VOLUME AND VAR **********************************/
    /*************************************************************************/
    const int rank = Teuchos::rank<int>(*comm);
    if ( rank==0 ) {
      std::stringstream nameVOL, nameVAR;
      nameVOL << "vol.txt";
      nameVAR << "var.txt";
      std::ofstream fileVOL, fileVAR;
      fileVOL.open(nameVOL.str());
      fileVAR.open(nameVAR.str());
      fileVOL << std::scientific << std::setprecision(15);
      fileVAR << std::scientific << std::setprecision(15);
      int size = var.size();
      for (int i = 0; i < size; ++i) {
        fileVOL << std::setw(25) << std::left << vol[i] << std::endl;
        fileVAR << std::setw(25) << std::left << var[i] << std::endl;
      }
      fileVOL.close();
      fileVAR.close();
    }

    // Get a summary from the time monitor.
    Teuchos::TimeMonitor::summarize();
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
