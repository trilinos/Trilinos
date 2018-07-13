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

/*! \file  example_09.cpp
    \brief Shows how to minimize volume subject to a constraint on
           compliance.
*/

#include "Teuchos_Comm.hpp"
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
#include "ROL_ConstraintFromObjective.hpp"
#include "ROL_OptimizationSolver.hpp"
#include "ROL_MonteCarloGenerator.hpp"
#include "ROL_TpetraTeuchosBatchManager.hpp"
#include "ROL_SingletonTeuchosBatchManager.hpp"

#include "../../TOOLS/pdeconstraint.hpp"
#include "../../TOOLS/linearpdeconstraint.hpp"
#include "../../TOOLS/pdeobjective.hpp"
#include "../../TOOLS/pdevector.hpp"
#include "../../TOOLS/integralconstraint.hpp"
#include "../../TOOLS/meshreader.hpp"
#include "obj_compliance.hpp"
#include "obj_volume.hpp"
#include "mesh_topo-opt.hpp"
#include "pde_elasticity.hpp"
#include "pde_filter.hpp"

typedef double RealT;

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
    RealT tol(1e-8);

    /*** Read in XML input ***/
    std::string filename = "input_ex09.xml";
    Teuchos::RCP<Teuchos::ParameterList> parlist = Teuchos::rcp( new Teuchos::ParameterList() );
    Teuchos::updateParametersFromXmlFile( filename, parlist.ptr() );

    // Retrieve parameters.
    const RealT cmpFactor    = parlist->sublist("Problem").get("Compliance Factor", 1.1);
    const RealT objFactor    = parlist->sublist("Problem").get("Objective Scaling", 1e-4);

    /*** Initialize main data structure. ***/
    int probDim = parlist->sublist("Problem").get("Problem Dimension",2);
    ROL::Ptr<MeshManager<RealT> > meshMgr;
    if (probDim == 2) {
      meshMgr = ROL::makePtr<MeshManager_TopoOpt<RealT> >(*parlist);
    } else if (probDim == 3) {
      meshMgr = ROL::makePtr<MeshReader<RealT> >(*parlist);
    }
    else {
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument,
        ">>> PDE-OPT/topo-opt/elasticity/example_09.cpp: Problem dim is not 2 or 3!");
    }
    // Initialize PDE describing elasticity equations.
    ROL::Ptr<PDE_Elasticity<RealT> > pde
      = ROL::makePtr<PDE_Elasticity<RealT> >(*parlist);
    ROL::Ptr<ROL::Constraint_SimOpt<RealT> > con
      = ROL::makePtr<PDE_Constraint<RealT> >(pde,meshMgr,serial_comm,*parlist,*outStream);
    // Initialize the filter PDE.
    ROL::Ptr<PDE_Filter<RealT> > pdeFilter
      = ROL::makePtr<PDE_Filter<RealT> >(*parlist);
    ROL::Ptr<ROL::Constraint_SimOpt<RealT> > conFilter
      = ROL::makePtr<Linear_PDE_Constraint<RealT> >(pdeFilter,meshMgr,serial_comm,*parlist,*outStream);
    // Cast the constraint and get the assembler.
    ROL::Ptr<PDE_Constraint<RealT> > pdecon
      = ROL::dynamicPtrCast<PDE_Constraint<RealT> >(con);
    ROL::Ptr<Assembler<RealT> > assembler = pdecon->getAssembler();
    con->setSolveParameters(*parlist);

    // Create vectors.
    ROL::Ptr<Tpetra::MultiVector<> > u_ptr, p_ptr, z_ptr, r_ptr;
    u_ptr = assembler->createStateVector();    u_ptr->putScalar(0.0);
    p_ptr = assembler->createStateVector();    p_ptr->putScalar(0.0);
    z_ptr = assembler->createControlVector();  z_ptr->putScalar(1.0);
    r_ptr = assembler->createResidualVector(); r_ptr->putScalar(0.0);
    ROL::Ptr<ROL::Vector<RealT> > up, pp, zp, rp;
    up = ROL::makePtr<PDE_PrimalSimVector<RealT> >(u_ptr,pde,assembler,*parlist);
    pp = ROL::makePtr<PDE_PrimalSimVector<RealT> >(p_ptr,pde,assembler,*parlist);
    zp = ROL::makePtr<PDE_PrimalOptVector<RealT> >(z_ptr,pde,assembler,*parlist);
    rp = ROL::makePtr<PDE_DualSimVector<RealT> >(r_ptr,pde,assembler,*parlist);

    // Build sampler.
    int nsamp  = parlist->sublist("Problem").get("Number of samples", 4);
    Teuchos::Array<RealT> loadMag
      = Teuchos::getArrayFromStringParameter<double>(parlist->sublist("Problem").sublist("Load"), "Magnitude");
    int nLoads = loadMag.size();
    int dim    = probDim;
    std::vector<ROL::Ptr<ROL::Distribution<RealT> > > distVec(dim*nLoads);
    for (int i = 0; i < nLoads; ++i) {
      std::stringstream sli;
      sli << "Stochastic Load " << i;
      // Magnitude distribution.
      Teuchos::ParameterList magList;
      magList.sublist("Distribution") = parlist->sublist("Problem").sublist(sli.str()).sublist("Magnitude");
      distVec[i*dim + 0] = ROL::DistributionFactory<RealT>(magList);
      // Polar angle distribution.
      Teuchos::ParameterList polList;
      polList.sublist("Distribution") = parlist->sublist("Problem").sublist(sli.str()).sublist("Polar Angle");
      distVec[i*dim + 1] = ROL::DistributionFactory<RealT>(polList);
      // Azimuth angle distribution.
      if (probDim == 3) {
        Teuchos::ParameterList aziList;
        aziList.sublist("Distribution") = parlist->sublist("Problem").sublist(sli.str()).sublist("Azimuth Angle");
        distVec[i*dim + 2] = ROL::DistributionFactory<RealT>(aziList);
      }
    }
    ROL::Ptr<ROL::BatchManager<RealT> > xbman
      = ROL::makePtr<ROL::TpetraTeuchosBatchManager<RealT> >(comm);
    ROL::Ptr<ROL::SampleGenerator<RealT> > xsampler
      = ROL::makePtr<ROL::MonteCarloGenerator<RealT> >(nsamp,distVec,xbman);
    ROL::Ptr<ROL::BatchManager<RealT> > cbman
      = ROL::makePtr<ROL::SingletonTeuchosBatchManager<RealT,int> >(comm);

    // Initialize "filtered" of "unfiltered" constraint.
    ROL::Ptr<ROL::Constraint_SimOpt<RealT> > pdeWithFilter;
    bool useFilter = parlist->sublist("Problem").get("Use Filter", true);
    if (useFilter) {
      bool useStorage = parlist->sublist("Problem").get("Use State Storage",true);
      pdeWithFilter
        = ROL::makePtr<ROL::CompositeConstraint_SimOpt<RealT>>(con, conFilter,
                       *rp, *rp, *up, *zp, *zp, useStorage);
    }
    else {
      pdeWithFilter = con;
    }
    pdeWithFilter->setSolveParameters(*parlist);

    // Initialize volume objective.
    ROL::Ptr<QoI<RealT> > qoi_vol
      = ROL::makePtr<QoI_Volume_TopoOpt<RealT> >(pde->getFE(),pde->getFieldHelper());
    ROL::Ptr<ROL::Objective<RealT> > vobj
      = ROL::makePtr<IntegralOptObjective<RealT> >(qoi_vol,assembler);

    // Initialize compliance inequality constraint.
    con->value(*rp, *up, *zp, tol);
    RealT objScaling = objFactor, rnorm2 = rp->dot(*rp);
    if (rnorm2 > 1e2*ROL::ROL_EPSILON<RealT>()) {
      objScaling /= rnorm2;
    }
    std::vector<ROL::Ptr<QoI<RealT> > > qoi_cmp(1,ROL::nullPtr);
    qoi_cmp[0]
      = ROL::makePtr<QoI_Compliance_TopoOpt<RealT> >(pde->getFE(), pde->getLoad(),
                                                     pde->getFieldHelper(), objScaling);
    ROL::Ptr<ROL::Objective_SimOpt<RealT> > cobj
      = ROL::makePtr<PDE_Objective<RealT> >(qoi_cmp,assembler);

    // Initialize reduced compliance objective function.
    bool storage     = parlist->sublist("Problem").get("Use state storage",true);
    std::string type = parlist->sublist("SOL").get("Stochastic Component Type","Risk Neutral");
    storage = (type == "Risk Neutral") ? false : storage;
    ROL::Ptr<ROL::SimController<RealT> > stateStore
      = ROL::makePtr<ROL::SimController<RealT>>();
    ROL::Ptr<ROL::Reduced_Objective_SimOpt<RealT> > cobjRed
      = ROL::makePtr<ROL::Reduced_Objective_SimOpt<RealT> >(cobj,
                     pdeWithFilter,stateStore,up,zp,pp,storage);

    // Create compliance constraint, multiplier and bounds
    RealT mycomp(0), comp(0);
    for (int i = 0; i < xsampler->numMySamples(); ++i) {
      cobjRed->setParameter(xsampler->getMyPoint(i));
      mycomp += xsampler->getMyWeight(i)*cobjRed->value(*zp,tol);
    }
    xsampler->sumAll(&mycomp,&comp,1);
    ROL::Ptr<ROL::ConstraintFromObjective<RealT> > icon
      = ROL::makePtr<ROL::ConstraintFromObjective<RealT> >(cobjRed);
    ROL::Ptr<std::vector<RealT> > imul_ptr, iup_ptr;
    ROL::Ptr<ROL::Vector<RealT> > imul, iup;
    imul = ROL::makePtr<ROL::SingletonVector<RealT> >(0);
    iup  = ROL::makePtr<ROL::SingletonVector<RealT> >(cmpFactor*comp);
    ROL::Ptr<ROL::BoundConstraint<RealT> > ibnd
      = ROL::makePtr<ROL::Bounds<RealT> >(*iup,false);

    // Initialize bound constraints.
    ROL::Ptr<Tpetra::MultiVector<> > lo_ptr, hi_ptr;
    lo_ptr = assembler->createControlVector(); lo_ptr->putScalar(0.0);
    hi_ptr = assembler->createControlVector(); hi_ptr->putScalar(1.0);
    ROL::Ptr<ROL::Vector<RealT> > lop, hip;
    lop = ROL::makePtr<PDE_PrimalOptVector<RealT> >(lo_ptr,pde,assembler);
    hip = ROL::makePtr<PDE_PrimalOptVector<RealT> >(hi_ptr,pde,assembler);
    ROL::Ptr<ROL::BoundConstraint<RealT> > bnd
      = ROL::makePtr<ROL::Bounds<RealT> >(lop,hip);

    // Build optimization problem.
    ROL::OptimizationProblem<RealT> optProb(vobj,zp,bnd,icon,imul,ibnd);
    optProb.setStochasticInequality(*parlist,xsampler,cbman);

    // Run derivative checks
    bool checkDeriv = parlist->sublist("Problem").get("Check derivatives",false);
    if ( checkDeriv ) {
      optProb.check(*outStream);
    }

    // Build optimization solver and solve.
    ROL::OptimizationSolver<RealT> optSolver(optProb,*parlist);
    Teuchos::Time algoTimer("Algorithm Time", true);
    optSolver.solve(*outStream);
    algoTimer.stop();
    *outStream << "Total optimization time = " << algoTimer.totalElapsedTime() << " seconds.\n";

    // Output.
    pdecon->printMeshData(*outStream);
    up->zero(); pp->zero();
    for (int i = 0; i < xsampler->numMySamples(); ++i) {
      con->setParameter(xsampler->getMyPoint(i));
      con->solve(*rp,*up,*zp,tol);
      pp->axpy(xsampler->getMyWeight(i),*up);
    }
    xsampler->sumAll(*pp,*up);
    pdecon->outputTpetraVector(u_ptr,"mean_state.txt");
    pdecon->outputTpetraVector(z_ptr,"density.txt");

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
