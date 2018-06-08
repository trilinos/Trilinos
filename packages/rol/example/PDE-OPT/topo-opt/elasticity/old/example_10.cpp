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

/*! \file  example_10.cpp
    \brief Shows how to minimize volume subject to a constraint on
           compliance or minimize compliance subject to a constraint
           on volume.
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
#include "mesh_wheel.hpp"
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
    std::string filename = "input_ex10.xml";
    Teuchos::RCP<Teuchos::ParameterList> parlist = Teuchos::rcp( new Teuchos::ParameterList() );
    Teuchos::updateParametersFromXmlFile( filename, parlist.ptr() );

    // Retrieve parameters.
    int probDim      = parlist->sublist("Problem").get("Problem Dimension",2);
    bool volMin      = parlist->sublist("Problem").get("Minimize Volume", true);
    RealT cmpScaling = parlist->sublist("Problem").get("Compliance Scaling", 1e-4);

    /*** Initialize main data structure. ***/
    ROL::Ptr<MeshManager<RealT> > meshMgr;
    if (probDim == 2) {
      meshMgr = ROL::makePtr<MeshManager_Wheel<RealT>>(*parlist);
    } else if (probDim == 3) {
      meshMgr = ROL::makePtr<MeshReader<RealT>>(*parlist);
    }
    else {
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument,
        ">>> PDE-OPT/topo-opt/elasticity/example_10.cpp: Problem dim is not 2 or 3!");
    }
    // Initialize PDE defining elasticity equations.
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
    con->setSolveParameters(*parlist);

    // Create vectors.
    ROL::Ptr<Tpetra::MultiVector<> > u_ptr, p_ptr, z_ptr, r_ptr;
    u_ptr = assembler->createStateVector();    u_ptr->putScalar(0.0);
    p_ptr = assembler->createStateVector();    p_ptr->putScalar(0.0);
    z_ptr = assembler->createControlVector();  z_ptr->putScalar(1.0);
    r_ptr = assembler->createResidualVector(); r_ptr->putScalar(0.0);
    ROL::Ptr<ROL::Vector<RealT> > up, pp, zp, rp;
    up = ROL::makePtr<PDE_PrimalSimVector<RealT>>(u_ptr,pde,assembler,*parlist);
    pp = ROL::makePtr<PDE_PrimalSimVector<RealT>>(p_ptr,pde,assembler,*parlist);
    zp = ROL::makePtr<PDE_PrimalOptVector<RealT>>(z_ptr,pde,assembler,*parlist);
    rp = ROL::makePtr<PDE_DualSimVector<RealT>>(r_ptr,pde,assembler,*parlist);

    // Build sampler.
    int nsamp    = parlist->sublist("Problem").get("Number of samples", 4);
    int stochDim = 0;
    std::vector<ROL::Ptr<ROL::Distribution<RealT> > > distVec;
    if (parlist->sublist("Problem").isSublist("Load")) {
      Teuchos::Array<RealT> loadMag
        = Teuchos::getArrayFromStringParameter<double>(parlist->sublist("Problem").sublist("Load"), "Magnitude");
      int nLoads = loadMag.size();
      stochDim   = probDim*nLoads;
      for (int i = 0; i < nLoads; ++i) {
        std::stringstream sli;
        sli << "Stochastic Load " << i;
        // Magnitude distribution.
        Teuchos::ParameterList magList;
        magList.sublist("Distribution") = parlist->sublist("Problem").sublist(sli.str()).sublist("Magnitude");
        distVec.push_back(ROL::DistributionFactory<RealT>(magList));
        // Polar angle distribution.
        Teuchos::ParameterList polList;
        polList.sublist("Distribution") = parlist->sublist("Problem").sublist(sli.str()).sublist("Polar Angle");
        distVec.push_back(ROL::DistributionFactory<RealT>(polList));
        // Azimuth angle distribution.
        if (probDim == 3) {
          Teuchos::ParameterList aziList;
          aziList.sublist("Distribution") = parlist->sublist("Problem").sublist(sli.str()).sublist("Azimuth Angle");
          distVec.push_back(ROL::DistributionFactory<RealT>(aziList));
        }
      }
    }
    if (parlist->sublist("Problem").isSublist("Traction")) {
      Teuchos::Array<RealT> tractionMag
        = Teuchos::getArrayFromStringParameter<double>(parlist->sublist("Problem").sublist("Traction"), "Magnitude");
      int nTractions = tractionMag.size();
      stochDim      += probDim*nTractions;
      for (int i = 0; i < nTractions; ++i) {
        std::stringstream sli;
        sli << "Stochastic Traction " << i;
        // Magnitude distribution.
        Teuchos::ParameterList magList;
        magList.sublist("Distribution") = parlist->sublist("Problem").sublist(sli.str()).sublist("Magnitude");
        distVec.push_back(ROL::DistributionFactory<RealT>(magList));
        // Polar angle distribution.
        Teuchos::ParameterList polList;
        polList.sublist("Distribution") = parlist->sublist("Problem").sublist(sli.str()).sublist("Polar Angle");
        distVec.push_back(ROL::DistributionFactory<RealT>(polList));
        // Azimuth angle distribution.
        if (probDim == 3) {
          Teuchos::ParameterList aziList;
          aziList.sublist("Distribution") = parlist->sublist("Problem").sublist(sli.str()).sublist("Azimuth Angle");
          distVec.push_back(ROL::DistributionFactory<RealT>(aziList));
        }
      }
    }
    ROL::Ptr<ROL::BatchManager<RealT> > bman
      = ROL::makePtr<ROL::TpetraTeuchosBatchManager<RealT>>(comm);
    ROL::Ptr<ROL::SampleGenerator<RealT> > sampler
      = ROL::makePtr<ROL::MonteCarloGenerator<RealT>>(nsamp,distVec,bman);

    // Initialize "filtered" of "unfiltered" constraint.
    ROL::Ptr<ROL::Constraint_SimOpt<RealT> > pdeWithFilter;
    bool useFilter = parlist->sublist("Problem").get("Use Filter", true);
    if (useFilter) {
      bool useStorage = parlist->sublist("Problem").get("Use State Storage",true);
      pdeWithFilter =
        ROL::makePtr<ROL::CompositeConstraint_SimOpt<RealT>>
          (con, conFilter, *rp, *rp, *up, *zp, *zp, useStorage);
    }
    else {
      pdeWithFilter = con;
    }
    pdeWithFilter->setSolveParameters(*parlist);

    // Compute compliance scaling
    con->value(*rp, *up, *zp, tol);
    RealT rnorm2 = rp->dot(*rp);
    if (rnorm2 > 1e2*ROL::ROL_EPSILON<RealT>()) {
      cmpScaling /= rnorm2;
    }

    // Initialize volume and compliance objective functions.
    ROL::Ptr<QoI<RealT> > qoi_vol
      = ROL::makePtr<QoI_Volume_TopoOpt<RealT>>(pde->getFE(),pde->getFieldHelper());
    ROL::Ptr<ROL::Objective<RealT> > obj_vol
      = ROL::makePtr<IntegralOptObjective<RealT>>(qoi_vol,assembler);
    ROL::Ptr<QoI<RealT> > qoi_com
      = ROL::makePtr<QoI_Compliance_TopoOpt<RealT>>(pde->getFE(), pde->getLoad(),// Volumetric Load
                                                    pde->getBdryFE(),pde->getBdryCellLocIds(),pde->getTraction(), // Traction
                                                    pde->getFieldHelper(),cmpScaling);
    ROL::Ptr<ROL::Objective_SimOpt<RealT> > obj_com
      = ROL::makePtr<PDE_Objective<RealT>>(qoi_com,assembler);

    // Initialize reduced compliance objective function.
    bool storage     = parlist->sublist("Problem").get("Use state storage",true);
    std::string type = parlist->sublist("SOL").get("Stochastic Component Type","Risk Neutral");
    storage = (type == "Risk Neutral") ? false : storage;
    ROL::Ptr<ROL::SimController<RealT> > stateStore
      = ROL::makePtr<ROL::SimController<RealT>>();
    ROL::Ptr<ROL::Reduced_Objective_SimOpt<RealT> > robj_com
      = ROL::makePtr<ROL::Reduced_Objective_SimOpt<RealT>>(obj_com,
                     pdeWithFilter,stateStore,up,zp,pp,storage);

    // Create objective, constraint, multiplier and bounds
    ROL::Ptr<ROL::Objective<RealT> > obj;
    ROL::Ptr<ROL::Constraint<RealT> > icon;
    ROL::Ptr<ROL::Vector<RealT> > iup;
    if (volMin) {
      obj  = obj_vol;
      icon = ROL::makePtr<ROL::ConstraintFromObjective<RealT>>(robj_com);
      // Set upper bound to average compliance for solid beam.
      RealT cmpFactor = parlist->sublist("Problem").get("Compliance Factor", 1.1);
      RealT mycomp(0), comp(0);
      for (int i = 0; i < sampler->numMySamples(); ++i) {
        robj_com->setParameter(sampler->getMyPoint(i));
        mycomp += sampler->getMyWeight(i)*robj_com->value(*zp,tol);
      }
      sampler->sumAll(&mycomp,&comp,1);
      iup  = ROL::makePtr<ROL::SingletonVector<RealT>>(cmpFactor*comp);
    } else {
      obj  = robj_com;
      icon = ROL::makePtr<ROL::ConstraintFromObjective<RealT>>(obj_vol);
      // Set upper bound to fraction of total volume.
      RealT domainWidth  = parlist->sublist("Geometry").get("Width", 2.0);
      RealT domainHeight = parlist->sublist("Geometry").get("Height", 1.0);
      RealT domainDepth  = parlist->sublist("Geometry").get("Depth", 1.0);
      RealT volFraction  = parlist->sublist("Problem").get("Volume Fraction", 0.4);
      RealT vol          = domainHeight*domainWidth*domainDepth;
      iup  = ROL::makePtr<ROL::SingletonVector<RealT>>(volFraction*vol);
    }
    ROL::Ptr<ROL::Vector<RealT> > imul
      = ROL::makePtr<ROL::SingletonVector<RealT>>(0);
    ROL::Ptr<ROL::BoundConstraint<RealT> > ibnd
      = ROL::makePtr<ROL::Bounds<RealT>>(*iup,false);

    // Initialize bound constraints.
    ROL::Ptr<Tpetra::MultiVector<> > ll_ptr, uu_ptr;
    ll_ptr = assembler->createControlVector(); ll_ptr->putScalar(0.0);
    uu_ptr = assembler->createControlVector(); uu_ptr->putScalar(1.0);
    ROL::Ptr<ROL::Vector<RealT> > llp, uup;
    llp = ROL::makePtr<PDE_PrimalOptVector<RealT>>(ll_ptr,pde,assembler);
    uup = ROL::makePtr<PDE_PrimalOptVector<RealT>>(uu_ptr,pde,assembler);
    ROL::Ptr<ROL::BoundConstraint<RealT> > bnd
      = ROL::makePtr<ROL::Bounds<RealT>>(llp,uup);

    // Build optimization problem.
    ROL::OptimizationProblem<RealT> optProb(obj,zp,bnd,icon,imul,ibnd);
    if (volMin) {
      ROL::Ptr<ROL::BatchManager<RealT> > cbman
        = ROL::makePtr<ROL::SingletonTeuchosBatchManager<RealT,int>>(comm);
      optProb.setStochasticInequality(*parlist,sampler,cbman);
    } else {
      optProb.setStochasticObjective(*parlist,sampler);
    }

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
    for (int i = 0; i < sampler->numMySamples(); ++i) {
      con->setParameter(sampler->getMyPoint(i));
      con->solve(*rp,*up,*zp,tol);
      pp->axpy(sampler->getMyWeight(i),*up);
    }
    sampler->sumAll(*pp,*up);
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
