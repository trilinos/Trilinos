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
    \brief Shows how to solve the stuctural topology optimization problem.

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

#include "ROL_OptimizationSolver.hpp"
#include "ROL_ScaledStdVector.hpp"
#include "ROL_Reduced_Objective_SimOpt.hpp"
#include "ROL_Bounds.hpp"
#include "ROL_CompositeConstraint_SimOpt.hpp"
#include "ROL_LinearCombinationObjective.hpp"

#include "../../TOOLS/pdeconstraint.hpp"
#include "../../TOOLS/linearpdeconstraint.hpp"
#include "../../TOOLS/pdeobjective.hpp"
#include "../../TOOLS/pdevector.hpp"
#include "../../TOOLS/integralconstraint.hpp"
#include "../../TOOLS/meshreader.hpp"
#include "src/obj_compliance.hpp"
#include "src/obj_volume.hpp"
#include "src/obj_phasefield.hpp"
#include "src/mesh_topo-opt.hpp"
#include "src/pde_elasticity.hpp"
#include "src/pde_filter.hpp"

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
    std::string filename = "input_ex01.xml";
    Teuchos::RCP<Teuchos::ParameterList> parlist = Teuchos::rcp( new Teuchos::ParameterList() );
    Teuchos::updateParametersFromXmlFile( filename, parlist.ptr() );

    // Retrieve parameters.
    const std::string example = parlist->sublist("Problem").get("Example", "Default");
    const RealT domainWidth   = parlist->sublist("Geometry").get("Width", 1.0);
    const RealT domainHeight  = parlist->sublist("Geometry").get("Height", 1.0);
    const RealT domainDepth   = parlist->sublist("Geometry").get("Depth", 1.0);
    const RealT volFraction   = parlist->sublist("Problem").get("Volume Fraction", 0.4);
    const RealT cmpFactor     = parlist->sublist("Problem").get("Compliance Factor", 1.1);
    RealT cmpScaling          = parlist->sublist("Problem").get("Compliance Scaling", 1e-4);
    int probDim               = parlist->sublist("Problem").get("Problem Dimension", 2);
    const std::string minType = parlist->sublist("Problem").get("Minimization Type", "Volume");
    const bool usePhaseField  = parlist->sublist("Problem").get("Use Phase Field", false);
    bool useFilter            = parlist->sublist("Problem").get("Use Filter", true);
    if (example == "2D Wheel"                   ||
        example == "2D Truss"                   ||
        example == "2D Cantilever with 1 Load"  ||
        example == "2D Cantilever with 3 Loads" ||
        example == "2D Beams"                   ||
        example == "2D Carrier Plate") {
      probDim = 2;
    }
    else if (example == "3D Cantilever") {
      probDim = 3;
    }
    if (usePhaseField) {
      useFilter = false;
    }
    *outStream << std::endl;
    *outStream << "Problem Data"          << std::endl;
    *outStream << "  Example:           " << example << std::endl;
    *outStream << "  Dimension:         " << probDim << std::endl;
    *outStream << "  Minimize Type:     " << minType << std::endl;
    *outStream << "  Use Phase Field:   " << usePhaseField << std::endl;
    if (!usePhaseField) {
      *outStream << "  SIMP Power:        "
                 << parlist->sublist("Problem").get("SIMP Power",3.0) << std::endl;
      *outStream << "  Use Filter:        " << useFilter << std::endl;
    }
    if (minType == "Volume") {
      *outStream << "  Compliance Factor: " << cmpFactor << std::endl;
    }
    else if (minType == "Compliance") {
      *outStream << "  Volume Fraction:   " << volFraction << std::endl;
    }
    *outStream << std::endl;

    /*** Initialize main data structure. ***/
    ROL::Ptr<MeshManager<RealT> > meshMgr;
    if (probDim == 2) {
      meshMgr = ROL::makePtr<MeshManager_TopoOpt<RealT>>(*parlist);
    } else if (probDim == 3) {
      meshMgr = ROL::makePtr<MeshReader<RealT>>(*parlist);
    }
    else {
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument,
        ">>> PDE-OPT/topo-opt/elasticity/example_01.cpp: Problem dim is not 2 or 3!");
    }
    // Initialize PDE describing elasticity equations.
    ROL::Ptr<PDE_Elasticity<RealT> > pde
      = ROL::makePtr<PDE_Elasticity<RealT>>(*parlist);
    ROL::Ptr<ROL::Constraint_SimOpt<RealT> > con
      = ROL::makePtr<PDE_Constraint<RealT>>(pde,meshMgr,comm,*parlist,*outStream);
    // Initialize the filter PDE.
    ROL::Ptr<PDE_Filter<RealT> > pdeFilter
      = ROL::makePtr<PDE_Filter<RealT>>(*parlist);
    ROL::Ptr<ROL::Constraint_SimOpt<RealT> > conFilter
      = ROL::makePtr<Linear_PDE_Constraint<RealT>>(pdeFilter,meshMgr,comm,*parlist,*outStream);
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
    //z_ptr = assembler->createControlVector();  z_ptr->putScalar((minType=="Compliance" ? volFraction : 1.0));
    r_ptr = assembler->createResidualVector(); r_ptr->putScalar(0.0);
    ROL::Ptr<ROL::Vector<RealT> > up, pp, zp, rp;
    up = ROL::makePtr<PDE_PrimalSimVector<RealT>>(u_ptr,pde,assembler,*parlist);
    pp = ROL::makePtr<PDE_PrimalSimVector<RealT>>(p_ptr,pde,assembler,*parlist);
    zp = ROL::makePtr<PDE_PrimalOptVector<RealT>>(z_ptr,pde,assembler,*parlist);
    rp = ROL::makePtr<PDE_DualSimVector<RealT>>(r_ptr,pde,assembler,*parlist);

    // Initialize "filtered" of "unfiltered" constraint.
    ROL::Ptr<ROL::Constraint_SimOpt<RealT> > pdeWithFilter;
    if (useFilter && !usePhaseField) {
      bool useStorage = parlist->sublist("Problem").get("Use State Storage",true);
      pdeWithFilter
        = ROL::makePtr<ROL::CompositeConstraint_SimOpt<RealT>>(con, conFilter,
                                                               *rp, *rp, *up, *zp, *zp,
                                                               useStorage);
    }
    else {
      pdeWithFilter = con;
    }
    pdeWithFilter->setSolveParameters(*parlist);

    // Initialize compliance objective function.
    con->value(*rp, *up, *zp, tol);
    RealT rnorm2 = rp->dot(*rp);
    if (rnorm2 > 1e2*ROL::ROL_EPSILON<RealT>()) {
      cmpScaling /= rnorm2;
    }
    ROL::Ptr<QoI<RealT>> qoi_com
      = ROL::makePtr<QoI_Compliance_TopoOpt<RealT>>(pde->getFE(),
                                                    pde->getLoad(),
                                                    pde->getBdryFE(),
                                                    pde->getBdryCellLocIds(),
                                                    pde->getTraction(),
                                                    pde->getFieldHelper(),
                                                    cmpScaling);
    ROL::Ptr<ROL::Objective_SimOpt<RealT>> obj_com
      = ROL::makePtr<PDE_Objective<RealT>>(qoi_com,assembler);

    // Initialize reduced compliance objective function.
    bool storage = parlist->sublist("Problem").get("Use state storage",true);
    ROL::Ptr<ROL::SimController<RealT> > stateStore
      = ROL::makePtr<ROL::SimController<RealT>>();
    ROL::Ptr<ROL::Reduced_Objective_SimOpt<RealT> > robj_com
      = ROL::makePtr<ROL::Reduced_Objective_SimOpt<RealT>>(obj_com,
                     pdeWithFilter,stateStore,up,zp,pp,storage);

    // Create objective, constraint, multiplier and bounds
    ROL::Ptr<ROL::Objective<RealT>>       obj;
    ROL::Ptr<ROL::Constraint<RealT>>      icon;
    ROL::Ptr<ROL::Vector<RealT>>          iup;
    ROL::Ptr<ROL::Vector<RealT>>          imul;
    ROL::Ptr<ROL::BoundConstraint<RealT>> ibnd;
    if (!usePhaseField) {
      // Build volume objective function.
      ROL::Ptr<QoI<RealT>> qoi_vol
        = ROL::makePtr<QoI_Volume_TopoOpt<RealT>>(pde->getFE(),
                                                  pde->getFieldHelper());
      ROL::Ptr<ROL::Objective<RealT>> obj_vol
        = ROL::makePtr<IntegralOptObjective<RealT>>(qoi_vol,assembler);
      if (minType == "Volume") {
        obj  = obj_vol;
        icon = ROL::makePtr<ROL::ConstraintFromObjective<RealT>>(robj_com);
        // Set upper bound to compliance for solid beam.
        RealT comp = robj_com->value(*zp,tol);
        iup  = ROL::makePtr<ROL::SingletonVector<RealT>>(cmpFactor*comp);
        imul = ROL::makePtr<ROL::SingletonVector<RealT>>(0);
        ibnd = ROL::makePtr<ROL::Bounds<RealT>>(*iup,false);
      }
      else if (minType == "Compliance") {
        obj  = robj_com;
        icon = ROL::makePtr<ROL::ConstraintFromObjective<RealT>>(obj_vol);
        // Set upper bound to fraction of total volume.
        RealT vol = domainHeight*domainWidth*domainDepth;
        iup  = ROL::makePtr<ROL::SingletonVector<RealT>>(volFraction*vol);
        imul = ROL::makePtr<ROL::SingletonVector<RealT>>(0);
        ibnd = ROL::makePtr<ROL::Bounds<RealT>>(*iup,false);
      }
      else if (minType == "Total") {
        std::vector<ROL::Ptr<ROL::Objective<RealT>>> obj_vec(2,ROL::nullPtr);
        obj_vec[0] = robj_com;
        obj_vec[1] = obj_vol;
        std::vector<RealT> weights(2,0.0);
        weights[0] = static_cast<RealT>(1);
        weights[1] = parlist->sublist("Problem").get("Volume Objective Scale",0.04096);
        obj  = ROL::makePtr<ROL::LinearCombinationObjective<RealT>>(weights,obj_vec);
        icon = ROL::nullPtr;
        iup  = ROL::nullPtr;
        imul = ROL::nullPtr;
        ibnd = ROL::nullPtr;
      }
      else {
        throw ROL::Exception::NotImplemented("Unknown minimization type!");
      }
    }
    else {
      // Build volume objective function.
      ROL::Ptr<QoI<RealT>> qoi_vol
        = ROL::makePtr<QoI_Volume_PhaseField<RealT>>(pde->getFE(),
                                                     pde->getFieldHelper());
      ROL::Ptr<IntegralOptObjective<RealT>> obj_vol
        = ROL::makePtr<IntegralOptObjective<RealT>>(qoi_vol,assembler);
      // Build Modica-Mortola Energy objective function.
      RealT penParam = parlist->sublist("Problem").get("Phase Field Penalty Parameter",1e-1);
      ROL::Ptr<QoI<RealT>> qoi_pfe
        = ROL::makePtr<QoI_ModicaMortolaEnergy_PhaseField<RealT>>(pde->getFE(),
                                                                  pde->getFieldHelper(),
                                                                  penParam);
      ROL::Ptr<IntegralOptObjective<RealT>> obj_pfe
        = ROL::makePtr<IntegralOptObjective<RealT>>(qoi_pfe,assembler);
      // Get weights for linear combination objective.
      std::vector<RealT> weights(1,0.0);
      weights[0] = parlist->sublist("Problem").get("Phase Field Energy Objective Scale",0.00064);
      std::vector<ROL::Ptr<ROL::Objective<RealT>>> obj_vec(1,ROL::nullPtr);
      obj_vec[0] = obj_pfe;
      if (minType == "Volume") {
        weights.push_back(parlist->sublist("Problem").get("Volume Objective Scale",0.04096));
        obj_vec.push_back(obj_vol);
        obj  = ROL::makePtr<ROL::LinearCombinationObjective<RealT>>(weights,obj_vec);
        icon = ROL::makePtr<ROL::ConstraintFromObjective<RealT>>(robj_com);
        // Set upper bound to compliance for solid beam.
        RealT comp = robj_com->value(*zp,tol);
        iup  = ROL::makePtr<ROL::SingletonVector<RealT>>(cmpFactor*comp);
        imul = ROL::makePtr<ROL::SingletonVector<RealT>>(0);
        ibnd = ROL::makePtr<ROL::Bounds<RealT>>(*iup,false);
      }
      else if (minType == "Compliance") {
        weights.push_back(static_cast<RealT>(1));
        obj_vec.push_back(robj_com);
        obj  = ROL::makePtr<ROL::LinearCombinationObjective<RealT>>(weights,obj_vec);
        icon = ROL::makePtr<ROL::ConstraintFromObjective<RealT>>(obj_vol);
        // Set upper bound to fraction of total volume.
        RealT vol = domainHeight*domainWidth*domainDepth;
        iup  = ROL::makePtr<ROL::SingletonVector<RealT>>(volFraction*vol);
        imul = ROL::makePtr<ROL::SingletonVector<RealT>>(0);
        ibnd = ROL::makePtr<ROL::Bounds<RealT>>(*iup,false);
      }
      else if (minType == "Total") {
        weights.push_back(static_cast<RealT>(1));
        obj_vec.push_back(robj_com);
        weights.push_back(parlist->sublist("Problem").get("Volume Objective Scale",0.04096));
        obj_vec.push_back(obj_vol);
        obj  = ROL::makePtr<ROL::LinearCombinationObjective<RealT>>(weights,obj_vec);
        icon = ROL::nullPtr;
        iup  = ROL::nullPtr;
        imul = ROL::nullPtr;
        ibnd = ROL::nullPtr;
      }
      else {
        throw ROL::Exception::NotImplemented("Unknown minimization type!");
      }
    }

    // Initialize bound constraints.
    ROL::Ptr<Tpetra::MultiVector<> > lo_ptr, hi_ptr;
    RealT lval = (usePhaseField ? -1.0 : 0.0), uval = 1.0;
    lo_ptr = assembler->createControlVector(); lo_ptr->putScalar(lval);
    hi_ptr = assembler->createControlVector(); hi_ptr->putScalar(uval);
    ROL::Ptr<ROL::Vector<RealT> > lop, hip;
    lop = ROL::makePtr<PDE_PrimalOptVector<RealT>>(lo_ptr,pde,assembler);
    hip = ROL::makePtr<PDE_PrimalOptVector<RealT>>(hi_ptr,pde,assembler);
    ROL::Ptr<ROL::BoundConstraint<RealT> > bnd
      = ROL::makePtr<ROL::Bounds<RealT>>(lop,hip);
    if (usePhaseField) {
      bnd->deactivate();
    }

    // Set up and solve.
    ROL::OptimizationProblem<RealT> problem(obj,zp,bnd,icon,imul,ibnd);
    bool derivCheck = parlist->sublist("Problem").get("Check derivatives",false);
    if (derivCheck) {
      problem.check(*outStream);
    }
    ROL::OptimizationSolver<RealT> solver(problem,*parlist);
    Teuchos::Time algoTimer("Algorithm Time", true);
    solver.solve(*outStream);
    algoTimer.stop();
    *outStream << "Total optimization time = " << algoTimer.totalElapsedTime() << " seconds.\n";

    // Output.
    pdecon->printMeshData(*outStream);
    con->solve(*rp,*up,*zp,tol);
    pdecon->outputTpetraVector(u_ptr,"state.txt");
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
