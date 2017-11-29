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

/*! \file  example_08.cpp
    \brief Shows how to minimize volume subject to a constraint on
           compliance.
*/

#include "Teuchos_Comm.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_oblackholestream.hpp"
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

#include "../../TOOLS/pdeconstraint.hpp"
#include "../../TOOLS/linearpdeconstraint.hpp"
#include "../../TOOLS/pdeobjective.hpp"
#include "../../TOOLS/pdevector.hpp"
#include "../../TOOLS/integralconstraint.hpp"
#include "../../TOOLS/meshreader.hpp"
#include "obj_topo-opt.hpp"
#include "mesh_topo-opt.hpp"
#include "pde_elasticity.hpp"
#include "pde_filter.hpp"

typedef double RealT;

int main(int argc, char *argv[]) {
  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  ROL::SharedPointer<std::ostream> outStream;
  Teuchos::oblackholestream bhs; // outputs nothing

  /*** Initialize communicator. ***/
  Teuchos::GlobalMPISession mpiSession (&argc, &argv, &bhs);
  ROL::SharedPointer<const Teuchos::Comm<int> > comm
    = Tpetra::DefaultPlatform::getDefaultPlatform().getComm();
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
    RealT tol(1e-8);

    /*** Read in XML input ***/
    std::string filename = "input_ex08.xml";
    ROL::SharedPointer<Teuchos::ParameterList> parlist = ROL::makeShared<Teuchos::ParameterList>();
    Teuchos::updateParametersFromXmlFile( filename, parlist.ptr() );

    // Retrieve parameters.
    const RealT cmpFactor    = parlist->sublist("Problem").get("Compliance Factor", 1.1);
    const RealT objFactor    = parlist->sublist("Problem").get("Objective Scaling", 1e-4);

    /*** Initialize main data structure. ***/
    int probDim = parlist->sublist("Problem").get("Problem Dimension",2);
    ROL::SharedPointer<MeshManager<RealT> > meshMgr;
    if (probDim == 2) {
      meshMgr = ROL::makeShared<MeshManager_TopoOpt<RealT>>(*parlist);
    } else if (probDim == 3) {
      meshMgr = ROL::makeShared<MeshReader<RealT>>(*parlist);
    }
    else {
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument,
        ">>> PDE-OPT/topo-opt/elasticity/example_08.cpp: Problem dim is not 2 or 3!");
    }
    // Initialize PDE describing elasticity equations.
    ROL::SharedPointer<PDE_Elasticity<RealT> > pde
      = ROL::makeShared<PDE_Elasticity<RealT>>(*parlist);
    ROL::SharedPointer<ROL::Constraint_SimOpt<RealT> > con
      = ROL::makeShared<PDE_Constraint<RealT>>(pde,meshMgr,comm,*parlist,*outStream);
    // Initialize the filter PDE.
    ROL::SharedPointer<PDE_Filter<RealT> > pdeFilter
      = ROL::makeShared<PDE_Filter<RealT>>(*parlist);
    ROL::SharedPointer<ROL::Constraint_SimOpt<RealT> > conFilter
      = ROL::makeShared<Linear_PDE_Constraint<RealT>>(pdeFilter,meshMgr,comm,*parlist,*outStream);
    // Cast the constraint and get the assembler.
    ROL::SharedPointer<PDE_Constraint<RealT> > pdecon
      = ROL::dynamicPointerCast<PDE_Constraint<RealT> >(con);
    ROL::SharedPointer<Assembler<RealT> > assembler = pdecon->getAssembler();
    con->setSolveParameters(*parlist);

    // Create state vector.
    ROL::SharedPointer<Tpetra::MultiVector<> > u_ptr, p_ptr, z_ptr, r_ptr;
    u_ptr = assembler->createStateVector();    u_ptr->putScalar(0.0);
    p_ptr = assembler->createStateVector();    p_ptr->putScalar(0.0);
    z_ptr = assembler->createControlVector();  z_ptr->putScalar(1.0);
    r_ptr = assembler->createResidualVector(); r_ptr->putScalar(0.0);
    ROL::SharedPointer<ROL::Vector<RealT> > up, pp, zp, rp;
    up = ROL::makeShared<PDE_PrimalSimVector<RealT>>(u_ptr,pde,assembler,*parlist);
    pp = ROL::makeShared<PDE_PrimalSimVector<RealT>>(p_ptr,pde,assembler,*parlist);
    zp = ROL::makeShared<PDE_PrimalOptVector<RealT>>(z_ptr,pde,assembler,*parlist);
    rp = ROL::makeShared<PDE_DualSimVector<RealT>>(r_ptr,pde,assembler,*parlist);

    // Initialize "filtered" of "unfiltered" constraint.
    ROL::SharedPointer<ROL::Constraint_SimOpt<RealT> > pdeWithFilter;
    bool useFilter = parlist->sublist("Problem").get("Use Filter", true);
    if (useFilter) {
      bool useStorage = parlist->sublist("Problem").get("Use State Storage",true);
      pdeWithFilter
        = ROL::makeShared<ROL::CompositeConstraint_SimOpt<RealT>>(con, conFilter,
                       *rp, *rp, *up, *zp, *zp, useStorage);
    }
    else {
      pdeWithFilter = con;
    }
    pdeWithFilter->setSolveParameters(*parlist);

    // Initialize volume objective.
    ROL::SharedPointer<QoI<RealT> > qoi_vol
      = ROL::makeShared<QoI_VolumeObj_TopoOpt<RealT>>(pde->getFE(),pde->getFieldHelper());
    ROL::SharedPointer<ROL::Objective<RealT> > vobj
      = ROL::makeShared<IntegralOptObjective<RealT>>(qoi_vol,assembler);

    // Initialize compliance inequality constraint.
    con->value(*rp, *up, *zp, tol);
    RealT objScaling = objFactor, rnorm2 = rp->dot(*rp);
    if (rnorm2 > 1e2*ROL::ROL_EPSILON<RealT>()) {
      objScaling /= rnorm2;
    }
    std::vector<ROL::SharedPointer<QoI<RealT> > > qoi_cmp(1,ROL::nullPointer);
    qoi_cmp[0]
      = ROL::makeShared<QoI_TopoOpt<RealT>>(pde->getFE(), pde->getLoad(),
                                            pde->getFieldHelper(), objScaling);
    ROL::SharedPointer<ROL::Objective_SimOpt<RealT> > cobj
      = ROL::makeShared<PDE_Objective<RealT>>(qoi_cmp,assembler);

    // Initialize reduced compliance objective function.
    bool storage = parlist->sublist("Problem").get("Use state storage",true);
    ROL::SharedPointer<ROL::SimController<RealT> > stateStore
      = ROL::makeShared<ROL::SimController<RealT>>();
    ROL::SharedPointer<ROL::Reduced_Objective_SimOpt<RealT> > cobjRed
      = ROL::makeShared<ROL::Reduced_Objective_SimOpt<RealT>>(cobj,
                     pdeWithFilter,stateStore,up,zp,pp,storage);

    // Create compliance constraint, multiplier and bounds
    RealT comp = cobjRed->value(*zp,tol);
    ROL::SharedPointer<ROL::Constraint<RealT> > icon
      = ROL::makeShared<ROL::ConstraintFromObjective<RealT>>(cobjRed);
    ROL::SharedPointer<std::vector<RealT> > imul_ptr, iup_ptr;
    ROL::SharedPointer<ROL::Vector<RealT> > imul, iup;
    imul = ROL::makeShared<ROL::SingletonVector<RealT>>(0);
    iup  = ROL::makeShared<ROL::SingletonVector<RealT>>(cmpFactor*comp);
    ROL::SharedPointer<ROL::BoundConstraint<RealT> > ibnd
      = ROL::makeShared<ROL::Bounds<RealT>>(*iup,false);

    // Initialize bound constraints.
    ROL::SharedPointer<Tpetra::MultiVector<> > lo_ptr, hi_ptr;
    lo_ptr = assembler->createControlVector(); lo_ptr->putScalar(0.0);
    hi_ptr = assembler->createControlVector(); hi_ptr->putScalar(1.0);
    ROL::SharedPointer<ROL::Vector<RealT> > lop, hip;
    lop = ROL::makeShared<PDE_PrimalOptVector<RealT>>(lo_ptr,pde,assembler);
    hip = ROL::makeShared<PDE_PrimalOptVector<RealT>>(hi_ptr,pde,assembler);
    ROL::SharedPointer<ROL::BoundConstraint<RealT> > bnd
      = ROL::makeShared<ROL::Bounds<RealT>>(lop,hip);

    // Build optimization problem.
    ROL::OptimizationProblem<RealT> optProb(vobj,zp,bnd,icon,imul,ibnd);

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
