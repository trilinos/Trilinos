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

#include "Tpetra_Core.hpp"
#include "Tpetra_Version.hpp"

#include <iostream>
#include <algorithm>

#include "ROL_Solver.hpp"
#include "ROL_Bounds.hpp"
#include "ROL_ConstraintFromObjective.hpp"
#include "ROL_LinearCombinationObjective.hpp"

#include "../../../TOOLS/pdevector.hpp"
#include "../../../TOOLS/pdeobjective.hpp"
#include "../../../TOOLS/meshreader.hpp"

#include "src/pde_elasticity.hpp"
#include "src/pde_filter.hpp"
#include "src/filtered_compliance_robj.hpp"
#include "src/volume_con.hpp"
#include "src/obj_volume.hpp"

typedef double RealT;

int main(int argc, char *argv[]) {
  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  ROL::Ptr<std::ostream> outStream;
  ROL::nullstream bhs; // outputs nothing

  /*** Initialize communicator. ***/
  Teuchos::GlobalMPISession mpiSession (&argc, &argv, &bhs);
  ROL::Ptr<const Teuchos::Comm<int>> comm
    = Tpetra::getDefaultComm();
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
    int probDim               = parlist->sublist("Problem").get("Problem Dimension", 2);
    const RealT volFraction   = parlist->sublist("Problem").get("Volume Fraction", 0.4);
    const RealT initDens      = volFraction;
    const std::string minType  = "Compliance";
    const std::string hessAppr = "None";
    if (probDim==3) {
      parlist->sublist("Mesh").set("File Name","meshfiles/truss3d.txt");
    }
    else {
      probDim = 2;
      parlist->sublist("Problem").set("Problem Dimension",probDim);
      parlist->sublist("Mesh").set("File Name","meshfiles/truss2d.txt");
    }
    parlist->sublist("Problem").set("Use Inexact Linear Solves",false);
    parlist->sublist("General").set("Inexact Objective Function",false);
    parlist->sublist("General").set("Inexact Gradient",false);

    /*** Initialize main data structure. ***/
    int nProcs = comm->getSize();
    if (nProcs == 1) nProcs = 0;
    ROL::Ptr<MeshManager<RealT>>
      meshMgr = ROL::makePtr<MeshReader<RealT>>(*parlist,nProcs);
    // Initialize PDE describing elasticity equations.
    ROL::Ptr<PDE_Elasticity<RealT>>
      pde = ROL::makePtr<PDE_Elasticity<RealT>>(*parlist);

    // Initialize the filter PDE.
    ROL::Ptr<PDE_Filter<RealT>> 
      pdeFilter = ROL::makePtr<PDE_Filter<RealT>>(*parlist);
    pde->setDensityFields(pdeFilter->getFields());

    // Initialize reduced compliance objective function.
    ROL::Ptr<TopOptFilteredComplianceObjective<RealT>>
      robj_com = ROL::makePtr<TopOptFilteredComplianceObjective<RealT>>(
                 pde,pdeFilter,meshMgr,comm,*parlist,*outStream);
    ROL::Ptr<Assembler<RealT>> assembler, assemblerFilter;
    assembler = robj_com->getAssembler();
    assemblerFilter = robj_com->getFilterAssembler();

    // Create vectors.
    ROL::Ptr<Tpetra::MultiVector<>> u_ptr, p_ptr, r_ptr, z_ptr;
    u_ptr = assembler->createStateVector();         u_ptr->putScalar(0.0);
    p_ptr = assembler->createStateVector();         p_ptr->putScalar(0.0);
    r_ptr = assembler->createResidualVector();      r_ptr->putScalar(0.0);
    z_ptr = assemblerFilter->createControlVector(); z_ptr->putScalar(1.0);
    ROL::Ptr<ROL::Vector<RealT>> up, pp, rp, zp, imul;
    up = ROL::makePtr<PDE_PrimalSimVector<RealT>>(u_ptr,pde,assembler,*parlist);
    pp = ROL::makePtr<PDE_PrimalSimVector<RealT>>(p_ptr,pde,assembler,*parlist);
    rp = ROL::makePtr<PDE_DualSimVector<RealT>>(r_ptr,pde,assembler,*parlist);
    zp = ROL::makePtr<PDE_PrimalOptVector<RealT>>(z_ptr,pdeFilter,assemblerFilter,*parlist);
    imul = ROL::makePtr<ROL::SingletonVector<RealT>>(0);

    // Build volume objective function.
    ROL::Ptr<QoI<RealT>>
      qoi_vol = ROL::makePtr<QoI_Volume_TopoOpt<RealT>>(pdeFilter->getDensityFE(),
                                                        volFraction);
    ROL::Ptr<TopOptVolumeConstraint<RealT>>
      con_vol = ROL::makePtr<TopOptVolumeConstraint<RealT>>(qoi_vol,assemblerFilter,zp);

    // Normalize compliance objective function
    zp->setScalar(initDens);
    RealT cs = ROL::dynamicPtrCast<TopOptFilteredComplianceObjective<RealT>>(robj_com)->normalize(*zp,tol);

    // Output problem details
    *outStream << std::endl;
    *outStream << "Problem Data"          << std::endl;
    *outStream << "  Dimension:         " << probDim << std::endl;
    *outStream << "  SIMP Power:        "
               << parlist->sublist("Problem").get("SIMP Power",3.0) << std::endl;
    *outStream << "  Volume Fraction:   " << volFraction << std::endl;
    *outStream << "  Initial Density:   " << initDens    << std::endl;
    *outStream << "  Compliance Scale:  " << cs  << std::endl;
    *outStream << std::endl;

    // Create objective, constraint, multiplier and bounds

    // Initialize bound constraints.
    RealT lval = 0.0, uval = 1.0;
    ROL::Ptr<ROL::Vector<RealT>> lop = zp->clone(); lop->setScalar(lval);
    ROL::Ptr<ROL::Vector<RealT>> hip = zp->clone(); hip->setScalar(uval);
    ROL::Ptr<ROL::BoundConstraint<RealT>>
      bnd = ROL::makePtr<ROL::Bounds<RealT>>(lop,hip);

    // Set up optimization problem.
    ROL::Ptr<ROL::Problem<RealT>>
      prob = ROL::makePtr<ROL::Problem<RealT>>(robj_com,zp);
    prob->addBoundConstraint(bnd);
    prob->addLinearConstraint("Volume",con_vol,imul);
    prob->setProjectionAlgorithm(*parlist);

    // Check derivatives.
    bool derivCheck = parlist->sublist("Problem").get("Check derivatives",false);
    if (derivCheck) prob->check(true,*outStream);

    // Solve optimization problem.
    ROL::Ptr<ROL::Solver<RealT>> solver;
    std::vector<RealT> vals(6);

    // SPG
    zp->setScalar(initDens);
    imul->zero();
    robj_com->reset();
    con_vol->reset();
    prob->finalize(false,true,*outStream);
    parlist->sublist("Step").set("Type","Spectral Gradient");
    Teuchos::Time spgTimer("SPG Time", true);
    solver = ROL::makePtr<ROL::Solver<RealT>>(prob,*parlist);
    solver->solve(*outStream);
    spgTimer.stop();
    *outStream << "Total optimization time = " << spgTimer.totalElapsedTime() << " seconds.\n";
    // Output.
    robj_com->summarize(*outStream);
    con_vol->summarize(*outStream);
    vals[0] = solver->getAlgorithmState()->value;

    // PQN
    zp->setScalar(initDens);
    imul->zero();
    robj_com->reset();
    con_vol->reset();
    parlist->sublist("Step").set("Type","Line Search");
    parlist->sublist("Step").sublist("Line Search").sublist("Descent Method").set("Type","Quasi-Newton Method");
    Teuchos::Time pqnTimer("PQN Time", true);
    solver = ROL::makePtr<ROL::Solver<RealT>>(prob,*parlist);
    solver->solve(*outStream);
    pqnTimer.stop();
    *outStream << "Total optimization time = " << pqnTimer.totalElapsedTime() << " seconds.\n";
    // Output.
    robj_com->summarize(*outStream);
    con_vol->summarize(*outStream);
    vals[1] = solver->getAlgorithmState()->value;

    // LMTR
    zp->setScalar(initDens);
    imul->zero();
    robj_com->reset();
    con_vol->reset();
    parlist->sublist("Step").set("Type","Trust Region");
    parlist->sublist("Step").sublist("Trust Region").set("Subproblem Model","Lin-More");
    Teuchos::Time lmtrTimer("LMTR Time", true);
    solver = ROL::makePtr<ROL::Solver<RealT>>(prob,*parlist);
    solver->solve(*outStream);
    lmtrTimer.stop();
    *outStream << "Total optimization time = " << lmtrTimer.totalElapsedTime() << " seconds.\n";
    // Output.
    robj_com->summarize(*outStream);
    con_vol->summarize(*outStream);
    vals[2] = solver->getAlgorithmState()->value;

    // TRSPG
    zp->setScalar(initDens);
    imul->zero();
    robj_com->reset();
    con_vol->reset();
    parlist->sublist("Step").set("Type","Trust Region");
    parlist->sublist("Step").sublist("Trust Region").set("Subproblem Model","SPG");
    Teuchos::Time trspgTimer("TRSPG Time", true);
    solver = ROL::makePtr<ROL::Solver<RealT>>(prob,*parlist);
    solver->solve(*outStream);
    trspgTimer.stop();
    *outStream << "Total optimization time = " << trspgTimer.totalElapsedTime() << " seconds.\n";
    // Output.
    robj_com->summarize(*outStream);
    con_vol->summarize(*outStream);
    vals[3] = solver->getAlgorithmState()->value;

    // AL-LMTR
    zp->setScalar(initDens);
    imul->zero();
    robj_com->reset();
    con_vol->reset();
    prob->edit();
    prob->finalize(true,true,*outStream);
    parlist->sublist("Step").set("Type","Augmented Lagrangian");
    parlist->sublist("Step").sublist("Augmented Lagrangian").set("Subproblem Step Type","Trust Region");
    parlist->sublist("Step").sublist("Trust Region").set("Subproblem Model","Lin-More");
    Teuchos::Time allmtrTimer("AL-LMTR Time", true);
    solver = ROL::makePtr<ROL::Solver<RealT>>(prob,*parlist);
    solver->solve(*outStream);
    allmtrTimer.stop();
    *outStream << "Total optimization time = " << allmtrTimer.totalElapsedTime() << " seconds.\n";
    // Output.
    robj_com->summarize(*outStream);
    con_vol->summarize(*outStream);
    vals[4] = solver->getAlgorithmState()->value;

    // AL-LMTR
    zp->setScalar(initDens);
    imul->zero();
    robj_com->reset();
    con_vol->reset();
    parlist->sublist("Step").set("Type","Augmented Lagrangian");
    parlist->sublist("Step").sublist("Augmented Lagrangian").set("Subproblem Step Type","Trust Region");
    parlist->sublist("Step").sublist("Trust Region").set("Subproblem Model","SPG");
    Teuchos::Time altrspgTimer("AL-TRSPG Time", true);
    solver = ROL::makePtr<ROL::Solver<RealT>>(prob,*parlist);
    solver->solve(*outStream);
    altrspgTimer.stop();
    *outStream << "Total optimization time = " << altrspgTimer.totalElapsedTime() << " seconds.\n";
    // Output.
    robj_com->summarize(*outStream);
    con_vol->summarize(*outStream);
    vals[5] = solver->getAlgorithmState()->value;

    RealT minval = std::min(vals[0],std::min(vals[1],std::min(vals[2],vals[3])));
    *outStream << std::scientific << std::setprecision(16);
    *outStream << "SPG      " << (vals[0]-minval)/minval << std::endl;
    *outStream << "PQN      " << (vals[1]-minval)/minval << std::endl;
    *outStream << "LMTR     " << (vals[2]-minval)/minval << std::endl;
    *outStream << "TRSPG    " << (vals[3]-minval)/minval << std::endl;
    *outStream << "AL-LMTR  " << (vals[4]-minval)/minval << std::endl;
    *outStream << "AL-TRSPG " << (vals[5]-minval)/minval << std::endl;
    *outStream << std::endl;

    robj_com->printToFile(*zp,*outStream);

    // Get a summary from the time monitor.
    Teuchos::TimeMonitor::summarize();
  }
  catch (std::logic_error& err) {
    *outStream << err.what() << "\n";
    errorFlag = -1000;
  }; // end try

  if (errorFlag != 0)
    std::cout << "End Result: TEST FAILED\n";
  else
    std::cout << "End Result: TEST PASSED\n";

  return 0;
}
