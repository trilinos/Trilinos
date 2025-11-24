// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  example_01.cpp
    \brief Shows how to solve the Poisson-Boltzmann problem.
*/

#include "Teuchos_Comm.hpp"
#include "ROL_GlobalMPISession.hpp"

#include "Tpetra_Core.hpp"
#include "Tpetra_Version.hpp"

#include <iostream>
#include <algorithm>

#include "ROL_Stream.hpp"
#include "ROL_Solver.hpp"
#include "ROL_Reduced_Objective_SimOpt.hpp"
#include "ROL_BoundConstraint_SimOpt.hpp"
#include "ROL_Bounds.hpp"

#include "../TOOLS/meshmanagerK.hpp"
#include "../TOOLS/pdeconstraintK.hpp"
#include "../TOOLS/pdeobjectiveK.hpp"
#include "../TOOLS/pdevectorK.hpp"
#include "pde_allen_cahnK.hpp"
#include "obj_allen_cahnK.hpp"

using RealT = double;
using DeviceT = Kokkos::HostSpace;

int main(int argc, char *argv[]) {
  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  ROL::Ptr<std::ostream> outStream;
  ROL::nullstream bhs; // outputs nothing

  /*** Initialize communicator. ***/
  ROL::GlobalMPISession mpiSession (&argc, &argv, &bhs);
  Kokkos::ScopeGuard kokkosScope (argc, argv);
  ROL::Ptr<const Teuchos::Comm<int> > comm
    = Tpetra::getDefaultComm();
  const int myRank = comm->getRank();
  if ((iprint > 0) && (myRank == 0))
    outStream = ROL::makePtrFromRef(std::cout);
  else
    outStream = ROL::makePtrFromRef(bhs);
  int errorFlag  = 0;

  // *** Example body.
  try {

    /*** Read in XML input ***/
    std::string filename = "input.xml";
    auto parlist = ROL::getParametersFromXmlFile(filename);

    /*** Initialize main data structure. ***/
    auto meshMgr = ROL::makePtr<MeshManager_Rectangle<RealT,DeviceT>>(*parlist);
    // Initialize PDE describe Poisson's equation
    auto pde = ROL::makePtr<PDE_Allen_Cahn<RealT,DeviceT>>(*parlist);
    auto con = ROL::makePtr<PDE_Constraint<RealT,DeviceT>>(pde,meshMgr,comm,*parlist,*outStream);
    auto assembler = con->getAssembler();
    assembler->printMeshData(*outStream);
    con->setSolveParameters(*parlist);

    /*************************************************************************/
    /***************** BUILD VECTORS *****************************************/
    /*************************************************************************/
    auto u_ptr = assembler->createStateVector();
    auto p_ptr = assembler->createStateVector();
    auto r_ptr = assembler->createResidualVector();
    auto z_ptr = assembler->createControlVector();
    u_ptr->randomize();  //u_ptr->putScalar(static_cast<RealT>(1));
    p_ptr->randomize();  //p_ptr->putScalar(static_cast<RealT>(1));
    r_ptr->randomize();  //r_ptr->putScalar(static_cast<RealT>(1));
    z_ptr->randomize();  //z_ptr->putScalar(static_cast<RealT>(1));
    auto up = ROL::makePtr<PDE_PrimalSimVector<RealT,DeviceT>>(u_ptr,pde,assembler,*parlist);
    auto pp = ROL::makePtr<PDE_PrimalSimVector<RealT,DeviceT>>(p_ptr,pde,assembler,*parlist);
    auto rp = ROL::makePtr<PDE_DualSimVector<RealT,DeviceT>>(r_ptr,pde,assembler,*parlist);
    auto zp = ROL::makePtr<PDE_PrimalOptVector<RealT,DeviceT>>(z_ptr,pde,assembler,*parlist);

    /*************************************************************************/
    /***************** BUILD COST FUNCTIONAL *********************************/
    /*************************************************************************/
    // Initialize quadratic objective function
    std::vector<ROL::Ptr<QoI<RealT,DeviceT>>> qoi_vec(2,ROL::nullPtr);
    qoi_vec[0] = ROL::makePtr<QoI_State_Cost_Allen_Cahn<RealT,DeviceT>>(pde->getFE());
    qoi_vec[1] = ROL::makePtr<QoI_Control_Cost_Allen_Cahn<RealT,DeviceT>>(pde->getFE(),
                                                                  pde->getBdryFE(),
                                                                  pde->getBdryCellLocIds());
    std::vector<RealT> weights(2);
    weights[0] = static_cast<RealT>(1);
    weights[1] = parlist->sublist("Problem").get("Control Penalty Parameter", 1e-4);
    // Build full-space objective
    auto obj = ROL::makePtr<PDE_Objective<RealT,DeviceT>>(qoi_vec,weights,assembler);
    // Build reduced-space objective
    bool storage = parlist->sublist("Problem").get("Use state storage",true);
    auto stateStore = ROL::makePtr<ROL::VectorController<RealT>>();
    auto robj = ROL::makePtr<ROL::Reduced_Objective_SimOpt<RealT>>(obj,con,
                                                                stateStore,
                                                                up,zp,pp,
                                                                storage);

    /*************************************************************************/
    /***************** BUILD BOUND CONSTRAINT ********************************/
    /*************************************************************************/
    // Control bounds
    auto zlo_ptr = assembler->createControlVector();
    auto zhi_ptr = assembler->createControlVector();
    RealT lo = parlist->sublist("Problem").get("Lower Bound",0.0);
    RealT hi = parlist->sublist("Problem").get("Upper Bound",1.0);
    zlo_ptr->putScalar(lo); zhi_ptr->putScalar(hi);
    auto zlop = ROL::makePtr<PDE_PrimalOptVector<RealT,DeviceT>>(zlo_ptr,pde,assembler);
    auto zhip = ROL::makePtr<PDE_PrimalOptVector<RealT,DeviceT>>(zhi_ptr,pde,assembler);
    auto zbnd = ROL::makePtr<ROL::Bounds<RealT>>(zlop,zhip);
    bool deactivate = parlist->sublist("Problem").get("Deactivate Bound Constraints",false);
    if (deactivate) zbnd->deactivate();
    // State bounds
    auto ulo_ptr = assembler->createStateVector();
    auto uhi_ptr = assembler->createStateVector();
    ulo_ptr->putScalar(ROL::ROL_NINF<RealT>()); uhi_ptr->putScalar(ROL::ROL_INF<RealT>());
    auto ulop = ROL::makePtr<PDE_PrimalSimVector<RealT,DeviceT>>(ulo_ptr,pde,assembler);
    auto uhip = ROL::makePtr<PDE_PrimalSimVector<RealT,DeviceT>>(uhi_ptr,pde,assembler);
    auto ubnd = ROL::makePtr<ROL::Bounds<RealT>>(ulop,uhip);
    ubnd->deactivate();

    // SimOpt bounds
    auto bnd = ROL::makePtr<ROL::BoundConstraint_SimOpt<RealT>>(ubnd,zbnd);

    // Create ROL SimOpt vectors
    auto du_ptr = assembler->createStateVector();
    auto dz_ptr = assembler->createControlVector();
    du_ptr->randomize(); //du_ptr->putScalar(static_cast<RealT>(0));
    dz_ptr->randomize(); //dz_ptr->putScalar(static_cast<RealT>(1));
    auto dup = ROL::makePtr<PDE_PrimalSimVector<RealT,DeviceT>>(du_ptr,pde,assembler,*parlist);
    auto dzp = ROL::makePtr<PDE_PrimalOptVector<RealT,DeviceT>>(dz_ptr,pde,assembler,*parlist);

    ROL::Vector_SimOpt<RealT> x(up,zp);
    ROL::Vector_SimOpt<RealT> d(dup,dzp);

    /*************************************************************************/
    /***************** RUN VECTOR AND DERIVATIVE CHECKS **********************/
    /*************************************************************************/
    bool checkDeriv = parlist->sublist("Problem").get("Check Derivatives",false);
    if ( checkDeriv ) {
      *outStream << "\n\nCheck Gradient of Full Objective Function\n";
      obj->checkGradient(x,d,true,*outStream);
      *outStream << "\n\nCheck Hessian of Full Objective Function\n";
      obj->checkHessVec(x,d,true,*outStream);
      *outStream << "\n\nCheck Full Jacobian of PDE Constraint\n";
      con->checkApplyJacobian(x,d,*rp,true,*outStream);
      *outStream << "\n\nCheck Jacobian_1 of PDE Constraint\n";
      con->checkApplyJacobian_1(*up,*zp,*dup,*rp,true,*outStream);
      *outStream << "\n\nCheck Jacobian_2 of PDE Constraint\n";
      con->checkApplyJacobian_2(*up,*zp,*dzp,*rp,true,*outStream);
      *outStream << "\n\nCheck Full Hessian of PDE Constraint\n";
      con->checkApplyAdjointHessian(x,*pp,d,x,true,*outStream);
      *outStream << "\n";
      con->checkAdjointConsistencyJacobian(*dup,d,x,true,*outStream);
      *outStream << "\n";
      con->checkInverseJacobian_1(*rp,*dup,*up,*zp,true,*outStream);
      *outStream << "\n";
      *outStream << "\n\nCheck Gradient of Reduced Objective Function\n";
      robj->checkGradient(*zp,*dzp,true,*outStream);
      *outStream << "\n\nCheck Hessian of Reduced Objective Function\n";
      robj->checkHessVec(*zp,*dzp,true,*outStream);
    }

    bool useFullSpace = parlist->sublist("Problem").get("Full space",false);
    ROL::Ptr<ROL::Problem<RealT>> optProb;
    if ( useFullSpace ) {
      optProb = ROL::makePtr<ROL::Problem<RealT>>(obj, makePtrFromRef(x));
      optProb->addBoundConstraint(bnd);
      optProb->addConstraint("PDE", con, pp);
    }
    else {
      optProb = ROL::makePtr<ROL::Problem<RealT>>(robj, zp);
      optProb->addBoundConstraint(zbnd);
    }
    optProb->finalize(false,true,*outStream);
    ROL::Solver<RealT> optSolver(optProb, *parlist);
    optSolver.solve(*outStream);

    // Output.
    RealT tol(1.e-8);
    con->solve(*rp,*up,*zp,tol);
    con->outputTpetraVector(u_ptr,"state.txt");
    con->outputTpetraVector(z_ptr,"control.txt");

    Teuchos::Array<RealT> res(1,0);
    con->value(*rp,*up,*zp,tol);
    r_ptr->norm2(res.view(0,1));
    *outStream << "Residual Norm: " << res[0] << std::endl;
    errorFlag += (res[0] > 1.e-6 ? 1 : 0);
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
