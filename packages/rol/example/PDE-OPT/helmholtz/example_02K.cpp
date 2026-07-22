// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  example_01.cpp
    \brief Shows how to solve the optimal control of Helmholtz problem.
*/

#include "Teuchos_Comm.hpp"
#include "Teuchos_Time.hpp"
#include "ROL_Stream.hpp"
#include "ROL_GlobalMPISession.hpp"

#include "Tpetra_Core.hpp"
#include "Tpetra_Version.hpp"

#include <iostream>
#include <algorithm>

#include "ROL_Solver.hpp"
#include "ROL_Reduced_Objective_SimOpt.hpp"

#include "../TOOLS/linearpdeconstraintK.hpp"
#include "../TOOLS/pdeobjectiveK.hpp"
#include "../TOOLS/pdevectorK.hpp"
#include "../TOOLS/meshmanagerK.hpp"
#include "../TOOLS/sysbuilder.hpp"

#include "pde_helmholtzK.hpp"
#include "obj_helmholtzK.hpp"

using RealT = double;
using DeviceT = Kokkos::HostSpace;

int main(int argc, char *argv[]) {
  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  ROL::Ptr<std::ostream> outStream;
  ROL::nullstream bhs; // outputs nothing

  /*** Initialize communicator. ***/
  ROL::GlobalMPISession mpiSession (&argc, &argv, &bhs);
  Kokkos::ScopeGuard kokkosScope(argc, argv);
  ROL::Ptr<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
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
    RealT tol(1e-8);// one(1);

    /*** Read in XML input ***/
    std::string filename = "input.xml";
    auto parlist = Teuchos::rcp( new Teuchos::ParameterList() );
    Teuchos::updateParametersFromXmlFile( filename, parlist.ptr() );

    /*** Initialize main data structure. ***/
    auto meshMgr = ROL::makePtr<MeshManager_Rectangle<RealT, DeviceT>>(*parlist);
    // Initialize PDE describing elasticity equations.
    auto pde = ROL::makePtr<PDE_Helmholtz<RealT, DeviceT>>(*parlist);
    ROL::Ptr<ROL::Constraint_SimOpt<RealT>> con
      = ROL::makePtr<Linear_PDE_Constraint<RealT, DeviceT>>(pde, meshMgr, comm, *parlist, *outStream);
    // Cast the constraint and get the assembler.
    ROL::Ptr<Linear_PDE_Constraint<RealT, DeviceT>> pdecon
      = ROL::dynamicPtrCast<Linear_PDE_Constraint<RealT, DeviceT>>(con);
    ROL::Ptr<Assembler<RealT, DeviceT>> assembler = pdecon->getAssembler();
    con->setSolveParameters(*parlist);

    // Create state vector.
    auto u_ptr = assembler->createStateVector();
    u_ptr->randomize();
    auto up = ROL::makePtr<PDE_PrimalSimVector<RealT, DeviceT>>(u_ptr,pde,assembler,*parlist);
    auto p_ptr = assembler->createStateVector();
    p_ptr->randomize();
    auto pp = ROL::makePtr<PDE_PrimalSimVector<RealT, DeviceT>>(p_ptr,pde,assembler,*parlist);
    // Create control vector.
    auto z_ptr = assembler->createControlVector();
    z_ptr->randomize();
    auto zp = ROL::makePtr<PDE_PrimalOptVector<RealT, DeviceT>>(z_ptr,pde,assembler,*parlist);
    // Create residual vector.
    auto r_ptr = assembler->createResidualVector();
    r_ptr->putScalar(0.0);
    auto rp = ROL::makePtr<PDE_DualSimVector<RealT, DeviceT>>(r_ptr,pde,assembler,*parlist);
    // Create state direction vector.
    auto du_ptr = assembler->createStateVector();
    du_ptr->randomize();
    //du_ptr->putScalar(0);
    auto dup = ROL::makePtr<PDE_PrimalSimVector<RealT, DeviceT>>(du_ptr,pde,assembler,*parlist);
    // Create control direction vector.
    auto dz_ptr = assembler->createControlVector();
    dz_ptr->randomize();
    //dz_ptr->putScalar(0);
    auto dzp = ROL::makePtr<PDE_PrimalOptVector<RealT, DeviceT>>(dz_ptr,pde,assembler,*parlist);
    // Create control test vector.
    auto rz_ptr = assembler->createControlVector();
    rz_ptr->randomize();
    auto rzp = ROL::makePtr<PDE_PrimalOptVector<RealT, DeviceT>>(rz_ptr,pde,assembler,*parlist);

    auto dualu_ptr = assembler->createStateVector();
    auto dualup = ROL::makePtr<PDE_DualSimVector<RealT, DeviceT>>(dualu_ptr,pde,assembler,*parlist);
    auto dualz_ptr = assembler->createControlVector();
    auto dualzp = ROL::makePtr<PDE_DualOptVector<RealT, DeviceT>>(dualz_ptr,pde,assembler,*parlist);

    // Create ROL SimOpt vectors.
    ROL::Vector_SimOpt<RealT> x(up,zp);
    ROL::Vector_SimOpt<RealT> d(dup,dzp);

    // Initialize compliance objective function.
    bool storage = parlist->sublist("Problem").get("Use Storage",true);
    std::vector<ROL::Ptr<QoI<RealT, DeviceT>>> qoi_vec(2, ROL::nullPtr);
    qoi_vec[0] = ROL::makePtr<QoI_Helmholtz_StateTracking<RealT, DeviceT>>(pde->getFE(),
                                                                           pde->getFieldHelper(),
                                                                           *parlist);
    qoi_vec[1] = ROL::makePtr<QoI_Helmholtz_ControlPenalty<RealT, DeviceT>>(pde->getFE(),
                                                                            pde->getFieldHelper(),
                                                                            *parlist);
    auto obj = ROL::makePtr<PDE_Objective<RealT, DeviceT>>(qoi_vec, assembler);
    auto robj = ROL::makePtr<ROL::Reduced_Objective_SimOpt<RealT>>(obj, con, up, zp, pp, storage, false);

    // Run derivative checks
    bool checkDeriv = parlist->sublist("Problem").get("Check derivatives",false);
    if ( checkDeriv ) {
      *outStream << "\n\nCheck Opt Vector\n";
      zp->checkVector(*dzp, *rzp, true, *outStream);

      std::vector<ROL::Ptr<ROL::Objective_SimOpt<RealT> > > obj_vec(2,ROL::nullPtr);
      obj_vec[0] = ROL::makePtr<IntegralObjective<RealT, DeviceT>>(qoi_vec[0], assembler);
      obj_vec[1] = ROL::makePtr<IntegralObjective<RealT, DeviceT>>(qoi_vec[1], assembler);

      *outStream << "\n\nCheck Gradient of State Objective Function\n";
      obj_vec[0]->checkGradient(x,d,true,*outStream);
      *outStream << "\n\nCheck Gradient_1 of State Objective Function\n";
      obj_vec[0]->checkGradient_1(*up,*zp,*dup,true,*outStream);
      *outStream << "\n\nCheck Gradient_2 of State Objective Function\n";
      obj_vec[0]->checkGradient_2(*up,*zp,*dzp,true,*outStream);
      *outStream << "\n\nCheck Hessian of State Objective Function\n";
      obj_vec[0]->checkHessVec(x,d,true,*outStream);
      *outStream << "\n\nCheck Hessian_11 of State Objective Function\n";
      obj_vec[0]->checkHessVec_11(*up,*zp,*dup,true,*outStream);
      *outStream << "\n\nCheck Hessian_12 of State Objective Function\n";
      obj_vec[0]->checkHessVec_12(*up,*zp,*dzp,true,*outStream);
      *outStream << "\n\nCheck Hessian_21 of State Objective Function\n";
      obj_vec[0]->checkHessVec_21(*up,*zp,*dup,true,*outStream);
      *outStream << "\n\nCheck Hessian_22 of State Objective Function\n";
      obj_vec[0]->checkHessVec_22(*up,*zp,*dzp,true,*outStream);

      *outStream << "\n\nCheck Gradient of Control Objective Function\n";
      obj_vec[1]->checkGradient(x,d,true,*outStream);
      *outStream << "\n\nCheck Gradient_1 of Control Objective Function\n";
      obj_vec[1]->checkGradient_1(*up,*zp,*dup,true,*outStream);
      *outStream << "\n\nCheck Gradient_2 of Control Objective Function\n";
      obj_vec[1]->checkGradient_2(*up,*zp,*dzp,true,*outStream);
      *outStream << "\n\nCheck Hessian of Control Objective Function\n";
      obj_vec[1]->checkHessVec(x,d,true,*outStream);
      *outStream << "\n\nCheck Hessian_11 of Control Objective Function\n";
      obj_vec[1]->checkHessVec_11(*up,*zp,*dup,true,*outStream);
      *outStream << "\n\nCheck Hessian_12 of State Objective Function\n";
      obj_vec[1]->checkHessVec_12(*up,*zp,*dzp,true,*outStream);
      *outStream << "\n\nCheck Hessian_21 of State Objective Function\n";
      obj_vec[1]->checkHessVec_21(*up,*zp,*dup,true,*outStream);
      *outStream << "\n\nCheck Hessian_22 of Control Objective Function\n";
      obj_vec[1]->checkHessVec_22(*up,*zp,*dzp,true,*outStream);

      *outStream << "\n\nCheck Gradient of Full Objective Function\n";
      obj->checkGradient(x,d,true,*outStream);
      *outStream << "\n\nCheck Gradient_1 of Full Objective Function\n";
      obj->checkGradient_1(*up,*zp,*dup,true,*outStream);
      *outStream << "\n\nCheck Gradient_2 of Full Objective Function\n";
      obj->checkGradient_2(*up,*zp,*dzp,true,*outStream);
      *outStream << "\n\nCheck Hessian of Full Objective Function\n";
      obj->checkHessVec(x,d,true,*outStream);
      *outStream << "\n\nCheck Hessian_11 of Full Objective Function\n";
      obj->checkHessVec_11(*up,*zp,*dup,true,*outStream);
      *outStream << "\n\nCheck Hessian_12 of State Objective Function\n";
      obj->checkHessVec_12(*up,*zp,*dzp,true,*outStream);
      *outStream << "\n\nCheck Hessian_21 of State Objective Function\n";
      obj->checkHessVec_21(*up,*zp,*dup,true,*outStream);
      *outStream << "\n\nCheck Hessian_22 of Full Objective Function\n";
      obj->checkHessVec_22(*up,*zp,*dzp,true,*outStream);

      *outStream << "\n\nCheck Full Jacobian of PDE Constraint\n";
      con->checkApplyJacobian(x,d,*rp,true,*outStream);
      *outStream << "\n\nCheck Jacobian_1 of PDE Constraint\n";
      con->checkApplyJacobian_1(*up,*zp,*dup,*rp,true,*outStream);
      *outStream << "\n\nCheck Jacobian_2 of PDE Constraint\n";
      con->checkApplyJacobian_2(*up,*zp,*dzp,*rp,true,*outStream);
      *outStream << "\n\nCheck Full Hessian of PDE Constraint\n";
      con->checkApplyAdjointHessian(x,*pp,d,x,true,*outStream);
      *outStream << "\n\nCheck Hessian_11 of PDE Constraint\n";
      con->checkApplyAdjointHessian_11(*up,*zp,*pp,*dup,*dualup,true,*outStream);
      *outStream << "\n\nCheck Hessian_21 of PDE Constraint\n";
      con->checkApplyAdjointHessian_21(*up,*zp,*pp,*dzp,*dualup,true,*outStream);
      *outStream << "\n\nCheck Hessian_12 of PDE Constraint\n";
      con->checkApplyAdjointHessian_12(*up,*zp,*pp,*dup,*dualzp,true,*outStream);
      *outStream << "\n\nCheck Hessian_22 of PDE Constraint\n";
      con->checkApplyAdjointHessian_22(*up,*zp,*pp,*dzp,*dualzp,true,*outStream);
      *outStream << "\n";
      con->checkAdjointConsistencyJacobian(*dup,d,x,true,*outStream);
      *outStream << "\n";
      con->checkInverseJacobian_1(*up,*up,*up,*zp,true,*outStream);
      *outStream << "\n";
      con->checkInverseAdjointJacobian_1(*up,*up,*up,*zp,true,*outStream);

      *outStream << "\n\nCheck Gradient of Reduced Objective Function\n";
      robj->checkGradient(*zp,*dzp,true,*outStream);
      *outStream << "\n\nCheck Hessian of Reduced Objective Function\n";
      robj->checkHessVec(*zp,*dzp,true,*outStream);
    }

    // Build KKT system.
    ROL::Ptr<Tpetra::CrsMatrix<>> A;
    ROL::Ptr<Tpetra::CrsMatrix<>> B;
    ROL::Ptr<Tpetra::CrsMatrix<>> M;
    ROL::Ptr<Tpetra::CrsMatrix<>> R;
    ROL::Ptr<Tpetra::MultiVector<>> simres;
    assembler->assemblePDEJacobian1(A, pde, u_ptr);
    assembler->assemblePDEJacobian2(B, pde, u_ptr);
    assembler->assemblePDERieszMap1(M, pde);
    assembler->assemblePDERieszMap2(R, pde);
    assembler->assemblePDEResidual(simres, pde, u_ptr, z_ptr);
    SysBuilder<RealT> kktbuilder;
    ROL::Ptr<Tpetra::CrsMatrix<>>   KKT;
    ROL::Ptr<Tpetra::MultiVector<>> BKKT;
    ROL::Ptr<Tpetra::MultiVector<>> XKKT;
    Teuchos::Time kktTimer("KKT Time", true);
    kktbuilder.buildMatrix(KKT, A, B, M, ROL::nullPtr, ROL::nullPtr, R);
    kktTimer.stop();
    *outStream << "KKT matrix assembly time = " << kktTimer.totalElapsedTime() << " seconds.\n";
    kktTimer.reset(); kktTimer.start();
    kktbuilder.buildVector(BKKT, simres, simres, simres);
    kktTimer.stop();
    *outStream << "KKT RHS assembly time = " << kktTimer.totalElapsedTime() << " seconds.\n";
    kktTimer.reset(); kktTimer.start();
    kktbuilder.buildVector(XKKT, simres, simres, simres);
    kktTimer.stop();
    *outStream << "KKT XVEC assembly time = " << kktTimer.totalElapsedTime() << " seconds.\n";
    //Tpetra::MatrixMarket::Writer< Tpetra::CrsMatrix<> > matWriter;
    //if (KKT != ROL::nullPtr) {
    //  matWriter.writeSparseFile("kkt.txt", KKT);
    //}
    Solver<RealT> kktsolver(parlist->sublist("Solver"));
    kktTimer.reset(); kktTimer.start();
    kktsolver.setA(KKT);
    kktTimer.stop();
    *outStream << "KKT factorization time = " << kktTimer.totalElapsedTime() << " seconds.\n";
    kktTimer.reset(); kktTimer.start();
    kktsolver.solve(XKKT, BKKT);
    kktTimer.stop();
    *outStream << "KKT solution time = " << kktTimer.totalElapsedTime() << " seconds.\n";

    // Schur matrix.
    Teuchos::Time schurTimer("Schur Time", true);
    ROL::Ptr<Tpetra::CrsMatrix<>> AAT =
      ROL::makePtr<Tpetra::CrsMatrix<>>(A->getRowMap(), 2*A->getGlobalMaxNumRowEntries());
    ROL::Ptr<Tpetra::CrsMatrix<>> BBT =
      ROL::makePtr<Tpetra::CrsMatrix<>>(B->getRowMap(), 2*B->getGlobalMaxNumRowEntries());
    ROL::Ptr<Tpetra::CrsMatrix<>> C;
    Tpetra::MatrixMatrix::Multiply(*A, false, *A, true, *AAT);
    Tpetra::MatrixMatrix::Multiply(*B, false, *B, true, *BBT);
    C = Tpetra::MatrixMatrix::add(1.0, false, *AAT, 1.0, false, *BBT);
    ROL::Ptr<Tpetra::MultiVector<>> Xschur;
    ROL::Ptr<Tpetra::MultiVector<>> Bschur;
    assembler->assemblePDEResidual(Xschur, pde, u_ptr, z_ptr);
    assembler->assemblePDEResidual(Bschur, pde, u_ptr, z_ptr);
    schurTimer.stop();
    *outStream << "Schur complement assembly time = " << schurTimer.totalElapsedTime() << " seconds.\n";
    Solver<RealT> schursolver(parlist->sublist("Solver"));
    schurTimer.reset(); schurTimer.start();
    schursolver.setA(C);
    schurTimer.stop();
    *outStream << "Schur complement factorization time = " << schurTimer.totalElapsedTime() << " seconds.\n";
    schurTimer.reset(); schurTimer.start();
    schursolver.solve(Xschur, Bschur);
    schurTimer.stop();
    *outStream << "Schur solution time = " << schurTimer.totalElapsedTime() << " seconds.\n";


    auto prob = ROL::makePtr<ROL::Problem<RealT>>(robj, zp);
    ROL::Solver<RealT> solver(prob, *parlist);
    Teuchos::Time algoTimer("Algorithm Time", true);
    solver.solve(*outStream);
    algoTimer.stop();
    *outStream << "Total optimization time = " << algoTimer.totalElapsedTime() << " seconds.\n";

    // Output.
    pdecon->printMeshData(*outStream);
    con->solve(*rp, *up, *zp, tol);
    pdecon->outputTpetraVector(u_ptr, "state.txt");
    pdecon->outputTpetraVector(z_ptr, "control.txt");

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
