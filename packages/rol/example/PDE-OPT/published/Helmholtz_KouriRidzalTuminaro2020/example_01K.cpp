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
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "ROL_GlobalMPISession.hpp"

#include "Tpetra_Core.hpp"
#include "Tpetra_Version.hpp"

#include <iostream>
#include <algorithm>

#include "ROL_Stream.hpp"
#include "ROL_Ptr.hpp"

#include "../../TOOLS/meshmanagerK.hpp"
#include "../../TOOLS/meshreaderK.hpp"
#include "../../TOOLS/assemblerK.hpp"

#include "pde_helmholtz_realK.hpp"
#include "pde_helmholtz_imagK.hpp"
#include "obj_helmholtzK.hpp"

using RealT = double;
using DeviceT = Kokkos::HostSpace;

int main(int argc, char *argv[]) {
  /*** Initialize communicator. ***/
  ROL::GlobalMPISession mpiSession (&argc, &argv);
  Kokkos::ScopeGuard kokkosScope (argc, argv);
  auto comm = Tpetra::getDefaultComm();

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  const int numProcs = (comm->getSize() > 1) ? comm->getSize() : 0;
  const int myRank = comm->getRank();
  auto outStream = ROL::makeStreamPtr( std::cout, (argc > 1) && (myRank==0) );

  int errorFlag  = 0;

  // *** Example body.
  try {
    //RealT tol(1e-8);// one(1);

    /*** Read in XML input ***/
    std::string filename = "input_ex01.xml";
    auto parlist = ROL::getParametersFromXmlFile(filename);
    int example = parlist->sublist("Problem").get("Example",1);

    /*** Initialize main data structure. ***/
    ROL::Ptr<MeshManager<RealT,DeviceT>> meshMgr;
    if (example==1)
      meshMgr = ROL::makePtr<MeshReader<RealT,DeviceT>>(*parlist, numProcs);
    else
      meshMgr = ROL::makePtr<MeshManager_Rectangle<RealT,DeviceT>>(*parlist);
    // Initialize PDE describing real components of the Helmholtz equation.
    auto pde_real = ROL::makePtr<PDE_Helmholtz_Real<RealT,DeviceT>>(*parlist);
    auto assembler_real = ROL::makePtr<Assembler<RealT,DeviceT>>(pde_real->getFields(),meshMgr,comm,*parlist,*outStream);
    assembler_real->setCellNodes(*pde_real);
    // Initialize PDE describing imaginary components of the Helmholtz equation.
    auto pde_imag = ROL::makePtr<PDE_Helmholtz_Imag<RealT,DeviceT>>(*parlist);
    auto assembler_imag = ROL::makePtr<Assembler<RealT,DeviceT>>(pde_imag->getFields(),meshMgr,comm,*parlist,*outStream);
    assembler_imag->setCellNodes(*pde_imag);
    // Initialize objective functions
    auto qoi_state_real = ROL::makePtr<QoI_Helmholtz_StateTracking<RealT,DeviceT>>(pde_real->getFE(), *parlist, 0);
    auto qoi_state_imag = ROL::makePtr<QoI_Helmholtz_StateTracking<RealT,DeviceT>>(pde_real->getFE(), *parlist, 1);
    auto qoi_ctrl       = ROL::makePtr<QoI_Helmholtz_ControlPenalty<RealT,DeviceT>>(pde_real->getFE(), *parlist);

    // Create state vector.
    auto u_ptr = assembler_real->createStateVector();    u_ptr->putScalar(0.0);
    auto z_ptr = assembler_real->createControlVector();  z_ptr->putScalar(0.0);
    auto r_ptr = assembler_real->createResidualVector(); r_ptr->putScalar(0.0);

    ROL::Ptr<Tpetra::CrsMatrix<>> A, B, L, M, C, R;
    ROL::Ptr<Tpetra::MultiVector<>> wr, wi;
    assembler_real->assemblePDEJacobian1(A,pde_real,u_ptr,z_ptr,ROL::nullPtr);
    assembler_imag->assemblePDEJacobian1(B,pde_imag,u_ptr,z_ptr,ROL::nullPtr);
    assembler_real->assemblePDEJacobian2(L,pde_real,u_ptr,z_ptr,ROL::nullPtr);
    assembler_real->assemblePDERieszMap2(M,pde_real);
    assembler_real->assembleQoIHessian11(C,qoi_state_real,u_ptr,z_ptr,ROL::nullPtr);
    assembler_real->assembleQoIHessian22(R,qoi_ctrl,u_ptr,z_ptr,ROL::nullPtr);
    assembler_real->assembleQoIGradient1(wr,qoi_state_real,u_ptr,z_ptr,ROL::nullPtr);
    assembler_real->assembleQoIGradient1(wi,qoi_state_imag,u_ptr,z_ptr,ROL::nullPtr);
    wr->scale(-1.0); wi->scale(-1.0);

    Tpetra::MatrixMarket::Writer<Tpetra::CrsMatrix<>> matWriter;
    matWriter.writeSparseFile("Amatrix.txt",A);
    matWriter.writeSparseFile("Bmatrix.txt",B);
    matWriter.writeSparseFile("Lmatrix.txt",L);
    matWriter.writeSparseFile("Mmatrix.txt",M);
    matWriter.writeSparseFile("Cmatrix.txt",C);
    matWriter.writeSparseFile("Rmatrix.txt",R);

    Tpetra::MatrixMarket::Writer<Tpetra::CrsMatrix<>> vecWriter;
    vecWriter.writeDenseFile("WRvector.txt", wr);
    vecWriter.writeDenseFile("WIvector.txt", wi);
    std::string mapfile = "map.txt" + filename;
    vecWriter.writeMapFile("map.txt", *wr->getMap());

    assembler_real->printMeshData(*outStream);

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
