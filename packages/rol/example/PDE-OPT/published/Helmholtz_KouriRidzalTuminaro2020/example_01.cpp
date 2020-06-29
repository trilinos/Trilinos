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
    \brief Shows how to solve the optimal control of Helmholtz problem.
*/

#include "Teuchos_Comm.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "Tpetra_Core.hpp"
#include "Tpetra_Version.hpp"

#include <iostream>
#include <algorithm>

#include "ROL_Stream.hpp"
#include "ROL_Ptr.hpp"

#include "../../TOOLS/meshmanager.hpp"
#include "../../TOOLS/meshreader.hpp"
#include "../../TOOLS/assembler.hpp"

#include "pde_helmholtz_real.hpp"
#include "pde_helmholtz_imag.hpp"
#include "obj_helmholtz.hpp"

typedef double RealT;

int main(int argc, char *argv[]) {
  /*** Initialize communicator. ***/
  Teuchos::GlobalMPISession mpiSession (&argc, &argv);
  ROL::Ptr<const Teuchos::Comm<int>> comm
    = Tpetra::getDefaultComm();

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  const int numProcs = (comm->getSize() > 1) ? comm->getSize() : 0;
  const int myRank = comm->getRank();
  ROL::Ptr<std::ostream> outStream = ROL::makeStreamPtr( std::cout, (argc > 1) && (myRank==0) );

  int errorFlag  = 0;

  // *** Example body.
  try {
    //RealT tol(1e-8);// one(1);

    /*** Read in XML input ***/
    std::string filename = "input_ex01.xml";
    Teuchos::RCP<Teuchos::ParameterList> parlist = Teuchos::rcp( new Teuchos::ParameterList() );
    Teuchos::updateParametersFromXmlFile( filename, parlist.ptr() );
    int example = parlist->sublist("Problem").get("Example",1);

    /*** Initialize main data structure. ***/
    ROL::Ptr<MeshManager<RealT>> meshMgr;
    if (example==1) {
      meshMgr = ROL::makePtr<MeshReader<RealT>>(*parlist, numProcs);
    }
    else {
      meshMgr = ROL::makePtr<MeshManager_Rectangle<RealT>>(*parlist);
    }
    // Initialize PDE describing real components of the Helmholtz equation.
    ROL::Ptr<PDE_Helmholtz_Real<RealT>> pde_real
      = ROL::makePtr<PDE_Helmholtz_Real<RealT>>(*parlist);
    ROL::Ptr<Assembler<RealT>> assembler_real
      = ROL::makePtr<Assembler<RealT>>(pde_real->getFields(),meshMgr,comm,*parlist,*outStream);
    assembler_real->setCellNodes(*pde_real);
    // Initialize PDE describing imaginary components of the Helmholtz equation.
    ROL::Ptr<PDE_Helmholtz_Imag<RealT>> pde_imag
      = ROL::makePtr<PDE_Helmholtz_Imag<RealT>>(*parlist);
    ROL::Ptr<Assembler<RealT>> assembler_imag
      = ROL::makePtr<Assembler<RealT>>(pde_imag->getFields(),meshMgr,comm,*parlist,*outStream);
    assembler_imag->setCellNodes(*pde_imag);
    // Initialize objective functions
    ROL::Ptr<QoI<RealT>> qoi_state_real, qoi_state_imag, qoi_ctrl;
    qoi_state_real = ROL::makePtr<QoI_Helmholtz_StateTracking<RealT>>(pde_real->getFE(),
                                                                      *parlist, 0);
    qoi_state_imag = ROL::makePtr<QoI_Helmholtz_StateTracking<RealT>>(pde_real->getFE(),
                                                                      *parlist, 1);
    qoi_ctrl       = ROL::makePtr<QoI_Helmholtz_ControlPenalty<RealT>>(pde_real->getFE(),
                                                                       *parlist);

    // Create state vector.
    ROL::Ptr<Tpetra::MultiVector<>> u_ptr, z_ptr, r_ptr;
    u_ptr = assembler_real->createStateVector();    u_ptr->putScalar(0.0);
    z_ptr = assembler_real->createControlVector();  z_ptr->putScalar(0.0);
    r_ptr = assembler_real->createResidualVector(); r_ptr->putScalar(0.0);

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

    Tpetra::MatrixMarket::Writer< Tpetra::CrsMatrix<>> vecWriter;
    vecWriter.writeDenseFile("WRvector.txt", wr);
    vecWriter.writeDenseFile("WIvector.txt", wi);
    std::string mapfile = "map.txt" + filename;
    vecWriter.writeMapFile("map.txt", *wr->getMap());

    assembler_real->printMeshData(*outStream);

/*
    // Initialize compliance objective function.
    bool storage = parlist->sublist("Problem").get("Use Storage",true);
    std::vector<ROL::Ptr<QoI<RealT> > > qoi_vec(2,ROL::nullPtr);
    qoi_vec[0] = ROL::makePtr<QoI_Helmholtz_StateTracking<RealT>>(pde->getFE(),
                                                                     pde->getFieldHelper(),
                                                                     *parlist);
    qoi_vec[1] = ROL::makePtr<QoI_Helmholtz_ControlPenalty<RealT>>(pde->getFE(),
                                                                      pde->getFieldHelper(),
                                                                      *parlist);
    ROL::Ptr<ROL::Objective_SimOpt<RealT> > obj
      = ROL::makePtr<PDE_Objective<RealT>>(qoi_vec,assembler);
    ROL::Ptr<ROL::Reduced_Objective_SimOpt<RealT> > robj
      = ROL::makePtr<ROL::Reduced_Objective_SimOpt<RealT>>(obj, con, up, zp, pp, storage, false);

    // Run derivative checks
    bool checkDeriv = parlist->sublist("Problem").get("Check derivatives",false);
    if ( checkDeriv ) {
      *outStream << "\n\nCheck Opt Vector\n";
      zp->checkVector(*dzp,*rzp,true,*outStream);

      std::vector<ROL::Ptr<ROL::Objective_SimOpt<RealT> > > obj_vec(2,ROL::nullPtr);
      obj_vec[0] = ROL::makePtr<IntegralObjective<RealT>>(qoi_vec[0],assembler);
      obj_vec[1] = ROL::makePtr<IntegralObjective<RealT>>(qoi_vec[1],assembler);

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

    // Output uncontrolled state.
    zp->zero();
    pdecon->printMeshData(*outStream);
    con->solve(*rp,*up,*zp,tol);
    pdecon->outputTpetraVector(u_ptr,"state_uncontrolled.txt");

    bool useFullSpace = parlist->sublist("Problem").get("Full space",false);
    ROL::Ptr<ROL::Algorithm<RealT> > algo;
    Teuchos::Time algoTimer("Algorithm Time", true);
    if ( useFullSpace ) {
      ROL::OptimizationProblem<RealT> optProb(obj, makePtrFromRef(x), con, rp);
      ROL::OptimizationSolver<RealT> optSolver(optProb, *parlist);
      optSolver.solve(*outStream);
    }
    else {
      ROL::Ptr<ROL::Step<RealT>>
        step = ROL::makePtr<ROL::TrustRegionStep<RealT>>(*parlist);
      ROL::Ptr<ROL::StatusTest<RealT>>
        status = ROL::makePtr<ROL::StatusTest<RealT>>(*parlist);
      ROL::Algorithm<RealT> algo(step,status,false);
      algo.run(*zp,*robj,true,*outStream);
    }
    algoTimer.stop();
    *outStream << "Total optimization time = " << algoTimer.totalElapsedTime() << " seconds.\n";

    // Output.
    pdecon->printMeshData(*outStream);
    con->solve(*rp,*up,*zp,tol);
    pdecon->outputTpetraVector(u_ptr,"state.txt");
    pdecon->outputTpetraVector(z_ptr,"control.txt");
*/
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
