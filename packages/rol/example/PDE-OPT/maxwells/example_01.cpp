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
    \brief Shows how to solve the optimal control of Maxwells problem.
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
#include "ROL_Reduced_Objective_SimOpt.hpp"

#include "../TOOLS/linearpdeconstraint.hpp"
#include "../TOOLS/pdeobjective.hpp"
#include "../TOOLS/pdevector.hpp"
#include "../TOOLS/meshmanager.hpp"

#include "pde_maxwells.hpp"
#include "obj_maxwells.hpp"

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
    RealT tol(1e-8);// one(1);

    /*** Read in XML input ***/
    std::string filename = "input.xml";
    Teuchos::RCP<Teuchos::ParameterList> parlist = Teuchos::rcp( new Teuchos::ParameterList() );
    Teuchos::updateParametersFromXmlFile( filename, parlist.ptr() );

    /*** Initialize main data structure. ***/
    ROL::Ptr<MeshManager<RealT> > meshMgr
      = ROL::makePtr<MeshManager_Brick<RealT>>(*parlist);
    // Initialize PDE describing elasticity equations.
    ROL::Ptr<PDE_Maxwells<RealT> > pde
      = ROL::makePtr<PDE_Maxwells<RealT>>(*parlist);
    ROL::Ptr<ROL::Constraint_SimOpt<RealT> > con
      = ROL::makePtr<Linear_PDE_Constraint<RealT>>(pde,meshMgr,comm,*parlist,*outStream);
    // Cast the constraint and get the assembler.
    ROL::Ptr<Linear_PDE_Constraint<RealT> > pdecon
      = ROL::dynamicPtrCast<Linear_PDE_Constraint<RealT> >(con);
    ROL::Ptr<Assembler<RealT> > assembler = pdecon->getAssembler();
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
    z_ptr->randomize();
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
    //dz_ptr->putScalar(0);
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

//    // Initialize compliance objective function.
//    bool storage = parlist->sublist("Problem").get("Use Storage",true);
//    std::vector<ROL::Ptr<QoI<RealT> > > qoi_vec(2,ROL::nullPtr);
//    qoi_vec[0] = ROL::makePtr<QoI_Maxwells_StateTracking<RealT>>(pde->getFE(),
//                                                                    pde->getFieldHelper(),
//                                                                    *parlist));
//    qoi_vec[1] = ROL::makePtr<QoI_Maxwells_ControlPenalty<RealT>>(pde->getFE(),
//                                                                     pde->getFieldHelper(),
//                                                                     *parlist));
//    ROL::Ptr<ROL::Objective_SimOpt<RealT> > obj
//      = ROL::makePtr<PDE_Objective<RealT>>(qoi_vec,assembler);
//    ROL::Ptr<ROL::Reduced_Objective_SimOpt<RealT> > robj
//      = ROL::makePtr<ROL::Reduced_Objective_SimOpt<RealT>>(obj, con, up, zp, pp, storage, false);

    // Run derivative checks
    bool checkDeriv = parlist->sublist("Problem").get("Check derivatives",false);
    if ( checkDeriv ) {
      *outStream << "\n\nCheck Opt Vector\n";
      zp->checkVector(*dzp,*rzp,true,*outStream);

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

//      std::vector<ROL::Ptr<ROL::Objective_SimOpt<RealT> > > obj_vec(2,ROL::nullPtr);
//      obj_vec[0] = ROL::makePtr<IntegralObjective<RealT>>(qoi_vec[0],assembler);
//      obj_vec[1] = ROL::makePtr<IntegralObjective<RealT>>(qoi_vec[1],assembler);
//
//      *outStream << "\n\nCheck Gradient of State Objective Function\n";
//      obj_vec[0]->checkGradient(x,d,true,*outStream);
//      *outStream << "\n\nCheck Gradient_1 of State Objective Function\n";
//      obj_vec[0]->checkGradient_1(*up,*zp,*dup,true,*outStream);
//      *outStream << "\n\nCheck Gradient_2 of State Objective Function\n";
//      obj_vec[0]->checkGradient_2(*up,*zp,*dzp,true,*outStream);
//      *outStream << "\n\nCheck Hessian of State Objective Function\n";
//      obj_vec[0]->checkHessVec(x,d,true,*outStream);
//      *outStream << "\n\nCheck Hessian_11 of State Objective Function\n";
//      obj_vec[0]->checkHessVec_11(*up,*zp,*dup,true,*outStream);
//      *outStream << "\n\nCheck Hessian_12 of State Objective Function\n";
//      obj_vec[0]->checkHessVec_12(*up,*zp,*dzp,true,*outStream);
//      *outStream << "\n\nCheck Hessian_21 of State Objective Function\n";
//      obj_vec[0]->checkHessVec_21(*up,*zp,*dup,true,*outStream);
//      *outStream << "\n\nCheck Hessian_22 of State Objective Function\n";
//      obj_vec[0]->checkHessVec_22(*up,*zp,*dzp,true,*outStream);
//
//      *outStream << "\n\nCheck Gradient of Control Objective Function\n";
//      obj_vec[1]->checkGradient(x,d,true,*outStream);
//      *outStream << "\n\nCheck Gradient_1 of Control Objective Function\n";
//      obj_vec[1]->checkGradient_1(*up,*zp,*dup,true,*outStream);
//      *outStream << "\n\nCheck Gradient_2 of Control Objective Function\n";
//      obj_vec[1]->checkGradient_2(*up,*zp,*dzp,true,*outStream);
//      *outStream << "\n\nCheck Hessian of Control Objective Function\n";
//      obj_vec[1]->checkHessVec(x,d,true,*outStream);
//      *outStream << "\n\nCheck Hessian_11 of Control Objective Function\n";
//      obj_vec[1]->checkHessVec_11(*up,*zp,*dup,true,*outStream);
//      *outStream << "\n\nCheck Hessian_12 of State Objective Function\n";
//      obj_vec[1]->checkHessVec_12(*up,*zp,*dzp,true,*outStream);
//      *outStream << "\n\nCheck Hessian_21 of State Objective Function\n";
//      obj_vec[1]->checkHessVec_21(*up,*zp,*dup,true,*outStream);
//      *outStream << "\n\nCheck Hessian_22 of Control Objective Function\n";
//      obj_vec[1]->checkHessVec_22(*up,*zp,*dzp,true,*outStream);
//
//      *outStream << "\n\nCheck Gradient of Full Objective Function\n";
//      obj->checkGradient(x,d,true,*outStream);
//      *outStream << "\n\nCheck Gradient_1 of Full Objective Function\n";
//      obj->checkGradient_1(*up,*zp,*dup,true,*outStream);
//      *outStream << "\n\nCheck Gradient_2 of Full Objective Function\n";
//      obj->checkGradient_2(*up,*zp,*dzp,true,*outStream);
//      *outStream << "\n\nCheck Hessian of Full Objective Function\n";
//      obj->checkHessVec(x,d,true,*outStream);
//      *outStream << "\n\nCheck Hessian_11 of Full Objective Function\n";
//      obj->checkHessVec_11(*up,*zp,*dup,true,*outStream);
//      *outStream << "\n\nCheck Hessian_12 of State Objective Function\n";
//      obj->checkHessVec_12(*up,*zp,*dzp,true,*outStream);
//      *outStream << "\n\nCheck Hessian_21 of State Objective Function\n";
//      obj->checkHessVec_21(*up,*zp,*dup,true,*outStream);
//      *outStream << "\n\nCheck Hessian_22 of Full Objective Function\n";
//      obj->checkHessVec_22(*up,*zp,*dzp,true,*outStream);
//
//      *outStream << "\n\nCheck Gradient of Reduced Objective Function\n";
//      robj->checkGradient(*zp,*dzp,true,*outStream);
//      *outStream << "\n\nCheck Hessian of Reduced Objective Function\n";
//      robj->checkHessVec(*zp,*dzp,true,*outStream);
    }

    // Output uncontrolled state.
    zp->zero();
    z_ptr->putScalar(static_cast<RealT>(1));
    pdecon->printMeshData(*outStream);
    con->solve(*rp,*up,*zp,tol);
    pdecon->outputTpetraVector(u_ptr,"state_uncontrolled.txt");

//    ROL::Algorithm<RealT> algo("Trust Region",*parlist,false);
//    Teuchos::Time algoTimer("Algorithm Time", true);
//    algo.run(*zp,*robj,true,*outStream);
//    algoTimer.stop();
//    *outStream << "Total optimization time = " << algoTimer.totalElapsedTime() << " seconds.\n";

    // Output.
    pdecon->printMeshData(*outStream);
    con->solve(*rp,*up,*zp,tol);
    pdecon->outputTpetraVector(u_ptr,"state.txt");
    pdecon->outputTpetraVector(z_ptr,"control.txt");
    assembler->serialPrintStateEdgeField(u_ptr,pde->getFieldHelper(),"stateCellCenter.txt",pde->getFE());

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
