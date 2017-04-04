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
    \brief Shows how to solve the Navier-Stokes control problem.
*/

#include "Teuchos_Comm.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "Tpetra_DefaultPlatform.hpp"
#include "Tpetra_Version.hpp"

#include <iostream>
#include <algorithm>
//#include <fenv.h>

#include "ROL_TpetraMultiVector.hpp"
#include "ROL_Algorithm.hpp"
#include "ROL_Reduced_Objective_SimOpt.hpp"

#include "../TOOLS/meshmanager.hpp"
#include "../TOOLS/pdeconstraint.hpp"
#include "../TOOLS/pdeobjective.hpp"
#include "../TOOLS/pdevector.hpp"
#include "pde_thermal-fluids.hpp"
#include "obj_thermal-fluids.hpp"

typedef double RealT;

int main(int argc, char *argv[]) {
//  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  Teuchos::RCP<std::ostream> outStream;
  Teuchos::oblackholestream bhs; // outputs nothing

  /*** Initialize communicator. ***/
  Teuchos::GlobalMPISession mpiSession (&argc, &argv, &bhs);
  Teuchos::RCP<const Teuchos::Comm<int> > comm
    = Tpetra::DefaultPlatform::getDefaultPlatform().getComm();
  const int myRank = comm->getRank();
  if ((iprint > 0) && (myRank == 0)) {
    outStream = Teuchos::rcp(&std::cout, false);
  }
  else {
    outStream = Teuchos::rcp(&bhs, false);
  }
  int errorFlag  = 0;

  // *** Example body.
  try {

    /*** Read in XML input ***/
    std::string filename = "input.xml";
    Teuchos::RCP<Teuchos::ParameterList> parlist = Teuchos::rcp( new Teuchos::ParameterList() );
    Teuchos::updateParametersFromXmlFile( filename, parlist.ptr() );

    /*** Initialize main data structure. ***/
    Teuchos::RCP<MeshManager<RealT> > meshMgr
      = Teuchos::rcp(new MeshManager_BackwardFacingStepChannel<RealT>(*parlist));
    // Initialize PDE describing Navier-Stokes equations.
    Teuchos::RCP<PDE_ThermalFluids<RealT> > pde
      = Teuchos::rcp(new PDE_ThermalFluids<RealT>(*parlist));
    Teuchos::RCP<ROL::EqualityConstraint_SimOpt<RealT> > con
      = Teuchos::rcp(new PDE_Constraint<RealT>(pde,meshMgr,comm,*parlist,*outStream));
    // Cast the constraint and get the assembler.
    Teuchos::RCP<PDE_Constraint<RealT> > pdecon
      = Teuchos::rcp_dynamic_cast<PDE_Constraint<RealT> >(con);
    Teuchos::RCP<Assembler<RealT> > assembler = pdecon->getAssembler();
    con->setSolveParameters(*parlist);
    pdecon->outputTpetraData();

    // Create state vector and set to zeroes
    Teuchos::RCP<Tpetra::MultiVector<> > u_rcp, p_rcp, du_rcp, r_rcp, z_rcp, dz_rcp;
    u_rcp  = assembler->createStateVector();     u_rcp->randomize();
    p_rcp  = assembler->createStateVector();     p_rcp->randomize();
    du_rcp = assembler->createStateVector();     du_rcp->randomize();
    r_rcp  = assembler->createResidualVector();  r_rcp->randomize();
    z_rcp  = assembler->createControlVector();   z_rcp->randomize();
    dz_rcp = assembler->createControlVector();   dz_rcp->randomize();
    Teuchos::RCP<ROL::Vector<RealT> > up, pp, dup, rp, zp, dzp;
    up  = Teuchos::rcp(new PDE_PrimalSimVector<RealT>(u_rcp,pde,assembler));
    pp  = Teuchos::rcp(new PDE_PrimalSimVector<RealT>(p_rcp,pde,assembler));
    dup = Teuchos::rcp(new PDE_PrimalSimVector<RealT>(du_rcp,pde,assembler));
    rp  = Teuchos::rcp(new PDE_DualSimVector<RealT>(r_rcp,pde,assembler));
    zp  = Teuchos::rcp(new PDE_PrimalOptVector<RealT>(z_rcp,pde,assembler));
    dzp = Teuchos::rcp(new PDE_PrimalOptVector<RealT>(dz_rcp,pde,assembler));
    // Create ROL SimOpt vectors
    ROL::Vector_SimOpt<RealT> x(up,zp);
    ROL::Vector_SimOpt<RealT> d(dup,dzp);

    // Initialize objective function.
    std::vector<Teuchos::RCP<QoI<RealT> > > qoi_vec(2,Teuchos::null);
    qoi_vec[0] = Teuchos::rcp(new QoI_State_ThermalFluids<RealT>(*parlist,
                                                                 pde->getVelocityFE(),
                                                                 pde->getPressureFE(),
                                                                 pde->getThermalFE(),
                                                                 pde->getFieldHelper()));
    qoi_vec[1] = Teuchos::rcp(new QoI_L2Penalty_ThermalFluids<RealT>(pde->getVelocityFE(),
                                                                     pde->getPressureFE(),
                                                                     pde->getThermalFE(),
                                                                     pde->getThermalBdryFE(),
                                                                     pde->getBdryCellLocIds(),
                                                                     pde->getFieldHelper()));
    Teuchos::RCP<StdObjective_ThermalFluids<RealT> > std_obj
      = Teuchos::rcp(new StdObjective_ThermalFluids<RealT>(*parlist));
    Teuchos::RCP<ROL::Objective_SimOpt<RealT> > obj
      = Teuchos::rcp(new PDE_Objective<RealT>(qoi_vec,std_obj,assembler));
    Teuchos::RCP<ROL::SimController<RealT> > stateStore
      = Teuchos::rcp(new ROL::SimController<RealT>());
    Teuchos::RCP<ROL::Reduced_Objective_SimOpt<RealT> > robj
      = Teuchos::rcp(new ROL::Reduced_Objective_SimOpt<RealT>(obj, con, stateStore, up, zp, pp, true, false));

    //up->zero();
    //zp->zero();
    //z_rcp->putScalar(1.e0);
    //dz_rcp->putScalar(0);

    // Run derivative checks
    bool checkDeriv = parlist->sublist("Problem").get("Check derivatives",false);
    if ( checkDeriv ) {
      *outStream << "Check Gradient of Full Objective Function" << std::endl;
      obj->checkGradient(x,d,true,*outStream);
      *outStream << std::endl << "Check Hessian of Full Objective Function" << std::endl;
      obj->checkHessVec(x,d,true,*outStream);
      *outStream << std::endl << "Check Jacobian of Constraint" << std::endl;
      con->checkApplyJacobian(x,d,*up,true,*outStream);
      *outStream << std::endl << "Check Jacobian_1 of Constraint" << std::endl;
      con->checkApplyJacobian_1(*up,*zp,*dup,*rp,true,*outStream);
      *outStream << std::endl << "Check Jacobian_2 of Constraint" << std::endl;
      con->checkApplyJacobian_2(*up,*zp,*dzp,*rp,true,*outStream);
      *outStream << std::endl << "Check Hessian of Constraint" << std::endl;
      con->checkApplyAdjointHessian(x,*dup,d,x,true,*outStream);
      *outStream << std::endl << "Check Hessian_11 of Constraint" << std::endl;
      con->checkApplyAdjointHessian_11(*up,*zp,*pp,*dup,*rp,true,*outStream);
      *outStream << std::endl << "Check Hessian_12 of Constraint" << std::endl;
      con->checkApplyAdjointHessian_12(*up,*zp,*pp,*dup,*dzp,true,*outStream);
      *outStream << std::endl << "Check Hessian_21 of Constraint" << std::endl;
      con->checkApplyAdjointHessian_21(*up,*zp,*pp,*dzp,*rp,true,*outStream);
      *outStream << std::endl << "Check Hessian_22 of Constraint" << std::endl;
      con->checkApplyAdjointHessian_22(*up,*zp,*pp,*dzp,*dzp,true,*outStream);

      *outStream << std::endl << "Check Adjoint Jacobian of Constraint" << std::endl;
      con->checkAdjointConsistencyJacobian(*dup,d,x,true,*outStream);
      *outStream << std::endl << "Check Adjoint Jacobian_1 of Constraint" << std::endl;
      con->checkAdjointConsistencyJacobian_1(*pp,*dup,*up,*zp,true,*outStream);
      *outStream << std::endl << "Check Adjoint Jacobian_2 of Constraint" << std::endl;
      con->checkAdjointConsistencyJacobian_2(*pp,*dzp,*up,*zp,true,*outStream);

      *outStream << std::endl << "Check Constraint Solve" << std::endl;
      con->checkSolve(*up,*zp,*rp,true,*outStream);
      *outStream << std::endl << "Check Inverse Jacobian_1 of Constraint" << std::endl;
      con->checkInverseJacobian_1(*rp,*dup,*up,*zp,true,*outStream);
      *outStream << std::endl << "Check Inverse Adjoint Jacobian_1 of Constraint" << std::endl;
      con->checkInverseAdjointJacobian_1(*rp,*pp,*up,*zp,true,*outStream);

      *outStream << std::endl << "Check Gradient of Reduced Objective Function" << std::endl;
      robj->checkGradient(*zp,*dzp,true,*outStream);
      *outStream << std::endl << "Check Hessian of Reduced Objective Function" << std::endl;
      robj->checkHessVec(*zp,*dzp,true,*outStream);
    }
    up->zero();
    zp->zero();

    RealT tol(1.e-8);
    con->solve(*rp,*up,*zp,tol);
    pdecon->outputTpetraVector(u_rcp,"state_uncontrolled.txt");

    bool useCompositeStep = parlist->sublist("Problem").get("Full space",false);
    Teuchos::RCP<ROL::Algorithm<RealT> > algo;
    if ( useCompositeStep ) {
      algo = Teuchos::rcp(new ROL::Algorithm<RealT>("Composite Step",*parlist,false));
      algo->run(x,*rp,*obj,*con,true,*outStream);
    }
    else {
      algo = Teuchos::rcp(new ROL::Algorithm<RealT>("Trust Region",*parlist,false));
      algo->run(*zp,*robj,true,*outStream);
      std::vector<RealT> param;
      stateStore->get(*up,param);
    }

    // Output.
    assembler->printMeshData(*outStream);
    Teuchos::Array<RealT> res(1,0);
    //con->solve(*rp,*up,*zp,tol);
    pdecon->outputTpetraVector(u_rcp,"state.txt");
    pdecon->outputTpetraVector(z_rcp,"control.txt");
    con->value(*rp,*up,*zp,tol);
    r_rcp->norm2(res.view(0,1));
    *outStream << "Residual Norm: " << res[0] << std::endl;
    errorFlag += (res[0] > 1.e-6 ? 1 : 0);
    //pdecon->outputTpetraData();

    Teuchos::RCP<ROL::Objective_SimOpt<RealT> > obj0
      = Teuchos::rcp(new IntegralObjective<RealT>(qoi_vec[0],assembler));
    RealT val = obj0->value(*up,*zp,tol);
    *outStream << "Vorticity Value: " << val << std::endl;
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
