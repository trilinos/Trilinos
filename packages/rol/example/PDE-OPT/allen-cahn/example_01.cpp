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
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "Tpetra_Core.hpp"
#include "Tpetra_Version.hpp"

#include <iostream>
#include <algorithm>

#include "ROL_Stream.hpp"
#include "ROL_Solver.hpp"
#include "ROL_Reduced_Objective_SimOpt.hpp"
#include "ROL_BoundConstraint_SimOpt.hpp"
#include "ROL_Bounds.hpp"

#include "../TOOLS/meshmanager.hpp"
#include "../TOOLS/pdeconstraint.hpp"
#include "../TOOLS/pdeobjective.hpp"
#include "../TOOLS/pdevector.hpp"
#include "pde_allen_cahn.hpp"
#include "obj_allen_cahn.hpp"

typedef double RealT;

int main(int argc, char *argv[]) {
  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  ROL::Ptr<std::ostream> outStream;
  ROL::nullstream bhs; // outputs nothing

  /*** Initialize communicator. ***/
  Teuchos::GlobalMPISession mpiSession (&argc, &argv, &bhs);
  ROL::Ptr<const Teuchos::Comm<int> > comm
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

    /*** Read in XML input ***/
    std::string filename = "input.xml";
    Teuchos::RCP<Teuchos::ParameterList> parlist = Teuchos::rcp( new Teuchos::ParameterList() );
    Teuchos::updateParametersFromXmlFile( filename, parlist.ptr() );

    /*** Initialize main data structure. ***/
    ROL::Ptr<MeshManager<RealT>>
      meshMgr = ROL::makePtr<MeshManager_Rectangle<RealT>>(*parlist);
    // Initialize PDE describe Poisson's equation
    ROL::Ptr<PDE_Allen_Cahn<RealT>>
      pde = ROL::makePtr<PDE_Allen_Cahn<RealT>>(*parlist);
    ROL::Ptr<PDE_Constraint<RealT>>
      con = ROL::makePtr<PDE_Constraint<RealT>>(pde,meshMgr,comm,*parlist,*outStream);
    const ROL::Ptr<Assembler<RealT>> assembler = con->getAssembler();
    assembler->printMeshData(*outStream);
    con->setSolveParameters(*parlist);

    /*************************************************************************/
    /***************** BUILD VECTORS *****************************************/
    /*************************************************************************/
    ROL::Ptr<Tpetra::MultiVector<>>  u_ptr = assembler->createStateVector();
    ROL::Ptr<Tpetra::MultiVector<>>  p_ptr = assembler->createStateVector();
    ROL::Ptr<Tpetra::MultiVector<>>  r_ptr = assembler->createResidualVector();
    ROL::Ptr<Tpetra::MultiVector<>>  z_ptr = assembler->createControlVector();
    ROL::Ptr<ROL::Vector<RealT>> up, pp, rp, zp;
    u_ptr->randomize();  //u_ptr->putScalar(static_cast<RealT>(1));
    p_ptr->randomize();  //p_ptr->putScalar(static_cast<RealT>(1));
    r_ptr->randomize();  //r_ptr->putScalar(static_cast<RealT>(1));
    z_ptr->randomize();  //z_ptr->putScalar(static_cast<RealT>(1));
    up = ROL::makePtr<PDE_PrimalSimVector<RealT>>(u_ptr,pde,assembler,*parlist);
    pp = ROL::makePtr<PDE_PrimalSimVector<RealT>>(p_ptr,pde,assembler,*parlist);
    rp = ROL::makePtr<PDE_DualSimVector<RealT>>(r_ptr,pde,assembler,*parlist);
    zp = ROL::makePtr<PDE_PrimalOptVector<RealT>>(z_ptr,pde,assembler,*parlist);

    /*************************************************************************/
    /***************** BUILD COST FUNCTIONAL *********************************/
    /*************************************************************************/
    // Initialize quadratic objective function
    std::vector<ROL::Ptr<QoI<RealT>>> qoi_vec(2,ROL::nullPtr);
    qoi_vec[0] = ROL::makePtr<QoI_State_Cost_Allen_Cahn<RealT>>(pde->getFE());
    qoi_vec[1] = ROL::makePtr<QoI_Control_Cost_Allen_Cahn<RealT>>(pde->getFE(),
                                                                  pde->getBdryFE(),
                                                                  pde->getBdryCellLocIds());
    std::vector<RealT> weights(2);
    weights[0] = static_cast<RealT>(1);
    weights[1] = parlist->sublist("Problem").get("Control Penalty Parameter", 1e-4);
    // Build full-space objective
    ROL::Ptr<PDE_Objective<RealT>>
      obj = ROL::makePtr<PDE_Objective<RealT>>(qoi_vec,weights,assembler);
    // Build reduced-space objective
    bool storage = parlist->sublist("Problem").get("Use state storage",true);
    ROL::Ptr<ROL::VectorController<RealT>>
      stateStore = ROL::makePtr<ROL::VectorController<RealT>>();
    ROL::Ptr<ROL::Reduced_Objective_SimOpt<RealT>>
      robj = ROL::makePtr<ROL::Reduced_Objective_SimOpt<RealT>>(obj,con,
                                                                stateStore,
                                                                up,zp,pp,
                                                                storage);

    /*************************************************************************/
    /***************** BUILD BOUND CONSTRAINT ********************************/
    /*************************************************************************/
    // Control bounds
    ROL::Ptr<Tpetra::MultiVector<>> zlo_ptr = assembler->createControlVector();
    ROL::Ptr<Tpetra::MultiVector<>> zhi_ptr = assembler->createControlVector();
    RealT lo = parlist->sublist("Problem").get("Lower Bound",0.0);
    RealT hi = parlist->sublist("Problem").get("Upper Bound",1.0);
    zlo_ptr->putScalar(lo); zhi_ptr->putScalar(hi);
    ROL::Ptr<ROL::Vector<RealT>> zlop, zhip;
    zlop = ROL::makePtr<PDE_PrimalOptVector<RealT>>(zlo_ptr,pde,assembler);
    zhip = ROL::makePtr<PDE_PrimalOptVector<RealT>>(zhi_ptr,pde,assembler);
    ROL::Ptr<ROL::BoundConstraint<RealT>>
      zbnd = ROL::makePtr<ROL::Bounds<RealT>>(zlop,zhip);
    bool deactivate = parlist->sublist("Problem").get("Deactivate Bound Constraints",false);
    if (deactivate) zbnd->deactivate();
    // State bounds
    ROL::Ptr<Tpetra::MultiVector<>> ulo_ptr = assembler->createStateVector();
    ROL::Ptr<Tpetra::MultiVector<>> uhi_ptr = assembler->createStateVector();
    ulo_ptr->putScalar(ROL::ROL_NINF<RealT>()); uhi_ptr->putScalar(ROL::ROL_INF<RealT>());
    ROL::Ptr<ROL::Vector<RealT>> ulop, uhip;
    ulop = ROL::makePtr<PDE_PrimalSimVector<RealT>>(ulo_ptr,pde,assembler);
    uhip = ROL::makePtr<PDE_PrimalSimVector<RealT>>(uhi_ptr,pde,assembler);
    ROL::Ptr<ROL::BoundConstraint<RealT>>
    ubnd = ROL::makePtr<ROL::Bounds<RealT>>(ulop,uhip);
    ubnd->deactivate();

    // SimOpt bounds
    ROL::Ptr<ROL::BoundConstraint<RealT>>
      bnd = ROL::makePtr<ROL::BoundConstraint_SimOpt<RealT>>(ubnd,zbnd);

    // Create ROL SimOpt vectors
    ROL::Ptr<Tpetra::MultiVector<>> du_ptr = assembler->createStateVector();
    ROL::Ptr<Tpetra::MultiVector<>> dz_ptr = assembler->createControlVector();
    du_ptr->randomize(); //du_ptr->putScalar(static_cast<RealT>(0));
    dz_ptr->randomize(); //dz_ptr->putScalar(static_cast<RealT>(1));
    ROL::Ptr<ROL::Vector<RealT>> dup, dzp;
    dup = ROL::makePtr<PDE_PrimalSimVector<RealT>>(du_ptr,pde,assembler,*parlist);
    dzp = ROL::makePtr<PDE_PrimalOptVector<RealT>>(dz_ptr,pde,assembler,*parlist);

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
