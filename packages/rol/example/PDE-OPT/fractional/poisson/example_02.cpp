// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  example_01.cpp
    \brief Shows how to solve the Poisson control problem.
*/

#include "Teuchos_Comm.hpp"
#include "ROL_Stream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "Tpetra_Core.hpp"
#include "Tpetra_Version.hpp"

#include <iostream>
#include <algorithm>

#include "ROL_Algorithm.hpp"
#include "ROL_Reduced_Objective_SimOpt.hpp"

#include "../../TOOLS/meshmanager.hpp"
#include "../../TOOLS/linearpdeconstraint.hpp"
#include "../../TOOLS/pdeobjective.hpp"
#include "../../TOOLS/pdevector.hpp"
#include "pde_fractional_poisson.hpp"
#include "obj_fractional_poisson.hpp"
#include "fractional_constraint.hpp"
#include "fractional_objective.hpp"

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

    // Set parameters for cylidner mesh
    RealT s     = parlist->sublist("Problem").get("Fractional Power",0.5);
    RealT gamma = (s == 0.5) ? static_cast<RealT>(1)
                  : static_cast<RealT>(3)/(static_cast<RealT>(2)*s) + static_cast<RealT>(1e-3);
    int NX      = parlist->sublist("Geometry").get("NX",20);
    int NY      = parlist->sublist("Geometry").get("NY",20);
    RealT NT    = static_cast<RealT>(NX*NY);
    RealT width = static_cast<RealT>(1) + std::log10(NT)/static_cast<RealT>(3);
    int NI      = static_cast<int>(width * std::sqrt(NT)); 
    parlist->sublist("Geometry").sublist("Cylinder").set("Grading Parameter", gamma);
    parlist->sublist("Geometry").sublist("Cylinder").set("NI",                NI);
    parlist->sublist("Geometry").sublist("Cylinder").set("Height",            width);

    *outStream << std::endl;
    *outStream << "Fractional Power:       " << s     << std::endl;
    *outStream << "Mesh Grading Parameter: " << gamma << std::endl;
    *outStream << "Cylinder Height:        " << width << std::endl;
    *outStream << "Number of Intervals:    " << NI    << std::endl;
    *outStream << std::endl;

    // Initialize 2D Poisson's equation
    ROL::Ptr<MeshManager<RealT> > meshMgr_local
      = ROL::makePtr<MeshManager_Rectangle<RealT>>(*parlist);
    ROL::Ptr<PDE_Fractional_Poisson_Local<RealT> > pde_local
      = ROL::makePtr<PDE_Fractional_Poisson_Local<RealT>>(*parlist);
    // Initialize 1D singular Poisson's equation
    ROL::Ptr<MeshManager<RealT> > meshMgr_cylinder
      = ROL::makePtr<MeshManager_Fractional_Cylinder<RealT>>(*parlist);
    ROL::Ptr<PDE_Fractional_Poisson_Cylinder<RealT> > pde_cylinder
      = ROL::makePtr<PDE_Fractional_Poisson_Cylinder<RealT>>(*parlist);
    // Build fractional constraint
    ROL::Ptr<ROL::Constraint_SimOpt<RealT> > con
      = ROL::makePtr<FractionalConstraint<RealT>>(pde_local,    meshMgr_local,    comm,
                                                     pde_cylinder, meshMgr_cylinder, comm,
                                                     *parlist, *outStream);
    ROL::Ptr<FractionalConstraint<RealT> > fracCon
      = ROL::dynamicPtrCast<FractionalConstraint<RealT> >(con);
    ROL::Ptr<Assembler<RealT> > assembler = fracCon->getLocalAssembler();

    // Build objective fuction
    std::vector<ROL::Ptr<QoI<RealT> > > qoi_vec(2,ROL::nullPtr);
    qoi_vec[0] = ROL::makePtr<QoI_L2Tracking_Fractional_Poisson<RealT>>(pde_local->getFE());
    qoi_vec[1] = ROL::makePtr<QoI_L2Penalty_Fractional_Poisson<RealT>>(pde_local->getFE());
    RealT stateCost   = parlist->sublist("Problem").get("State Cost",1e0);
    RealT controlCost = parlist->sublist("Problem").get("Control Cost",1e0);
    std::vector<RealT> wts = {stateCost, controlCost};
    ROL::Ptr<ROL::Objective_SimOpt<RealT> > pdeobj
      = ROL::makePtr<PDE_Objective<RealT>>(qoi_vec,wts,assembler);
    ROL::Ptr<ROL::Objective_SimOpt<RealT> > obj
      = ROL::makePtr<FractionalObjective<RealT>>(pdeobj);

    // Build vectors
    ROL::Ptr<Tpetra::MultiVector<> > stateVec = assembler->createStateVector();
    ROL::Ptr<Tpetra::MultiVector<> > u_ptr    = ROL::makePtr<Tpetra::MultiVector<>>(stateVec->getMap(),NI+1);
    ROL::Ptr<Tpetra::MultiVector<> > p_ptr    = ROL::makePtr<Tpetra::MultiVector<>>(stateVec->getMap(),NI+1);
    ROL::Ptr<Tpetra::MultiVector<> > du_ptr   = ROL::makePtr<Tpetra::MultiVector<>>(stateVec->getMap(),NI+1);
    u_ptr->randomize();  //u_ptr->putScalar(static_cast<RealT>(1));
    p_ptr->randomize();  //p_ptr->putScalar(static_cast<RealT>(1));
    du_ptr->randomize(); //du_ptr->putScalar(static_cast<RealT>(0));
    ROL::Ptr<ROL::Vector<RealT> > up  = ROL::makePtr<ROL::TpetraMultiVector<RealT>>(u_ptr);
    ROL::Ptr<ROL::Vector<RealT> > pp  = ROL::makePtr<ROL::TpetraMultiVector<RealT>>(p_ptr);
    ROL::Ptr<ROL::Vector<RealT> > dup = ROL::makePtr<ROL::TpetraMultiVector<RealT>>(du_ptr);
    // Create residual vectors
    ROL::Ptr<Tpetra::MultiVector<> > conVec = assembler->createResidualVector();
    ROL::Ptr<Tpetra::MultiVector<> > r_ptr  = ROL::makePtr<Tpetra::MultiVector<>>(conVec->getMap(),NI+1);
    r_ptr->randomize(); //r_ptr->putScalar(static_cast<RealT>(1));
    ROL::Ptr<ROL::Vector<RealT> > rp = ROL::makePtr<ROL::TpetraMultiVector<RealT>>(r_ptr);
    // Create control vector and set to ones
    ROL::Ptr<Tpetra::MultiVector<> > z_ptr  = assembler->createControlVector();
    ROL::Ptr<Tpetra::MultiVector<> > dz_ptr = assembler->createControlVector();
    z_ptr->randomize();  //z_ptr->putScalar(static_cast<RealT>(1));
    dz_ptr->randomize(); //dz_ptr->putScalar(static_cast<RealT>(0));
    ROL::Ptr<ROL::Vector<RealT> > zp  = ROL::makePtr<ROL::TpetraMultiVector<RealT>>(z_ptr);
    ROL::Ptr<ROL::Vector<RealT> > dzp = ROL::makePtr<ROL::TpetraMultiVector<RealT>>(dz_ptr);
    // Create ROL SimOpt vectors
    ROL::Vector_SimOpt<RealT> x(up,zp);
    ROL::Vector_SimOpt<RealT> d(dup,dzp);

    // Build reduced objective function
    bool storage = parlist->sublist("Problem").get("Use State and Adjoint Storage",true);
    ROL::Ptr<ROL::Reduced_Objective_SimOpt<RealT> > objReduced
      = ROL::makePtr<ROL::Reduced_Objective_SimOpt<RealT>>(obj, con, up, zp, pp, storage, false);

    // Check derivatives.
    bool checkDeriv = parlist->sublist("Problem").get("Check Derivatives",false);
    if ( checkDeriv ) {
      *outStream << "\n\nCheck Gradient_1 of Full Objective Function\n";
      obj->checkGradient_1(*up,*zp,*dup,true,*outStream);
      *outStream << "\n\nCheck Gradient_2 of Full Objective Function\n";
      obj->checkGradient_2(*up,*zp,*dzp,true,*outStream);
      *outStream << "\n\nCheck Gradient of Full Objective Function\n";
      obj->checkGradient(x,d,true,*outStream);
      *outStream << "\n\nCheck Hessian_11 of Full Objective Function\n";
      obj->checkHessVec_11(*up,*zp,*dup,true,*outStream);
      *outStream << "\n\nCheck Hessian_12 of Full Objective Function\n";
      obj->checkHessVec_12(*up,*zp,*dzp,true,*outStream);
      *outStream << "\n\nCheck Hessian_21 of Full Objective Function\n";
      obj->checkHessVec_21(*up,*zp,*dup,true,*outStream);
      *outStream << "\n\nCheck Hessian_22 of Full Objective Function\n";
      obj->checkHessVec_22(*up,*zp,*dzp,true,*outStream);
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
      *outStream << "\n\nCheck Hessian_11 of PDE Constraint\n";
      con->checkApplyAdjointHessian_11(*up,*zp,*pp,*dup,*dup,true,*outStream);
      *outStream << "\n\nCheck Hessian_21 of PDE Constraint\n";
      con->checkApplyAdjointHessian_21(*up,*zp,*pp,*dzp,*dup,true,*outStream);
      *outStream << "\n\nCheck Hessian_12 of PDE Constraint\n";
      con->checkApplyAdjointHessian_12(*up,*zp,*pp,*dup,*dzp,true,*outStream);
      *outStream << "\n\nCheck Hessian_22 of PDE Constraint\n";
      con->checkApplyAdjointHessian_22(*up,*zp,*pp,*dzp,*dzp,true,*outStream);
      *outStream << "\n";
      con->checkAdjointConsistencyJacobian(*dup,d,x,true,*outStream);
      *outStream << "\n";
      con->checkInverseJacobian_1(*up,*up,*up,*zp,true,*outStream);
      *outStream << "\n";
      *outStream << "\n\nCheck Gradient of Reduced Objective Function\n";
      objReduced->checkGradient(*zp,*dzp,true,*outStream);
      *outStream << "\n\nCheck Hessian of Reduced Objective Function\n";
      objReduced->checkHessVec(*zp,*dzp,true,*outStream);
    }

    // Run optimization
    ROL::Ptr<ROL::Algorithm<RealT> > algo
      = ROL::makePtr<ROL::Algorithm<RealT>>("Trust Region",*parlist,false);
    zp->zero();
    algo->run(*zp,*objReduced,true,*outStream);

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
