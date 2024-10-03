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

//#include "ROL_Algorithm.hpp"
//#include "ROL_Reduced_Objective_SimOpt.hpp"
#include "ROL_KrylovFactory.hpp"
//#include "ROL_BelosKrylov.hpp"

#include "../../TOOLS/meshmanager.hpp"
#include "../../TOOLS/linearpdeconstraint.hpp"
//#include "../../TOOLS/pdeobjective.hpp"
#include "../../TOOLS/pdevector.hpp"
#include "pde_fractional_poisson.hpp"
//#include "obj_fractional_poisson.hpp"
#include "fractional_operator.hpp"
#include "fractional_vector.hpp"

typedef double RealT;

template<class Real>
void checkDerivatives(const ROL::Ptr<ROL::Constraint_SimOpt<Real> > & con,
                      std::ostream & outStream = std::cout) {
  // Cast the constraint and get the assembler.
  ROL::Ptr<Linear_PDE_Constraint<Real> > pdecon
    = ROL::dynamicPtrCast<Linear_PDE_Constraint<Real> >(con);
  ROL::Ptr<Assembler<Real> > assembler = pdecon->getAssembler();
  ROL::Ptr<PDE<Real> > pde = pdecon->getPDE();
  // Create state vector and set to zeroes
  ROL::Ptr<Tpetra::MultiVector<> > u_ptr = assembler->createStateVector();
  u_ptr->randomize();
  ROL::Ptr<ROL::Vector<Real> > up
    = ROL::makePtr<PDE_PrimalSimVector<Real>>(u_ptr,pde,assembler);
  // Create control vector and set to ones
  ROL::Ptr<Tpetra::MultiVector<> > z_ptr = assembler->createControlVector();
  z_ptr->randomize();
  ROL::Ptr<ROL::Vector<Real> > zp
    = ROL::makePtr<PDE_PrimalOptVector<Real>>(z_ptr,pde,assembler);
  // Create residual vector and set to zeros
  ROL::Ptr<Tpetra::MultiVector<> > r_ptr = assembler->createResidualVector();
  r_ptr->putScalar(0.0);
  ROL::Ptr<ROL::Vector<Real> > rp
    = ROL::makePtr<PDE_DualSimVector<Real>>(r_ptr,pde,assembler);
  // Create multiplier vector and set to zeros
  ROL::Ptr<Tpetra::MultiVector<> > p_ptr = assembler->createStateVector();
  p_ptr->putScalar(0.0);
  ROL::Ptr<ROL::Vector<Real> > pp
    = ROL::makePtr<PDE_PrimalSimVector<Real>>(p_ptr,pde,assembler);
  // Create state direction vector and set to random
  ROL::Ptr<Tpetra::MultiVector<> > du_ptr = assembler->createStateVector();
  du_ptr->randomize();
  ROL::Ptr<ROL::Vector<Real> > dup
    = ROL::makePtr<PDE_PrimalSimVector<Real>>(du_ptr,pde,assembler);
  // Create control direction vector and set to random
  ROL::Ptr<Tpetra::MultiVector<> > dz_ptr = assembler->createControlVector();
  dz_ptr->randomize();
  ROL::Ptr<ROL::Vector<Real> > dzp
    = ROL::makePtr<PDE_PrimalOptVector<Real>>(dz_ptr,pde,assembler);
  // Create ROL SimOpt vectors
  ROL::Vector_SimOpt<Real> x(up,zp);
  ROL::Vector_SimOpt<Real> d(dup,dzp);
  // Check derivatives.
  con->checkApplyJacobian(x,d,*up,true,outStream);
  con->checkApplyAdjointHessian(x,*dup,d,x,true,outStream);
  con->checkAdjointConsistencyJacobian(*dup,d,x,true,outStream);
  con->checkInverseJacobian_1(*up,*up,*up,*zp,true,outStream);
  con->checkInverseAdjointJacobian_1(*up,*up,*up,*zp,true,outStream);
}

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

    /*** Initialize main data structure. ***/
    // Initialize PDE describe Poisson's equation
    ROL::Ptr<MeshManager<RealT> > meshMgr_local
      = ROL::makePtr<MeshManager_Rectangle<RealT>>(*parlist);
    ROL::Ptr<PDE_Fractional_Poisson_Local<RealT> > pde_local
      = ROL::makePtr<PDE_Fractional_Poisson_Local<RealT>>(*parlist);
    ROL::Ptr<ROL::Constraint_SimOpt<RealT> > con_local
      = ROL::makePtr<Linear_PDE_Constraint<RealT>>(pde_local,meshMgr_local,comm,*parlist,*outStream);
    // Initialize PDE describe Poisson's equation
    ROL::Ptr<MeshManager<RealT> > meshMgr_cylinder
      = ROL::makePtr<MeshManager_Fractional_Cylinder<RealT>>(*parlist);
    ROL::Ptr<PDE_Fractional_Poisson_Cylinder<RealT> > pde_cylinder
      = ROL::makePtr<PDE_Fractional_Poisson_Cylinder<RealT>>(*parlist);
    ROL::Ptr<ROL::Constraint_SimOpt<RealT> > con_cylinder
      = ROL::makePtr<Linear_PDE_Constraint<RealT>>(pde_cylinder,meshMgr_cylinder,comm,*parlist,*outStream);

    // Check derivatives.
    checkDerivatives<RealT>(con_cylinder, *outStream);
    checkDerivatives<RealT>(con_local, *outStream);

    // Build fractional problem.
    ROL::Ptr<ROL::LinearOperator<RealT> > A
      = ROL::makePtr<FractionalOperator<RealT>>(pde_local,meshMgr_local,pde_cylinder,meshMgr_cylinder,comm,*parlist,*outStream);
    ROL::Ptr<ROL::LinearOperator<RealT> > M
      = ROL::makePtr<FractionalPreconditioner<RealT>>(pde_local,meshMgr_local,pde_cylinder,meshMgr_cylinder,comm,*parlist,*outStream);
    FractionalVector<RealT> F(pde_local,meshMgr_local,pde_cylinder,meshMgr_cylinder,comm,*parlist,*outStream);
    const ROL::Ptr<const ROL::Vector<RealT> > Fvec = F.get();

    // Test preconditioner
    RealT tol(1e-8);
    ROL::Ptr<ROL::Vector<RealT> > PFvec = Fvec->clone();
    ROL::Ptr<ROL::Vector<RealT> > PPFvec = Fvec->clone();
    M->applyInverse(*PFvec,*Fvec,tol);
    M->apply(*PPFvec,*PFvec,tol);
    PPFvec->axpy(static_cast<RealT>(-1),*Fvec);
    *outStream << "Preconditioner error: " << PPFvec->norm() << std::endl;

    // Build CG object.
    //ROL::Ptr<ROL::Krylov<RealT> > solver = ROL::makePtr<ROL::BelosKrylov<RealT>>(*parlist);
    ROL::Ptr<ROL::Krylov<RealT> > solver = ROL::KrylovFactory<RealT>(*parlist);

    // Run CG
    int iter(0), flag(0);
    ROL::Ptr<ROL::Vector<RealT> > Xvec = Fvec->clone();
    solver->run(*Xvec,*A,*Fvec,*M,iter,flag);
    *outStream << "GMRES Iterations: " << iter << "  GMRES Flag: " << flag << std::endl;

    ROL::Ptr<ROL::Vector<RealT> > Rvec = Fvec->clone();
    A->apply(*Rvec,*Xvec,tol);
    Rvec->axpy(static_cast<RealT>(-1),*Fvec);
    *outStream << "GMRES Residual: " << Rvec->norm() << std::endl;

    // Print solution
    ROL::Ptr<Linear_PDE_Constraint<RealT> > pdecon
      = ROL::dynamicPtrCast<Linear_PDE_Constraint<RealT> >(con_local);
    pdecon->printMeshData(*outStream);
    ROL::Ptr<const Tpetra::MultiVector<> > Xptr
      = ROL::dynamicPtrCast<const ROL::TpetraMultiVector<RealT> >(Xvec)->getVector();
    Teuchos::Array<size_t> col(1,0);
    pdecon->outputTpetraVector(Xptr->subView(col()),"state.txt");

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
