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
    \brief Shows how to solve the Poisson control problem.
*/

#include "Teuchos_Comm.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "Tpetra_DefaultPlatform.hpp"
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
void checkDerivatives(const ROL::SharedPointer<ROL::Constraint_SimOpt<Real> > & con,
                      std::ostream & outStream = std::cout) {
  // Cast the constraint and get the assembler.
  ROL::SharedPointer<Linear_PDE_Constraint<Real> > pdecon
    = ROL::dynamicPointerCast<Linear_PDE_Constraint<Real> >(con);
  ROL::SharedPointer<Assembler<Real> > assembler = pdecon->getAssembler();
  ROL::SharedPointer<PDE<Real> > pde = pdecon->getPDE();
  // Create state vector and set to zeroes
  ROL::SharedPointer<Tpetra::MultiVector<> > u_ptr = assembler->createStateVector();
  u_ptr->randomize();
  ROL::SharedPointer<ROL::Vector<Real> > up
    = ROL::makeShared<PDE_PrimalSimVector<Real>>(u_ptr,pde,assembler);
  // Create control vector and set to ones
  ROL::SharedPointer<Tpetra::MultiVector<> > z_ptr = assembler->createControlVector();
  z_ptr->randomize();
  ROL::SharedPointer<ROL::Vector<Real> > zp
    = ROL::makeShared<PDE_PrimalOptVector<Real>>(z_ptr,pde,assembler);
  // Create residual vector and set to zeros
  ROL::SharedPointer<Tpetra::MultiVector<> > r_ptr = assembler->createResidualVector();
  r_ptr->putScalar(0.0);
  ROL::SharedPointer<ROL::Vector<Real> > rp
    = ROL::makeShared<PDE_DualSimVector<Real>>(r_ptr,pde,assembler);
  // Create multiplier vector and set to zeros
  ROL::SharedPointer<Tpetra::MultiVector<> > p_ptr = assembler->createStateVector();
  p_ptr->putScalar(0.0);
  ROL::SharedPointer<ROL::Vector<Real> > pp
    = ROL::makeShared<PDE_PrimalSimVector<Real>>(p_ptr,pde,assembler);
  // Create state direction vector and set to random
  ROL::SharedPointer<Tpetra::MultiVector<> > du_ptr = assembler->createStateVector();
  du_ptr->randomize();
  ROL::SharedPointer<ROL::Vector<Real> > dup
    = ROL::makeShared<PDE_PrimalSimVector<Real>>(du_ptr,pde,assembler);
  // Create control direction vector and set to random
  ROL::SharedPointer<Tpetra::MultiVector<> > dz_ptr = assembler->createControlVector();
  dz_ptr->randomize();
  ROL::SharedPointer<ROL::Vector<Real> > dzp
    = ROL::makeShared<PDE_PrimalOptVector<Real>>(dz_ptr,pde,assembler);
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
  ROL::SharedPointer<std::ostream> outStream;
  Teuchos::oblackholestream bhs; // outputs nothing

  /*** Initialize communicator. ***/
  Teuchos::GlobalMPISession mpiSession (&argc, &argv, &bhs);
  ROL::SharedPointer<const Teuchos::Comm<int> > comm
    = Tpetra::DefaultPlatform::getDefaultPlatform().getComm();
  const int myRank = comm->getRank();
  if ((iprint > 0) && (myRank == 0)) {
    outStream = ROL::makeSharedFromRef(std::cout);
  }
  else {
    outStream = ROL::makeSharedFromRef(bhs);
  }
  int errorFlag  = 0;

  // *** Example body.
  try {

    /*** Read in XML input ***/
    std::string filename = "input.xml";
    ROL::SharedPointer<Teuchos::ParameterList> parlist = ROL::makeShared<Teuchos::ParameterList>();
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
    ROL::SharedPointer<MeshManager<RealT> > meshMgr_local
      = ROL::makeShared<MeshManager_Rectangle<RealT>>(*parlist);
    ROL::SharedPointer<PDE_Fractional_Poisson_Local<RealT> > pde_local
      = ROL::makeShared<PDE_Fractional_Poisson_Local<RealT>>(*parlist);
    ROL::SharedPointer<ROL::Constraint_SimOpt<RealT> > con_local
      = ROL::makeShared<Linear_PDE_Constraint<RealT>>(pde_local,meshMgr_local,comm,*parlist,*outStream);
    // Initialize PDE describe Poisson's equation
    ROL::SharedPointer<MeshManager<RealT> > meshMgr_cylinder
      = ROL::makeShared<MeshManager_Fractional_Cylinder<RealT>>(*parlist);
    ROL::SharedPointer<PDE_Fractional_Poisson_Cylinder<RealT> > pde_cylinder
      = ROL::makeShared<PDE_Fractional_Poisson_Cylinder<RealT>>(*parlist);
    ROL::SharedPointer<ROL::Constraint_SimOpt<RealT> > con_cylinder
      = ROL::makeShared<Linear_PDE_Constraint<RealT>>(pde_cylinder,meshMgr_cylinder,comm,*parlist,*outStream);

    // Check derivatives.
    checkDerivatives<RealT>(con_cylinder, *outStream);
    checkDerivatives<RealT>(con_local, *outStream);

    // Build fractional problem.
    ROL::SharedPointer<ROL::LinearOperator<RealT> > A
      = ROL::makeShared<FractionalOperator<RealT>>(pde_local,meshMgr_local,pde_cylinder,meshMgr_cylinder,comm,*parlist,*outStream);
    ROL::SharedPointer<ROL::LinearOperator<RealT> > M
      = ROL::makeShared<FractionalPreconditioner<RealT>>(pde_local,meshMgr_local,pde_cylinder,meshMgr_cylinder,comm,*parlist,*outStream);
    FractionalVector<RealT> F(pde_local,meshMgr_local,pde_cylinder,meshMgr_cylinder,comm,*parlist,*outStream);
    const ROL::SharedPointer<const ROL::Vector<RealT> > Fvec = F.get();

    // Test preconditioner
    RealT tol(1e-8);
    ROL::SharedPointer<ROL::Vector<RealT> > PFvec = Fvec->clone();
    ROL::SharedPointer<ROL::Vector<RealT> > PPFvec = Fvec->clone();
    M->applyInverse(*PFvec,*Fvec,tol);
    M->apply(*PPFvec,*PFvec,tol);
    PPFvec->axpy(static_cast<RealT>(-1),*Fvec);
    *outStream << "Preconditioner error: " << PPFvec->norm() << std::endl;

    // Build CG object.
    //ROL::SharedPointer<ROL::Krylov<RealT> > solver = ROL::makeShared<ROL::BelosKrylov<RealT>>(*parlist);
    ROL::SharedPointer<ROL::Krylov<RealT> > solver = ROL::KrylovFactory<RealT>(*parlist);

    // Run CG
    int iter(0), flag(0);
    ROL::SharedPointer<ROL::Vector<RealT> > Xvec = Fvec->clone();
    solver->run(*Xvec,*A,*Fvec,*M,iter,flag);
    *outStream << "GMRES Iterations: " << iter << "  GMRES Flag: " << flag << std::endl;

    ROL::SharedPointer<ROL::Vector<RealT> > Rvec = Fvec->clone();
    A->apply(*Rvec,*Xvec,tol);
    Rvec->axpy(static_cast<RealT>(-1),*Fvec);
    *outStream << "GMRES Residual: " << Rvec->norm() << std::endl;

    // Print solution
    ROL::SharedPointer<Linear_PDE_Constraint<RealT> > pdecon
      = ROL::dynamicPointerCast<Linear_PDE_Constraint<RealT> >(con_local);
    pdecon->printMeshData(*outStream);
    ROL::SharedPointer<const Tpetra::MultiVector<> > Xptr
      = ROL::dynamicPointerCast<const ROL::TpetraMultiVector<RealT> >(Xvec)->getVector();
    Teuchos::Array<size_t> col(1,0);
    pdecon->outputTpetraVector(Xptr->subView(col()),"state.txt");

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
