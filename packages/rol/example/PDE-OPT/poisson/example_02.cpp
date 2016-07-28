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

/*! \file  example_02.cpp
    \brief Shows how to solve the Poisson problem.
*/

#include "Teuchos_Comm.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "Tpetra_DefaultPlatform.hpp"
#include "Tpetra_Version.hpp"

#include <iostream>
#include <algorithm>

#include "../TOOLS/meshmanager.hpp"
#include "../TOOLS/pdefem.hpp"
#include "../TOOLS/pdeconstraint.hpp"
#include "pde_poisson.hpp"

typedef double RealT;

int main(int argc, char *argv[]) {
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
      = Teuchos::rcp(new MeshManager_Rectangle<RealT>(*parlist));
    Teuchos::RCP<PDE<RealT> > pde
      = Teuchos::rcp(new PDE_Poisson<RealT>(*parlist));
    Teuchos::RCP<PDE_FEM<RealT> > fem
      = Teuchos::rcp(new PDE_FEM<RealT>(pde,meshMgr,comm,*parlist,*outStream));
    fem->printMeshData(*outStream);
    Teuchos::RCP<ROL::EqualityConstraint_SimOpt<RealT> > con
      = Teuchos::rcp(new PDE_Constraint<RealT>(fem));

    // Create state vector and set to zeroes
    Teuchos::RCP<Tpetra::MultiVector<> > u_rcp = fem->createStateVector();
    u_rcp->putScalar(0.0);
    Teuchos::RCP<ROL::Vector<RealT> > up
      = Teuchos::rcp(new ROL::TpetraMultiVector<RealT>(u_rcp));
    // Create control vector and set to ones
    Teuchos::RCP<Tpetra::MultiVector<> > z_rcp = fem->createControlVector();
    z_rcp->putScalar(1.0);
    Teuchos::RCP<ROL::Vector<RealT> > zp
      = Teuchos::rcp(new ROL::TpetraMultiVector<RealT>(z_rcp));
    // Create residual vector and set to zeros
    Teuchos::RCP<Tpetra::MultiVector<> > r_rcp = fem->createResidualVector();
    r_rcp->putScalar(0.0);
    Teuchos::RCP<ROL::Vector<RealT> > rp
      = Teuchos::rcp(new ROL::TpetraMultiVector<RealT>(r_rcp));
    // Create state direction vector and set to random
    Teuchos::RCP<Tpetra::MultiVector<> > du_rcp = fem->createStateVector();
    du_rcp->randomize();
    Teuchos::RCP<ROL::Vector<RealT> > dup
      = Teuchos::rcp(new ROL::TpetraMultiVector<RealT>(du_rcp));
    // Create control direction vector and set to random
    Teuchos::RCP<Tpetra::MultiVector<> > dz_rcp = fem->createControlVector();
    dz_rcp->randomize();
    Teuchos::RCP<ROL::Vector<RealT> > dzp
      = Teuchos::rcp(new ROL::TpetraMultiVector<RealT>(dz_rcp));
    // Create ROL SimOpt vectors
    ROL::Vector_SimOpt<RealT> x(up,zp);
    ROL::Vector_SimOpt<RealT> d(dup,dzp);

    // Run derivative checks
    con->checkApplyJacobian(x,d,*up,true,*outStream);
    con->checkApplyAdjointHessian(x,*dup,d,x,true,*outStream);
    con->checkAdjointConsistencyJacobian(*dup,d,x,true,*outStream);
    con->checkInverseJacobian_1(*up,*up,*up,*zp,true,*outStream);
    con->checkInverseAdjointJacobian_1(*up,*up,*up,*zp,true,*outStream);

    fem->assembleResidual(*u_rcp,*z_rcp);
    fem->assembleJacobian1(*u_rcp,*z_rcp);
    r_rcp->scale(-1.0,*(fem->getResidual()));
    fem->solve(u_rcp,r_rcp);
    fem->outputTpetraVector(u_rcp,"solution.txt");

    Teuchos::Array<RealT> res(1,0);
    fem->assembleResidual(*u_rcp,*z_rcp);
    fem->getResidual()->norm2(res.view(0,1));
    *outStream << "Residual Norm: " << res[0] << std::endl;
    errorFlag += (res[0] > 1.e-6 ? 1 : 0);
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
