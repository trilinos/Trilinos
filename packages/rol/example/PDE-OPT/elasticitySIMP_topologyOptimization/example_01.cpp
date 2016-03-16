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
    \brief Shows how to solve the mother problem of PDE-constrained optimization:
*/

#include "Teuchos_Comm.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "Tpetra_DefaultPlatform.hpp"
#include "Tpetra_Version.hpp"

#include "ROL_Algorithm.hpp"
#include "ROL_TrustRegionStep.hpp"
#include "ROL_CompositeStep.hpp"
#include "ROL_Reduced_Objective_SimOpt.hpp"
#include "ROL_TpetraMultiVector.hpp"

#include <iostream>
#include <algorithm>

#include "data.hpp"
#include "objective.hpp"
#include "constraint.hpp"

typedef double RealT;

int main(int argc, char *argv[]) {

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint = argc - 1;
  Teuchos::RCP<std::ostream> outStream;
  Teuchos::oblackholestream bhs; // outputs nothing

  /*** Initialize communicator. ***/
  Teuchos::GlobalMPISession mpiSession (&argc, &argv, &bhs);
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::DefaultPlatform::getDefaultPlatform().getComm();
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
    Teuchos::RCP<ElasticitySIMPOperators<RealT> > data = Teuchos::rcp(new ElasticitySIMPOperators<RealT>(comm, parlist, outStream));
    /*** Build vectors and dress them up as ROL vectors. ***/
    Teuchos::RCP<const Tpetra::Map<> > vecmap_u = data->getDomainMapA();
    Teuchos::RCP<const Tpetra::Map<> > vecmap_z = data->getCellMap();
    Teuchos::RCP<Tpetra::MultiVector<> > u_rcp = Teuchos::rcp(new Tpetra::MultiVector<>(vecmap_u, 1, true));
    Teuchos::RCP<Tpetra::MultiVector<> > z_rcp = Teuchos::rcp(new Tpetra::MultiVector<>(vecmap_z, 1, true));
    Teuchos::RCP<Tpetra::MultiVector<> > du_rcp = Teuchos::rcp(new Tpetra::MultiVector<>(vecmap_u, 1, true));
    Teuchos::RCP<Tpetra::MultiVector<> > dw_rcp = Teuchos::rcp(new Tpetra::MultiVector<>(vecmap_u, 1, true));
    Teuchos::RCP<Tpetra::MultiVector<> > dz_rcp = Teuchos::rcp(new Tpetra::MultiVector<>(vecmap_z, 1, true));
    Teuchos::RCP<Tpetra::MultiVector<> > dz2_rcp = Teuchos::rcp(new Tpetra::MultiVector<>(vecmap_z, 1, true));
    // Set all values to 1 in u, z.
    u_rcp->putScalar(1.0);
    z_rcp->putScalar(1.0);

   //test     
   /*data->updateMaterialDensity (z_rcp);
    Teuchos::RCP<Tpetra::MultiVector<RealT> > rhs = Teuchos::rcp(new Tpetra::MultiVector<> (data->getVecF()->getMap(), 1, true));
    data->ApplyMatAToVec(rhs, u_rcp);
    data->outputTpetraVector(rhs, "KU0.txt");
    data->ApplyInverseJacobian1ToVec(u_rcp, rhs, false);
    data->outputTpetraVector(u_rcp, "KKU0.txt");
    
    data->ApplyJacobian1ToVec(rhs, u_rcp);
    data->outputTpetraVector(rhs, "KU1.txt");
    data->ApplyInverseJacobian1ToVec(u_rcp, rhs, false);
    data->outputTpetraVector(u_rcp, "KKU1.txt");
  */
    u_rcp->putScalar(1.0);
    z_rcp->putScalar(1.0);
    // Randomize d vectors.
    du_rcp->randomize();
    dw_rcp->randomize();
    dz_rcp->randomize();
    dz2_rcp->randomize();
    // Create ROL::TpetraMultiVectors.
    Teuchos::RCP<ROL::Vector<RealT> > up = Teuchos::rcp(new ROL::TpetraMultiVector<RealT>(u_rcp));
    Teuchos::RCP<ROL::Vector<RealT> > zp = Teuchos::rcp(new ROL::TpetraMultiVector<RealT>(z_rcp));
    Teuchos::RCP<ROL::Vector<RealT> > dup = Teuchos::rcp(new ROL::TpetraMultiVector<RealT>(du_rcp));
    Teuchos::RCP<ROL::Vector<RealT> > dwp = Teuchos::rcp(new ROL::TpetraMultiVector<RealT>(dw_rcp));
    Teuchos::RCP<ROL::Vector<RealT> > dzp = Teuchos::rcp(new ROL::TpetraMultiVector<RealT>(dz_rcp));
    Teuchos::RCP<ROL::Vector<RealT> > dz2p = Teuchos::rcp(new ROL::TpetraMultiVector<RealT>(dz2_rcp));
    // Create ROL SimOpt vectors.
    ROL::Vector_SimOpt<RealT> x(up,zp);
    ROL::Vector_SimOpt<RealT> d(dup,dzp);
    ROL::Vector_SimOpt<RealT> d2(dwp,dz2p);

    /*** Build objective function, constraint and reduced objective function. ***/
    Teuchos::RCP<Objective_PDEOPT_ElasticitySIMP<RealT> > obj = Teuchos::rcp(new Objective_PDEOPT_ElasticitySIMP<RealT>(data, parlist));
    Teuchos::RCP<EqualityConstraint_PDEOPT_ElasticitySIMP<RealT> > con = Teuchos::rcp(new EqualityConstraint_PDEOPT_ElasticitySIMP<RealT>(data, parlist));
    Teuchos::RCP<ROL::Reduced_Objective_SimOpt<RealT> > objReduced = Teuchos::rcp(new ROL::Reduced_Objective_SimOpt<RealT>(obj, con, up, dwp));
    /*** Check functional interface. ***/
    *outStream << "Checking Objective:" << "\n";
    obj->checkGradient(x,d,true,*outStream);
    obj->checkHessVec(x,d,true,*outStream);
    obj->checkHessSym(x,d,d2,true,*outStream);
    *outStream << "Checking Constraint:" << "\n";
    con->checkAdjointConsistencyJacobian(*dup,d,x,true,*outStream);
    con->checkAdjointConsistencyJacobian_1(*dwp, *dup, *up, *zp, true, *outStream);
    con->checkAdjointConsistencyJacobian_2(*dwp, *dzp, *up, *zp, true, *outStream);
    con->checkInverseJacobian_1(*up,*up,*up,*zp,true,*outStream);
    con->checkInverseAdjointJacobian_1(*up,*up,*up,*zp,true,*outStream);
    con->checkApplyJacobian(x,d,*up,true,*outStream);
    con->checkApplyAdjointHessian(x,*dup,d,x,true,*outStream);
    *outStream << "Checking Reduced Objective:" << "\n";
    objReduced->checkGradient(*zp,*dzp,true,*outStream);
    objReduced->checkHessVec(*zp,*dzp,true,*outStream);
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
