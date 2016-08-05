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

/*! \file  example_04.cpp
    \brief Shows how to solve the Stefan_Boltzmann problem in new pde framework.
*/

#include "Teuchos_Comm.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "Tpetra_DefaultPlatform.hpp"
#include "Tpetra_Version.hpp"

#include <iostream>
#include <algorithm>

#include "ROL_TpetraMultiVector.hpp"
#include "ROL_Algorithm.hpp"

//#include "../TOOLS/meshmanager.hpp"
#include "mesh_stefan_boltzmann.hpp"
#include "../TOOLS/pdeconstraint.hpp"
#include "../TOOLS/pdeobjective.hpp"
#include "pde_stefan_boltzmann.hpp"
#include "obj_stefan_boltzmann.hpp"

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
      = Teuchos::rcp(new MeshManager_Stefan_Boltzmann<RealT>(*parlist));
    std::cout << "Created mesh manager" << std::endl;
    // Initialize PDE describe Stefan Boltzmann's equation
    Teuchos::RCP<PDE_Stefan_Boltzmann<RealT> > pde
      = Teuchos::rcp(new PDE_Stefan_Boltzmann<RealT>(*parlist));
    std::cout << "Created PDE" << std::endl;
    Teuchos::RCP<PDE_Constraint<RealT> > con
      = Teuchos::rcp(new PDE_Constraint<RealT>(pde,meshMgr,comm,*parlist,*outStream));
    std::cout << "Created constraint" << std::endl;
    con->getAssembler()->printMeshData(*outStream);
    // Initialize quadratic objective function
    std::vector<Teuchos::RCP<QoI<RealT> > > qoi_vec(2,Teuchos::null);
    qoi_vec[0] = Teuchos::rcp(new QoI_L2Tracking_Stefan_Boltzmann<RealT>(pde->getFE_VOL()));
    qoi_vec[1] = Teuchos::rcp(new QoI_L2Penalty_Stefan_Boltzmann <RealT>(pde->getFE_VOL()));
    Teuchos::RCP<StdObjective_Stefan_Boltzmann<RealT> > std_obj
      = Teuchos::rcp(new StdObjective_Stefan_Boltzmann<RealT>(*parlist));
    Teuchos::RCP<PDE_Objective<RealT> > obj
      = Teuchos::rcp(new PDE_Objective<RealT>(qoi_vec,std_obj,con->getAssembler()));

    // Create state vector and set to zeroes
    Teuchos::RCP<Tpetra::MultiVector<> > u_rcp = con->getAssembler()->createStateVector();
    u_rcp->randomize();
    Teuchos::RCP<ROL::Vector<RealT> > up
      = Teuchos::rcp(new ROL::TpetraMultiVector<RealT>(u_rcp));
    // Create control vector and set to ones
    Teuchos::RCP<Tpetra::MultiVector<> > z_rcp = con->getAssembler()->createControlVector();
    z_rcp->putScalar(1.0);
    Teuchos::RCP<ROL::Vector<RealT> > zp
      = Teuchos::rcp(new ROL::TpetraMultiVector<RealT>(z_rcp));
    // Create residual vector and set to zeros
    Teuchos::RCP<Tpetra::MultiVector<> > r_rcp = con->getAssembler()->createResidualVector();
    r_rcp->putScalar(0.0);
    Teuchos::RCP<ROL::Vector<RealT> > rp
      = Teuchos::rcp(new ROL::TpetraMultiVector<RealT>(r_rcp));
    // Create state direction vector and set to random
    Teuchos::RCP<Tpetra::MultiVector<> > du_rcp = con->getAssembler()->createStateVector();
    du_rcp->randomize();
    Teuchos::RCP<ROL::Vector<RealT> > dup
      = Teuchos::rcp(new ROL::TpetraMultiVector<RealT>(du_rcp));
    // Create control direction vector and set to random
    Teuchos::RCP<Tpetra::MultiVector<> > dz_rcp = con->getAssembler()->createControlVector();
    dz_rcp->putScalar(0.0);
    //dz_rcp->randomize();
    Teuchos::RCP<ROL::Vector<RealT> > dzp
      = Teuchos::rcp(new ROL::TpetraMultiVector<RealT>(dz_rcp));
    // Create ROL SimOpt vectors
    ROL::Vector_SimOpt<RealT> x(up,zp);
    ROL::Vector_SimOpt<RealT> d(dup,dzp);

    // Run derivative checks
    obj->checkGradient(x,d,true,*outStream);
    obj->checkHessVec(x,d,true,*outStream);
    con->checkApplyJacobian(x,d,*up,true,*outStream);
    con->checkApplyAdjointHessian(x,*dup,d,x,true,*outStream);
    con->checkAdjointConsistencyJacobian(*dup,d,x,true,*outStream);
    con->checkInverseJacobian_1(*up,*up,*up,*zp,true,*outStream);
    con->checkInverseAdjointJacobian_1(*up,*up,*up,*zp,true,*outStream);

    ROL::Algorithm<RealT> algo("Composite Step",*parlist,false);
    algo.run(x,*rp,*obj,*con,true,*outStream);

    RealT tol(1.e-8);
    con->solve(*rp,*up,*zp,tol);
    con->getAssembler()->outputTpetraVector(u_rcp,"state.txt");
    con->getAssembler()->outputTpetraVector(z_rcp,"control.txt");

    Teuchos::Array<RealT> res(1,0);
    con->value(*rp,*up,*zp,tol);
    r_rcp->norm2(res.view(0,1));
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
