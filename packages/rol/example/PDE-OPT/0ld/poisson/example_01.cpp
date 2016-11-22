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
#include "ROL_BoundConstraint_SimOpt.hpp"
#include "ROL_TpetraMultiVector.hpp"

#include <iostream>
#include <algorithm>

#include "data.hpp"
#include "objective.hpp"
#include "constraint.hpp"

typedef double RealT;

int main(int argc, char *argv[]) {

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
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

    RealT z_lo_bound = parlist->sublist("Problem").get("Control lower bound", ROL::ROL_NINF<RealT>());
    RealT z_up_bound = parlist->sublist("Problem").get("Control upper bound",  ROL::ROL_INF<RealT>());
    RealT u_lo_bound = parlist->sublist("Problem").get("State lower bound", ROL::ROL_NINF<RealT>());
    RealT u_up_bound = parlist->sublist("Problem").get("State upper bound",  ROL::ROL_INF<RealT>());

    /*** Initialize main data structure. ***/
    Teuchos::RCP<PoissonData<RealT> > data = Teuchos::rcp(new PoissonData<RealT>(comm, parlist, outStream));

    /*** Build vectors and dress them up as ROL vectors. ***/
    Teuchos::RCP<const Tpetra::Map<> > vecmap_u = data->getMatA()->getDomainMap();
    Teuchos::RCP<const Tpetra::Map<> > vecmap_z = data->getMatB()->getDomainMap();
    Teuchos::RCP<const Tpetra::Map<> > vecmap_c = data->getMatA()->getRangeMap();
    Teuchos::RCP<Tpetra::MultiVector<> > u_rcp = Teuchos::rcp(new Tpetra::MultiVector<>(vecmap_u, 1, true));
    Teuchos::RCP<Tpetra::MultiVector<> > p_rcp = Teuchos::rcp(new Tpetra::MultiVector<>(vecmap_u, 1, true));
    Teuchos::RCP<Tpetra::MultiVector<> > u_lo_rcp = Teuchos::rcp(new Tpetra::MultiVector<>(vecmap_u, 1, true));
    Teuchos::RCP<Tpetra::MultiVector<> > u_up_rcp = Teuchos::rcp(new Tpetra::MultiVector<>(vecmap_u, 1, true));
    Teuchos::RCP<Tpetra::MultiVector<> > z_rcp = Teuchos::rcp(new Tpetra::MultiVector<>(vecmap_z, 1, true));
    Teuchos::RCP<Tpetra::MultiVector<> > z_lo_rcp = Teuchos::rcp(new Tpetra::MultiVector<>(vecmap_z, 1, true));
    Teuchos::RCP<Tpetra::MultiVector<> > z_up_rcp = Teuchos::rcp(new Tpetra::MultiVector<>(vecmap_z, 1, true));
    Teuchos::RCP<Tpetra::MultiVector<> > c_rcp = Teuchos::rcp(new Tpetra::MultiVector<>(vecmap_c, 1, true));
    Teuchos::RCP<Tpetra::MultiVector<> > du_rcp = Teuchos::rcp(new Tpetra::MultiVector<>(vecmap_u, 1, true));
    Teuchos::RCP<Tpetra::MultiVector<> > dz_rcp = Teuchos::rcp(new Tpetra::MultiVector<>(vecmap_z, 1, true));
    // Set all values to 1 in u, z and c.
    u_rcp->putScalar(1.0);
    p_rcp->putScalar(1.0);
    u_lo_rcp->putScalar(u_lo_bound);
    u_up_rcp->putScalar(u_up_bound);
    z_rcp->putScalar(1.0);
    z_lo_rcp->putScalar(z_lo_bound);
    z_up_rcp->putScalar(z_up_bound);
    c_rcp->putScalar(1.0);
    // Randomize d vectors.
    du_rcp->randomize();
    dz_rcp->randomize();
    // Create ROL::TpetraMultiVectors.
    Teuchos::RCP<ROL::Vector<RealT> > up = Teuchos::rcp(new ROL::TpetraMultiVector<RealT>(u_rcp));
    Teuchos::RCP<ROL::Vector<RealT> > pp = Teuchos::rcp(new ROL::TpetraMultiVector<RealT>(p_rcp));
    Teuchos::RCP<ROL::Vector<RealT> > u_lo_p = Teuchos::rcp(new ROL::TpetraMultiVector<RealT>(u_lo_rcp));
    Teuchos::RCP<ROL::Vector<RealT> > u_up_p = Teuchos::rcp(new ROL::TpetraMultiVector<RealT>(u_up_rcp));
    Teuchos::RCP<ROL::Vector<RealT> > zp = Teuchos::rcp(new ROL::TpetraMultiVector<RealT>(z_rcp));
    Teuchos::RCP<ROL::Vector<RealT> > z_lo_p = Teuchos::rcp(new ROL::TpetraMultiVector<RealT>(z_lo_rcp));
    Teuchos::RCP<ROL::Vector<RealT> > z_up_p = Teuchos::rcp(new ROL::TpetraMultiVector<RealT>(z_up_rcp));
    Teuchos::RCP<ROL::Vector<RealT> > cp = Teuchos::rcp(new ROL::TpetraMultiVector<RealT>(c_rcp));
    Teuchos::RCP<ROL::Vector<RealT> > dup = Teuchos::rcp(new ROL::TpetraMultiVector<RealT>(du_rcp));
    Teuchos::RCP<ROL::Vector<RealT> > dzp = Teuchos::rcp(new ROL::TpetraMultiVector<RealT>(dz_rcp));
    // Create ROL SimOpt vectors.
    ROL::Vector_SimOpt<RealT> x(up,zp);
    ROL::Vector_SimOpt<RealT> d(dup,dzp);

    /*** Build objective function, equality constraint, reduced objective function and bound constraint. ***/
    Teuchos::RCP<ROL::Objective_SimOpt<RealT> > obj =
      Teuchos::rcp(new Objective_PDEOPT_Poisson<RealT>(data, parlist));
    Teuchos::RCP<ROL::EqualityConstraint_SimOpt<RealT> > con =
      Teuchos::rcp(new EqualityConstraint_PDEOPT_Poisson<RealT>(data, parlist));
    Teuchos::RCP<ROL::Objective<RealT> > objReduced =
      Teuchos::rcp(new ROL::Reduced_Objective_SimOpt<RealT>(obj, con, up, pp));
    Teuchos::RCP<ROL::BoundConstraint<RealT> > bcon_control =
      Teuchos::rcp(new ROL::BoundConstraint<RealT>(z_lo_p, z_up_p));
    Teuchos::RCP<ROL::BoundConstraint<RealT> > bcon_state =
      Teuchos::rcp(new ROL::BoundConstraint<RealT>(u_lo_p, u_up_p));
    // Initialize SimOpt bound constraint.
    Teuchos::RCP<ROL::BoundConstraint<RealT> > bcon_simopt =
      Teuchos::rcp(new ROL::BoundConstraint_SimOpt<RealT>(bcon_state,bcon_control));

    /*** Check functional interface. ***/
    obj->checkGradient(x,d,true,*outStream);
    obj->checkHessVec(x,d,true,*outStream);
    con->checkApplyJacobian(x,d,*up,true,*outStream);
    con->checkApplyAdjointHessian(x,*dup,d,x,true,*outStream);
    con->checkAdjointConsistencyJacobian(*dup,d,x,true,*outStream);
    con->checkInverseJacobian_1(*up,*up,*up,*zp,true,*outStream);
    con->checkInverseAdjointJacobian_1(*up,*up,*up,*zp,true,*outStream);
    objReduced->checkGradient(*zp,*dzp,true,*outStream);
    objReduced->checkHessVec(*zp,*dzp,true,*outStream);

    /*** Solve optimization problem. ***/

    ROL::Algorithm<RealT> algo_tr("Trust Region",*parlist,false);
    zp->zero(); // set zero initial guess
    algo_tr.run(*zp, *objReduced, *bcon_control, true, *outStream);
    RealT tol = 1e-8;
    con->solve(*cp, *up, *zp, tol);

    data->outputTpetraVector(u_rcp, "state.txt");
    data->outputTpetraVector(z_rcp, "control.txt");

    ROL::MoreauYosidaPenalty<RealT> obj_my(obj, bcon_simopt, x, 10.0);
    ROL::Algorithm<RealT> algo_my("Moreau-Yosida Penalty", *parlist, false);
    x.zero(); // set zero initial guess
    std::vector<std::string> algo_output;
    algo_output = algo_my.run(x, *cp, obj_my, *con, *bcon_simopt, true, *outStream);
    for (unsigned i=0; i<algo_output.size(); ++i) {
      *outStream << algo_output[i];
    }
    
    //ROL::Algorithm<RealT> algo_cs("Composite Step",*parlist,false);
    //x.zero(); // set zero initial guess
    //algo_cs.run(x, *cp, *obj, *con, true, *outStream);

    data->outputTpetraVector(u_rcp, "stateFS.txt");
    data->outputTpetraVector(z_rcp, "controlFS.txt");

    *outStream << std::endl << "|| u_approx - u_analytic ||_L2 = "
               << data->computeStateError(u_rcp) << std::endl;

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
