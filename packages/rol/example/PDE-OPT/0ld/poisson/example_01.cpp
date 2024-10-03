// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  example_01.cpp
    \brief Shows how to solve the mother problem of PDE-constrained optimization:
*/

#include "Teuchos_Comm.hpp"
#include "ROL_Stream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "Tpetra_Core.hpp"
#include "Tpetra_Version.hpp"

#include "ROL_Algorithm.hpp"
#include "ROL_ConstraintStatusTest.hpp"
#include "ROL_Bounds.hpp"
#include "ROL_TrustRegionStep.hpp"
#include "ROL_MoreauYosidaPenaltyStep.hpp"
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
  ROL::Ptr<std::ostream> outStream;
  ROL::nullstream bhs; // outputs nothing

  /*** Initialize communicator. ***/
  Teuchos::GlobalMPISession mpiSession (&argc, &argv, &bhs);
  ROL::Ptr<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
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

    RealT z_lo_bound = parlist->sublist("Problem").get("Control lower bound", ROL::ROL_NINF<RealT>());
    RealT z_up_bound = parlist->sublist("Problem").get("Control upper bound",  ROL::ROL_INF<RealT>());
    RealT u_lo_bound = parlist->sublist("Problem").get("State lower bound", ROL::ROL_NINF<RealT>());
    RealT u_up_bound = parlist->sublist("Problem").get("State upper bound",  ROL::ROL_INF<RealT>());

    /*** Initialize main data structure. ***/
    ROL::Ptr<PoissonData<RealT> > data = ROL::makePtr<PoissonData<RealT>>(comm, parlist, outStream);

    /*** Build vectors and dress them up as ROL vectors. ***/
    ROL::Ptr<const Tpetra::Map<> > vecmap_u = data->getMatA()->getDomainMap();
    ROL::Ptr<const Tpetra::Map<> > vecmap_z = data->getMatB()->getDomainMap();
    ROL::Ptr<const Tpetra::Map<> > vecmap_c = data->getMatA()->getRangeMap();
    ROL::Ptr<Tpetra::MultiVector<> > u_ptr = ROL::makePtr<Tpetra::MultiVector<>>(vecmap_u, 1, true);
    ROL::Ptr<Tpetra::MultiVector<> > p_ptr = ROL::makePtr<Tpetra::MultiVector<>>(vecmap_u, 1, true);
    ROL::Ptr<Tpetra::MultiVector<> > u_lo_ptr = ROL::makePtr<Tpetra::MultiVector<>>(vecmap_u, 1, true);
    ROL::Ptr<Tpetra::MultiVector<> > u_up_ptr = ROL::makePtr<Tpetra::MultiVector<>>(vecmap_u, 1, true);
    ROL::Ptr<Tpetra::MultiVector<> > z_ptr = ROL::makePtr<Tpetra::MultiVector<>>(vecmap_z, 1, true);
    ROL::Ptr<Tpetra::MultiVector<> > z_lo_ptr = ROL::makePtr<Tpetra::MultiVector<>>(vecmap_z, 1, true);
    ROL::Ptr<Tpetra::MultiVector<> > z_up_ptr = ROL::makePtr<Tpetra::MultiVector<>>(vecmap_z, 1, true);
    ROL::Ptr<Tpetra::MultiVector<> > c_ptr = ROL::makePtr<Tpetra::MultiVector<>>(vecmap_c, 1, true);
    ROL::Ptr<Tpetra::MultiVector<> > du_ptr = ROL::makePtr<Tpetra::MultiVector<>>(vecmap_u, 1, true);
    ROL::Ptr<Tpetra::MultiVector<> > dz_ptr = ROL::makePtr<Tpetra::MultiVector<>>(vecmap_z, 1, true);
    // Set all values to 1 in u, z and c.
    u_ptr->putScalar(1.0);
    p_ptr->putScalar(1.0);
    u_lo_ptr->putScalar(u_lo_bound);
    u_up_ptr->putScalar(u_up_bound);
    z_ptr->putScalar(1.0);
    z_lo_ptr->putScalar(z_lo_bound);
    z_up_ptr->putScalar(z_up_bound);
    c_ptr->putScalar(1.0);
    // Randomize d vectors.
    du_ptr->randomize();
    dz_ptr->randomize();
    // Create ROL::TpetraMultiVectors.
    ROL::Ptr<ROL::Vector<RealT> > up = ROL::makePtr<ROL::TpetraMultiVector<RealT>>(u_ptr);
    ROL::Ptr<ROL::Vector<RealT> > pp = ROL::makePtr<ROL::TpetraMultiVector<RealT>>(p_ptr);
    ROL::Ptr<ROL::Vector<RealT> > u_lo_p = ROL::makePtr<ROL::TpetraMultiVector<RealT>>(u_lo_ptr);
    ROL::Ptr<ROL::Vector<RealT> > u_up_p = ROL::makePtr<ROL::TpetraMultiVector<RealT>>(u_up_ptr);
    ROL::Ptr<ROL::Vector<RealT> > zp = ROL::makePtr<ROL::TpetraMultiVector<RealT>>(z_ptr);
    ROL::Ptr<ROL::Vector<RealT> > z_lo_p = ROL::makePtr<ROL::TpetraMultiVector<RealT>>(z_lo_ptr);
    ROL::Ptr<ROL::Vector<RealT> > z_up_p = ROL::makePtr<ROL::TpetraMultiVector<RealT>>(z_up_ptr);
    ROL::Ptr<ROL::Vector<RealT> > cp = ROL::makePtr<ROL::TpetraMultiVector<RealT>>(c_ptr);
    ROL::Ptr<ROL::Vector<RealT> > dup = ROL::makePtr<ROL::TpetraMultiVector<RealT>>(du_ptr);
    ROL::Ptr<ROL::Vector<RealT> > dzp = ROL::makePtr<ROL::TpetraMultiVector<RealT>>(dz_ptr);
    // Create ROL SimOpt vectors.
    ROL::Vector_SimOpt<RealT> x(up,zp);
    ROL::Vector_SimOpt<RealT> d(dup,dzp);

    /*** Build objective function, equality constraint, reduced objective function and bound constraint. ***/
    ROL::Ptr<ROL::Objective_SimOpt<RealT> > obj =
      ROL::makePtr<Objective_PDEOPT_Poisson<RealT>>(data, parlist);
    ROL::Ptr<ROL::Constraint_SimOpt<RealT> > con =
      ROL::makePtr<EqualityConstraint_PDEOPT_Poisson<RealT>>(data, parlist);
    ROL::Ptr<ROL::Objective<RealT> > objReduced =
      ROL::makePtr<ROL::Reduced_Objective_SimOpt<RealT>>(obj, con, up, zp, pp);
    ROL::Ptr<ROL::BoundConstraint<RealT> > bcon_control =
      ROL::makePtr<ROL::Bounds<RealT>>(z_lo_p, z_up_p);
    ROL::Ptr<ROL::BoundConstraint<RealT> > bcon_state =
      ROL::makePtr<ROL::Bounds<RealT>>(u_lo_p, u_up_p);
    // Initialize SimOpt bound constraint.
    ROL::Ptr<ROL::BoundConstraint<RealT> > bcon_simopt =
      ROL::makePtr<ROL::BoundConstraint_SimOpt<RealT>>(bcon_state,bcon_control);

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

    ROL::Ptr<ROL::Step<RealT>>         step_tr = ROL::makePtr<ROL::TrustRegionStep<RealT>>(*parlist);
    ROL::Ptr<ROL::StatusTest<RealT>> status_tr = ROL::makePtr<ROL::StatusTest<RealT>>(*parlist);
    ROL::Algorithm<RealT> algo_tr(step_tr,status_tr,false);
    zp->zero(); // set zero initial guess
    algo_tr.run(*zp, *objReduced, *bcon_control, true, *outStream);
    RealT tol = 1e-8;
    con->solve(*cp, *up, *zp, tol);

    data->outputTpetraVector(u_ptr, "state.txt");
    data->outputTpetraVector(z_ptr, "control.txt");

    ROL::Ptr<ROL::Step<RealT>>         step_my = ROL::makePtr<ROL::MoreauYosidaPenaltyStep<RealT>>(*parlist);
    ROL::Ptr<ROL::StatusTest<RealT>> status_my = ROL::makePtr<ROL::ConstraintStatusTest<RealT>>(*parlist);
    ROL::Algorithm<RealT> algo_my(step_my,status_my,false);
    ROL::MoreauYosidaPenalty<RealT> obj_my(obj, bcon_simopt, x, *parlist);
    x.zero(); // set zero initial guess
    std::vector<std::string> algo_output;
    algo_output = algo_my.run(x, *cp, obj_my, *con, *bcon_simopt, true, *outStream);
    for (unsigned i=0; i<algo_output.size(); ++i) {
      *outStream << algo_output[i];
    }
    
    //ROL::Algorithm<RealT> algo_cs("Composite Step",*parlist,false);
    //x.zero(); // set zero initial guess
    //algo_cs.run(x, *cp, *obj, *con, true, *outStream);

    data->outputTpetraVector(u_ptr, "stateFS.txt");
    data->outputTpetraVector(z_ptr, "controlFS.txt");

    *outStream << std::endl << "|| u_approx - u_analytic ||_L2 = "
               << data->computeStateError(u_ptr) << std::endl;

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
