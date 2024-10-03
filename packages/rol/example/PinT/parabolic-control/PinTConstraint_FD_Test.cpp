// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <iostream>
#include <iomanip>
#include <random>
#include <utility>

#include "Teuchos_GlobalMPISession.hpp"

#include "ROL_Stream.hpp"
#include "ROL_ParameterList.hpp"
#include "ROL_Bounds.hpp"
#include "ROL_RandomVector.hpp"
#include "ROL_Vector_SimOpt.hpp"
#include "ROL_PinTConstraint.hpp"
#include "ROL_PinTVectorCommunication_StdVector.hpp"

#include "dynamicConstraint.hpp"

#define CHECK_ASSERT(expr) \
    {bool test = expr; \
    if(not test) { \
      std::stringstream ss; \
      ss << myRank << ". FAILED - Assertion failed on line " << __LINE__ << ": " #expr " " << std::endl; \
      throw std::logic_error(ss.str()); \
    } \
    else if(myRank==0) { \
      std::stringstream ss; \
      ss << myRank << ". Assertion passed on line " << __LINE__ << ": " #expr " " << std::endl; \
      *outStream << ss.str() << std::endl; \
    }}

#define CHECK_EQUALITY(expr1,expr2) \
    if(expr1!=expr2) { \
      std::stringstream ss; \
      ss << myRank << ". FAILED - Equality assertion failed on line " << __LINE__ << std::endl; \
      ss << myRank << ".  " << expr1 << " != " << expr2 << std::endl; \
      throw std::logic_error(ss.str()); \
    } else if(myRank==0) \
    *outStream << myRank << ".  CHECK_EQUALITY line " << __LINE__ << " (passed): " << expr1 << " == " << expr2 << std::endl; \

using RealT = double;
using size_type = std::vector<RealT>::size_type;

void run_test(MPI_Comm comm, const ROL::Ptr<std::ostream> & outStream);

int main( int argc, char* argv[] ) 
{
  using ROL::Ptr;
  using ROL::makePtr;
  using ROL::makePtrFromRef;

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);  
  int myRank = Teuchos::GlobalMPISession::getRank();
//  int numProcs = Teuchos::GlobalMPISession::getNProc();

//  int iprint     = argc - 1;
  auto outStream = ROL::makeStreamPtr( std::cout, argc > 1 && myRank==0 );

  int errorFlag  = 0;

  try {

    // because of the splitting this actually runs two jobs at once
    run_test(MPI_COMM_WORLD, outStream);

  }
  catch (std::logic_error& err) {
    *outStream << err.what() << "\n";
    errorFlag = -1000;
  }; // end try

  if(myRank==0) {
    if (errorFlag != 0)
      std::cout << "End Result: TEST FAILED\n";
    else
      std::cout << "End Result: TEST PASSED\n";
  }

  return 0;
}

void run_test(MPI_Comm comm, const ROL::Ptr<std::ostream> & outStream)
{
  using ROL::Ptr;
  using ROL::makePtr;
  using ROL::makePtrFromRef;

  using RealT             = double;
  using size_type         = std::vector<RealT>::size_type;
//  using ValidateFunction  = ROL::ValidateFunction<RealT>;
//  using Bounds            = ROL::Bounds<RealT>;
//  using PartitionedVector = ROL::PartitionedVector<RealT>;

  int numRanks = -1;
  int myRank = -1;

  MPI_Comm_size(comm, &numRanks);
  MPI_Comm_rank(comm, &myRank);

  *outStream << "Proc " << myRank << "/" << numRanks << std::endl;

  ROL::Ptr<const ROL::PinTCommunicators> communicators = ROL::makePtr<ROL::PinTCommunicators>(comm,1);
  ROL::Ptr<const ROL::PinTVectorCommunication<RealT>> vectorComm = ROL::makePtr<ROL::PinTVectorCommunication_StdVector<RealT>>();

  // Parse input parameter list
  ROL::Ptr<ROL::ParameterList> pl = ROL::getParametersFromXmlFile("input_ex01.xml");
  // bool derivCheck = pl->get("Derivative Check",         true); // Check derivatives.
  uint nx         = pl->get("Spatial Discretization",     64); // Set spatial discretization.
  uint nt         = pl->get("Temporal Discretization",   128); // Set temporal discretization.
  RealT T         = pl->get("End Time",                  1.0); // Set end time.
  RealT dt        = T/(static_cast<RealT>(nt));

  // Initialize objective function.
  ROL::Ptr<ROL::DynamicConstraint<RealT>> dyn_con
    = ROL::makePtr<Constraint_ParabolicControl<RealT>>(*pl);

  // Create control vectors.
  ROL::Ptr<ROL::StdVector<RealT>>         zk = ROL::makePtr<ROL::StdVector<RealT>>(ROL::makePtr<std::vector<RealT>>(nx+2));

  // Create initial state vector.
  ROL::Ptr<ROL::StdVector<RealT>> u0 = ROL::makePtr<ROL::StdVector<RealT>>(ROL::makePtr<std::vector<RealT>>(nx,0.0));

  ROL::Ptr< ROL::PinTVector<RealT>> state;
  ROL::Ptr< ROL::PinTVector<RealT>> control;

  state        = ROL::buildStatePinTVector<RealT>(   communicators, vectorComm, nt,     u0); // for Euler, Crank-Nicolson, stencil = [-1,0]
  control      = ROL::buildControlPinTVector<RealT>( communicators, vectorComm, nt,     zk); // time discontinous, stencil = [0]

  int numLocalSteps = control->numOwnedSteps();
  auto timeStamp = ROL::makePtr<std::vector<ROL::TimeStamp<RealT>>>(numLocalSteps);
  for( size_type k=0; k<timeStamp->size(); ++k ) {
    timeStamp->at(k).t.resize(2);
    timeStamp->at(k).t.at(0) = k*dt;
    timeStamp->at(k).t.at(1) = (k+1)*dt;
  }

  ROL::PinTConstraint<RealT> pint_constraint(dyn_con,u0,timeStamp);
  pint_constraint.applyMultigrid(2,communicators,vectorComm);

  double tol = 1e-10;

  auto c   = state->clone();
  auto jv  = c->clone();
  auto v_u = state->clone();
  auto w_u = state->clone();
  auto v_z = control->clone();
  ROL::RandomizeVector<RealT>(*state);
  ROL::RandomizeVector<RealT>(*control);
  ROL::RandomizeVector<RealT>(*v_u);
  ROL::RandomizeVector<RealT>(*w_u);
  ROL::RandomizeVector<RealT>(*v_z);

  pint_constraint.setGlobalScale(1.0);

  // check the solve
  //////////////////////////////////////////////////////////////////////
  {
    if(myRank==0)
      *outStream << "Checking solve" << std::endl;

    double solveNorm = pint_constraint.checkSolve(*state,*control,*c,true,*outStream);

    CHECK_ASSERT(solveNorm < tol);
  }

  // check the Jacobian_1
  /////////////////////////////////////////////////////////////////////////////
  {
    if(myRank==0)
      *outStream << "Checking apply Jacobian 1" << std::endl;

    auto errors = pint_constraint.checkApplyJacobian_1(*state,*control,*v_u,*jv,true,*outStream);
    CHECK_ASSERT(errors[6][3]/errors[6][1] < 1e-6);
  }

  // check the Jacobian_2
  /////////////////////////////////////////////////////////////////////////////
  {
    if(myRank==0)
      *outStream << "Checking apply Jacobian 2" << std::endl;

    auto errors = pint_constraint.checkApplyJacobian_2(*state,*control,*v_z,*jv,true,*outStream);
    CHECK_ASSERT(errors[6][3]/errors[6][1] < 1e-6);
  }

  // check the Adjoint Jacobian_1
  /////////////////////////////////////////////////////////////////////////////
  {
    if(myRank==0)
      *outStream << "Checking apply Adjoint Jacobian 1" << std::endl;

    auto error = pint_constraint.checkAdjointConsistencyJacobian_1(*w_u,*v_u,*state,*control,true,*outStream);
    CHECK_ASSERT(error<1e-8);
  }

  // check the Adjoint Jacobian_2
  /////////////////////////////////////////////////////////////////////////////
  {
    if(myRank==0)
      *outStream << "Checking apply Adjoint Jacobian 2" << std::endl;

    auto error = pint_constraint.checkAdjointConsistencyJacobian_2(*w_u,*v_z,*state,*control,true,*outStream);
    CHECK_ASSERT(error<1e-8);
  }

  // check inverse Jacobian_1
  /////////////////////////////////////////////////////////////////////////////
  {
    RealT inv_1     = pint_constraint.checkInverseJacobian_1(*jv,*v_u,*state,*control,true,*outStream);

    CHECK_ASSERT(inv_1 < tol);
  }

  // check inverse adjoint Jacobian_1
  /////////////////////////////////////////////////////////////////////////////
  {
    RealT adj_inv_1 = pint_constraint.checkInverseAdjointJacobian_1(*jv,*v_u,*state,*control,true,*outStream);

    CHECK_ASSERT(adj_inv_1 < tol);
  }

  auto x   = makePtr<ROL::Vector_SimOpt<RealT>>(state,control);
  auto v_1 = x->clone();
  auto v_2 = state->clone();
  auto r_1 = v_1->clone();
  auto r_2 = v_2->clone();
  ROL::RandomizeVector<RealT>(*r_1);
  ROL::RandomizeVector<RealT>(*r_2);

  // check local forward and adjoint operators (invertTimeStep* and apply*Jacobian_1_leveled_approx)
  {
    using ROL::PinTVector;

    auto v  = state->clone();
    auto fv  = state->clone();
    auto lfv  = state->clone();   // local fv
    auto afv  = state->clone();
    auto jv = state->clone();
    auto ajv = state->clone();
    PinTVector<RealT> & pint_v = dynamic_cast<PinTVector<RealT>&>(*v);
    PinTVector<RealT> & pint_jv = dynamic_cast<PinTVector<RealT>&>(*jv);
    PinTVector<RealT> & pint_ajv = dynamic_cast<PinTVector<RealT>&>(*ajv);

    ROL::RandomizeVector<RealT>(*v);
    // v->setScalar(3.3);
    // state->setScalar(1.1);
    // control->setScalar(1.1);

    double tol = 1e-12;

    std::stringstream ss;

    // compute jacobian action
    pint_constraint.applyJacobian_1(*jv,
                                    *v,
                                    *state,
                                    *control,tol);

    // compute jacobian action
    pint_constraint.applyAdjointJacobian_1(*ajv,
                                           *v,
                                           *state,
                                           *control,tol);

    // fix up right hand side to give the right solution for the sub domain solvers
    int numSteps = pint_jv.numOwnedSteps();

    // for the adjoint call
    pint_v.boundaryExchangeRightToLeft();

    if(myRank+1<numRanks) {
      // the last rank is handled special to include the last constraint
      pint_jv.getVectorPtr(2*(numSteps-1)+1)->axpy(1.0,*pint_v.getVectorPtr(2*(numSteps-1))); 

      pint_ajv.getVectorPtr(2*(numSteps-1))->axpy(1.0,*pint_v.getVectorPtr(2*(numSteps-1)+1));
    }

    int level = 0;
    pint_constraint.invertTimeStepJacobian(dynamic_cast<PinTVector<RealT>&>(*fv),
                                           dynamic_cast<const PinTVector<RealT>&>(*jv),
                                           dynamic_cast<const PinTVector<RealT>&>(*state),
                                           dynamic_cast<const PinTVector<RealT>&>(*control),
                                           tol,
                                           level);

    pint_constraint.applyJacobian_1_leveled_approx(*lfv,
                                                   *fv,
                                                   *state,
                                                   *control,
                                                   tol,
                                                   level);

    RealT v_norm   = v->norm();
    RealT jv_norm  = jv->norm();

    // check forward time step jacobian
    {
      fv->axpy(-1.0,*v);
      RealT fv_norm  = fv->norm();
  
      if(myRank==0)
        *outStream << "Testing vector norm = " << v_norm  << " (error = " << fv_norm/v_norm << ")" << std::endl;

      // check norms
      CHECK_ASSERT(fv_norm/v_norm <= 1e-12);

      if(myRank==0)
        *outStream << "Forward block Jacobi subdomain inversion PASSED with proc " 
                   << " (relative error = " << fv_norm / v_norm  << ")" << std::endl;

      lfv->axpy(-1.0,*jv);
      RealT lfv_norm  = lfv->norm();
      CHECK_ASSERT(lfv_norm/jv_norm <= 1e-13);
    }

    pint_constraint.invertAdjointTimeStepJacobian(dynamic_cast<PinTVector<RealT>&>(*afv),
                                                  dynamic_cast<const PinTVector<RealT>&>(*ajv),
                                                  dynamic_cast<const PinTVector<RealT>&>(*state),
                                                  dynamic_cast<const PinTVector<RealT>&>(*control),
                                                  tol,
                                                  level);

    pint_constraint.applyAdjointJacobian_1_leveled_approx(*lfv,
                                                          *afv,
                                                          *state,
                                                          *control,
                                                          tol,
                                                          level);

    RealT ajv_norm  = ajv->norm();
    
    // check backward time step jacobian
    {
      afv->axpy(-1.0,*v);
      RealT afv_norm = afv->norm();

      afv->print(ss);

      if(myRank==0)
        *outStream << "Testing vector norm = " << v_norm  << " (error = " << afv_norm/v_norm << ")" << std::endl;

      // std::cout << "====================Q===============\n"
      //           << "PROCESSOR " << myRank << "\n"
      //           << ss.str() 
      //           << "\n==================Q=================\n" << std::endl;

      CHECK_ASSERT(afv_norm/v_norm <= 1e-12);
  
      if(myRank==0)
        *outStream << "Adjoint block Jacobi subdomain inversion PASSED with proc " 
                   << " (relative error = " << afv_norm / v_norm  << ")\n" << std::endl;
      
      lfv->axpy(-1.0,*ajv);
      RealT lfv_norm  = lfv->norm();

      CHECK_ASSERT(lfv_norm/ajv_norm <= 1e-13);

      if(myRank==0)
        *outStream << "Adjoint block Jacobi subdomain apply PASSED with proc " 
                   << " (relative error = " << lfv_norm / ajv_norm  << ")\n" << std::endl;
    }

  }
}
