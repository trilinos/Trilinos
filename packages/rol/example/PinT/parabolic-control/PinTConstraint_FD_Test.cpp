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
    if(expr) { \
      std::stringstream ss; \
      ss << "Assertion failed on line " << __LINE__ << std::endl; \
      throw std::logic_error(ss.str()); \
    }

#define CHECK_EQUALITY(expr1,expr2) \
    if(expr1!=expr2) { \
      std::stringstream ss; \
      ss << myRank << ". Equality assertion failed on line " << __LINE__ << std::endl; \
      ss << myRank << ".  " << expr1 << " != " << expr2 << std::endl; \
      throw std::logic_error(ss.str()); \
    } \
    std::cout << myRank << ".  CHECK_EQUALITY line " << __LINE__ << " (passed): " << expr1 << " == " << expr2 << std::endl; \

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
  auto outStream = ROL::makeStreamPtr( std::cout, argc > 1 );

  int errorFlag  = 0;

  try {

    // because of the splitting this actually runs two jobs at once
    run_test(MPI_COMM_WORLD, outStream);

  }
  catch (std::logic_error err) {
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
  using Bounds            = ROL::Bounds<RealT>;
  using PartitionedVector = ROL::PartitionedVector<RealT>;

  int numRanks = -1;
  int myRank = -1;

  MPI_Comm_size(comm, &numRanks);
  MPI_Comm_rank(comm, &myRank);

  *outStream << "Proc " << myRank << "/" << numRanks << std::endl;

  ROL::Ptr<const ROL::PinTCommunicators> communicators = ROL::makePtr<ROL::PinTCommunicators>(comm,1);
  ROL::Ptr<const ROL::PinTVectorCommunication<RealT>> vectorComm = ROL::makePtr<ROL::PinTVectorCommunication_StdVector<RealT>>();

  // Parse input parameter list
  ROL::Ptr<ROL::ParameterList> pl = ROL::getParametersFromXmlFile("input_ex01.xml");
  bool derivCheck = pl->get("Derivative Check",         true); // Check derivatives.
  uint nx         = pl->get("Spatial Discretization",     64); // Set spatial discretization.
  uint nt         = pl->get("Temporal Discretization",   100); // Set temporal discretization.
  RealT T         = pl->get("End Time",                  1.0); // Set end time.
  RealT dt        = T/(static_cast<RealT>(nt)-1.0);

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

  double tol = 1e-10;

  // check the pint constraint
  {
    auto c   = state->clone();
    auto jv  = c->clone();
    auto v_u = state->clone();
    ROL::RandomizeVector<RealT>(*state);
    ROL::RandomizeVector<RealT>(*control);
    ROL::RandomizeVector<RealT>(*v_u);

    pint_constraint.setGlobalScale(1.037);
    pint_constraint.checkSolve(*state,*control,*c,true,*outStream);
    pint_constraint.checkApplyJacobian_1(*state,*control,*v_u,*jv,true,*outStream);
    RealT inv_1     = pint_constraint.checkInverseJacobian_1(*jv,*v_u,*state,*control,true,*outStream);

    CHECK_ASSERT(inv_1 > tol);

    return;

    RealT adj_inv_1 = pint_constraint.checkInverseAdjointJacobian_1(*jv,*v_u,*state,*control,true,*outStream);
    if(adj_inv_1 > tol) {
      std::stringstream ss;
      ss << "Adjoint Jacobian inverse inversion FAILED: error = " << adj_inv_1  << std::endl;
      throw std::logic_error(ss.str());
    }
    
  }

  auto x   = makePtr<ROL::Vector_SimOpt<RealT>>(state,control);
  auto v_1 = x->clone();
  auto v_2 = state->clone();
  auto r_1 = v_1->clone();
  auto r_2 = v_2->clone();
  ROL::RandomizeVector<RealT>(*r_1);
  ROL::RandomizeVector<RealT>(*r_2);

  // check parallel communication for invertTimeStep methods
  {
    using ROL::PinTVector;

    auto v  = state->clone();
    auto fv  = state->clone();
    auto afv  = state->clone();
    auto jv = state->clone();
    auto ajv = state->clone();
    PinTVector<RealT> & pint_v = dynamic_cast<PinTVector<RealT>&>(*v);
    PinTVector<RealT> & pint_jv = dynamic_cast<PinTVector<RealT>&>(*jv);
    PinTVector<RealT> & pint_ajv = dynamic_cast<PinTVector<RealT>&>(*ajv);

    ROL::RandomizeVector<RealT>(*v);

    double tol = 1e-12;

    // make sure we are globally consistent
    state->boundaryExchange();
    control->boundaryExchange();
    pint_v.boundaryExchange();

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

    pint_jv.boundaryExchange();

    // fix up right hand side to give the right solution for the sub domain solvers
    if(myRank!=0) {
      pint_jv.getVectorPtr(-1)->axpy(1.0,*pint_v.getRemoteBufferPtr(-1)); 

      pint_ajv.getRemoteBufferPtr(-1)->set(*pint_v.getVectorPtr(-1));
    }

    pint_ajv.boundaryExchangeSumInto();

    int level = 0;
    pint_constraint.invertTimeStepJacobian(dynamic_cast<PinTVector<RealT>&>(*fv),
                                           dynamic_cast<const PinTVector<RealT>&>(*jv),
                                           dynamic_cast<const PinTVector<RealT>&>(*state),
                                           dynamic_cast<const PinTVector<RealT>&>(*control),
                                           tol,
                                           level);

    pint_constraint.invertAdjointTimeStepJacobian(dynamic_cast<PinTVector<RealT>&>(*afv),
                                                  dynamic_cast<const PinTVector<RealT>&>(*ajv),
                                                  dynamic_cast<const PinTVector<RealT>&>(*state),
                                                  dynamic_cast<const PinTVector<RealT>&>(*control),
                                                  tol,
                                                  level);

    RealT v_norm   = v->norm();

    // check forward time step jacobian
    {
      fv->axpy(-1.0,*v);
  
      RealT fv_norm  = fv->norm();
  
      if(myRank==0) {
        *outStream << "Testing vector norm = " << v_norm << std::endl;
      }
  
      // check norms
      if(fv_norm/v_norm > 1e-13) {
        std::stringstream ss;
        ss << "Forward block Jacobi subdomain inversion FAILED with proc " 
           << " (relative error = " << fv_norm / v_norm  << ")" << std::endl;
        throw std::logic_error(ss.str());
      }
      else {
        if(myRank==0)
          *outStream << "Forward block Jacobi subdomain inversion PASSED with proc " 
                     << " (relative error = " << fv_norm / v_norm  << ")" << std::endl;
      }
    }
    
    // check backward time step jacobian
    {
      afv->axpy(-1.0,*v);
      RealT afv_norm = afv->norm();
  
      // check norms (adjoint)
      if(afv_norm/v_norm > 1e-13) {
        std::stringstream ss;
        ss << "Adjoint block Jacobi subdomain inversion FAILED with proc " 
           << " (relative error = " << afv_norm / v_norm  << ")" << std::endl;
        throw std::logic_error(ss.str());
      }
      else {
        if(myRank==0)
          *outStream << "Adjoint block Jacobi subdomain inversion PASSED with proc " 
                     << " (relative error = " << afv_norm / v_norm  << ")" << std::endl;
      }
    }

  }
}
