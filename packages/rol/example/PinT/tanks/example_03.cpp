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

#include "Tanks_DynamicConstraint.hpp"
#include "Tanks_SerialConstraint.hpp"
#include "Tanks_ConstraintCheck.hpp"

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
  using State             = Tanks::StateVector<RealT>;
  using Control           = Tanks::ControlVector<RealT>;

  int numRanks = -1;
  int myRank = -1;

  MPI_Comm_size(comm, &numRanks);
  MPI_Comm_rank(comm, &myRank);

  *outStream << "Proc " << myRank << "/" << numRanks << std::endl;

  ROL::Ptr<const ROL::PinTCommunicators> communicators = ROL::makePtr<ROL::PinTCommunicators>(comm,1);

  ROL::Ptr< ROL::PinTVector<RealT>> state;
  ROL::Ptr< ROL::PinTVector<RealT>> control;

  auto  pl_ptr = ROL::getParametersFromXmlFile("tank-parameters.xml");
  auto& pl     = *pl_ptr;
  auto  con    = Tanks::DynamicConstraint<RealT>::create(pl);
  auto  height = pl.get("Height of Tank",              10.0  );
  auto  Qin00  = pl.get("Corner Inflow",               100.0 );
  auto  h_init = pl.get("Initial Fluid Level",         2.0   );
  auto  nrows  = static_cast<size_type>( pl.get("Number of Rows"   ,3) );
  auto  ncols  = static_cast<size_type>( pl.get("Number of Columns",3) );
  auto  Nt     = static_cast<size_type>( pl.get("Number of Time Steps",100) );
  auto  T       = pl.get("Total Time", 20.0);

  RealT dt = T/Nt;

  ROL::Ptr<ROL::Vector<RealT>> initial_cond;
  {
    // control
    auto  z      = Control::create( pl, "Control (z)"     );    
    auto  vz     = z->clone( "Control direction (vz)"     );
    auto  z_lo   = z->clone( "Control Lower Bound (z_lo)" );
    auto  z_bnd  = makePtr<Bounds>( *z_lo );
    z_lo->zero();

    // State
    auto u_new     = State::create( pl, "New state (u_new)"   );
    auto u_old     = u_new->clone( "Old state (u_old)"        );
    auto u_initial = u_new->clone( "Initial conditions"       );
    auto u_new_lo  = u_new->clone( "State lower bound (u_lo)" );
    auto u_new_up  = u_new->clone( "State upper bound (u_up)" );
    auto u         = PartitionedVector::create( { u_old,    u_new    } );
    auto u_lo      = PartitionedVector::create( { u_new_lo, u_new_lo } );
    auto u_up      = PartitionedVector::create( { u_new_up, u_new_up } );
  
    u_lo->zero();
    u_up->setScalar( height );
    auto u_bnd = makePtr<Bounds>(u_new_lo,u_new_up);

    (*z)(0,0) = Qin00;

    for( size_type i=0; i<nrows; ++i ) {
      for( size_type j=0; j<ncols; ++j ) {
        u_old->h(i,j) = h_init;
        u_initial->h(i,j) = h_init;
      }
    }

    state        = ROL::buildStatePinTVector<RealT>(   communicators, Nt,     u_old);
    control      = ROL::buildControlPinTVector<RealT>( communicators, Nt,         z);

    initial_cond = u_initial;
    state->getVectorPtr(-1)->set(*u_initial);   // set the initial condition
  }

  auto timeStamp = ROL::makePtr<std::vector<ROL::TimeStamp<RealT>>>(state->numOwnedSteps());
  for( size_type k=0; k<timeStamp->size(); ++k ) {
    timeStamp->at(k).t.resize(2);
    timeStamp->at(k).t.at(0) = k*dt;
    timeStamp->at(k).t.at(1) = (k+1)*dt;
  }

  ROL::PinTConstraint<RealT> pint_constraint(con,initial_cond,timeStamp);

  double tol = 1e-10;

  // check the pint constraint
  {
    auto c   = state->clone();
    auto jv  = c->clone();
    auto v_u = state->clone();
    ROL::RandomizeVector<RealT>(*state);
    ROL::RandomizeVector<RealT>(*control);
    ROL::RandomizeVector<RealT>(*v_u);

    pint_constraint.setGlobalScale(1.0);
    pint_constraint.checkSolve(*state,*control,*c,true,*outStream);
    pint_constraint.checkApplyJacobian_1(*state,*control,*v_u,*jv,true,*outStream);
    RealT inv_1     = pint_constraint.checkInverseJacobian_1(*jv,*v_u,*state,*control,true,*outStream);
    RealT adj_inv_1 = pint_constraint.checkInverseAdjointJacobian_1(*jv,*v_u,*state,*control,true,*outStream);

    if(inv_1 > tol) {
      std::stringstream ss;
      ss << "Forward Jacobian inverse inversion FAILED: error = " << inv_1  << std::endl;
      throw std::logic_error(ss.str());
    }

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
