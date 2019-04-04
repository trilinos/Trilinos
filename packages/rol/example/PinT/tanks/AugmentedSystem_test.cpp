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

#include "ROL_Ptr.hpp"

#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

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
void run_test_kkt(MPI_Comm comm, const ROL::Ptr<std::ostream> & outStream);

int main( int argc, char* argv[] ) 
{
  using ROL::Ptr;
  using ROL::makePtr;
  using ROL::makePtrFromRef;

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);  
  int myRank = Teuchos::GlobalMPISession::getRank();

  auto outStream = ROL::makeStreamPtr( std::cout, argc > 1 );

  int errorFlag  = 0;

  try {

    run_test(MPI_COMM_WORLD, outStream);
    run_test_kkt(MPI_COMM_WORLD, outStream);

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
  auto  T      = pl.get("Total Time", 20.0);

  auto  maxlvs = pl.get("MGRIT Levels", 5);
  auto  sweeps = pl.get("MGRIT Sweeps", 1);
  auto  omega  = pl.get("MGRIT Relax", 2.0/3.0);

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


  // build the parallel in time constraint from the user constraint
  Ptr<ROL::PinTConstraint<RealT>> pint_con = makePtr<ROL::PinTConstraint<RealT>>(con,initial_cond,timeStamp);
  pint_con->applyMultigrid(maxlvs);
  pint_con->setSweeps(sweeps);
  pint_con->setRelaxation(omega);

  // check the pint constraint
  {
    auto c   = state->clone();
    auto jv  = c->clone();
    auto v_u = state->clone();
    ROL::RandomizeVector<RealT>(*state);
    ROL::RandomizeVector<RealT>(*control);
    ROL::RandomizeVector<RealT>(*v_u);

    pint_con->checkSolve(*state,*control,*c,true,*outStream);
    pint_con->checkApplyJacobian_1(*state,*control,*v_u,*jv,true,*outStream);
  }

  auto x   = makePtr<ROL::Vector_SimOpt<RealT>>(state,control);
  auto v_1 = x->clone();
  auto v_2 = state->clone();
  auto r_1 = v_1->clone();
  auto r_2 = v_2->clone();
  ROL::RandomizeVector<RealT>(*r_1);
  ROL::RandomizeVector<RealT>(*r_2);

  int numSolves = 1;
  std::clock_t timer_total = 0;
  for(int i=0;i<numSolves;i++) {
    v_1->zero();
    v_2->zero();

    double tol = 1e-6;

    auto b1 = r_1->clone();
    b1->set(*r_1);
    auto b2 = r_2->clone();
    b2->set(*r_2);

    std::clock_t timer_start = std::clock();
    auto res = pint_con->solveAugmentedSystem(*v_1,*v_2,*b1,*b2,*x,tol);
    timer_total += std::clock()-timer_start;

    if(myRank==0) {
      for(size_t k=0;k<res.size();k++)
        *outStream << "  " << k << ". " << res[k] << " (rel = " << res[k]/res[0] << ")" << std::endl;
    }
  }
 
  if(myRank==0)
    *outStream << numSolves << " Solves in " << (timer_total)/(RealT) CLOCKS_PER_SEC << " seconds." << std::endl;

  for(int chkProc=0;chkProc<numRanks;++chkProc) {
    using ROL::PinTVector;

    auto v  = state->clone();
    auto fv  = state->clone();
    auto afv  = state->clone();
    auto jv = state->clone();
    auto ajv = state->clone();

    ROL::RandomizeVector<RealT>(*v);

    double tol = 1e-12;

    // zero all contribution off the processor being checked
    if(myRank!=chkProc)
      v->scale(0.0);

    // make sure we are globally consistent
    state->boundaryExchange();
    control->boundaryExchange();
    dynamic_cast<PinTVector<RealT>&>(*v).boundaryExchange();

    // compute jacobian action
    pint_con->applyJacobian_1(*jv,
                              *v,
                              *state,
                              *control,tol);

    // compute jacobian action
    pint_con->applyAdjointJacobian_1(*ajv,
                                     *v,
                                     *state,
                                     *control,tol);

    // zero all contribution off the processor being checked
    if(myRank!=chkProc) {
      jv->zero();
      ajv->zero();
    }

    int level = 0;
    pint_con->invertTimeStepJacobian(dynamic_cast<PinTVector<RealT>&>(*fv),
                                     dynamic_cast<const PinTVector<RealT>&>(*jv),
                                     dynamic_cast<const PinTVector<RealT>&>(*state),
                                     dynamic_cast<const PinTVector<RealT>&>(*control),
                                     tol,
                                     level);

    pint_con->invertAdjointTimeStepJacobian(dynamic_cast<PinTVector<RealT>&>(*afv),
                                     dynamic_cast<const PinTVector<RealT>&>(*ajv),
                                     dynamic_cast<const PinTVector<RealT>&>(*state),
                                     dynamic_cast<const PinTVector<RealT>&>(*control),
                                     tol,
                                     level);

    fv->axpy(-1.0,*v);
    afv->axpy(-1.0,*v);

    RealT fv_norm  = fv->norm();
    RealT afv_norm = afv->norm();
    RealT v_norm   = v->norm();

    if(myRank==0) {
      *outStream << "Testing vector norm = " << v_norm << std::endl;
    }

    // check norms
    if(fv_norm/v_norm > 1e-13) {
      std::stringstream ss;
      ss << "Forward block Jacobi inversion FAILED with proc " << chkProc 
         << " (relative error = " << fv_norm / v_norm  << ")" << std::endl;
      throw std::logic_error(ss.str());
    }
    else {
       *outStream << "Forward block Jacobi inversion PASSED with proc " << chkProc 
                  << " (relative error = " << fv_norm / v_norm  << ")" << std::endl;
    }

    // check norms (adjoint)
    if(afv_norm/v_norm > 1e-13) {
      std::stringstream ss;
      *outStream << "Adjoint block Jacobi inversion FAILED with proc " << chkProc 
                 << " (relative error = " << afv_norm / v_norm  << ")" << std::endl;
      throw std::logic_error(ss.str());
    }
    else {
      *outStream << "Adjoint block Jacobi inversion PASSED with proc " << chkProc 
                 << " (relative error = " << afv_norm / v_norm  << ")" << std::endl;
    }
  }
}

void run_test_kkt(MPI_Comm comm, const ROL::Ptr<std::ostream> & outStream)
{
  using ROL::Ptr;
  using ROL::makePtr;
  using ROL::makePtrFromRef;

  using RealT             = double;
  using size_type         = std::vector<RealT>::size_type;
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
  auto  T      = pl.get("Total Time", 20.0);

  auto  maxlvs = pl.get("MGRIT Levels",           5);
  auto  sweeps = pl.get("MGRIT Sweeps",           3);
  auto  omega  = pl.get("MGRIT Relaxation", 2.0/3.0);
  auto  globalScale  = pl.get("Constraint Scale", 1.0e3);

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


  // build the parallel in time constraint from the user constraint
  Ptr<ROL::PinTConstraint<RealT>> pint_con = makePtr<ROL::PinTConstraint<RealT>>(con,initial_cond,timeStamp);
  pint_con->applyMultigrid(maxlvs);
  pint_con->setSweeps(sweeps);
  pint_con->setRelaxation(omega);
  pint_con->setGlobalScale(globalScale);

  double tol = 1e-12;

  // make sure we are globally consistent
  state->boundaryExchange();
  control->boundaryExchange();

  Ptr<ROL::Vector<RealT>> kkt_vector = makePtr<PartitionedVector>({state->clone(),control->clone(),state->clone()});
  auto kkt_x_in  = kkt_vector->clone();
  auto kkt_b     = kkt_vector->clone();

  ROL::RandomizeVector(*kkt_x_in);
  kkt_b->zero();
  pint_con->applyAugmentedKKT(*kkt_b,*kkt_x_in,*state,*control,tol);

  // check exact
  {
    auto kkt_x_out = kkt_vector->clone();
    kkt_x_out->zero();

    // apply_invkkt(*kkt_x_out,*kkt_b,*state,*control);
    pint_con->applyAugmentedInverseKKT(*kkt_x_out,*kkt_b,*state,*control,tol);

    kkt_x_out->axpy(-1.0,*kkt_x_in);

    RealT norm = kkt_x_out->norm();
    if(myRank==0)
      (*outStream) << "NORM EXACT= " << norm << std::endl;
  }

  if(myRank==0)
    (*outStream) << std::endl;

  // check jacobi
  {
    auto kkt_x_out = kkt_vector->clone();
    auto kkt_diff = kkt_vector->clone();
    auto kkt_err = kkt_vector->clone();
    auto kkt_res = kkt_vector->clone();
    kkt_x_out->zero();

    for(int i=0;i<sweeps;i++) {
      pint_con->applyAugmentedKKT(*kkt_res,*kkt_x_out,*state,*control,tol);
      kkt_res->scale(-1.0);
      kkt_res->axpy(1.0,*kkt_b);
    
      pint_con->applyAugmentedInverseKKT(*kkt_diff,*kkt_res,*state,*control,tol,true);

      kkt_x_out->axpy(omega,*kkt_diff);

      kkt_err->set(*kkt_x_out);
      kkt_err->axpy(-1.0,*kkt_x_in);
      RealT norm = kkt_err->norm() / kkt_b->norm();
      if(myRank==0)
        (*outStream) << "NORM JACOBI= " << norm << std::endl;
    }
  }

  if(myRank==0)
    (*outStream) << std::endl;

  // check MG
  MPI_Barrier(MPI_COMM_WORLD);
  double t0 = MPI_Wtime();
  {
    auto kkt_x_out = kkt_vector->clone();
    auto kkt_diff = kkt_vector->clone();
    auto kkt_err = kkt_vector->clone();
    auto kkt_res = kkt_vector->clone();
    kkt_x_out->zero();

    pint_con->applyAugmentedKKT(*kkt_res,*kkt_x_out,*state,*control,tol);
    kkt_res->scale(-1.0);
    kkt_res->axpy(1.0,*kkt_b);

    RealT res0 = kkt_res->norm();

    if(myRank==0)
      (*outStream) << "Multigrid initial res = " << res0 << std::endl;

    for(int i=0;i<100;i++) {
      pint_con->applyMultigridAugmentedKKT(*kkt_diff,*kkt_res,*state,*control,tol);

      kkt_x_out->axpy(omega,*kkt_diff);

      kkt_err->set(*kkt_x_out);
      kkt_err->axpy(-1.0,*kkt_x_in);
   
      pint_con->applyAugmentedKKT(*kkt_res,*kkt_x_out,*state,*control,tol);
      kkt_res->scale(-1.0);
      kkt_res->axpy(1.0,*kkt_b);

      RealT res = kkt_res->norm();
      if(myRank==0)
        (*outStream) << " " << i+1 << ". " << res/res0 << std::endl;

      // check the residual
      if(res/res0 < 1e-9) {
        // how many iterations
        if(myRank==0)
          (*outStream) << "\nIterations = " << i+1 << std::endl;
        
        break;
      }
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);
  double tf = MPI_Wtime();

  if(myRank==0)
    (*outStream) << "\nMG Time = " << tf-t0 << std::endl;
}


