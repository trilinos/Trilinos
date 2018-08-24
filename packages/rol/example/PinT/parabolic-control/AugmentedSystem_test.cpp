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
#include "Teuchos_StackedTimer.hpp"

#include "ROL_Stream.hpp"
#include "ROL_ParameterList.hpp"
#include "ROL_Bounds.hpp"
#include "ROL_RandomVector.hpp"
#include "ROL_Vector_SimOpt.hpp"
#include "ROL_PinTConstraint.hpp"
#include "ROL_GMRES.hpp"

#include "dynamicConstraint.hpp"
#include "dynamicObjective.hpp"

using RealT = double;
using size_type = std::vector<RealT>::size_type;

void run_test(MPI_Comm comm, const ROL::Ptr<std::ostream> & outStream);
void run_test_simple(MPI_Comm comm, const ROL::Ptr<std::ostream> & outStream);
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

    // run_test(MPI_COMM_WORLD, outStream);
    // run_test_simple(MPI_COMM_WORLD, outStream);
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

class KKTOperator : public ROL::LinearOperator<RealT> {
public:
  ROL::Ptr<ROL::PinTConstraint<RealT>> pint_con;
  ROL::Ptr<ROL::Vector<RealT>> state;
  ROL::Ptr<ROL::Vector<RealT>> control;

  void apply( ROL::Vector<RealT> &Hv, const ROL::Vector<RealT> &v, RealT &tol ) const override
  {
    pint_con->applyAugmentedKKT(Hv,v,*state,*control,tol);
  }
};

class MGRITKKTOperator : public ROL::LinearOperator<RealT> {
public:
  ROL::Ptr<ROL::PinTConstraint<RealT>> pint_con;
  ROL::Ptr<ROL::Vector<RealT>> state;
  ROL::Ptr<ROL::Vector<RealT>> control;

  void apply( ROL::Vector<RealT> &Hv, const ROL::Vector<RealT> &v, RealT &tol ) const override
  {
    assert(false);
    pint_con->applyMultigridAugmentedKKT(Hv,v,*state,*control,tol);
  }
  void applyInverse( ROL::Vector<RealT> &Hv, const ROL::Vector<RealT> &v, RealT &tol ) const override
  {
    pint_con->applyMultigridAugmentedKKT(Hv,v,*state,*control,tol);
  }
};

void run_test_kkt(MPI_Comm comm, const ROL::Ptr<std::ostream> & outStream)
{
  using ROL::Ptr;
  using ROL::makePtr;
  using ROL::makePtrFromRef;

  using RealT             = double;
  using size_type         = std::vector<RealT>::size_type;
  using Bounds            = ROL::Bounds<RealT>;
  using PartitionedVector = ROL::PartitionedVector<RealT>;

  auto timer = Teuchos::TimeMonitor::getStackedTimer();

  int numRanks = -1;
  int myRank = -1;

  MPI_Comm_size(comm, &numRanks);
  MPI_Comm_rank(comm, &myRank);

  *outStream << "Proc " << myRank << "/" << numRanks << std::endl;

  ROL::Ptr<const ROL::PinTCommunicators> communicators = ROL::makePtr<ROL::PinTCommunicators>(comm,1);

  // Parse input parameter list
  ROL::Ptr<ROL::ParameterList> pl = ROL::getParametersFromXmlFile("input_ex01.xml");
  bool derivCheck = pl->get("Derivative Check",         true); // Check derivatives.
  uint nx         = pl->get("Spatial Discretization",     64); // Set spatial discretization.
  uint nt         = pl->get("Temporal Discretization",   100); // Set temporal discretization.
  RealT T         = pl->get("End Time",                  1.0); // Set end time.
  RealT dt        = T/(static_cast<RealT>(nt)-1.0);

  // Initialize objective function.
  ROL::Ptr<ROL::DynamicObjective<RealT>> dyn_obj
    = ROL::makePtr<Objective_ParabolicControl<RealT>>(*pl);
  ROL::Ptr<ROL::DynamicConstraint<RealT>> dyn_con
    = ROL::makePtr<Constraint_ParabolicControl<RealT>>(*pl);

  // Create control vectors.
  ROL::Ptr<ROL::StdVector<RealT>>         zk = ROL::makePtr<ROL::StdVector<RealT>>(ROL::makePtr<std::vector<RealT>>(nx+2));

  // Create initial state vector.
  ROL::Ptr<ROL::StdVector<RealT>> u0 = ROL::makePtr<ROL::StdVector<RealT>>(ROL::makePtr<std::vector<RealT>>(nx,0.0));

  ROL::Ptr< ROL::PinTVector<RealT>> state;
  ROL::Ptr< ROL::PinTVector<RealT>> control;

  state        = ROL::buildStatePinTVector<RealT>(   communicators, nt,     u0); // for Euler, Crank-Nicolson, stencil = [-1,0]
  control      = ROL::buildControlPinTVector<RealT>( communicators, nt,     zk); // time discontinous, stencil = [0]

  // Construct reduced dynamic objective
  auto timeStamp = ROL::makePtr<std::vector<ROL::TimeStamp<RealT>>>(state->numOwnedSteps());
  for( uint k=0; k<timeStamp->size(); ++k ) {
    timeStamp->at(k).t.resize(2);
    timeStamp->at(k).t.at(0) = k*dt;
    timeStamp->at(k).t.at(1) = (k+1)*dt;
  }

  // Add MGRIT parameter list stuff
  int sweeps    = pl->get("MGRIT Sweeps", 1); 
  RealT omega   = pl->get("MGRIT Relaxation",2.0/3.0);
  int numLevels = pl->get("MGRIT Levels",3);
  double relTol = pl->get("MGRIT Krylov Relative Tolerance",1e-4);
  
  if(myRank==0) {
    (*outStream) << "Sweeps = " << sweeps    << std::endl;
    (*outStream) << "Omega = "  << omega     << std::endl;
    (*outStream) << "Levels = " << numLevels << std::endl;
  }

  // build the parallel in time constraint from the user constraint
  Ptr<ROL::PinTConstraint<RealT>> pint_con = makePtr<ROL::PinTConstraint<RealT>>(dyn_con,u0,timeStamp);
  // pint_con->setGlobalScale(1.0);
  pint_con->applyMultigrid(numLevels);
  // pint_con->applyMultigrid(2);
  pint_con->setSweeps(sweeps);
  pint_con->setRelaxation(omega);

  double tol = 1e-12;

  // make sure we are globally consistent
  state->boundaryExchange();
  control->boundaryExchange();

  Ptr<ROL::Vector<RealT>> kkt_vector = makePtr<PartitionedVector>({state->clone(),control->clone(),state->clone()});
  auto kkt_x_in  = kkt_vector->clone();
  auto kkt_b     = kkt_vector->clone();

  ROL::RandomizeVector(*kkt_x_in);
  kkt_b->zero();
  // apply_kkt(*kkt_b,*kkt_x_in,*state,*control);
  pint_con->applyAugmentedKKT(*kkt_b,*kkt_x_in,*state,*control,tol); // b = A * x

  // check exact: Update so this works right (TODO)
  {
    auto kkt_x_out = kkt_vector->clone();
    kkt_x_out->zero();

    pint_con->applyAugmentedInverseKKT(*kkt_x_out,*kkt_b,*state,*control,tol); // apply Wathen style preconditioner serially

    kkt_x_out->axpy(-1.0,*kkt_x_in);

    RealT norm = kkt_x_out->norm();
    if(myRank==0)
      (*outStream) << "NORM EXACT= " << norm << std::endl;
  }

  if(myRank==0)
    (*outStream) << std::endl;

  // check jacobi is reducing the norm 
  {
    auto kkt_x_out = kkt_vector->clone();
    auto kkt_diff = kkt_vector->clone();
    auto kkt_err = kkt_vector->clone();
    auto kkt_res = kkt_vector->clone();
    kkt_x_out->zero();

    for(int i=0;i<sweeps;i++) {
      pint_con->applyAugmentedKKT(*kkt_res,*kkt_x_out,*state,*control,tol);
      kkt_res->scale(-1.0);
      kkt_res->axpy(1.0,*kkt_b);                                                       // r = b - A*x_i
    
      pint_con->applyAugmentedInverseKKT(*kkt_diff,*kkt_res,*state,*control,tol,true); // dx_i = M^{-1} r

      kkt_x_out->axpy(omega,*kkt_diff);                                                // x_{i+1} = x_i + omega*dx_i

      kkt_err->set(*kkt_x_out);
      kkt_err->axpy(-1.0,*kkt_x_in);                                                   // err = x_{i+1} - x
      RealT norm = kkt_err->norm() / kkt_b->norm();
      if(myRank==0)
        (*outStream) << "NORM JACOBI= " << norm << std::endl;
    }
  }

  if(myRank==0)
    (*outStream) << std::endl;

  {
    RealT res0 = kkt_b->norm();
    kkt_b->scale(1.0/res0);
    res0 = kkt_b->norm();

    assert(std::fabs(res0-1.0) <= 1.0e-14);

    auto kkt_x_out = kkt_vector->clone();
    kkt_x_out->zero();

    KKTOperator kktOperator;
    kktOperator.pint_con = pint_con;
    kktOperator.state = state;
    kktOperator.control = control;

    MGRITKKTOperator mgOperator;
    mgOperator.pint_con = pint_con;
    mgOperator.state = state;
    mgOperator.control = control;

    Teuchos::ParameterList parlist;
    parlist.set("Absolute Tolerance",1.e-1);
    parlist.set("Relative Tolerance",relTol);
 
    ROL::GMRES<RealT> krylov(parlist); // TODO: Do Belos

    int flag = 0, iter = 0;

    MPI_Barrier(MPI_COMM_WORLD);
    double t0 = MPI_Wtime();

    pint_con->clearTimePerLevel();

    timer->start("krylov");
    RealT finalTol = krylov.run(*kkt_x_out,kktOperator,*kkt_b,mgOperator,iter,flag);
    timer->stop("krylov");

    MPI_Barrier(MPI_COMM_WORLD);
    double tf = MPI_Wtime();

    std::stringstream ss;
    timer->report(ss);

    if(myRank==0) {
      (*outStream) << "Krylov Iteration = " << iter << " " << (finalTol / res0) << " " << tf-t0 << std::endl;

      const std::vector<RealT> & timePerLevel = pint_con->getTimePerLevel();
      (*outStream) << "Time per level = ";
      for(RealT time : timePerLevel) 
        (*outStream) << time << " ";
      (*outStream) << std::endl;

      (*outStream) << std::endl;
      (*outStream) << ss.str();
    }
  }
}


