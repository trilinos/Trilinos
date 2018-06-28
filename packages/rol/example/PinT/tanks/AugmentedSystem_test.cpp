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

#include "TankConstraint.hpp"

#include "ROL_Constraint_PinTSimOpt.hpp"

using RealT = double;
using size_type = std::vector<RealT>::size_type;

int main( int argc, char* argv[] ) 
{
  using ROL::Ptr;
  using ROL::makePtr;
  using ROL::makePtrFromRef;

  using StateVector   = TankStateVector<RealT>;
  using ControlVector = TankControlVector<RealT>;

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);  
  int myrank = Teuchos::GlobalMPISession::getRank();
  int numProcs = Teuchos::GlobalMPISession::getNProc();

  int iprint     = argc - 1;
  ROL::Ptr<std::ostream> outStream;
  Teuchos::oblackholestream bhs; // outputs nothing
  if (iprint > 0)
    outStream = ROL::makePtrFromRef(std::cout);
  else
    outStream = ROL::makePtrFromRef(bhs);

  int errorFlag  = 0;

  try {

    int spatial_procs = 1;
    ROL::Ptr<const ROL::PinTCommunicators> communicators = ROL::makePtr<ROL::PinTCommunicators>(MPI_COMM_WORLD,spatial_procs);

    auto tank_parameters = makePtr<Teuchos::ParameterList>();

    auto  pl_ptr  = ROL::getParametersFromXmlFile("tank-parameters.xml");
    auto& pl      = *pl_ptr;

    auto nrows    = static_cast<size_type>( pl.get("Number of Rows",3) );
    auto ncols    = static_cast<size_type>( pl.get("Number of Columns",3) );
    auto numSteps = static_cast<size_type>( pl.get("Number of Time Steps",100) );

    if(myrank==0)
      *outStream << "Rows = " << nrows << ", Cols = " << ncols << ", Time steps = " << numSteps << std::endl;

    auto u    = makePtr<StateVector>( nrows, ncols, "New state (u)" );
    auto z    = makePtr<ControlVector>( nrows, ncols, "Control (z)" );    

    auto tankState = makePtr<TankState<RealT>>(pl);
    auto con = makePtr<TankConstraint<RealT>>(tankState,pl);

    auto state   = ROL::buildPinTVector<RealT>(communicators,numSteps,{-1,0}, u);
    auto control = ROL::buildPinTVector<RealT>(communicators,numSteps,   {0}, z);

    // build the parallel in time constraint from the user constraint
    Ptr<ROL::Constraint_PinTSimOpt<RealT>> pint_con = makePtr<ROL::Constraint_PinTSimOpt<RealT>>(con);

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
    double tol = 1e-12;
    std::clock_t timer_total = 0;
    for(int i=0;i<numSolves;i++) {
      v_1->zero();
      v_2->zero();

      auto b1 = r_1->clone();
      b1->set(*r_1);
      auto b2 = r_2->clone();
      b2->set(*r_2);

      std::clock_t timer_start = std::clock();
      auto res = pint_con->solveAugmentedSystem(*v_1,*v_2,*b1,*b2,*x,tol);
      timer_total += std::clock()-timer_start;

      if(myrank==0) {
        for(size_t k=0;k<res.size();k++)
          *outStream << "  " << k << ". " << res[k] << " (rel = " << res[k]/res[0] << ")" << std::endl;
      }
    }
 
    if(myrank==0)
      *outStream << numSolves << " Solves in " << (timer_total)/(RealT) CLOCKS_PER_SEC << " seconds." << std::endl;

    for(int chkProc=0;chkProc<numProcs;++chkProc) {
      using ROL::PinTVector;

      auto v  = state->clone();
      auto fv  = state->clone();
      auto afv  = state->clone();
      auto jv = state->clone();
      auto ajv = state->clone();

      ROL::RandomizeVector<RealT>(*v);

      double tol = 1e-12;

      // zero all contribution off the processor being checked
      if(myrank!=chkProc)
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
      if(myrank!=chkProc) {
        jv->zero();
        ajv->zero();
      }

      pint_con->invertTimeStepJacobian(dynamic_cast<PinTVector<RealT>&>(*fv),
                                       dynamic_cast<const PinTVector<RealT>&>(*jv),
                                       dynamic_cast<const PinTVector<RealT>&>(*state),
                                       dynamic_cast<const PinTVector<RealT>&>(*control),
                                       tol);

      pint_con->invertAdjointTimeStepJacobian(dynamic_cast<PinTVector<RealT>&>(*afv),
                                       dynamic_cast<const PinTVector<RealT>&>(*ajv),
                                       dynamic_cast<const PinTVector<RealT>&>(*state),
                                       dynamic_cast<const PinTVector<RealT>&>(*control),
                                       tol);

      fv->axpy(-1.0,*v);
      afv->axpy(-1.0,*v);

      RealT fv_norm  = fv->norm();
      RealT afv_norm = afv->norm();
      RealT v_norm   = v->norm();

      if(myrank==0) {
        *outStream << "Testing vector norm = " << v_norm << std::endl;
      }

      // check norms
      if(fv_norm/v_norm > 1e-13) {
        errorFlag = 1000;
        if(myrank==0)
          *outStream << "Forward block Jacobi inversion FAILED with proc " << chkProc 
                     << " (relative error = " << fv_norm / v_norm  << ")" << std::endl;
      }
      else {
        if(myrank==0)
          *outStream << "Forward block Jacobi inversion PASSED with proc " << chkProc 
                     << " (relative error = " << fv_norm / v_norm  << ")" << std::endl;
      }

      // check norms (adjoint)
      if(afv_norm/v_norm > 1e-13) {
        errorFlag = 1000;
        if(myrank==0)
          *outStream << "Adjoint block Jacobi inversion FAILED with proc " << chkProc 
                     << " (relative error = " << afv_norm / v_norm  << ")" << std::endl;
      }
      else {
        if(myrank==0)
          *outStream << "Adjoint block Jacobi inversion PASSED with proc " << chkProc 
                     << " (relative error = " << afv_norm / v_norm  << ")" << std::endl;
      }
    }
  }
  catch (std::logic_error err) {
    *outStream << err.what() << "\n";
    errorFlag = -1000;
  }; // end try

  if(myrank==0) {
    if (errorFlag != 0)
      std::cout << "End Result: TEST FAILED\n";
    else
      std::cout << "End Result: TEST PASSED\n";
  }

  return 0;
}
