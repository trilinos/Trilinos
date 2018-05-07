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

    std::string tank_xml("tank-parameters.xml");
    Teuchos::updateParametersFromXmlFile(tank_xml, tank_parameters.ptr());
    auto& pl = *tank_parameters;

    auto nrows    = static_cast<size_type>( pl.get("Number of Rows",3) );
    auto ncols    = static_cast<size_type>( pl.get("Number of Columns",3) );
    auto numSteps = static_cast<size_type>( pl.get("Number of Time Steps",100) );

    auto u    = makePtr<StateVector>( nrows, ncols, "New state (u)" );
    auto z    = makePtr<ControlVector>( nrows, ncols, "Control (z)" );    

    auto tankState = makePtr<TankState<RealT>>(pl);
    auto con = makePtr<TankConstraint<RealT>>(tankState,pl);

    auto state   = ROL::buildPinTVector<RealT>(communicators,numSteps,{-1,0}, u);
    auto control = ROL::buildPinTVector<RealT>(communicators,numSteps,   {0}, z);


    // build the parallel in time constraint from the user constraint
    Ptr<ROL::Constraint_PinTSimOpt<RealT>> pint_con = makePtr<ROL::Constraint_PinTSimOpt<RealT>>(con);

    auto x   = makePtr<ROL::Vector_SimOpt<RealT>>(state,control);
    auto v_1 = x->clone();
    auto v_2 = state->clone();
    auto r_1 = v_1->clone();
    auto r_2 = v_2->clone();

    ROL::RandomizeVector<RealT>(*r_1);
    ROL::RandomizeVector<RealT>(*r_2);

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

    int numSolves = 5;
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

      if(myrank==0)
        *outStream << i+1 << ". res ("<< res.size() <<") = " << res[res.size()-3] << std::endl;
    }
 
    if(myrank==0)
      *outStream << numSolves << " Solves in " << (timer_total)/(RealT) CLOCKS_PER_SEC << " seconds." << std::endl;
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
