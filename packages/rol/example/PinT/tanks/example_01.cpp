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

#include "Teuchos_GlobalMPISession.hpp"

#include "ROL_Stream.hpp"
#include "ROL_ParameterList.hpp"
#include "ROL_Bounds.hpp"
#include "ROL_RandomVector.hpp"
//#include "ROL_PinTVector.hpp"

#include "TankConstraint.hpp"
#include "TankVector.hpp"
#include "LowerBandedMatrix.hpp"

#include <iostream>


using RealT = double;
using size_type = std::vector<RealT>::size_type;


int main( int argc, char* argv[] ) {

  using ROL::Ptr;
  using ROL::makePtr;
  using ROL::makePtrFromRef;

  using Bounds                      = ROL::Bounds<RealT>;
  using PartitionedVector           = ROL::PartitionedVector<RealT>;
 
  using StateVector   = TankStateVector<RealT>;
  using ControlVector = TankControlVector<RealT>;

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);  

  int iprint     = argc - 1;
  Ptr<std::ostream> outStream;
  ROL::nullstream bhs; // outputs nothing
  if (iprint > 0)
    outStream = makePtrFromRef(std::cout);
  else
    outStream = makePtrFromRef(bhs);

  int errorFlag  = 0;

  // *** Example body.

  try {   

    std::string tank_xml("tank-parameters.xml");

    auto tank_parameters = ROL::getParametersFromXmlFile(tank_xml);
    auto& pl = *tank_parameters;

    auto tankState = makePtr<TankState<RealT>>(pl);
    auto con = makePtr<TankConstraint<RealT>>(tankState,pl);

    auto height = pl.get("Height of Tank",              10.0  );
    auto Qin00  = pl.get("Corner Inflow",               100.0 );
    auto h_init = pl.get("Initial Fluid Level",         2.0   );
//    auto T      = pl.get("Total Time",                  20.0  );
//    auto nt     = pl.get("Number of Time Steps",        100   );

    auto nrows  = static_cast<size_type>( pl.get("Number of Rows",3) );
    auto ncols  = static_cast<size_type>( pl.get("Number of Columns",3) );

    auto z      = makePtr<ControlVector>( nrows, ncols, "Control (z)" );    
    auto vz     = z->clone( "Control direction (vz)"       );
    auto z_lo   = z->clone( "Control Lower Bound (z_lo)"   );
    z_lo->zero();

    auto z_bnd  = makePtr<Bounds>( *z_lo );

    // State
    auto u_new    = makePtr<StateVector>( nrows, ncols, "New state (u_new)" );
    auto u_old    = u_new->clone( "Old state (u_old)"        );
    auto u_new_lo = u_new->clone( "State lower bound (u_lo)" );
    auto u_new_up = u_new->clone( "State upper bound (u_up)" );

    auto u    = PartitionedVector::create( { u_old, u_new } );
    auto u_lo = PartitionedVector::create( { u_new_lo, u_new_lo } );
    auto u_up = PartitionedVector::create( { u_new_up, u_new_up } );

    u_lo->zero();
    u_up->setScalar( height );
    
    // State direction 
    auto vu_new = u_new->clone( "New state direction (vu_new)" );
    auto vu_old = u_old->clone( "Old state direction (vu_old)" );
    auto vu     = PartitionedVector::create( { vu_old, vu_new } );

    auto c = u_new->clone( "State residual (c)"              );
    auto x = ROL::makePtr<ROL::Vector_SimOpt<RealT>>( u,  z  );
    auto v = ROL::makePtr<ROL::Vector_SimOpt<RealT>>( vu, vz );

    (*z)(0,0) = Qin00;

    for( size_type i=0; i<nrows; ++i )
      for( size_type j=0; j<ncols; ++j )
        u_old->h(i,j) = h_init;

    auto u_new_bnd = makePtr<Bounds>( u_new_lo, u_new_up );
    RandomizeFeasibleVector( *u_new, *u_new_bnd );
    RandomizeVector( *c ) ;
    RandomizeVector( *v );

    auto u_old_bnd = u_new_bnd; 
    auto u_bnd = makePtr<Bounds>( u_lo, u_up ); 

    con->print_tankstate_parameters( *outStream );
    con->checkApplyJacobian( *x, *v, *c, true, *outStream );
    con->checkAdjointConsistencyJacobian( *c, *v, *x, true, *outStream );
    con->checkInverseJacobian_1_new( *c, *u_new, *u_old, *z, *vu_new, true, *outStream );
    con->checkInverseAdjointJacobian_1_new( *c, *u_new, *u_old, *z, *vu_new, true, *outStream );

    RandomizeVector( *c ) ;
    con->checkSolve(*u, *z, *c, true, *outStream ); 

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
