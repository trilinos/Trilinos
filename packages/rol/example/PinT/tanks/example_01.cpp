// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_GlobalMPISession.hpp"
#include "ROL_Stream.hpp"
#include "ROL_ParameterList.hpp"
#include "ROL_SerialConstraint.hpp"
#include "ROL_Bounds.hpp"
#include "ROL_RandomVector.hpp"
#include "ROL_Vector_SimOpt.hpp"
#include "Tanks_DynamicConstraint.hpp"
#include "Tanks_SerialConstraint.hpp"
#include "Tanks_ConstraintCheck.hpp"

#include <iostream>

int main( int argc, char* argv[] ) {
 
  using ROL::Ptr;
  using ROL::makePtr;

  using RealT             = double;
  using size_type         = std::vector<RealT>::size_type;
  using ValidateFunction  = ROL::ValidateFunction<RealT>;
  using Bounds            = ROL::Bounds<RealT>;
  using PartitionedVector = ROL::PartitionedVector<RealT>;
  using TimeStamp         = ROL::TimeStamp<RealT>;
  using State             = Tanks::StateVector<RealT>;
  using Control           = Tanks::ControlVector<RealT>;

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);  

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  auto outStream = ROL::makeStreamPtr( std::cout, argc > 1 );
  int errorFlag  = 0;

  try {     // *** Example body.

    auto  pl_ptr  = ROL::getParametersFromXmlFile("tank-parameters.xml");
    auto& pl      = *pl_ptr;
    auto  dyn_con = Tanks::DynamicConstraint<RealT>::create(pl);
    auto  height  = pl.get("Height of Tank",              10.0  );
    auto  Qin00   = pl.get("Corner Inflow",               100.0 );
    auto  h_init  = pl.get("Initial Fluid Level",         2.0   );
    auto  nrows   = static_cast<size_type>( pl.get("Number of Rows"   ,3) );
    auto  ncols   = static_cast<size_type>( pl.get("Number of Columns",3) );
    auto  Nt      = static_cast<size_type>( pl.get("Number of Time Stamps",100) );
    auto  T       = pl.get("Total Time", 20.0);

    RealT dt = T/Nt;

    auto  z       = Control::create( pl, "Control (z)"     );    
    auto  vz      = z->clone( "Control direction (vz)"     );
    auto  z_lo    = z->clone( "Control Lower Bound (z_lo)" );
    auto  z_bnd   = makePtr<Bounds>( *z_lo );
    z_lo->zero();

    // State
    auto u_new    = State::create( pl, "New state (u_new)"   );
    auto u_old    = u_new->clone( "Old state (u_old)"        );
    auto u_new_lo = u_new->clone( "State lower bound (u_lo)" );
    auto u_new_up = u_new->clone( "State upper bound (u_up)" );
    auto u        = PartitionedVector::create( { u_old,    u_new    } );
    auto u_lo     = PartitionedVector::create( { u_new_lo, u_new_lo } );
    auto u_up     = PartitionedVector::create( { u_new_up, u_new_up } );

    u_lo->zero();
    u_up->setScalar( height );
    auto u_bnd = makePtr<Bounds>(u_new_lo,u_new_up);

    (*z)(0,0) = Qin00;

    for( size_type i=0; i<nrows; ++i )
      for( size_type j=0; j<ncols; ++j )
        u_old->h(i,j) = h_init;

    // Check the Tanks::DynamicConstraint methods
    ValidateFunction validator( 1, 13, 20, 11, true, *outStream);
    Tanks::check( *dyn_con, *u_bnd, *z_bnd, validator );

    // Check the Tanks::SerialConstraint methods
//    Tanks::SerialConstraint<RealT> serial_con(pl);
    auto timeStamp = makePtr<std::vector<TimeStamp>>(Nt);
    for( size_type k=0; k<Nt; ++k ) {
      timeStamp->at(k).t.resize(2);
      timeStamp->at(k).t.at(0) = k*dt;
      timeStamp->at(k).t.at(1) = (k+1)*dt;
    }

    ROL::SerialConstraint<RealT> serial_con( dyn_con, *u_old, timeStamp );
    Tanks::check( serial_con, *u_bnd, *z_bnd, *outStream );
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
