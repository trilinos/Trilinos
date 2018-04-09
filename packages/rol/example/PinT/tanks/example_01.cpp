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

#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "TankConstraint.hpp"
#include "TankVector.hpp"
#include "LowerBandedMatrix.hpp"

#include <iostream>


using RealT = double;
using size_type = std::vector<RealT>::size_type;

int main( int argc, char* argv[] ) {
  
  using std::vector;

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);  

  int iprint     = argc - 1;
  ROL::Ptr<std::ostream> outStream;
  Teuchos::oblackholestream bhs; // outputs nothing
  if (iprint > 0)
    outStream = ROL::makePtrFromRef(std::cout);
  else
    outStream = ROL::makePtrFromRef(bhs);

  int errorFlag  = 0;

  auto print_vector = [outStream]( const vector<RealT>& x) { 
    *outStream << "\n";
    for( auto e : x ) *outStream << e << " ";
    *outStream << "\n";
  };
  // *** Example body.

//  try {   

    auto tank_parameters = ROL::makePtr<Teuchos::ParameterList>();
    std::string tank_xml("tank-parameters.xml");
    Teuchos::updateParametersFromXmlFile(tank_xml, tank_parameters.ptr());
    auto& pl = *tank_parameters;

    TankConstraint<RealT> con(pl);

    auto Qin00 = pl.get("Corner Inflow",100.0);
    auto hinit = pl.get("Initial Fluid Level", 2.0);
    auto nrows  = static_cast<size_type>( pl.get("Number of Rows",3) );
    auto ncols  = static_cast<size_type>( pl.get("Number of Columns",3) );

    TankControlVector<RealT> z(nrows,ncols,"z");    
    TankStateVector<RealT>   un(nrows,ncols,"u_new");
    TankStateVector<RealT>   uo(nrows,ncols,"u_old");
    TankStateVector<RealT>   c(nrows,ncols,"c");

    z(0,0) = Qin00;

    for( size_type i=0; i<nrows; ++i ) {
      for( size_type j=0; j<ncols; ++j ) {
        uo.h(i,j) = hinit;
      }
    }


    uo.print(*outStream);

    RealT tol = 0;
    con.solve( c, uo, un, z, tol );
 
    un.print(*outStream);
    c.print(*outStream);

    con.value( c, uo, un, z, tol );
    c.print(*outStream);


//  }
//  catch (std::logic_error err) {
//    *outStream << err.what() << "\n";
//    errorFlag = -1000;
//  }; // end try

  if (errorFlag != 0)
    std::cout << "End Result: TEST FAILED\n";
  else
    std::cout << "End Result: TEST PASSED\n";

  return 0;
}
