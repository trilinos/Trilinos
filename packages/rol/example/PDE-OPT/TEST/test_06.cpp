// @HEADER
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

/*! \file test_06.cpp
    \brief Unit test for Automatic Differentiation evalutators
*/

#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"


#include "../TOOLS/template_tools.hpp"

int main(int argc, char *argv[] ) {

  using RealT = double;

  using Teuchos::RCP; using Teuchos::rcp; 
  Teuchos::oblackholestream bhs;

  RCP<std::ostream> os;
  if(argc>1)   os = rcp(&std::cout,false);
  else         os = rcp(&bhs,false);
  
  int errorFlag = 0;

  RealT x = 2.0;

  using DFad = Sacado::Fad::DFad<RealT>;
  DFad x_fad(1,0,x);
  DFad f_fad = x_fad*x_fad;

  errorFlag += static_cast<int>(AD::has_dx<RealT>::value);
  errorFlag += static_cast<int>(!AD::has_dx<DFad>::value);

  *os << "Does type RealT = double have a method .dx()?... " << AD::has_dx<RealT>::value << std::endl;
  *os << "Does type DFad<RealT> have a method .dx()?...... " << AD::has_dx<DFad>::value << std::endl;

  auto xval  = AD::value(x);
  auto xfval = AD::value(x_fad);
  auto ffval = AD::value(f_fad);

  auto xder  = AD::derivative(x,0);
  auto xfder = AD::derivative(x_fad,0);
  auto ffder = AD::derivative(f_fad,0);

  *os << "RealT x(2.0);"             << std::endl;
  *os << "DFad x_fad(1,x); "         << std::endl;
  *os << "DFad f_fad = x_fad*x_fad " << std::endl;  

  *os << "AD::value(x)     = " << xval  << std::endl;
  *os << "AD::value(x_fad) = " << xfval << std::endl;
  *os << "AD::value(f_fad) = " << ffval << std::endl;

  errorFlag += (xval  != 2);
  errorFlag += (xfval != 2);
  errorFlag += (ffval != 4);

  *os << "AD::derivative(x,0)      = "  << xder  << std::endl;
  *os << "AD::derivative(x_fad,0)  = "  << xfder << std::endl;
  *os << "AD::derivative(f_fad,0)  = "  << ffder << std::endl;

  errorFlag += (xder  != 0);
  errorFlag += (xfder != 1);
  errorFlag += (ffder != 4);

  if (errorFlag != 0)
    std::cout << "End Result: TEST FAILED\n";
  else
    std::cout << "End Result: TEST PASSED\n";

  return 0;

}
