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


/*! \file  trust_region_algorithm.cpp
    \brief Test solving an unconstrained problem using a Trust Region Method
*/

#include "ROL2.hpp"
#include "ROL2_TypeU_TestProblems.hpp"

int main( int argc, char *argv[] ) {

  using RealT = double;

  auto os_ptr = ROL2::makeStreamPtr(std::cout, argc);
  auto& os = *os_ptr;

  int errorFlag  = 0;

  auto errtol = ROL2::default_tolerance<RealT>();

  try {

    // Number of elements in StdVector objects
    int dim = 10; 

    // Create optimization vector and set initial guess to all 1's
    ROL2::StdVector<RealT> x(dim);
    x.setScalar(1.0);

    // gradient vector
    auto g_ptr = x.clone();
    auto& g = *g_ptr;

    // Create Zakharov objective
    auto k_ptr = x.clone();
    auto& k_data =  ROL2::StdVector<RealT>::getData(*k_ptr);
    for( int i=0; i<10; ++i ) k_data[i] = 1+i;
    auto obj = ROL2::TypeU::Objective_Zakharov<RealT>(k_ptr);

    ROL2::ParameterList parlist;

    auto algo = ROL2::TypeU::TrustRegionAlgorithm<RealT>( parlist ); 

     

  }
  catch (std::logic_error& err) {
    os << err.what() << "\n";
    errorFlag = -1000;
  }; // end try

  if (errorFlag != 0) std::cout << "End Result: TEST FAILED\n";
  else                std::cout << "End Result: TEST PASSED\n";

  return 0;

}

