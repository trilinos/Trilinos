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
// Questions: Alejandro Mota (amota@sandia.gov)
//
// ************************************************************************
// @HEADER

#include <gtest/gtest.h>
#include "ROL_Algorithm.hpp"
#include "ROL_LineSearchStep.hpp"
#include "ROL_StatusTest.hpp"
#include <Intrepid2_MiniTensor_FunctionSet.h>
#include "ROL_MiniTensor_Function.hpp"

typedef double Real;
using Real = double;

int
main(int ac, char * av[])
{
  Kokkos::initialize();

  // Disables elapsed time and output by default.
  ::testing::GTEST_FLAG(print_time) = false;

  ::testing::InitGoogleTest(&ac, av);

  auto const
  retval = RUN_ALL_TESTS();

  Kokkos::finalize();

  return retval;
}

TEST(MiniTensor_ROL, Paraboloid)
{
  bool const
  print_output = ::testing::GTEST_FLAG(print_time);

  // outputs nothing
  Teuchos::oblackholestream
  bhs;

  std::ostream &
  os = (print_output == true) ? std::cout : bhs;

  constexpr Intrepid2::Index
  DIM{2};

  ROL::MiniTensor_Objective<Intrepid2::Paraboloid, Real, DIM>
  obj;

  // Set parameters.
  Teuchos::ParameterList
  parlist;

  parlist.sublist("Step").sublist("Line Search").sublist("Descent Method")
      .set("Type", "Newton-Krylov");

  parlist.sublist("Status Test").set("Gradient Tolerance", 10.e-12);
  parlist.sublist("Status Test").set("Step Tolerance", 1.0e-14);
  parlist.sublist("Status Test").set("Iteration Limit", 128);

  // Define algorithm.
  ROL::Algorithm<Real>
  algo("Line Search", parlist);

  // Set Initial Guess
  Intrepid2::Vector<Real, DIM>
  xval(Intrepid2::RANDOM);

  ROL::MiniTensorVector<Real, DIM> x(xval);

  // Run Algorithm
  algo.run(x, obj, true, os);

  Intrepid2::Vector<Real, DIM> const
  sol = ROL::MTfromROL<Real, DIM>(x);

  os << "Solution : " << sol << '\n';

  Real const
  epsilon{Intrepid2::machine_epsilon<Real>()};

  Real const
  error = Intrepid2::norm(sol);

  ASSERT_LE(error, epsilon);
}

