/*
// @HEADER
// ************************************************************************
//             FEI: Finite Element Interface to Linear Solvers
//                  Copyright (2005) Sandia Corporation.
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
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
// Questions? Contact Alan Williams (william@sandia.gov) 
//
// ************************************************************************
// @HEADER
*/


#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>

#include <fei_MatrixTraits.hpp>
#include <fei_MatrixTraits_FillableMat.hpp>

#include <fei_iostream.hpp>

#include <cmath>
#include <limits>

TEUCHOS_UNIT_TEST(MatrixTraits, FillableMat_1)
{
  fei::FillableMat fm;

  fm.putCoef(0, 0, 0.0);
  fm.putCoef(1, 1, 1.1);
  fm.putCoef(1, 2, 1.2);
  fm.putCoef(2, 2, 2.2);

  int numRows = 0;
  fei::MatrixTraits<fei::FillableMat>::getNumLocalRows(&fm, numRows);
  TEUCHOS_TEST_EQUALITY(numRows, 3, out, success);

  std::vector<int> indices1(2);
  std::vector<double> coefs1(2);

  fei::MatrixTraits<fei::FillableMat>::copyOutRow(&fm, 1, 2, &coefs1[0], &indices1[0]);

  TEUCHOS_TEST_EQUALITY(indices1[0], 1, out, success);
  TEUCHOS_TEST_EQUALITY(indices1[1], 2, out, success);

  const double eps = std::numeric_limits<double>::epsilon();

  TEUCHOS_TEST_EQUALITY(std::abs(coefs1[0] - 1.1) < eps, true, out, success);
  TEUCHOS_TEST_EQUALITY(std::abs(coefs1[1] - 1.2) < eps, true, out, success);
}

