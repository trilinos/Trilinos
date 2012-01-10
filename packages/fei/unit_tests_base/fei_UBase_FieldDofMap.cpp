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

#include <fei_iostream.hpp>
#include <fei_FieldDofMap.hpp>

#include <vector>
#include <cmath>

namespace {

TEUCHOS_UNIT_TEST(FieldDofMap, test1)
{
  fei::FieldDofMap<int> fdmap;

  fdmap.add_field(0, 1);
  fdmap.add_field(2, 2);
  fdmap.add_field(5, 5);

  int dof_id = fdmap.get_dof_id(0, 0);
  int dof_id_correct = fei::UNKNOWN;

  TEUCHOS_TEST_EQUALITY(dof_id, dof_id_correct, out, success);

  dof_id = fdmap.get_dof_id(2, 1);
  dof_id_correct = fei::UNKNOWN + 2;

  TEUCHOS_TEST_EQUALITY(dof_id, dof_id_correct, out, success);

  dof_id = fdmap.get_dof_id(5, 3);
  dof_id_correct = fei::UNKNOWN + 6;

  TEUCHOS_TEST_EQUALITY(dof_id, dof_id_correct, out, success);

  TEUCHOS_TEST_THROW(fdmap.get_dof_id(0, 1), std::runtime_error, out, success);
}

}//namespace <anonymous>

