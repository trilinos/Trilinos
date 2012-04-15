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
#include <fei_ParameterSet.hpp>
#include <Teuchos_ParameterList.hpp>
#include <fei_Trilinos_Helpers.hpp>

#include <vector>
#include <cmath>

#undef fei_file
#define fei_file "fei_unit_Params.cpp"
#include "fei_ErrMacros.hpp"

TEUCHOS_UNIT_TEST(Params, Params_test1)
{
  fei::ParameterSet paramset;
  Teuchos::ParameterList paramlist;

  paramset.add(fei::Param("int1", 1));
  paramset.add(fei::Param("double1",1.0));
  paramset.add(fei::Param("str1", "1.0"));

  Trilinos_Helpers::copy_parameterset(paramset, paramlist);

  TEUCHOS_TEST_EQUALITY(paramlist.get<int>("int1"), 1, out, success);
  TEUCHOS_TEST_EQUALITY((paramlist.get<double>("double1") -1.0) < 1.e-14, true, out, success);
  TEUCHOS_TEST_EQUALITY(paramlist.get<std::string>("str1"), std::string("1.0"), out, success);
}

TEUCHOS_UNIT_TEST(Params, Params_test2)
{
  fei::ParameterSet paramset;
  Teuchos::ParameterList paramlist;

  paramlist.set("int1", 1);
  paramlist.set("double1",1.0);
  paramlist.set("str1", "1.0");

  Trilinos_Helpers::copy_parameterlist(paramlist, paramset);

  int i1;
  int err = paramset.getIntParamValue("int1", i1);
  double d1;
  err += paramset.getDoubleParamValue("double1", d1);
  std::string str1;
  err += paramset.getStringParamValue("str1", str1);

  TEUCHOS_TEST_EQUALITY(err, 0, out, success);
}

