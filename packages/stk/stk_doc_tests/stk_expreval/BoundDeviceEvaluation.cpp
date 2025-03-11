// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//

#include <Kokkos_Core.hpp>
#include <gtest/gtest.h>
#include <stk_ngp_test/ngp_test.hpp>
#include <stk_expreval/Evaluator.hpp>
#include <fstream>
#include <iostream>
#include <limits>
#include <iomanip>
#include <cmath>
#include <memory>

namespace {

//-BEGIN

using ViewInt1DHostType = Kokkos::View<int*, Kokkos::LayoutRight, Kokkos::HostSpace>;

double perform_device_evaluation(const std::string& expression)
{
  stk::expreval::Eval eval(expression);
  eval.parse();

  //For device-side variable binding and evaluation, need to generate a unique index for each variable.
  int yIndex = eval.get_variable_index("y");
  int zIndex = eval.get_variable_index("z");

  //create ParsedEval that holds all necessary info for device
  auto & parsedEval = eval.get_parsed_eval();

  //evaluate the expression on device
  double result = 0.0;
  Kokkos::parallel_reduce(stk::ngp::DeviceRangePolicy(0, 1),
    KOKKOS_LAMBDA (const int& /*i*/, double& localResult) {

      //device data that will be bound to expression variables
      double yDeviceValue = 3.0;
      double zDeviceValue = 4.0;

      //create DeviceVariableMap, which will contain device-side data
      stk::expreval::DeviceVariableMap<> deviceVariableMap(parsedEval);

      //bind variable values via the DeviceVariableMap
      deviceVariableMap.bind(yIndex, yDeviceValue, 1, 1);
      deviceVariableMap.bind(zIndex, zDeviceValue, 1, 1);

      localResult = parsedEval.evaluate(deviceVariableMap);

    }, result);

  return result;
}

TEST(DeviceEvaluation, bindScalar)
{
  double result = perform_device_evaluation("x=5; y=y+x; y+z");
  EXPECT_DOUBLE_EQ(result, 12);
}
//-END

} // namespace <unnamed>
