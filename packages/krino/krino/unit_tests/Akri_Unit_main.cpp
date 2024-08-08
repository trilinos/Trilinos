// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <gtest/gtest.h>
#include <stk_util/environment/Env.hpp>
#include <stk_util/diag/WriterRegistry.hpp>
#include <stk_unit_test_utils/ParallelGtestOutput.hpp>

#include <Akri_DiagWriter.hpp>
#include <Akri_RegisterProduct.hpp>
#include <Akri_Startup.hpp>

int main(int argc, char **argv) {
  sierra::Env::set_input_file_required(false);

  testing::InitGoogleTest(&argc, argv);

  krino::Startup startup__(argc, argv);
    
  stk::unit_test_util::create_parallel_output(sierra::Env::parallel_rank());

  return RUN_ALL_TESTS();
}
