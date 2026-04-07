// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <gtest/gtest.h>
#include <stk_util/environment/Env.hpp>
#include <stk_util/environment/EnvData.hpp>
#include <stk_util/environment/OutputLog.hpp>
#include <stk_unit_test_utils/ParallelGtestOutput.hpp>
#include <stk_unit_test_utils/CommandLineArgs.hpp>
#include <stk_util/diag/WriterRegistry.hpp>
#include <Akri_DiagWriter.hpp>
#include <mpi.h>

int main(int argc, char **argv)
{
    sierra::Env::set_mpi_communicator( stk::initialize(&argc, &argv) );

    sierra::Diag::registerWriter("krinolog", krinolog, krino::theDiagWriterParser());
    stk::EnvData::instance().m_outputP0 = &sierra::out();
    //stk::bind_output_streams("outfile=krinolog.txt out>outfile+pout dout>out pout>null"); // necessary for krinolog to work, otherwise you will get segfault
    stk::bind_output_streams("dout>out out>pout pout>null"); // necessary for krinolog to work, otherwise you will get segfault
    stk::unit_test_util::create_parallel_output_with_comm(sierra::Env::parallel_rank(), stk::EnvData::instance().parallel_comm());

    stk::unit_test_util::GlobalCommandLineArguments::self().set_values(argc, argv);

    testing::InitGoogleTest(&argc, argv);
    int returnVal = RUN_ALL_TESTS();

    stk::finalize();

    return returnVal;
}
