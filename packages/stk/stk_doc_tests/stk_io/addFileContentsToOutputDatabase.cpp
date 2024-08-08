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

#include <gtest/gtest.h>                // for AssertHelper, EXPECT_STREQ, etc
#include <stddef.h>                     // for size_t
#include <unistd.h>                     // for unlink
#include <ostream>                      // for basic_ostream, operator<<, etc
#include <stk_io/StkMeshIoBroker.hpp>   // for StkMeshIoBroker
#include <string>                       // for string, basic_string, etc
#include <vector>                       // for vector
#include "Ioss_DBUsage.h"               // for DatabaseUsage::READ_MODEL
#include "Ioss_Property.h"              // for Property
#include "Ioss_Region.h"                // for Region
#include "stk_io/DatabasePurpose.hpp"   // for DatabasePurpose::READ_MESH, etc
#include "stk_util/parallel/Parallel.hpp"
namespace Ioss { class DatabaseIO; }

namespace {

TEST(StkMeshIoBrokerHowTo, addFileContentsToOutputDatabase)
{
  std::string filename = "information_records.e";
  MPI_Comm communicator = MPI_COMM_WORLD;
  int numProcs = stk::parallel_machine_size(communicator);
  if (numProcs != 1) {
    return;
  }

  //-BEGIN
  // ============================================================
  //+ SETUP
  std::string input_file = "application_input_file.i";
  std::string info1("This is the first line of the input file.");
  std::string info2("This is the second line of the input file. "
                    "It is longer than 80 characters, so it should be wrapped.");
  std::string info3("This is the third line of the input file.");
  std::string info4("This is the fourth and last line of the input file.");

  std::string additional_info_record = "This is an info record added explicitly,"
                                       " not from the input file.";
  {
    std::ofstream my_file(input_file.c_str());
    my_file << info1 <<"\n" << info2 <<"\n" << info3 <<"\n" << info4 <<"\n";
  }

  {
    // ============================================================
    //+ EXAMPLE
    stk::io::StkMeshIoBroker stkIo(communicator);
    size_t ifh = stkIo.add_mesh_database("9x9x9|shell:xyzXYZ", "generated", stk::io::READ_MESH);
    stkIo.set_active_mesh(ifh);
    stkIo.create_input_mesh();
    stkIo.populate_bulk_data();

    // Output...
    size_t fh = stkIo.create_output_mesh(filename,
                                         stk::io::WRITE_RESULTS);
    Ioss::Region *io_reg = stkIo.get_output_ioss_region(fh).get();

    //+ Add the data from the file "application_input_file.i"
    //+    as information records on this file.
    io_reg->property_add(Ioss::Property("input_file_name",input_file)); /*@\label{io:info:file}*/

    //+ Add the data from the "additional_info_record" vector as
    //+    information records on this file.
    io_reg->add_information_record(additional_info_record); /*@\label{io:info:vector}*/

    stkIo.write_output_mesh(fh);
    // ... Verification deleted
    //-END
  }

  // ============================================================
  // VERIFICATION

  {
    stk::io::StkMeshIoBroker stkIo(communicator);
    // Verify output mesh contains the data in
    // 'input_file' as information records...  Note that
    // the output mesh will contain all element blocks; however, the
    // non-shell element block will have zero elements.  This is due
    // to the subset_selector subsetting the entities and not the
    // parts...
    size_t index = stkIo.add_mesh_database(filename, stk::io::READ_MESH);
    stkIo.set_active_mesh(index);
    stkIo.create_input_mesh();

    const std::vector<std::string> &info_records = stkIo.get_info_records();
    // First line of info records is host information (node name,
    // os version) (2) Next record is the file name of the input file
    // data that follows (1) File contains 4 records; 1 is longer than
    // 80 characters, so it wraps (4+1) Next line is the
    // "additional_info_record" added above (1)
    // Last records are the IOSS configuration summary (35).
    size_t expected_min_info_record_count = 1 + (4+1) + 1 + 1;
    EXPECT_TRUE(expected_min_info_record_count <= info_records.size());

    EXPECT_STREQ(input_file.c_str(), info_records[1].c_str());

    EXPECT_STREQ(info1.c_str(), info_records[2].c_str());
    EXPECT_STREQ(info2.substr(0,79).c_str(), info_records[3].substr(0,79).c_str());
    EXPECT_STREQ(info2.substr(79).c_str(), info_records[4].c_str());
    EXPECT_STREQ(info3.c_str(), info_records[5].c_str());
    EXPECT_STREQ(info4.c_str(), info_records[6].c_str());
    EXPECT_STREQ(additional_info_record.c_str(), info_records[7].c_str());

    //      unlink(filename.c_str());
    unlink(input_file.c_str());
  }
}

}
//TODO:
//. Check only on processor 0
