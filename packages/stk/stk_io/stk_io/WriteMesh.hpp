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

#ifndef PACKAGES_STK_STK_IO_STK_IO_WRITEMESH_HPP_
#define PACKAGES_STK_STK_IO_STK_IO_WRITEMESH_HPP_

// #######################  Start Clang Header Tool Managed Headers ########################
// clang-format off
#include <stk_io/DatabasePurpose.hpp>  // for WRITE_RESULTS, DatabasePurpose
#include <string>                      // for string
#include <vector>
namespace stk { namespace io { class StkMeshIoBroker; } }
namespace stk { namespace mesh { class BulkData; } }
namespace stk { namespace mesh { class Selector; } }
// clang-format on
// #######################   End Clang Header Tool Managed Headers  ########################

namespace stk
{
namespace io
{

void write_mesh(const std::string &filename,
                stk::mesh::BulkData &bulkData,
                stk::io::DatabasePurpose databasePurpose = stk::io::WRITE_RESULTS);

void write_mesh(const std::string &filename,
                stk::io::StkMeshIoBroker& ioBroker,
                stk::io::DatabasePurpose databasePurpose = stk::io::WRITE_RESULTS);

void write_mesh_with_canonical_name(const std::string &filename,
                                    stk::mesh::BulkData &bulkData,
                                    stk::io::DatabasePurpose databasePurpose = stk::io::WRITE_RESULTS);

void write_mesh_with_large_ids(const std::string &filename,
                stk::mesh::BulkData &bulkData,
                stk::io::DatabasePurpose databasePurpose = stk::io::WRITE_RESULTS);

void write_mesh_subset(const std::string &filename,
                stk::mesh::BulkData &bulkData,
                stk::mesh::Selector& subsetSelector,
                stk::io::DatabasePurpose databasePurpose = stk::io::WRITE_RESULTS);

void write_mesh_with_fields(const std::string& filename,
                            stk::mesh::BulkData &bulkData,
                            int step=0, double time=0.0,
                            stk::io::DatabasePurpose databasePurpose = stk::io::WRITE_RESULTS);

void write_mesh_with_specified_fields(const std::string& filename, stk::mesh::BulkData &bulkData,
							const std::vector<std::string>& fieldNames,
							int step, double time, stk::io::DatabasePurpose databasePurpose);

void write_mesh_with_fields(const std::string& filename,
                            StkMeshIoBroker &broker,
                            int step=0, double time=0.0,
                            stk::io::DatabasePurpose databasePurpose = stk::io::WRITE_RESULTS);

void write_mesh_with_specified_fields(const std::string& filename,
                            StkMeshIoBroker &broker,
							const std::vector<std::string>& fieldNames,
                            int step=0, double time=0.0,
                            stk::io::DatabasePurpose databasePurpose = stk::io::WRITE_RESULTS);

void set_64bit_properties(stk::io::StkMeshIoBroker &broker);

void write_mesh_with_large_ids_and_fields(const std::string& filename,
                                          stk::mesh::BulkData &bulkData,
                                          int step=0, double time=0.0,
                                          stk::io::DatabasePurpose databasePurpose = stk::io::WRITE_RESULTS);

} // namespace unit_test_util
} // namespace stk

#endif /* PACKAGES_STK_STK_IO_STK_IO_WRITEMESH_HPP_ */
