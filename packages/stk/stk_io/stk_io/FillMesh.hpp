// Copyright (c) 2013, Sandia Corporation.
 // Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 // the U.S. Government retains certain rights in this software.
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
 //     * Neither the name of Sandia Corporation nor the names of its
 //       contributors may be used to endorse or promote products derived
 //       from this software without specific prior written permission.
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

#ifndef PACKAGES_STK_STK_IO_STK_IO_FILLMESH_HPP_
#define PACKAGES_STK_STK_IO_STK_IO_FILLMESH_HPP_

#include <string>
#include "stk_io/DatabasePurpose.hpp"  // for DatabasePurpose::READ_MESH

namespace stk { namespace mesh { class BulkData; }}
namespace stk { namespace io { class StkMeshIoBroker; } }

namespace stk
{
namespace io
{

void fill_mesh(const std::string &meshSpec, stk::mesh::BulkData &bulkData);
void fill_mesh(const std::string &meshSpec, stk::mesh::BulkData &bulkData, stk::io::StkMeshIoBroker &stkIo);
void fill_mesh_with_auto_decomp(const std::string &meshSpec, stk::mesh::BulkData &bulkData);
void fill_mesh_with_auto_decomp(const std::string &meshSpec, stk::mesh::BulkData &bulkData, stk::io::StkMeshIoBroker &stkIo);
void fill_mesh_preexisting(stk::io::StkMeshIoBroker & stkIo, const std::string& meshSpec,
                           stk::mesh::BulkData& bulkData, stk::io::DatabasePurpose purpose = stk::io::READ_MESH);
void fill_mesh_save_step_info(const std::string& inFile, stk::mesh::BulkData& inBulk, int &numSteps, double &maxTime);

} // namespace unit_test_util
} // namespace stk

#endif /* PACKAGES_STK_STK_IO_STK_IO_FILLMESH_HPP_ */
