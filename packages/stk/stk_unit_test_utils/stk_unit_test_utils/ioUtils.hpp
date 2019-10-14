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

#ifndef ioUtils_hpp
#define ioUtils_hpp

#include <string>
#include <vector>

#include <stk_io/DatabasePurpose.hpp> // for stk::io::DatabasePurpose
#include <stk_io/FillMesh.hpp>
#include <stk_io/WriteMesh.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_util/parallel/Parallel.hpp>  // for ParallelMachine, etc
namespace stk { namespace mesh { class BulkData; }}
namespace stk { namespace mesh { class FieldBase; }}
namespace stk { namespace io { class StkMeshIoBroker; } }

namespace stk
{
namespace unit_test_util
{

// Example:  meshSizeSpec = "2x2x1"
void generated_mesh_to_file_in_serial(const std::string &meshSizeSpec, const std::string &fileName);

void read_from_serial_file_and_decompose(const std::string& fileName, stk::mesh::BulkData &mesh, const std::string &decompositionMethod);

// This avoid the GeneratedMesh limitation on the z-dimension >= num_processors,
// and allows cases like 2x2x1 with cyclic decomposition on four processors.
void generate_mesh_from_serial_spec_and_load_in_parallel_with_auto_decomp(const std::string &meshSizeSpec, stk::mesh::BulkData & mesh, const std::string &decompositionMethod);

class FieldValueSetter
{
public:
    virtual ~FieldValueSetter() {}
    virtual void populate_field(stk::mesh::BulkData &bulk, stk::mesh::FieldBase* field, const unsigned step, const double time) const = 0;
};

class IdAndTimeFieldValueSetter : public FieldValueSetter
{
public:
    virtual void populate_field(stk::mesh::BulkData &bulk, stk::mesh::FieldBase* field, const unsigned step, const double time) const;
};

void generated_mesh_with_transient_data_to_file_in_serial(const std::string &meshSizeSpec,
                                                          const std::string &fileName,
                                                          const std::string& fieldName,
                                                          const std::string& globalVariableName,
                                                          const std::vector<double>& timeSteps,
                                                          const FieldValueSetter &fieldValueSetter);

} // namespace unit_test_util
} // namespace stk

#endif
