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

#ifndef GENERATED_MESH_TO_FILE_H_
#define GENERATED_MESH_TO_FILE_H_

// #######################  Start Clang Header Tool Managed Headers ########################
// clang-format off
#include <stddef.h>                             // for size_t
#include <stk_io/StkMeshIoBroker.hpp>           // for StkMeshIoBroker
#include <stk_mesh/base/BulkData.hpp>           // for BulkData, etc
#include <stk_mesh/base/CoordinateSystems.hpp>  // for Cartesian
#include <stk_mesh/base/Field.hpp>              // for Field
#include <stk_mesh/base/MetaData.hpp>           // for MetaData
#include <string>                               // for string
#include <vector>                               // for vector
#include "stk_topology/topology.hpp"            // for topology, etc
#include "stk_util/parallel/Parallel.hpp"       // for ParallelMachine
#include "stk_unit_test_utils/BuildMesh.hpp"
namespace stk { namespace unit_test_util { class FieldValueSetter; } }
// clang-format on
// #######################   End Clang Header Tool Managed Headers  ########################

namespace stk
{
namespace unit_test_util
{

class GeneratedMeshToFile
{
public:
  GeneratedMeshToFile(stk::ParallelMachine comm, stk::mesh::BulkData::AutomaticAuraOption auraOption);

  void setup_mesh(const std::string &meshSizeSpec,
                  const std::string &outputFileName,
                  bool useBigIds = true);
  void write_mesh();

protected:
  std::shared_ptr<stk::mesh::BulkData> bulkPtr;
  stk::mesh::BulkData& bulk;
  stk::mesh::MetaData& meta;
  stk::io::StkMeshIoBroker broker;
  size_t outputFileIndex = 0;

private:
  GeneratedMeshToFile();
};

class GeneratedMeshToFileWithTransientFields : public GeneratedMeshToFile
{
public:
  GeneratedMeshToFileWithTransientFields(stk::ParallelMachine comm,
                                         stk::mesh::BulkData::AutomaticAuraOption auraOption,
                                         const std::string& fieldBaseName,
                                         stk::topology::rank_t rank);

  virtual ~GeneratedMeshToFileWithTransientFields() = default;

  void write_mesh_with_field(const std::vector<double>& timeSteps,
                             const FieldValueSetter &fieldValueSetter,
                             const std::string& globalVariableName);

protected:
  stk::topology::rank_t fieldRank;
  stk::mesh::Field<double> &scalarField;
  stk::mesh::Field<double> &vectorField;

private:
  GeneratedMeshToFileWithTransientFields();
};

namespace simple_fields {

class STK_DEPRECATED_MSG("Please use the non-simple_fields-namespaced version of this class instead")
GeneratedMeshToFile
{
public:
  GeneratedMeshToFile(stk::ParallelMachine comm, stk::mesh::BulkData::AutomaticAuraOption auraOption);

  void setup_mesh(const std::string &meshSizeSpec,
                  const std::string &outputFileName,
                  bool useBigIds = true);
  void write_mesh();

protected:
  std::shared_ptr<stk::mesh::BulkData> bulkPtr;
  stk::mesh::BulkData& bulk;
  stk::mesh::MetaData& meta;
  stk::io::StkMeshIoBroker broker;
  size_t outputFileIndex = 0;

private:
  GeneratedMeshToFile();
};

class STK_DEPRECATED_MSG("Please use the non-simple_fields-namespaced version of this class instead")
GeneratedMeshToFileWithTransientFields : public GeneratedMeshToFile
{
public:
  GeneratedMeshToFileWithTransientFields(stk::ParallelMachine comm,
                                         stk::mesh::BulkData::AutomaticAuraOption auraOption,
                                         const std::string& fieldBaseName,
                                         stk::topology::rank_t rank);

  virtual ~GeneratedMeshToFileWithTransientFields() = default;

  void write_mesh_with_field(const std::vector<double>& timeSteps,
                             const FieldValueSetter &fieldValueSetter,
                             const std::string& globalVariableName);

protected:
  stk::topology::rank_t fieldRank;
  stk::mesh::Field<double> &scalarField;
  stk::mesh::Field<double> &vectorField;

private:
  GeneratedMeshToFileWithTransientFields();
};

} // namespace simple_fields

}
}

#endif
