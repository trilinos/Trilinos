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
#include <memory>

#include "stk_io/DatabasePurpose.hpp" // for stk::io::DatabasePurpose
#include "stk_io/FillMesh.hpp"
#include "stk_io/WriteMesh.hpp"
#include "stk_mesh/base/Types.hpp"
#include "stk_util/parallel/Parallel.hpp"  // for ParallelMachine, etc
#include "stk_mesh/base/BulkData.hpp"
#include "stk_mesh/base/MetaData.hpp"
#include "stk_io/StkMeshIoBroker.hpp"

namespace stk
{
namespace unit_test_util
{
void initialize_stk_io_for_text_mesh();

// Example:  meshSizeSpec = "2x2x1"
void generated_mesh_to_file_in_serial(const std::string& meshSizeSpec, const std::string& fileName);
void text_mesh_to_file_in_serial(const std::string& meshDesc, const std::string& fileName);

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

class MeshFromFile
{
public:
  MeshFromFile(const MPI_Comm& c);

  void fill_from_serial(const std::string& fileName);
  void fill_from_parallel(const std::string& baseName);
  bool is_empty() const { return m_empty; }

private:
  MPI_Comm m_comm;
  bool m_empty;

public:
  std::shared_ptr<stk::mesh::BulkData> bulk;
  stk::mesh::MetaData& meta;
  stk::io::StkMeshIoBroker broker;
};

class TransientVerifier
{
public:
  TransientVerifier(const MPI_Comm& c);

  void verify_num_transient_fields(const MeshFromFile& mesh, unsigned expectedNumFields) const;
  void verify_time_steps(const MeshFromFile& mesh, const std::vector<double>& expectedTimeSteps) const;
  void verify_global_variables_at_each_time_step(MeshFromFile& mesh,
                                                 const std::string& globalVariableName,
                                                 const std::vector<double>& expectedTimeSteps) const;
  void verify_sideset_orientation(const MeshFromFile& mesh,
                                  int expectedProc,
                                  const stk::mesh::EntityId expectedId,
                                  const stk::mesh::ConnectivityOrdinal expectedOrdinal) const;
  void compare_entity_rank_names(const MeshFromFile& meshA, const MeshFromFile& meshB) const;
  void verify_transient_field_names(const MeshFromFile& mesh, const std::string& fieldBaseName) const;
  void verify_transient_fields(MeshFromFile& mesh) const;
  void verify_decomp(MeshFromFile& mesh, const stk::mesh::EntityIdProcVec& expectedDecomp) const;

private:
  void verify_global_variable_names(const MeshFromFile& mesh, const std::string& baseName) const;
  void verify_global_double(const MeshFromFile& mesh, const std::string& variable, double goldValue) const;
  void verify_global_int(const MeshFromFile& mesh, const std::string& variable, int goldValue) const;
  void verify_global_real_vec(const MeshFromFile& mesh, const std::string& variable, double goldValue) const;
  void verify_transient_field_name(stk::mesh::FieldBase* field, const std::string& fieldName) const;
  void verify_transient_field_values(const stk::mesh::BulkData& bulk, stk::mesh::FieldBase* field, double timeStep) const;

  MPI_Comm m_comm;
  const double m_epsilon;
};

void generated_mesh_with_transient_data_to_file_in_serial(const std::string &meshSizeSpec,
                                                          const std::string &fileName,
                                                          const std::string& fieldName,
                                                          stk::topology::rank_t fieldRank,
                                                          const std::string& globalVariableName,
                                                          const std::vector<double>& timeSteps,
                                                          const FieldValueSetter &fieldValueSetter);

void read_from_serial_file_and_decompose(const std::string& fileName, stk::mesh::BulkData &mesh,
                                         const std::string &decompositionMethod);


namespace simple_fields {

STK_DEPRECATED_MSG("Please use the non-simple_fields-namespaced version of this function instead")
void generated_mesh_to_file_in_serial(const std::string& meshSizeSpec, const std::string& fileName);

STK_DEPRECATED_MSG("Please use the non-simple_fields-namespaced version of this function instead")
void text_mesh_to_file_in_serial(const std::string& meshDesc, const std::string& fileName);

STK_DEPRECATED_MSG("Please use the non-simple_fields-namespaced version of this function instead")
void generate_mesh_from_serial_spec_and_load_in_parallel_with_auto_decomp(const std::string &meshSizeSpec, stk::mesh::BulkData & mesh, const std::string &decompositionMethod);

class MeshFromFile
{
public:
  MeshFromFile(const MPI_Comm& c);

  void fill_from_serial(const std::string& fileName);
  void fill_from_parallel(const std::string& baseName);
  bool is_empty() const { return m_empty; }

private:
  MPI_Comm m_comm;
  bool m_empty;

public:
  std::shared_ptr<stk::mesh::BulkData> bulk;
  stk::mesh::MetaData& meta;
  stk::io::StkMeshIoBroker broker;
};

class STK_DEPRECATED_MSG("Please use the non-simple_fields-namespaced version of this class instead")
TransientVerifier
{
public:
  TransientVerifier(const MPI_Comm& c);

  void verify_num_transient_fields(const MeshFromFile& mesh, unsigned expectedNumFields) const;
  void verify_time_steps(const MeshFromFile& mesh, const std::vector<double>& expectedTimeSteps) const;
  void verify_global_variables_at_each_time_step(MeshFromFile& mesh,
                                                 const std::string& globalVariableName,
                                                 const std::vector<double>& expectedTimeSteps) const;
  void verify_sideset_orientation(const MeshFromFile& mesh,
                                  int expectedProc,
                                  const stk::mesh::EntityId expectedId,
                                  const stk::mesh::ConnectivityOrdinal expectedOrdinal) const;
  void compare_entity_rank_names(const MeshFromFile& meshA, const MeshFromFile& meshB) const;
  void verify_transient_field_names(const MeshFromFile& mesh, const std::string& fieldBaseName) const;
  void verify_transient_fields(MeshFromFile& mesh) const;
  void verify_decomp(MeshFromFile& mesh, const stk::mesh::EntityIdProcVec& expectedDecomp) const;

private:
  void verify_global_variable_names(const MeshFromFile& mesh, const std::string& baseName) const;
  void verify_global_double(const MeshFromFile& mesh, const std::string& variable, double goldValue) const;
  void verify_global_int(const MeshFromFile& mesh, const std::string& variable, int goldValue) const;
  void verify_global_real_vec(const MeshFromFile& mesh, const std::string& variable, double goldValue) const;
  void verify_transient_field_name(stk::mesh::FieldBase* field, const std::string& fieldName) const;
  void verify_transient_field_values(const stk::mesh::BulkData& bulk, stk::mesh::FieldBase* field, double timeStep) const;

  MPI_Comm m_comm;
  const double m_epsilon;
};

STK_DEPRECATED_MSG("Please use the non-simple_fields-namespaced version of this function instead")
void generated_mesh_with_transient_data_to_file_in_serial(const std::string &meshSizeSpec,
                                                          const std::string &fileName,
                                                          const std::string& fieldName,
                                                          stk::topology::rank_t fieldRank,
                                                          const std::string& globalVariableName,
                                                          const std::vector<double>& timeSteps,
                                                          const FieldValueSetter &fieldValueSetter);

STK_DEPRECATED_MSG("Please use the non-simple_fields-namespaced version of this function instead")
void read_from_serial_file_and_decompose(const std::string& fileName, stk::mesh::BulkData &mesh,
                                         const std::string &decompositionMethod);

} // namespace simple_fields

} // namespace unit_test_util
} // namespace stk

#endif
