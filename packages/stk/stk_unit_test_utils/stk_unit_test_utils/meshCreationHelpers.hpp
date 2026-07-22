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

#ifndef stk_unit_test_utils_meshcreationhelpers_hpp
#define stk_unit_test_utils_meshcreationhelpers_hpp

#include <stk_util/stk_config.h>
#include <stk_util/parallel/Parallel.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_unit_test_utils/FieldEvaluator.hpp>
#include <Ioss_Field.h>
#include <string>
#include <memory>

namespace stk
{
namespace unit_test_util
{
void create_mesh_without_time_steps(const std::string & filename, MPI_Comm communicator);
void create_mesh_with__field_1__field_2__field_3(const std::string & filename, MPI_Comm communicator);

void write_hex20_mesh_file(MPI_Comm comm, size_t nx, size_t ny, size_t nz,
                           const std::string& exoFileName);

std::shared_ptr<stk::mesh::BulkData>
create_mesh_with_field(const stk::ParallelMachine comm,
                       const std::string& fieldName = "field_1",
                       const double fieldValue = 0.0,
                       const stk::mesh::EntityRank fieldRank = stk::topology::NODE_RANK,
                       const Ioss::Field::RoleType fieldRole = Ioss::Field::TRANSIENT,
                       const std::string generatedMeshString = "generated:1x1x1|nodeset:xyz|sideset:xyz",
                       const int numCopies = 1,
                       const int spatialDim = 3);

inline
std::shared_ptr<stk::mesh::BulkData>
create_mesh_with_field(const stk::ParallelMachine comm,
                       const stk::mesh::EntityRank fieldRank)
{
  return create_mesh_with_field(comm, "field_1", 0.0, fieldRank);
}

inline
std::shared_ptr<stk::mesh::BulkData>
create_mesh_with_field(const stk::ParallelMachine comm,
                       const stk::mesh::EntityRank fieldRank,
                       const int numCopies)
{
  return create_mesh_with_field(comm, "field_1", 0.0, fieldRank,
                                Ioss::Field::TRANSIENT,
                                "generated:1x1x1|nodeset:xyz|sideset:xyz",
                                numCopies);
}

inline
std::shared_ptr<stk::mesh::BulkData>
create_mesh_with_field(const stk::ParallelMachine comm,
                       const stk::mesh::EntityRank fieldRank,
                       const std::string generatedMeshString,
                       const int numCopies = 1)
{
  return create_mesh_with_field(comm, "field_1", 0.0, fieldRank, 
                                Ioss::Field::TRANSIENT, generatedMeshString, numCopies);
}

std::shared_ptr<stk::mesh::BulkData>
create_mesh_with_field(const stk::ParallelMachine comm,
                       const std::string& fieldName,
                       const FieldEvaluator& fieldEval,
                       const stk::mesh::EntityRank fieldRank = stk::topology::NODE_RANK,
                       const Ioss::Field::RoleType fieldRole = Ioss::Field::TRANSIENT,
                       const std::string generatedMeshString = "generated:1x1x1|nodeset:xyz|sideset:xyz",
                       const int numCopies = 1,
                       const int spatialDim = 3);

inline
std::shared_ptr<stk::mesh::BulkData>
create_mesh_with_field_2d(const stk::ParallelMachine comm,
                       const std::string& fieldName = "field_1",
                       const double fieldValue = 0.0,
                       const stk::mesh::EntityRank fieldRank = stk::topology::NODE_RANK,
                       const Ioss::Field::RoleType fieldRole = Ioss::Field::TRANSIENT,
                       const std::string generatedMeshString = "textmesh:0,1,QUAD_4_2D,1,2,3,4 |coordinates: 0,0, 1,0, 1,1, 0,1 |dimension:2",
                       const int numCopies = 1)
{
  return create_mesh_with_field(comm, fieldName, fieldValue, fieldRank, 
                                      fieldRole, generatedMeshString, numCopies, 2);
}

inline
std::shared_ptr<stk::mesh::BulkData>
create_mesh_with_field_2d(const stk::ParallelMachine comm,
                       const stk::mesh::EntityRank fieldRank,
                       const std::string generatedMeshString,
                       const int numCopies = 1)
{
  return create_mesh_with_field_2d(comm, "field_1", 0.0, fieldRank, 
                                Ioss::Field::TRANSIENT, generatedMeshString, numCopies);
}

void scale_to_unit_bbox(stk::mesh::BulkData& mesh);

void add_field_to_mesh(stk::mesh::BulkData& mesh,
                       const std::string& fieldName,
                       const FieldEvaluator& fieldEval,
                       const stk::mesh::EntityRank fieldRank = stk::topology::NODE_RANK,
                       const Ioss::Field::RoleType fieldRole = Ioss::Field::TRANSIENT);
}
}

#endif

