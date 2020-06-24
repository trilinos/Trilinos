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

#include <gtest/gtest.h>
#include <stddef.h>
#include <stk_io/StkMeshIoBroker.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_topology/topology.hpp>
#include <string>
#include "stk_mesh/base/Types.hpp"
#include "stk_mesh/base/GetEntities.hpp"
#include "stk_mesh/base/Field.hpp"
#include "stk_transfer/copy_by_id/TransferCopyById.hpp"
#include "stk_transfer/copy_by_id/TransferCopyByIdStkMeshAdapter.hpp"
#include "stk_transfer/copy_by_id/SearchByIdGeometric.hpp"

namespace
{
//BEGIN
TEST(StkMeshHowTo, useCopyTransfer)
{
    MPI_Comm communicator = MPI_COMM_WORLD;
    if (stk::parallel_machine_size(communicator) != 2) { return; }
    // syntax creates a 2x2x4 'bar' of hex-8 elements
    const std::string generatedMeshSpecification = "generated:2x2x4";

    double init_vals = std::numeric_limits<double>::max();
    typedef stk::mesh::Field<double> ScalarField;

    // Set up MeshA
    stk::io::StkMeshIoBroker meshReaderA(communicator);
    meshReaderA.add_mesh_database(generatedMeshSpecification, stk::io::READ_MESH);
    meshReaderA.create_input_mesh();
    stk::mesh::MetaData & metaA = meshReaderA.meta_data();
    ScalarField & scalarFieldNodeA = metaA.declare_field<ScalarField>(stk::topology::NODE_RANK, "Node Scalar Field");
    stk::mesh::put_field_on_mesh(scalarFieldNodeA, metaA.universal_part(), &init_vals);
    meshReaderA.populate_bulk_data();
    stk::mesh::BulkData & meshA = meshReaderA.bulk_data();

    // Set up MeshB
    stk::io::StkMeshIoBroker meshReaderB(communicator);
    meshReaderB.add_mesh_database(generatedMeshSpecification, stk::io::READ_MESH);
    meshReaderB.create_input_mesh();
    stk::mesh::MetaData & metaB = meshReaderB.meta_data();
    ScalarField & scalarFieldNodeB = metaB.declare_field<ScalarField>(stk::topology::NODE_RANK, "Node Scalar Field");
    stk::mesh::put_field_on_mesh(scalarFieldNodeB, metaB.universal_part(), &init_vals);
    meshReaderB.populate_bulk_data();
    stk::mesh::BulkData & meshB = meshReaderB.bulk_data();

    // Change parallel decomposition of MeshB
    stk::mesh::EntityVector locally_owned_elements;
    stk::mesh::get_selected_entities(metaB.locally_owned_part(),
                                     meshB.buckets(stk::topology::ELEMENT_RANK),
                                     locally_owned_elements);
    stk::mesh::EntityVector locally_owned_nodes;
    stk::mesh::get_selected_entities(metaB.locally_owned_part(),
                                     meshB.buckets(stk::topology::NODE_RANK),
                                     locally_owned_nodes);

    int myProc = meshA.parallel_rank();
    int otherProc = 1-myProc;

    std::vector<std::pair<stk::mesh::Entity, int> > entityProcPairs;
    entityProcPairs.reserve(locally_owned_elements.size()+locally_owned_nodes.size());
    for (size_t i=0 ; i<locally_owned_elements.size() ; ++i) {
        entityProcPairs.push_back(std::make_pair(locally_owned_elements[i], otherProc));
    }
    for (size_t i=0 ; i<locally_owned_nodes.size() ; ++i) {
        entityProcPairs.push_back(std::make_pair(locally_owned_nodes[i], otherProc));
    }
    meshB.change_entity_owner(entityProcPairs);

    // Fill nodal field on meshA with identifier of each node.
    {
        const stk::mesh::BucketVector & entityBuckets =
                meshA.get_buckets(stk::topology::NODE_RANK,metaA.locally_owned_part());
        for (size_t bucketIndex = 0; bucketIndex < entityBuckets.size(); ++bucketIndex) {
            stk::mesh::Bucket & entityBucket = * entityBuckets[bucketIndex];
            for (size_t entityIndex = 0; entityIndex < entityBucket.size(); ++entityIndex) {
                stk::mesh::Entity entity = entityBucket[entityIndex];
                double * scalar = stk::mesh::field_data(scalarFieldNodeA, entity);
                *scalar = static_cast<double>(meshA.identifier(entity));
            }
        }
    }

    // Set up CopyTransfer
    stk::mesh::EntityVector entitiesA;
    stk::mesh::get_selected_entities(metaA.locally_owned_part(),
                                     meshA.buckets(stk::topology::NODE_RANK),
                                     entitiesA);
    std::vector<stk::mesh::FieldBase*> fieldsA;
    fieldsA.push_back(&scalarFieldNodeA);
    stk::transfer::TransferCopyByIdStkMeshAdapter transferMeshA(meshA,entitiesA,fieldsA);

    stk::mesh::EntityVector entitiesB;
    stk::mesh::get_selected_entities(metaB.locally_owned_part(),
                                     meshB.buckets(stk::topology::NODE_RANK),
                                     entitiesB);
    std::vector<stk::mesh::FieldBase*> fieldsB;
    fieldsB.push_back(&scalarFieldNodeB);
    stk::transfer::TransferCopyByIdStkMeshAdapter transferMeshB(meshB,entitiesB,fieldsB);

    stk::transfer::SearchByIdGeometric copySearch;

    stk::transfer::TransferCopyById copyTransfer(copySearch,transferMeshA,transferMeshB);
    copyTransfer.initialize();

    // Apply CopyTransfer
    copyTransfer.apply();

    // Verify nodal fields on meshB are correct
    {
        const double tolerance = 1.0e-8;
        const stk::mesh::BucketVector & entityBuckets =
                meshB.get_buckets(stk::topology::NODE_RANK,metaB.locally_owned_part());
        for (size_t bucketIndex = 0; bucketIndex < entityBuckets.size(); ++bucketIndex) {
            stk::mesh::Bucket & entityBucket = * entityBuckets[bucketIndex];
            for (size_t entityIndex = 0; entityIndex < entityBucket.size(); ++entityIndex) {
                stk::mesh::Entity entity = entityBucket[entityIndex];
                double * scalar = stk::mesh::field_data(scalarFieldNodeB, entity);
                EXPECT_NEAR( static_cast<double>(meshB.identifier(entity)), *scalar, tolerance);
            }
        }
    }

}
//END

} // namespace
