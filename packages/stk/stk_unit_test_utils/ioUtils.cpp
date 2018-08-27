// #######################  Start Clang Header Tool Managed Headers ########################
// clang-format off
#include "ioUtils.hpp"
#include <stddef.h>                                  // for size_t
#include <unistd.h>                                  // for unlink
#include <stk_io/StkIoUtils.hpp>
#include <string>                                    // for allocator, etc
#include <vector>                                    // for vector
#include "GeneratedMeshToFile.hpp"
#include "Ioss_Property.h"                           // for Property
#include "mpi.h"                                     // for MPI_COMM_SELF, etc
#include "stk_io/DatabasePurpose.hpp"
#include "stk_io/StkMeshIoBroker.hpp"                // for StkMeshIoBroker
#include "stk_mesh/base/BulkData.hpp"                // for BulkData, etc
#include "stk_mesh/base/BulkDataInlinedMethods.hpp"
#include "stk_mesh/base/Entity.hpp"                  // for Entity
#include "stk_mesh/base/FieldBase.hpp"               // for FieldBase, etc
#include "stk_mesh/base/GetEntities.hpp"             // for get_entities
#include "stk_mesh/base/Types.hpp"                   // for EntityRank, etc
#include "stk_topology/topology.hpp"                 // for topology, etc
#include "stk_util/parallel/Parallel.hpp"
// clang-format on
// #######################   End Clang Header Tool Managed Headers  ########################

namespace stk
{
namespace unit_test_util
{

void generated_mesh_to_file_in_serial(const std::string &meshSizeSpec, const std::string &fileName)
{
    if (stk::parallel_machine_rank(MPI_COMM_WORLD) == 0)
    {
        GeneratedMeshToFile gMesh(MPI_COMM_SELF, stk::mesh::BulkData::NO_AUTO_AURA);

        gMesh.setup_mesh(meshSizeSpec, fileName);
        gMesh.write_mesh();
    }
}

void IdAndTimeFieldValueSetter::populate_field(stk::mesh::BulkData &bulk, stk::mesh::FieldBase* field, const unsigned step, const double time) const
{
    stk::mesh::EntityRank fieldRank = field->entity_rank();

    std::vector<stk::mesh::Entity> entities;
    stk::mesh::get_entities(bulk, fieldRank, entities);

    stk::mesh::FieldVector allTransientFields = stk::io::get_transient_fields(bulk.mesh_meta_data());

    for(stk::mesh::FieldBase * transientField : allTransientFields)
    {
        for(size_t i = 0; i < entities.size(); i++)
        {
            unsigned numEntriesPerEntity = stk::mesh::field_scalars_per_entity(*transientField, entities[i]);
            double value = 100.0 * time + static_cast<double>(bulk.identifier(entities[i]));
            double *data = static_cast<double*> (stk::mesh::field_data(*transientField, entities[i]));
            for(unsigned j=0; j<numEntriesPerEntity; j++)
                data[j] = value + j;
        }
    }
}

void generated_mesh_with_transient_data_to_file_in_serial(const std::string &meshSizeSpec,
                                                          const std::string &fileName,
                                                          const std::string& fieldName,
                                                          const std::string& globalVariableName,
                                                          const std::vector<double>& timeSteps,
                                                          const FieldValueSetter &fieldValueSetter)
{
    if (stk::parallel_machine_rank(MPI_COMM_WORLD) == 0)
    {
        GeneratedMeshToFileWithTransientFields gMesh(MPI_COMM_SELF, stk::mesh::BulkData::NO_AUTO_AURA, fieldName, stk::topology::NODE_RANK);

        gMesh.setup_mesh(meshSizeSpec, fileName);
        gMesh.write_mesh_with_field(timeSteps, fieldValueSetter, globalVariableName);
    }
}

void read_from_serial_file_and_decompose(const std::string& fileName, stk::mesh::BulkData &mesh, const std::string &decompositionMethod)
{
    stk::io::StkMeshIoBroker broker;
    broker.set_bulk_data(mesh);
    broker.property_add(Ioss::Property("DECOMPOSITION_METHOD", decompositionMethod));
    broker.add_mesh_database(fileName, stk::io::READ_MESH);
    broker.create_input_mesh();
    broker.populate_bulk_data();
}

void generate_mesh_from_serial_spec_and_load_in_parallel_with_auto_decomp(const std::string &meshSizeSpec, stk::mesh::BulkData &mesh, const std::string &decompositionMethod)
{
    // meshSizeSpec should NOT include generated:, just "2x2x1" for example.
    // decomposition methods: "linear", "rcb", "rib", "hsfc", "block", "cyclic", "random", "kway", "geom_kway", "metis_sfc"
    const std::string tempFilename = "exodus_" + meshSizeSpec + ".e";
    generated_mesh_to_file_in_serial(meshSizeSpec,tempFilename);

    read_from_serial_file_and_decompose(tempFilename, mesh, decompositionMethod);
    unlink(tempFilename.c_str());
}


} // namespace unit_test_util
} // namespace stk

