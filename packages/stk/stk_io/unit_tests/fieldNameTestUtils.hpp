#include <mpi.h>
#include <stk_io/MeshReadWriteUtils.hpp>
#include <Ioss_SubSystem.h>

inline stk::mesh::MetaData& generateMetaData(stk::io::MeshData &stkIo)
{
    const std::string exodusFileName = "generated:1x1x1";
    stkIo.open_mesh_database(exodusFileName);
    stkIo.create_input_mesh();
    return stkIo.meta_data();
}

inline bool fieldWithNameChangedIsOutput(stk::io::MeshData &stkIo, MPI_Comm communicator, const std::string &goldFieldName)
{
    const std::string outputFileName = "resultsOutput.exo";
    stkIo.create_output_mesh(outputFileName);
    stkIo.define_output_fields();

    Ioss::Region *outputRegion = stkIo.output_io_region().get();
    Ioss::NodeBlock *nodeBlockAssociatedWithDisplacementField = outputRegion->get_node_blocks()[0];
    Ioss::NameList fieldNames;
    nodeBlockAssociatedWithDisplacementField->field_describe(Ioss::Field::TRANSIENT, &fieldNames);

    unlink(outputFileName.c_str());

    return (goldFieldName == fieldNames[0]);
}
