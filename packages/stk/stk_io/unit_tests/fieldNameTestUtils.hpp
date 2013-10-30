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

inline bool fieldWithNameChangedIsOutput(stk::io::MeshData &stkIo, MPI_Comm communicator, const size_t resultsOutputIndex, const std::string &goldFieldName)
{
    double dummyTime = 0;
    stkIo.process_output_request(dummyTime, resultsOutputIndex);
    Ioss::Region *outputRegion = stkIo.get_output_io_region(resultsOutputIndex).get();
    Ioss::NodeBlock *nodeBlockAssociatedWithDisplacementField = outputRegion->get_node_blocks()[0];
    Ioss::NameList fieldNames;
    nodeBlockAssociatedWithDisplacementField->field_describe(Ioss::Field::TRANSIENT, &fieldNames);

    return (goldFieldName == fieldNames[0]);
}
