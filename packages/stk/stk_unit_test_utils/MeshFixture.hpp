#include <mpi.h>
#include <gtest/gtest.h>
#include <stk_mesh/base/MetaData.hpp>   // for MetaData, put_field
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_unit_test_utils/ioUtils.hpp>
#include <stk_util/environment/ReportHandler.hpp>

class MeshFixture : public ::testing::Test
{
protected:
    MeshFixture()
    : communicator(MPI_COMM_WORLD), metaData(), bulkData(nullptr)
    {

    }

    virtual ~MeshFixture()
    {
        delete bulkData;
    }

    void setup_mesh(const std::string &meshSpecification, stk::mesh::BulkData::AutomaticAuraOption auraOption)
    {
        bulkData = new stk::mesh::BulkData(metaData, communicator, auraOption);
        stk::unit_test_util::fill_mesh_using_stk_io(meshSpecification, *bulkData, communicator);
    }

    MPI_Comm comm()
    {
        return communicator;
    }

    stk::mesh::MetaData& meta()
    {
        return metaData;
    }

    stk::mesh::BulkData& bulk()
    {
        ThrowRequireMsg(bulkData!=nullptr, "Unit test error. Trying to get bulk data before it has been initialized.");
        return *bulkData;
    }

private:
    MPI_Comm communicator;
    stk::mesh::MetaData metaData;
    stk::mesh::BulkData *bulkData;
};
