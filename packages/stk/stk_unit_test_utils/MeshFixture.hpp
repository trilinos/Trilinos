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
    : comm(MPI_COMM_WORLD), meta(), bulk(nullptr)
    {

    }

    virtual ~MeshFixture()
    {
        delete bulk;
    }

    void setup_mesh(const std::string &meshSpecification, stk::mesh::BulkData::AutomaticAuraOption auraOption)
    {
        bulk = new stk::mesh::BulkData(meta, comm, auraOption);
        stk::unit_test_util::fill_mesh_using_stk_io(meshSpecification, *bulk, comm);
    }

    stk::mesh::BulkData& bulkData()
    {
        ThrowRequireMsg(bulk!=nullptr, "Unit test error. Trying to get bulk data before it has been initialized.");
        return *bulk;
    }

    MPI_Comm comm;
    stk::mesh::MetaData meta;
    stk::mesh::BulkData *bulk;
};
