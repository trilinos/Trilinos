#ifndef UNITTEST_MESHFIXTURE_HPP
#define UNITTEST_MESHFIXTURE_HPP

#include <mpi.h>
#include <gtest/gtest.h>
#include <stk_mesh/base/MetaData.hpp>   // for MetaData, put_field
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_unit_test_utils/ioUtils.hpp>
#include <stk_util/environment/ReportHandler.hpp>

namespace stk
{
namespace unit_test_util
{

class MeshFixture : public ::testing::Test
{
protected:
    MeshFixture()
    : communicator(MPI_COMM_WORLD), metaData(nullptr), bulkData(nullptr)
    {
        allocate_meta();
    }

    MeshFixture(unsigned spatial_dim)
    : communicator(MPI_COMM_WORLD), metaData(nullptr), bulkData(nullptr)
    {
        allocate_meta(spatial_dim);
    }

    virtual ~MeshFixture()
    {
        delete bulkData;
        delete metaData;
    }

    void setup_empty_mesh(stk::mesh::BulkData::AutomaticAuraOption auraOption)
    {
        allocate_bulk(auraOption);
    }
    void setup_mesh(const std::string &meshSpecification, stk::mesh::BulkData::AutomaticAuraOption auraOption)
    {
        allocate_bulk(auraOption);
        stk::unit_test_util::fill_mesh_using_stk_io(meshSpecification, *bulkData);
    }

    void setup_mesh_with_cyclic_decomp(const std::string &meshSpecification, stk::mesh::BulkData::AutomaticAuraOption auraOption)
    {
        allocate_bulk(auraOption);
        stk::unit_test_util::generate_mesh_from_serial_spec_and_load_in_parallel_with_auto_decomp(meshSpecification,*bulkData,"cyclic");
    }

    MPI_Comm get_comm()
    {
        return communicator;
    }

    void initialize_mesh()
    {
        delete bulkData;
        delete metaData;

        bulkData = nullptr;
        metaData = nullptr;
    }

    virtual stk::mesh::MetaData& get_meta()
    {
        ThrowRequireMsg(metaData!=nullptr, "Unit test error. Trying to get meta data before it has been initialized.");
        return *metaData;
    }

    stk::mesh::BulkData& get_bulk()
    {
        ThrowRequireMsg(bulkData!=nullptr, "Unit test error. Trying to get bulk data before it has been initialized.");
        return *bulkData;
    }

    virtual void allocate_bulk(stk::mesh::BulkData::AutomaticAuraOption auraOption)
    {
        if(nullptr == metaData)
            allocate_meta();

        bulkData = new stk::mesh::BulkData(get_meta(), communicator, auraOption);
    }

    virtual void allocate_meta(unsigned spatial_dim = 3)
    {
        ThrowRequireMsg(metaData==nullptr, "Unit test error. Trying to reset non NULL meta data.");
        metaData = new stk::mesh::MetaData(spatial_dim);
    }

    void set_bulk(stk::mesh::BulkData *inBulkData)
    {
        ThrowRequireMsg(bulkData==nullptr, "Unit test error. Trying to reset non NULL bulk data.");
        bulkData = inBulkData;
    }

    void delete_meta()
    {
        ThrowRequireMsg(bulkData==nullptr, "Unit test error. Trying to delete meta with non NULL bulk data.");
        delete metaData;
        metaData = nullptr;
    }

private:
    MPI_Comm communicator;
    stk::mesh::MetaData *metaData = nullptr;
    stk::mesh::BulkData *bulkData = nullptr;
};

class MeshTestFixture : public MeshFixture
{
protected:
    void run_test_on_num_procs(int numProcs, stk::mesh::BulkData::AutomaticAuraOption auraOption)
    {
        if(stk::parallel_machine_size(get_comm()) == numProcs)
        {
            run_test(auraOption);
        }
    }

    void run_test_on_num_procs_or_less(int numProcs, stk::mesh::BulkData::AutomaticAuraOption auraOption)
    {
        if(stk::parallel_machine_size(get_comm()) <= numProcs)
        {
            run_test(auraOption);
        }
    }

    virtual void run_test(stk::mesh::BulkData::AutomaticAuraOption auraOption) = 0;
};

}}

#endif

