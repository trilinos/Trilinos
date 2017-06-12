#ifndef UNITTEST_MESHFILEFIXTURE_HPP
#define UNITTEST_MESHFILEFIXTURE_HPP

#include <mpi.h>

#include <gtest/gtest.h>

#include <stk_io/StkMeshIoBroker.hpp>
#include <stk_mesh/base/MetaData.hpp>   // for MetaData, put_field
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_unit_test_utils/ioUtils.hpp>
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_util/environment/ReportHandler.hpp>

namespace stk
{

namespace unit_test_util
{

class MeshFileFixture : public MeshFixture
{
protected:
    MeshFileFixture()
    : stkIo(get_comm())
    {
        allocate_bulk(stk::mesh::BulkData::NO_AUTO_AURA);
    }
    void read_mesh(const std::string &fileToRead)
    {
        stk::io::fill_mesh_preexisting(stkIo, fileToRead, get_bulk());
    }
    void TearDown()
    {
        unlink(filename.c_str());
    }
    const std::string filename = "filename.exo";
    stk::io::StkMeshIoBroker stkIo;
};

}}

#endif

