#ifndef STKMESHFROMGENERATEDMESH_H_
#define STKMESHFROMGENERATEDMESH_H_

#include <stk_mesh/base/BulkData.hpp>
#include <unit_tests/BulkDataTester.hpp>
#include <stk_unit_test_utils/ioUtils.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_io/StkMeshIoBroker.hpp>
#include <mpi.h>

namespace unitTestUtils
{

namespace exampleMeshes
{


class StkMeshCreator
{

public:
    StkMeshCreator(const std::string& generatedMeshSpec, MPI_Comm communicator) :
    m_stkMeshMetaData(NULL), m_stkMeshBulkData(NULL)
    {
        const int spatialDim = 3;
        m_stkMeshMetaData = new stk::mesh::MetaData(spatialDim);
        m_stkMeshBulkData = new stk::mesh::unit_test::BulkDataTester(*m_stkMeshMetaData, communicator);

        readExodusFileIntoStkMesh(generatedMeshSpec, *m_stkMeshBulkData, communicator);
    }

    ~StkMeshCreator()
    {
        delete m_stkMeshBulkData;
        delete m_stkMeshMetaData;
    }

    stk::mesh::MetaData* getMetaData() { return m_stkMeshMetaData; }
    stk::mesh::unit_test::BulkDataTester* getBulkData() { return m_stkMeshBulkData; }

private:

    void readExodusFileIntoStkMesh(const std::string& generatedMeshSpecification, stk::mesh::BulkData &stkMeshBulkData, MPI_Comm communicator)
    {
        stk::unit_test_util::fill_mesh_using_stk_io(generatedMeshSpecification, stkMeshBulkData, communicator);
    }

private:
    stk::mesh::MetaData *m_stkMeshMetaData;
    stk::mesh::unit_test::BulkDataTester *m_stkMeshBulkData;
};

} // namespace exampleMeshes
} // namespace unitTestUtils

#endif
