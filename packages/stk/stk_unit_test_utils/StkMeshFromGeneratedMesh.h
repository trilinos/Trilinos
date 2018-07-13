#ifndef STKMESHFROMGENERATEDMESH_H_
#define STKMESHFROMGENERATEDMESH_H_

#include <stk_mesh/base/BulkData.hpp>
#include <stk_unit_test_utils/BulkDataTester.hpp>
#include <stk_unit_test_utils/ioUtils.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_io/StkMeshIoBroker.hpp>
#include "mpi.h"

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
        m_stkMeshBulkData = new stk::unit_test_util::BulkDataTester(*m_stkMeshMetaData, communicator);

        readExodusFileIntoStkMesh(generatedMeshSpec, *m_stkMeshBulkData, communicator);
    }

    ~StkMeshCreator()
    {
        delete m_stkMeshBulkData;
        delete m_stkMeshMetaData;
    }

    stk::mesh::MetaData* getMetaData() { return m_stkMeshMetaData; }
    stk::unit_test_util::BulkDataTester* getBulkData() { return m_stkMeshBulkData; }

private:

    void readExodusFileIntoStkMesh(const std::string& generatedMeshSpecification, stk::mesh::BulkData &stkMeshBulkData, MPI_Comm communicator)
    {
        stk::io::fill_mesh(generatedMeshSpecification, stkMeshBulkData);
    }

private:
    stk::mesh::MetaData *m_stkMeshMetaData;
    stk::unit_test_util::BulkDataTester *m_stkMeshBulkData;
};

} // namespace exampleMeshes
} // namespace unitTestUtils

#endif
