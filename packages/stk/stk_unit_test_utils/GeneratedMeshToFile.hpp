#ifndef GENERATED_MESH_TO_FILE_H_
#define GENERATED_MESH_TO_FILE_H_

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/CoordinateSystems.hpp>
#include <stk_io/StkMeshIoBroker.hpp>

#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldBase.hpp>
#include <mpi.h>

#include "ioUtils.hpp"

namespace stk
{
namespace unit_test_util
{

class GeneratedMeshToFile
{
public:
    GeneratedMeshToFile(stk::ParallelMachine comm, stk::mesh::BulkData::AutomaticAuraOption auraOption);

    void setup_mesh(const std::string &meshSizeSpec, const std::string &outputFileName);
    void write_mesh();

protected:
    stk::mesh::MetaData meta;
    stk::mesh::BulkData bulk;
    stk::io::StkMeshIoBroker broker;
    size_t outputFileIndex = 0;

private:
    GeneratedMeshToFile();
};

class GeneratedMeshToFileWithTransientFields : public GeneratedMeshToFile
{
public:
    GeneratedMeshToFileWithTransientFields(stk::ParallelMachine comm,
                                          stk::mesh::BulkData::AutomaticAuraOption auraOption,
                                          const std::string& fieldBaseName,
                                          stk::topology::rank_t rank);

    virtual ~GeneratedMeshToFileWithTransientFields() {};

    void write_mesh_with_field(const std::vector<double>& timeSteps,
                               const FieldValueSetter &fieldValueSetter,
                               const std::string& globalVariableName);

protected:
    stk::topology::rank_t fieldRank;
    stk::mesh::Field<double> &scalarField;
    stk::mesh::Field<double, stk::mesh::Cartesian> &vectorField;

private:
    GeneratedMeshToFileWithTransientFields();
};

}
}

#endif
