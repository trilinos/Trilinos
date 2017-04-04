#ifndef ioUtils_hpp
#define ioUtils_hpp

#include <string>
#include <vector>

#include <stk_mesh/base/Types.hpp>
#include <stk_util/parallel/Parallel.hpp>  // for ParallelMachine, etc
#include <stk_io/DatabasePurpose.hpp> // for stk::io::DatabasePurpose
#include <stk_io/FillMesh.hpp>
#include <stk_io/WriteMesh.hpp>
namespace stk { namespace mesh { class BulkData; }}
namespace stk { namespace mesh { class FieldBase; }}
namespace stk { namespace io { class StkMeshIoBroker; } }

namespace stk
{
namespace unit_test_util
{

// Example:  meshSizeSpec = "2x2x1"
void generated_mesh_to_file_in_serial(const std::string &meshSizeSpec, const std::string &fileName);

void read_from_serial_file_and_decompose(const std::string& fileName, stk::mesh::BulkData &mesh, const std::string &decompositionMethod);

// This avoid the GeneratedMesh limitation on the z-dimension >= num_processors,
// and allows cases like 2x2x1 with cyclic decomposition on four processors.
void generate_mesh_from_serial_spec_and_load_in_parallel_with_auto_decomp(const std::string &meshSizeSpec, stk::mesh::BulkData & mesh, const std::string &decompositionMethod);

class FieldValueSetter
{
public:
    virtual ~FieldValueSetter() {}
    virtual void populate_field(stk::mesh::BulkData &bulk, stk::mesh::FieldBase* field, const unsigned step, const double time) const = 0;
};

class IdAndTimeFieldValueSetter : public FieldValueSetter
{
public:
    virtual void populate_field(stk::mesh::BulkData &bulk, stk::mesh::FieldBase* field, const unsigned step, const double time) const;
};

void generated_mesh_with_transient_data_to_file_in_serial(const std::string &meshSizeSpec,
                                                          const std::string &fileName,
                                                          const std::string& fieldName,
                                                          const std::string& globalVariableName,
                                                          const std::vector<double>& timeSteps,
                                                          const FieldValueSetter &fieldValueSetter);

} // namespace unit_test_util
} // namespace stk

#endif
