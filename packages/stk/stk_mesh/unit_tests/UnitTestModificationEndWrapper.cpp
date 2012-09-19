
#include <unit_tests/UnitTestModificationEndWrapper.hpp>

namespace stk {
namespace mesh {

bool UnitTestModificationEndWrapper::wrap(stk::mesh::BulkData& mesh, bool generate_aura)
{
  return mesh.internal_modification_end(generate_aura);
}

} // namespace mesh

namespace unit_test {

bool modification_end_wrapper(stk::mesh::BulkData& mesh, bool generate_aura)
{
  return stk::mesh::UnitTestModificationEndWrapper::wrap(mesh, generate_aura);
}

} // namespace unit_test
} // namespace stk
