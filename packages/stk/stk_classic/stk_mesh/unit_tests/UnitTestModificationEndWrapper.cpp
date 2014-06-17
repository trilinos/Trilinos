
#include <unit_tests/UnitTestModificationEndWrapper.hpp>

namespace stk_classic {
namespace mesh {

bool UnitTestModificationEndWrapper::wrap(stk_classic::mesh::BulkData& mesh, bool generate_aura)
{
  return mesh.internal_modification_end(generate_aura);
}

} // namespace mesh

namespace unit_test {

bool modification_end_wrapper(stk_classic::mesh::BulkData& mesh, bool generate_aura)
{
  return stk_classic::mesh::UnitTestModificationEndWrapper::wrap(mesh, generate_aura);
}

} // namespace unit_test
} // namespace stk_classic
