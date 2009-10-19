#ifndef stk_mesh_diag_EntityKey_hpp
#define stk_mesh_diag_EntityKey_hpp

#include <stk_util/diag/Writer.hpp>
#include <stk_mesh/base/EntityKey.hpp>


namespace stk {
namespace mesh {

inline
stk::diag::Writer &operator<<(stk::diag::Writer &dout, const EntityKey &entity_key)  {
  return dout << entity_type(entity_key) << ":" << entity_id(entity_key);
}

} // namespace stk
} // namespace mesh

#endif // stk_mesh_diag_EntityKey_hpp
