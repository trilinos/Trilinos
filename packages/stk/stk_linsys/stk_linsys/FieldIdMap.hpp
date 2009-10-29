#ifndef stk_linsys_FieldIdMap_hpp
#define stk_linsys_FieldIdMap_hpp

#include <map>
#include <stk_mesh/base/Field.hpp>

namespace stk {
namespace linsys {

/** Mappings from stk::mesh::Field objects to integer ids used by fei objects.
 */
typedef std::map<const stk::mesh::FieldBase*,int> FieldIdMap;

}//namespace linsys
}//namespace stk

#endif
