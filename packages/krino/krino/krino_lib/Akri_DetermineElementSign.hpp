#ifndef KRINO_KRINO_KRINO_LIB_AKRI_DETERMINEELEMENTSIGN_HPP_
#define KRINO_KRINO_KRINO_LIB_AKRI_DETERMINEELEMENTSIGN_HPP_

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/Selector.hpp>

namespace krino { class Surface; }
namespace krino { class FieldRef; }

namespace krino {
typedef std::map<stk::mesh::Entity, std::vector<int8_t>> ElementToSignsMap;

ElementToSignsMap determine_element_signs(const stk::mesh::BulkData & mesh,
    const FieldRef coordsField,
    const stk::mesh::Selector & elementSelector,
    const std::vector<stk::mesh::Selector> & perSurfaceElementSelector,
    const std::vector<const Surface*> & surfaces);
}

#endif /* KRINO_KRINO_KRINO_LIB_AKRI_DETERMINEELEMENTSIGN_HPP_ */
