#ifndef stk_search_util_PrintEntityProc_hpp
#define stk_search_util_PrintEntityProc_hpp

#include <string>
#include <vector>

#include <stk_mesh/base/Comm.hpp>
#include <stk_util/diag/Writer_fwd.hpp>
#include <stk_search/CoarseSearch.hpp>
#include <stk_search/IdentProc.hpp>

#include <stk_search_util/stk_mesh/CreateBoundingBox.hpp>


namespace stk {
namespace mesh {
class Entity;
}

namespace search_util {
/**
 * Used to output the results of a coarse or direct search to
 * verify which entity contains another entity.
 */
void print_entity_map(stk::diag::Writer &writer,
                      const std::vector<std::pair<stk::mesh::Entity*,
                      stk::mesh::Entity*> >& entity_map,
                      const std::string & relation);

/**
 * Used to output a sharing or ghosting vector in human readable
 * form.
 * The "action" argument will typically be "Share " or "Ghost "
 * The "to_from" argument will typically be
 * - for sharing " with "
 * - for ghosting " from " or " to "
 *
 * Decodes the entity key and prints as entity type, entity id
 *
 * Example output: "Share NODE 37 with processor 12"
 */
void print_entity_proc_map(stk::diag::Writer &writer,
                           const std::vector<stk::mesh::EntityProc>& entity_proc,
                           const std::string &action,
                           const std::string &to_from);

/**
 * Used to output a sharing or ghosting vector in human readable
 * form using the owner processor.
 * The "action" argument will typically be "Share " or "Ghost "
 * The "to_from" argument will typically be
 * - for sharing " with "
 * - for ghosting " from " or " to "
 *
 * Decodes the entity key and prints as entity type, entity id
 *
 * Example output: "Share NODE 37 with processor 12"
 */
void print_entity_proc_map(stk::diag::Writer &writer,
                           const std::vector<stk::mesh::Entity*>& entity_proc,
                           const std::string &action,
                           const std::string &to_from);

/**
 * Used to output the results of a relation vector in human
 * readable form.  This function cannot be used as the default
 * output of an IdentProcRelation since it is using knowledge that
 * what is really being stored in the IdentProc is stk::mesh
 * entity keys.
 */
void print_stk_mesh_relation_map(stk::diag::Writer &writer,
                                 IdentProcRelation relation);
} // namespace search_util
} // namespace stk
#endif
