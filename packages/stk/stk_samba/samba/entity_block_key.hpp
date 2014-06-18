#ifndef SAMBA_SAMBA_ENTITY_SET_KEY_HPP
#define SAMBA_SAMBA_ENTITY_SET_KEY_HPP

#include <samba/utility.hpp>

namespace samba {

/**
 * An ordinal for application-requested set of entities.
 *
 * In older mesh implementations, this was called a "part ordinal"
 */
SAMBA_MAKE_PRIMITIVE_WITH_HASHABLE_TAG(entity_block_key,entity_block,uint32_t,3);

} // namespace samba

SAMBA_IS_PRIMITIVE(samba::entity_block_key)

#endif //SAMBA_SAMBA_ENTITY_SET_KEY_HPP
