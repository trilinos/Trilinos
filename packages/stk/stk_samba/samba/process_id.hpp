#ifndef SAMBA_SAMBA_PROCESS_ID_HPP
#define SAMBA_SAMBA_PROCESS_ID_HPP

#include <samba/utility.hpp>

namespace samba {

/**
 * Represents a process identifier. In MPI, this would be the comm rank.
 * This is a primitive data type.
 *
 * Process-ids are limited to a certain number of bits so that they can
 * be incoporated into entity_keys without taking up too much space in the key.
 */
SAMBA_MAKE_BIT_RESTRICTED_PRIMITIVE(process_id, process_id, int32_t, 24);

} // namespace samba

SAMBA_IS_PRIMITIVE(samba::process_id)

#endif //SAMBA_SAMBA_PROCESS_ID_HPP
