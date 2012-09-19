#ifndef SAMBA_SAMBA_PARTITION_DESCRIPTOR_HPP
#define SAMBA_SAMBA_PARTITION_DESCRIPTOR_HPP

#include <samba/utility.hpp>

namespace samba {

/**
 * A unique identifier for a partition. Identifiers for partitions *can*
 * change when the mesh is modified.
 */
SAMBA_MAKE_BIT_RESTRICTED_PRIMITIVE(partition_id, partition_id, uint32_t, 24);

} // namespace samba

SAMBA_IS_PRIMITIVE(samba::partition_id)

#endif //SAMBA_SAMBA_PARTITION_DESCRIPTOR_HPP
