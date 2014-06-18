#ifndef SAMBA_SAMBA_TYPES_HPP
#define SAMBA_SAMBA_TYPES_HPP

#include <samba/partition_index.hpp>
#include <samba/entity_key.hpp>
#include <samba/entity_part.hpp>
#include <samba/connectivity_kind.hpp>
#include <samba/connectivity_map.hpp>
#include <samba/connectivity_ordinal.hpp>
#include <samba/connectivity_orientation.hpp>
#include <samba/spatial_dimension.hpp>

#include <boost/range.hpp>

#include <vector>

namespace samba {

//*****************************************************************************
//part vector and range
//*****************************************************************************
typedef std::vector<entity_part>  entity_part_vector;
typedef boost::sub_range<const entity_part_vector> entity_part_range;

//*****************************************************************************
//entity_key_interval
//*****************************************************************************
typedef interval<entity_key> entity_key_interval;

typedef const entity_key*                                                         entity_key_iterator;
typedef entity_key*                                                               non_const_entity_key_iterator;
typedef std::pair<entity_key_iterator, entity_key_iterator>                       entity_key_range;
typedef std::pair<non_const_entity_key_iterator, non_const_entity_key_iterator>   non_const_entity_key_range;

typedef const partition_index* partition_index_iterator;
typedef partition_index*       non_const_partition_index_iterator;

typedef std::pair<partition_index_iterator, partition_index_iterator>                        partition_index_range;
typedef std::pair<non_const_partition_index_iterator, non_const_partition_index_iterator>    non_const_partition_index_range;


typedef const connectivity_ordinal* ordinal_iterator;
typedef connectivity_ordinal*       non_const_ordinal_iterator;

typedef std::pair<ordinal_iterator, ordinal_iterator>                      ordinal_range;
typedef std::pair<non_const_ordinal_iterator, non_const_ordinal_iterator>  non_const_ordinal_range;


typedef const connectivity_orientation* orientation_iterator;
typedef connectivity_orientation*       non_const_orientation_iterator;

typedef std::pair<orientation_iterator, orientation_iterator>                      orientation_range;
typedef std::pair<non_const_orientation_iterator, non_const_orientation_iterator>  non_const_orientation_range;

} //namespace samba

#endif //SAMBA_SAMBA_TYPES_HPP
