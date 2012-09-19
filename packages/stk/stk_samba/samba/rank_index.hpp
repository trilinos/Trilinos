#ifndef SAMBA_SAMBA_RANK_INDEX_HPP
#define SAMBA_SAMBA_RANK_INDEX_HPP

#include <samba/utility.hpp>
#include <samba/entity_rank.hpp>

namespace samba {

SAMBA_MAKE_PRIMITIVE(node_index,    node_index,    uint32_t);
SAMBA_MAKE_PRIMITIVE(edge_index,    edge_index,    uint32_t);
SAMBA_MAKE_PRIMITIVE(face_index,    face_index,    uint32_t);
SAMBA_MAKE_PRIMITIVE(element_index, element_index, uint32_t);

// Meta function for getting the entity_rank type associated with a rank index
template <typename RankIndex>
struct index_type_to_entity_rank;

template <>
struct index_type_to_entity_rank<node_index>
{
  typedef entity_rank::node_type type;
};

template <>
struct index_type_to_entity_rank<edge_index>
{
  typedef entity_rank::edge_type type;
};

template <>
struct index_type_to_entity_rank<face_index>
{
  typedef entity_rank::face_type type;
};

template <>
struct index_type_to_entity_rank<element_index>
{
  typedef entity_rank::element_type type;
};

// Meta function for getting the rank index type associated with an entity_rank
template <typename EntityRank>
struct entity_rank_to_index_type;

template <>
struct entity_rank_to_index_type<entity_rank::node_type>
{
  typedef node_index type;
};

template <>
struct entity_rank_to_index_type<entity_rank::edge_type>
{
  typedef edge_index type;
};

template <>
struct entity_rank_to_index_type<entity_rank::face_type>
{
  typedef face_index type;
};

template <>
struct entity_rank_to_index_type<entity_rank::element_type>
{
  typedef element_index type;
};

}

SAMBA_IS_PRIMITIVE(samba::node_index);
SAMBA_IS_PRIMITIVE(samba::edge_index);
SAMBA_IS_PRIMITIVE(samba::face_index);
SAMBA_IS_PRIMITIVE(samba::element_index);

#endif //SAMBA_SAMBA_RANK_INDEX_HPP
