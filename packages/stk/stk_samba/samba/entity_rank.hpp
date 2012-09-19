#ifndef SAMBA_SAMBA_ENTITY_RANK_HPP
#define SAMBA_SAMBA_ENTITY_RANK_HPP

#include <samba/utility.hpp>

namespace samba {

struct entity_rank
{
  struct tag
  {
    typedef int value_type;
    static const int value = 0; //sort first by entity_rank
    friend inline std::ostream& operator<<(std::ostream& out,tag)
    { return out << "entity_rank"; }
  };

  typedef uint8_t value_type;

  static const int num_bits = 8;

  template <entity_rank::value_type Rank>
  struct rank_type
  {
    typedef entity_rank::value_type value_type;
    static const value_type value = Rank;
  };

  typedef rank_type<0> node_type;
  typedef rank_type<1> edge_type;
  typedef rank_type<2> face_type;
  typedef rank_type<3> element_type;
  typedef rank_type<4> invalid_type;

  static const entity_rank node()
  { static entity_rank d = {node_type::value}; return d; }

  static const entity_rank edge()
  { static entity_rank d = {edge_type::value}; return d; }

  static const entity_rank face()
  { static entity_rank d = {face_type::value}; return d; }

  static const entity_rank element()
  { static entity_rank d = {element_type::value}; return d; }

  static const entity_rank invalid()
  { static entity_rank d = {invalid_type::value}; return d; }

  //create a new entity_rank
  static const entity_rank create(value_type v)
  { entity_rank d = {v}; return d; }

  value_type operator()() const { return m_value; }

  SAMBA_PRIMITIVE_COMPARISONS(entity_rank,value_type)
  SAMBA_ARITHMETIC_OPERATORS(entity_rank,value_type)

  entity_rank & operator=(entity_rank::value_type v)
  { m_value = v; return *this; }

  template <typename Archive>
  void serialize(Archive &ar, const unsigned version)
  { ar & m_value; }

  value_type m_value;
};

} //namespace samba

SAMBA_IS_PRIMITIVE(samba::entity_rank)

#include <samba/entity_rank/entity_rank.tcc>

#endif //SAMBA_SAMBA_ENTITY_RANK_HPP
