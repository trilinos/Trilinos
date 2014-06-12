#ifndef SAMBA_SAMBA_CONNECTIVITY_ORIENTATION_HPP
#define SAMBA_SAMBA_CONNECTIVITY_ORIENTATION_HPP

#include <samba/utility.hpp>

#include <samba/connectivity_ordinal.hpp>

namespace samba {

struct connectivity_orientation
{
  typedef uint8_t value_type;

  static const int num_bits = 8;

  static const connectivity_orientation invalid()
  { static connectivity_orientation d = {integer_max<value_type>::value}; return d; }

  //create a new connectivity_orientation
  static const connectivity_orientation create(value_type v)
  { connectivity_orientation d = {v}; return d; }

  value_type operator()() const { return m_value; }

  SAMBA_PRIMITIVE_COMPARISONS(connectivity_orientation,value_type)

  connectivity_orientation & operator=(connectivity_orientation::value_type v)
  { m_value = v; return *this; }


  //for adjacent entities the orientation is actually a connectivity_ordinal
  //representing the side id of the adjacent entity side.
  connectivity_ordinal ordinal() const
  { return connectivity_ordinal::create(m_value); }

  template <typename Archive>
  void serialize(Archive &ar, const unsigned version)
  { ar & m_value; }

  value_type m_value;
};

} //namespace samba


SAMBA_IS_PRIMITIVE(samba::connectivity_orientation)

#endif //SAMBA_SAMBA_CONNECTIVITY_ORIENTATION_HPP

