#ifndef SAMBA_SAMBA_CONNECTIVITY_KIND_HPP
#define SAMBA_SAMBA_CONNECTIVITY_KIND_HPP

#include <samba/utility.hpp>

namespace samba {

struct connectivity_kind
{
  struct tag
  {
    typedef int value_type;
    friend inline std::ostream& operator<<(std::ostream& out,tag)
    { return out << "connectivity_kind"; }
  };

  typedef uint8_t value_type;

  static const int num_bits = 8;

  template <connectivity_kind::value_type Kind>
  struct kind_type
  {
    static const connectivity_kind::value_type value = Kind;
  };

  typedef kind_type<0> fixed_type;
  typedef kind_type<1> dynamic_type;
  typedef kind_type<2> invalid_type;

  static const connectivity_kind fixed()
  { static connectivity_kind d = {fixed_type::value}; return d; }

  static const connectivity_kind dynamic()
  { static connectivity_kind d = {dynamic_type::value}; return d; }

  static const connectivity_kind invalid()
  { static connectivity_kind d = {invalid_type::value}; return d; }

  //create a new connectivity_kind
  static const connectivity_kind create(value_type v)
  { connectivity_kind d = {v}; return d; }

  value_type operator()() const { return m_value; }

  SAMBA_PRIMITIVE_COMPARISONS(connectivity_kind,value_type)

  connectivity_kind & operator=(connectivity_kind::value_type v)
  { m_value = v; return *this; }

  template <typename Archive>
  void serialize(Archive &ar, const unsigned version)
  { ar & m_value; }

  value_type m_value;
};

} //namespace samba

SAMBA_IS_PRIMITIVE(samba::connectivity_kind)

#include <samba/connectivity_kind/connectivity_kind.tcc>

#endif //SAMBA_SAMBA_CONNECTIVITY_KIND_HPP

