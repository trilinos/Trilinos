#ifndef SAMBA_SAMBA_ENTITY_STATE_HPP
#define SAMBA_SAMBA_ENTITY_STATE_HPP

#include <samba/utility.hpp>

namespace samba {

/**
 * entity_state is a primitive with a small number of well-defined
 * values that represent the state an entity can be in:
 *   - universe: A way to refer to everything
 *   - modified: Going away?
 *   - owned:    Owned on this process, can intersect with shared
 *   - shared:   Shared with another process, can intersect with ghosted and owned, can modify
 *   - ghosted:  Not owned on this process, can intersect with shared
 */
struct entity_state
{
  struct tag
  {
    typedef int value_type;
    static const int value = 2;
    friend inline std::ostream& operator<<(std::ostream& out,tag)
    { return out << "entity_state"; }
  };

  typedef uint8_t value_type;

  template <entity_state::value_type State>
  struct state_type
  {
    typedef entity_state::value_type value_type;
    static const value_type value = State;
    typedef entity_state type;
    operator type() const
    {
      type t = {value};
      return t;
    }
  };

  typedef state_type<0> universe_type;
  typedef state_type<1> modified_type;
  typedef state_type<2> owned_type;
  typedef state_type<3> shared_type;
  typedef state_type<4> ghosted_type;
  typedef state_type<5> invalid_type;

  static const entity_state universe()
  { static entity_state d = {universe_type::value}; return d; }

  static const entity_state modified()
  { static entity_state d = {modified_type::value}; return d; }

  static const entity_state owned()
  { static entity_state d = {owned_type::value}; return d; }

  static const entity_state shared()
  { static entity_state d = {shared_type::value}; return d; }

  static const entity_state ghosted()
  { static entity_state d = {ghosted_type::value}; return d; }

  static const entity_state invalid()
  { static entity_state d = {invalid_type::value}; return d; }

  static const entity_state create(value_type v)
  { entity_state d = {v}; return d; }

  value_type operator()() const { return m_value; }

  SAMBA_PRIMITIVE_COMPARISONS(entity_state,value_type)
  SAMBA_ARITHMETIC_OPERATORS(entity_state,value_type)

  entity_state & operator=(entity_state::value_type v)
  { m_value = v; return *this; }

  template <typename Archive>
  void serialize(Archive &ar, const unsigned version)
  { ar & m_value; }

  value_type m_value;
};

} //namespace samba

SAMBA_IS_PRIMITIVE(samba::entity_state)

#include <samba/entity_state/entity_state.tcc>

#endif //SAMBA_SAMBA_ENTITY_STATE_HPP
