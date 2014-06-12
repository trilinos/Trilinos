#ifndef SAMBA_SAMBA_ENTITY_STATE_ENTITY_STATE_TCC
#define SAMBA_SAMBA_ENTITY_STATE_ENTITY_STATE_TCC

#include <samba/detail/apply_uint_8_functor.hpp>
#include <samba/entity_state/types.hpp>

namespace samba {


namespace detail {
struct entity_state_output_impl
{
  std::ostream& m_out;

  entity_state_output_impl(std::ostream& out)
    : m_out(out)
  {}

  typedef void result_type;

  template <entity_state::value_type State>
  void operator()(entity_state::state_type<State> r) const
  { m_out << r; }
};
} //namespace detail


inline std::ostream& operator<<(std::ostream& out, entity_state s)
{
  detail::apply_functor< entity_state::state_type
                        ,detail::entity_state_output_impl
                       > apply;

  detail::entity_state_output_impl f(out);
  out << "{" << entity_state::tag() << ":";

  apply(f,s());

  out << "}";
  return out;
}


} //namespace samba

#endif //SAMBA_SAMBA_ENTITY_STATE_ENTITY_STATE_TCC

