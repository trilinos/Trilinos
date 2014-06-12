#ifndef SAMBA_SAMBA_ENTITY_RANK_ENTITY_RANK_TCC
#define SAMBA_SAMBA_ENTITY_RANK_ENTITY_RANK_TCC

#include <samba/entity_rank/types.hpp>
#include <samba/detail/apply_uint_8_functor.hpp>

namespace samba {

namespace detail {
struct entity_rank_output_impl
{
  std::ostream& m_out;

  entity_rank_output_impl(std::ostream& out)
    : m_out(out)
  {}

  typedef void result_type;

  template <entity_rank::value_type Rank>
  void operator()(entity_rank::rank_type<Rank> r) const
  { m_out << r; }
};
} //namespace detail


inline std::ostream& operator<<(std::ostream& out, entity_rank r)
{
  detail::entity_rank_output_impl f(out);

  detail::apply_functor< entity_rank::rank_type
                        ,detail::entity_rank_output_impl
                       > apply;

  out << "{" << entity_rank::tag() << ":";

  apply(f,r());

  out << "}";
  return out;
}

} //namespace samba

#endif //SAMBA_SAMBA_ENTITY_RANK_ENTITY_RANK_TCC

