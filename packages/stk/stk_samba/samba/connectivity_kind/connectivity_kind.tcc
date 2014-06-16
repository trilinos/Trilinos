#ifndef SAMBA_SAMBA_CONNECTIVITY_KIND_CONNECTIVITY_KIND_TCC
#define SAMBA_SAMBA_CONNECTIVITY_KIND_CONNECTIVITY_KIND_TCC

#include <samba/detail/apply_uint_8_functor.hpp>
#include <samba/connectivity_kind/types.hpp>

namespace samba {

namespace detail {
struct connectivity_kind_output_impl
{
  std::ostream& m_out;

  connectivity_kind_output_impl(std::ostream& out)
    : m_out(out)
  {}

  typedef void result_type;

  template <connectivity_kind::value_type Kind>
  void operator()(connectivity_kind::kind_type<Kind> k) const
  { m_out << k; }
};
} //namespace detail

inline std::ostream& operator<<(std::ostream& out, connectivity_kind rk)
{
  detail::apply_functor< connectivity_kind::kind_type
                        ,detail::connectivity_kind_output_impl
                       > apply;

  detail::connectivity_kind_output_impl f(out);
  out << "{" << connectivity_kind::tag() << ":";

  apply(f,rk());

  out << "}";
  return out;
}

} //namespace samba

#endif //SAMBA_SAMBA_CONNECTIVITY_KIND_CONNECTIVITY_KIND_TCC
