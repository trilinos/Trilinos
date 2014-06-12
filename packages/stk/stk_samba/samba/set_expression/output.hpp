#ifndef SAMBA_SAMBA_SET_EXPRESSION_OUTPUT_HPP
#define SAMBA_SAMBA_SET_EXPRESSION_OUTPUT_HPP

#include <samba/set_expression/set_expression_impl.hpp>

#include <ostream>

namespace samba {

namespace detail {

struct print_set_expression_impl :
  boost::static_visitor<std::ostream &>
{
  std::ostream & m_out;

  print_set_expression_impl(std::ostream & out)
    : m_out(out)
  {}

  std::ostream& operator() (entity_part d) const
  { return m_out << d; }

  template <typename Operation>
  std::ostream& operator() (detail::binary_op_<Operation> const& expr) const
  {
    m_out << "(";
    boost::apply_visitor( print_set_expression_impl(m_out), expr.m_lhs );
    m_out << Operation();
    boost::apply_visitor( print_set_expression_impl(m_out), expr.m_rhs );
    m_out << ")";
    return m_out;
  }
};

} //namespace detail

inline
std::ostream& operator << (std::ostream &out, set_expression const& expr)
{ return boost::apply_visitor(detail::print_set_expression_impl(out),expr()); }

} //namespace samba

#endif //SAMBA_SAMBA_SET_EXPRESSION_OUTPUT_HPP


