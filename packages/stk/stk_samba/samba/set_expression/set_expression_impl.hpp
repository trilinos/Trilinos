#ifndef SAMBA_SAMBA_SET_EXPRESSION_SET_EXPRESSION_IMPL_HPP
#define SAMBA_SAMBA_SET_EXPRESSION_SET_EXPRESSION_IMPL_HPP

#include <samba/types.hpp>

#include <boost/type_traits.hpp>
#include <boost/variant.hpp>
#include <boost/foreach.hpp>

#ifdef SAMBA_ENABLE_PARALLEL
#include <boost/serialization/variant.hpp>
#endif

#include <ostream>

namespace samba { namespace detail {

struct union_
{
  friend inline std::ostream& operator<<(std::ostream& out, union_)
  { return out << "|"; }
};

struct intersect_
{
  friend inline std::ostream& operator<<(std::ostream& out, intersect_)
  { return out << "&"; }
};

struct difference_
{
  friend inline std::ostream& operator<<(std::ostream& out, difference_)
  { return out << "-"; }
};

template <typename Op>
struct binary_op_;

typedef
  boost::variant<
     entity_part
   , boost::recursive_wrapper< binary_op_<union_> >
   , boost::recursive_wrapper< binary_op_<intersect_> >
   , boost::recursive_wrapper< binary_op_<difference_> >
  > set_expression_impl;

template <typename Op>
struct binary_op_
{
  set_expression_impl m_lhs;
  set_expression_impl m_rhs;

  binary_op_()
    : m_lhs()
    , m_rhs()
  {}

  binary_op_(set_expression_impl const& lhs, set_expression_impl const& rhs)
    : m_lhs(lhs)
    , m_rhs(rhs)
  {}

  template <typename Archive>
  void serialize (Archive & ar, const unsigned version)
  {
    ar & m_lhs;
    ar & m_rhs;
  }
};


struct contains_impl :
  boost::static_visitor<bool>
{
  entity_part_range m_partition;

  contains_impl(entity_part_range partition)
    : m_partition(partition)
  {}

  bool operator() (entity_part t) const
  {
    BOOST_FOREACH(entity_part p, m_partition) {
      if ( equivalent(p,t) ) return true;
    }
    return false;
  }

  bool operator() (binary_op_<union_> const& set) const
  { return boost::apply_visitor(*this, set.m_lhs) || boost::apply_visitor(*this, set.m_rhs); }

  bool operator() (binary_op_<intersect_> const& set) const
  { return boost::apply_visitor(*this, set.m_lhs) && boost::apply_visitor(*this, set.m_rhs); }

  bool operator() (binary_op_<difference_> const& set) const
  { return boost::apply_visitor(*this, set.m_lhs) && !boost::apply_visitor(*this, set.m_rhs); }

};

}} //namespace samba::detail

namespace boost {

template <typename Op>
struct has_nothrow_copy<samba::detail::binary_op_<Op> > : mpl::true_ {};

template <typename Op>
struct has_nothrow_constructor<samba::detail::binary_op_<Op> > : mpl::true_ {};

} // namespace boost

namespace samba {

struct set_expression
{
  typedef detail::set_expression_impl value_type;

  //default
  set_expression()
    : m_value(entity_state::invalid())
  {}

  //from m_value
  set_expression( value_type const& p)
    : m_value(p)
  {}

  //from base terminal type
  template <typename T>
  set_expression(  T t
                 , typename boost::enable_if< detail::convertable_to_entity_part<T>,int>::type=0
                )
    : m_value(entity_part(t))
  {}

  //from binary_op
  template <typename Op>
  set_expression(  detail::binary_op_<Op> const& t )
    : m_value(t)
  {}

  static set_expression intersect_(set_expression const& l, set_expression const& r)
  { return set_expression(detail::binary_op_<detail::intersect_>(l.m_value,r.m_value)); }

  static set_expression union_(set_expression const& l, set_expression const& r)
  { return set_expression(detail::binary_op_<detail::union_>(l.m_value,r.m_value)); }

  static set_expression difference_(set_expression const& l, set_expression const& r)
  { return set_expression(detail::binary_op_<detail::difference_>(l.m_value,r.m_value)); }

  value_type const& operator()() const { return m_value; }

  template <typename Partition>
  bool contains( Partition partition) const;

  private:
    value_type m_value;
};

inline bool contains(set_expression const& set, entity_part_range partition )
{ return boost::apply_visitor(detail::contains_impl(partition),set()); }

template <typename Partition>
inline
bool set_expression::contains(Partition partition) const
{ return samba::contains(*this, partition.parts()); }

} //namespace samba

SAMBA_HAS_NO_THROW_COPY(samba::set_expression)
SAMBA_HAS_NO_THROW_CONSTRUCTOR(samba::set_expression)

#endif //SAMBA_SAMBA_SET_EXPRESSION_SET_EXPRESSION_IMPL_HPP
