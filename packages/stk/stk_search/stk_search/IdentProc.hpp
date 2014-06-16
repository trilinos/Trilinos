/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef STK_SEARCH_IDENTPROC_HPP
#define STK_SEARCH_IDENTPROC_HPP

#include <cstdlib>
#include <iostream>
#include <utility>

#if STK_SEARCH_HAVE_MODERN_COMPILER
#include <boost/type_traits.hpp>
#endif

namespace stk { namespace search {

template <typename Ident, typename Proc = int>
class IdentProc
{
public:
  typedef Ident ident_type;
  typedef Proc proc_type;

  typedef IdentProc<ident_type,proc_type> self_type;

  IdentProc()
    : m_value()
  {}

  IdentProc( ident_type const& i , proc_type const& p )
    : m_value(p,i)
  {}


  ident_type const& id() const {return m_value.second; }
  proc_type const& proc() const {return m_value.first; }

  void set_id(ident_type const& x_id) { m_value.second = x_id; }
  void set_proc(proc_type const& x_proc) { m_value.first = x_proc; }

  bool operator==(self_type const& rhs) const { return m_value == rhs.m_value; }
  bool operator!=(self_type const& rhs) const { return m_value != rhs.m_value; }
  bool operator< (self_type const& rhs) const { return m_value < rhs.m_value; }
  bool operator> (self_type const& rhs) const { return m_value > rhs.m_value; }
  bool operator<=(self_type const& rhs) const { return m_value <= rhs.m_value; }
  bool operator>=(self_type const& rhs) const { return m_value >= rhs.m_value; }

  friend std::ostream& operator<<(std::ostream& out, IdentProc<ident_type,proc_type> const& ip)
  {
    out << "{id:" << ip.id() << ",proc:" << ip.proc() << "}";
    return out;
  }

private:
  std::pair<proc_type,ident_type> m_value;
};


// If you have a modern compiler, some lack of semantic safety (intent) can be
// caught at compilation time.

template <typename T>
struct get_proc
{
#if STK_SEARCH_HAVE_MODERN_COMPILER
  typedef boost::false_type supported;
#endif
  int operator()(T const& id) const
  {
    std::cerr << "get_proc::operator()(..) called on unsupported type." << std::endl;
    std::abort();
    return -1;
  }
};

template <typename T>
struct get_proc<std::pair<T, int> >
{
#if STK_SEARCH_HAVE_MODERN_COMPILER
  typedef boost::false_type supported;
#endif
  int operator()(std::pair<T, int> const& id) const
  {
    std::cerr << "get_proc::operator()(..) called on unsupported type." << std::endl;
    std::abort();
    return -1;
  }
};

template <typename Ident, typename Proc>
struct get_proc< stk::search::IdentProc<Ident,Proc> >
{
#if STK_SEARCH_HAVE_MODERN_COMPILER
  typedef boost::true_type supported;
#endif
  int operator()(stk::search::IdentProc<Ident,Proc> const& id) const
  {
    return id.proc();
  }
};

}} // namespace stk::search

#endif
