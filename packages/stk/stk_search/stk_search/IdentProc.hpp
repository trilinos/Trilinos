/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef STK_SEARCH_IDENTPROC_HPP
#define STK_SEARCH_IDENTPROC_HPP

#include <iostream>
#include <utility>
#include <boost/type_traits.hpp>

namespace stk { namespace search {

#if 1
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
#else
template <typename Ident, typename Proc = int>
class IdentProc
{
public:
  typedef Ident ident_type;
  typedef Proc proc_type;

  typedef IdentProc<ident_type,proc_type> self_type;

  IdentProc()
    : m_proc(), m_id()
  {}

  IdentProc( ident_type const& i , proc_type const& p )
    : m_proc(p), m_id(i)
  {}

  ident_type const& id() const {return m_id; }
  proc_type const& proc() const {return m_proc; }

  void set_id(ident_type const& x_id) { m_id = x_id; }
  void set_proc(proc_type const& x_proc) { m_proc = x_proc; }

  bool operator==(self_type const& rhs) const
  {
    return m_proc == rhs.m_proc && m_id == rhs.m_id;
  }
  bool operator!=(self_type const& rhs) const
  {
    return m_proc != rhs.m_proc || m_id != rhs.m_id;
  }
  bool operator< (self_type const& rhs) const
  {
    return m_proc < rhs.m_proc || (m_proc == rhs.m_proc && m_id < rhs.m_id);
  }
  bool operator> (self_type const& rhs) const
  {
    return m_proc > rhs.m_proc || (m_proc == rhs.m_proc && m_id > rhs.m_id);
  }
  bool operator<=(self_type const& rhs) const
  {
    return m_proc <= rhs.m_proc || (m_proc == rhs.m_proc && m_id <= rhs.m_id);
  }
  bool operator>=(self_type const& rhs) const
  {
    return m_proc >= rhs.m_proc || (m_proc == rhs.m_proc && m_id >= rhs.m_id);
  }

  friend std::ostream& operator<<(std::ostream& out, IdentProc<ident_type,proc_type> const& ip)
  {
    out << "{id:" << ip.id() << ",proc:" << ip.proc() << "}";
    return out;
  }

private:
  proc_type  m_proc;
  ident_type m_id;
  bool padding;
};
#endif


template <typename Ident>
struct get_proc
{
  typedef boost::false_type supported;
};

#if 0
// Expresses what SHOULD NOT BE DONE.
template <typename T>
struct get_proc<std::pair<T, int> >
{
  typedef boost::false_type supported;
  int operator()(std::pair<T, int> const& id) const
  {
    std::cerr << "get_proc::operator()(..) called on unsupported type." << std::endl;
    std::abort();
    return -1;
  }
};
#endif

template <typename Ident, typename Proc>
struct get_proc< stk::search::IdentProc<Ident,Proc> >
{
  typedef boost::true_type supported;
  int operator()(stk::search::IdentProc<Ident,Proc> const& id) const
  {
    return id.proc();
  }
};

}} // namespace stk::search

#endif
