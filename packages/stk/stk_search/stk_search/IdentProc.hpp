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

}} // namespace stk::search

#endif
