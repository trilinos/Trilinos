// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
// 
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 

#ifndef STK_SEARCH_IDENTPROC_HPP
#define STK_SEARCH_IDENTPROC_HPP

#include <cstdlib>
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


// If you have a modern compiler, some lack of semantic safety (intent) can be
// caught at compilation time.

template <typename T>
struct get_proc
{
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
  int operator()(stk::search::IdentProc<Ident,Proc> const& id) const
  {
    return id.proc();
  }
};

}} // namespace stk::search

#endif
