// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef percept_Edge_hpp
#define percept_Edge_hpp

  namespace percept
  {

    template<typename IdType>
    class MyEdge
    {
    public:
      MyEdge(IdType i0, IdType i1) : m_i0(i0<i1?i0:i1), m_i1(i0<i1?i1:i0) {}

      bool operator == ( const MyEdge & e ) const
      { return m_i0 == e.m_i0 && m_i1 == e.m_i1; }

      bool operator != ( const MyEdge & e ) const
      { return not operator==(e); }


      bool operator < ( const MyEdge & e ) const 
      {
        if (operator==(e)) return false;
        if (m_i0 == e.m_i0) return m_i1 < e.m_i1;
        return m_i0 < e.m_i0;
      }
      IdType getId0() const { return m_i0; }
      IdType getId1() const { return m_i1; }

    private:
      IdType m_i0;
      IdType m_i1;
    };

  }//namespace percept

#endif
