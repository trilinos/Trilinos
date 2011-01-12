#ifndef stk_percept_Edge_hpp
#define stk_percept_Edge_hpp

namespace stk
{
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
}//namespace stk

#endif
