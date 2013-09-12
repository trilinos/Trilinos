#ifndef stk_adapt_SDCEntityType_hpp
#define stk_adapt_SDCEntityType_hpp

#include <iostream>
#include <stdexcept>
#include <string>
#include <sstream>
#include <vector>

#include <stk_mesh/base/Entity.hpp>

#include <stk_percept/stk_mesh.hpp>

#include <stk_util/environment/CPUTime.hpp>

#include <stk_percept/NoMallocArray.hpp>
#include <stk_percept/PerceptMesh.hpp>
#include <stk_percept/Util.hpp>

#include <stk_percept/PerceptBoostArray.hpp>

#include <stk_adapt/SubDimCell.hpp>

namespace stk {
  namespace adapt {

    // FIXME - this is a hack to avoid a valgrind error that will go away after
    //  migrating to the current stk_migration branch - srk 3/15/13
    extern bool s_compare_using_entity_impl;

    // tuple of node id and the owning element for a node on a sub-dimensional entity (like a face or edge),
    //    and other data...
    enum SubDimCellDataEnum {
      SDC_DATA_GLOBAL_NODE_IDS,
      SDC_DATA_OWNING_ELEMENT_KEY,
      SDC_DATA_OWNING_ELEMENT_ORDINAL,  // the ordinal of the face/edge owning this sub-dim entity
      SDC_DATA_SPACING
    };

    typedef stk::mesh::Entity SDCEntityType;

    class MyEntityLess {
    public:
      stk::percept::PerceptMesh *m_eMesh;
      MyEntityLess(stk::percept::PerceptMesh *eMesh=0) : m_eMesh(eMesh) {}

      /** \brief  Comparison operator */
      bool operator()(const stk::mesh::Entity lhs, const stk::mesh::Entity rhs) const
      {
        if (s_compare_using_entity_impl)
          {
            return lhs.local_offset() < rhs.local_offset();
          }
        else
          {
            const stk::mesh::EntityKey lhs_key = m_eMesh->is_valid(lhs) ? m_eMesh->key(lhs) : stk::mesh::EntityKey();
            const stk::mesh::EntityKey rhs_key = m_eMesh->is_valid(rhs) ? m_eMesh->key(rhs) : stk::mesh::EntityKey();
            return lhs_key < rhs_key;
          }
      }

    }; //class MyEntityLess

    class MyEntityEqual
    {
    public:
      const stk::percept::PerceptMesh *m_eMesh;
      MyEntityEqual(const stk::percept::PerceptMesh *eMesh=0) : m_eMesh(eMesh) {}

      bool operator()(const stk::mesh::Entity lhs, const stk::mesh::Entity rhs) const
      {
        if (s_compare_using_entity_impl)
          {
            return lhs.local_offset() == rhs.local_offset();
          }
        else
          {
            const stk::mesh::EntityKey lhs_key = m_eMesh->is_valid(lhs) ? m_eMesh->key(lhs) : stk::mesh::EntityKey();
            const stk::mesh::EntityKey rhs_key = m_eMesh->is_valid(rhs) ? m_eMesh->key(rhs) : stk::mesh::EntityKey();
            return lhs_key == rhs_key;
          }
      }
    };

    struct CompareSDCEntityType {
      stk::percept::PerceptMesh *m_eMesh;
      CompareSDCEntityType(stk::percept::PerceptMesh *eMesh=0) : m_eMesh(eMesh) {}
      bool operator() (SDCEntityType i, SDCEntityType j) {
        //return (m_eMesh.identifier(i) < m_eMesh.identifier(j));
        MyEntityLess el = MyEntityLess(m_eMesh);
        return el(i,j);
      }
    };

    template<class T, std::size_t N=4>
    class MySDCHashCode
    {
    public:

      typedef  stk::percept::NoMallocArray<T, N> sdc_type;

      stk::percept::PerceptMesh* m_eMesh;
      MySDCHashCode(stk::percept::PerceptMesh* eMesh=0) : m_eMesh(eMesh) {}

      int operator()(sdc_type& sdc)
      {
        size_t sum = 0;
        for ( typename sdc_type::iterator i = sdc.begin(); i != sdc.end(); i++)
          {

            if (s_compare_using_entity_impl)
              {
                //sum += size_t(i->local_offset());
                sum += size_t(i->local_offset());
              }
            else
              {
                //sum += static_cast<size_t>((*i).identifier());
                sum += static_cast<size_t>(m_eMesh->identifier(*i));
              }
          }
        return sum;
      }
    };

    template<class T, std::size_t N=4, class CompareClass = SubDimCellCompare<T>, class HC = MySDCHashCode<T, N> >
    class MySubDimCell : public SubDimCell<T, N, CompareClass, HC>
    {
    public:
      typedef SubDimCell<T, N, CompareClass, HC> base_type;

      stk::percept::PerceptMesh& m_eMesh;
      MySubDimCell(stk::percept::PerceptMesh& eMesh) : base_type(), m_eMesh(eMesh) {
        base_type::m_HashCode = HC(&eMesh);
        base_type::m_CompareClass = CompareClass(&eMesh);
      }
      MySubDimCell(stk::percept::PerceptMesh& eMesh, unsigned num_ids) : base_type(num_ids), m_eMesh(eMesh)
      {
        base_type::m_HashCode = HC(&eMesh);
        base_type::m_CompareClass = CompareClass(&eMesh);
      }

    };


    template<>
    struct my_fast_hash<SDCEntityType, 2> : public std::unary_function< MySubDimCell<SDCEntityType, 2, CompareSDCEntityType>, std::size_t>
    {
      typedef MySubDimCell<SDCEntityType, 2, CompareSDCEntityType> _Tp ;

      inline std::size_t
      operator()(const _Tp& x) const
      {
        return x.getHash();
      }

    };


    template<>
    struct my_fast_hash<SDCEntityType, 4> : public std::unary_function< MySubDimCell<SDCEntityType, 4, CompareSDCEntityType>, std::size_t>
    {
      typedef MySubDimCell<SDCEntityType, 4, CompareSDCEntityType> _Tp ;

      inline std::size_t
      operator()(const _Tp& x) const
      {
        return x.getHash();
      }

    };


    template<>
    struct my_fast_equal_to<SDCEntityType, 2> :  public std::binary_function< MySubDimCell<SDCEntityType, 2, CompareSDCEntityType>,
                                                                              MySubDimCell<SDCEntityType, 2, CompareSDCEntityType>, bool>
    {
      typedef MySubDimCell<SDCEntityType, 2, CompareSDCEntityType> _Tp ;
      inline bool
      operator()(const _Tp& x, const _Tp& y) const
      {
        if (x.getHash() != y.getHash()) return false;
        if (x.size() != y.size()) return false;
        _Tp::const_iterator ix = x.begin();
        _Tp::const_iterator iy = y.begin();
        MyEntityEqual ee = MyEntityEqual(&(x.m_eMesh));
        for (; ix != x.end(); ix++, iy++)
          {
            if (!ee(*ix, *iy)) return false;
          }
        return true;
      }
    };

    template<>
    struct my_fast_equal_to<SDCEntityType, 4> :  public std::binary_function< MySubDimCell<SDCEntityType, 4, CompareSDCEntityType>,
                                                                              MySubDimCell<SDCEntityType, 4, CompareSDCEntityType>, bool>
    {
      typedef MySubDimCell<SDCEntityType, 4, CompareSDCEntityType> _Tp ;
      inline bool
      operator()(const _Tp& x, const _Tp& y) const
      {
        if (x.getHash() != y.getHash()) return false;
        if (x.size() != y.size()) return false;
        _Tp::const_iterator ix = x.begin();
        _Tp::const_iterator iy = y.begin();
        MyEntityEqual ee = MyEntityEqual(&(x.m_eMesh));
        for (; ix != x.end(); ix++, iy++)
          {
            if (!ee(*ix, *iy)) return false;
          }
        return true;
      }
    };

  }
}
#endif
