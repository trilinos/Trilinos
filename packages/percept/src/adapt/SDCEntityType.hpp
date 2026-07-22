// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef adapt_SDCEntityType_hpp
#define adapt_SDCEntityType_hpp

#include <iostream>
#include <stdexcept>
#include <string>
#include <sstream>
#include <vector>

#include <stk_mesh/base/Entity.hpp>

#include <percept/stk_mesh.hpp>

#include <stk_util/environment/CPUTime.hpp>

#include <percept/NoMallocArray.hpp>
#include <percept/PerceptMesh.hpp>
#include <percept/Util.hpp>

#include <percept/PerceptBoostArray.hpp>

#include <adapt/SubDimCell.hpp>

  namespace percept {

    // FIXME - this is a hack to avoid a valgrind error that will go away after
    //  migrating to the current stk_migration branch - srk 3/15/13
    extern bool s_compare_using_entity_impl;

    // tuple of node id and the owning element for a node on a sub-dimensional entity (like a face or edge),
    //    and other data...
    enum SubDimCellDataEnum {
      SDC_DATA_GLOBAL_NODE_IDS,
      SDC_DATA_OWNING_ELEMENT_KEY,
      SDC_DATA_OWNING_SUBDIM_RANK,  // the rank (2,1,...) of the face/edge owning this sub-dim entity
      SDC_DATA_OWNING_SUBDIM_ORDINAL,  // the ordinal of the face/edge owning this sub-dim entity
      SDC_DATA_SPACING
    };

    typedef stk::mesh::Entity SDCEntityType;

    class MyEntityLess {
    public:
      percept::PerceptMesh *m_eMesh;
      MyEntityLess(percept::PerceptMesh *eMesh=0) : m_eMesh(eMesh) {}

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
      const percept::PerceptMesh *m_eMesh;
      MyEntityEqual(const percept::PerceptMesh *eMesh=0) : m_eMesh(eMesh) {}

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
      percept::PerceptMesh *m_eMesh;
      CompareSDCEntityType(percept::PerceptMesh *eMesh=0) : m_eMesh(eMesh) {}
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

      typedef  percept::NoMallocArray<T, N> sdc_type;

      percept::PerceptMesh* m_eMesh;
      MySDCHashCode(percept::PerceptMesh* eMesh=0) : m_eMesh(eMesh) {}

      int operator()(sdc_type& sdc)
      {
        size_t sum = 0;
        const typename sdc_type::const_iterator i_begin = sdc.begin();
        const typename sdc_type::const_iterator i_end = sdc.end();

        if (s_compare_using_entity_impl)
          {
            for (typename sdc_type::const_iterator i = i_begin; i != i_end; ++i)
              {
                sum += size_t(i->local_offset());
              }
          }
        else
          {
            for ( typename sdc_type::const_iterator i = i_begin; i != i_end; ++i)
              {
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

      percept::PerceptMesh* m_eMesh;
      MySubDimCell(percept::PerceptMesh* eMesh) : base_type(), m_eMesh(eMesh) {
        base_type::m_HashCode = HC(eMesh);
        base_type::m_CompareClass = CompareClass(eMesh);
      }
      MySubDimCell(percept::PerceptMesh* eMesh, unsigned num_ids) : base_type(num_ids), m_eMesh(eMesh)
      {
        base_type::m_HashCode = HC(eMesh);
        base_type::m_CompareClass = CompareClass(eMesh);
      }
      MySubDimCell()
      {
          m_eMesh = 0;
      }

      MySubDimCell& operator=(const MySubDimCell& from);
    };


    template<>
    struct my_fast_hash<SDCEntityType, 2> : public std::function< std::size_t(MySubDimCell<SDCEntityType, 2, CompareSDCEntityType>)>
    {
      typedef MySubDimCell<SDCEntityType, 2, CompareSDCEntityType> _Tp ;

      inline std::size_t
      operator()(const _Tp& x) const
      {
        return x.getHash();
      }

    };


    template<>
    struct my_fast_hash<SDCEntityType, 4> : public std::function< std::size_t(MySubDimCell<SDCEntityType, 4, CompareSDCEntityType>)>
    {
      typedef MySubDimCell<SDCEntityType, 4, CompareSDCEntityType> _Tp ;

      inline std::size_t
      operator()(const _Tp& x) const
      {
        return x.getHash();
      }

    };


    template<>
    struct my_fast_equal_to<SDCEntityType, 2> :  public std::function<bool(MySubDimCell<SDCEntityType, 2, CompareSDCEntityType>,
										  MySubDimCell<SDCEntityType, 2, CompareSDCEntityType>)>
    {
      typedef MySubDimCell<SDCEntityType, 2, CompareSDCEntityType> _Tp ;
      inline bool
      operator()(const _Tp& x, const _Tp& y) const
      {
        if (x.getHash() != y.getHash()) return false;
        if (x.size() != y.size()) return false;
        _Tp::const_iterator ix = x.begin();
        _Tp::const_iterator iy = y.begin();
        MyEntityEqual ee = MyEntityEqual((x.m_eMesh));
        for (; ix != x.end(); ix++, iy++)
          {
            if (!ee(*ix, *iy)) return false;
          }
        return true;
      }
    };

    template<>
    struct my_fast_equal_to<SDCEntityType, 4> :  public std::function< bool(MySubDimCell<SDCEntityType, 4, CompareSDCEntityType>,
									    MySubDimCell<SDCEntityType, 4, CompareSDCEntityType>)>
    {
      typedef MySubDimCell<SDCEntityType, 4, CompareSDCEntityType> _Tp ;
      inline bool
      operator()(const _Tp& x, const _Tp& y) const
      {
        if (x.getHash() != y.getHash()) return false;
        if (x.size() != y.size()) return false;
        _Tp::const_iterator ix = x.begin();
        _Tp::const_iterator iy = y.begin();
        MyEntityEqual ee = MyEntityEqual((x.m_eMesh));
        for (; ix != x.end(); ix++, iy++)
          {
            if (!ee(*ix, *iy)) return false;
          }
        return true;
      }
    };

  }

#endif
