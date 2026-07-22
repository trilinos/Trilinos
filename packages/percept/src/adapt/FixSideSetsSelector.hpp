// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef adapt_FixSideSetsSelector_hpp
#define adapt_FixSideSetsSelector_hpp

#include <adapt/Refiner.hpp>


namespace percept {

  // for the refine pass
  class FixSideSetsSelectorRefine : public RefinerSelector
  {

  public:
    PerceptMesh& m_eMesh;
    SetOfEntities m_node_set;
    FixSideSetsSelectorRefine(PerceptMesh& eMesh) : m_eMesh(eMesh)
    {
      m_node_set.clear();
    }

    // incoming is the list of all new elements
    void add_elements(vector<stk::mesh::Entity>::iterator beg, vector<stk::mesh::Entity>::iterator end)
    {
      PerceptMesh& eMesh = m_eMesh;

      for (auto elem = beg; elem != end; ++elem)
        {
          if (!eMesh.is_valid(*elem))
            continue;
          if (eMesh.entity_rank(*elem) != eMesh.element_rank())
            continue;
          unsigned nnode= eMesh.get_bulk_data()->num_nodes(*elem);
          stk::mesh::Entity const *elem_nodes = eMesh.get_bulk_data()->begin_nodes(*elem);

          for (unsigned ii=0; ii < nnode; ii++)
            {
              m_node_set.insert(elem_nodes[ii]);
            }
        }
    }

    virtual bool use_batch_filter() override { return true; }

    virtual void batch_filter(stk::mesh::EntityRank /*rank*/, std::vector<stk::mesh::Entity>& /*elements*/) override { throw std::runtime_error("not impl"); }

    virtual void batch_filter(stk::mesh::EntityRank /*rank*/, SetOfEntities& sides) override
    {
      bool tmp = true;
      if (tmp) return;
      SetOfEntities sides_copy = sides;
      size_t orig_size = sides.size();
      sides.clear();
      for (auto side : sides_copy)
        {
          unsigned nnode= m_eMesh.get_bulk_data()->num_nodes(side);
          stk::mesh::Entity const *side_nodes = m_eMesh.get_bulk_data()->begin_nodes(side);
          for (unsigned ii=0; ii < nnode; ii++)
            {
              if (m_node_set.find(side_nodes[ii]) != m_node_set.end())
                {
                  sides.insert(side);
                  break;
                }
            }
        }
      size_t new_size = sides.size();
      stk::all_reduce( m_eMesh.parallel(), stk::ReduceSum<1>( &new_size ) );
      stk::all_reduce( m_eMesh.parallel(), stk::ReduceSum<1>( &orig_size ) );

      if (m_eMesh.get_rank()==0)
        {
          std::cerr << m_eMesh.rank() << "ref cpu= " << m_eMesh.cpu_time() << " new= " << new_size << " old= " << orig_size << " = " << 100.0*double(new_size)/double(std::max(size_t(1),orig_size)) << "%" << std::endl;
        }
    }

  };

  // for the unrefine pass
  class FixSideSetsSelectorUnrefine : public RefinerSelector
  {

  public:
    PerceptMesh& m_eMesh;
    SetOfEntities m_node_set;
    FixSideSetsSelectorUnrefine(PerceptMesh& eMesh) : m_eMesh(eMesh)
    {
      m_node_set.clear();
    }

    // incoming is the list of all new elements
    void add_elements(SetOfEntities& elem_set)
    {
      PerceptMesh& eMesh = m_eMesh;

      for (auto elem : elem_set)
        {
          if (!eMesh.is_valid(elem))
            continue;

          if (eMesh.entity_rank(elem) != eMesh.element_rank())
            continue;
          unsigned nnode= eMesh.get_bulk_data()->num_nodes(elem);
          stk::mesh::Entity const *elem_nodes = eMesh.get_bulk_data()->begin_nodes(elem);

          for (unsigned ii=0; ii < nnode; ii++)
            {
              m_node_set.insert(elem_nodes[ii]);
            }
        }
    }

    virtual bool use_batch_filter() override { return true; }

    virtual void batch_filter(stk::mesh::EntityRank /*rank*/, std::vector<stk::mesh::Entity>& /*elements*/) override { throw std::runtime_error("still not impl"); }

    virtual void batch_filter(stk::mesh::EntityRank /*rank*/, SetOfEntities& sides) override
    {
      bool tmp = true;
      if (tmp) return;
      SetOfEntities sides_copy = sides;
      size_t orig_size = sides.size();
      sides.clear();
      for (auto side : sides_copy)
        {
          unsigned nnode= m_eMesh.get_bulk_data()->num_nodes(side);
          stk::mesh::Entity const *side_nodes = m_eMesh.get_bulk_data()->begin_nodes(side);
          for (unsigned ii=0; ii < nnode; ii++)
            {
              if (m_node_set.find(side_nodes[ii]) != m_node_set.end())
                {
                  sides.insert(side);
                  break;
                }
            }
        }
      size_t new_size = sides.size();
      stk::all_reduce( m_eMesh.parallel(), stk::ReduceSum<1>( &new_size ) );
      stk::all_reduce( m_eMesh.parallel(), stk::ReduceSum<1>( &orig_size ) );

      if (m_eMesh.get_rank()==0)
        {
          std::cerr << m_eMesh.rank() << "unref cpu= " << m_eMesh.cpu_time() << " new= " << new_size << " old= " << orig_size << " = " << 100.0*double(new_size)/double(std::max(size_t(1),orig_size)) << "%" << std::endl;
        }
    }

  };

}

#endif
