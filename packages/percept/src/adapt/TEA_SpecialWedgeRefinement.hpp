// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef adapt_TEA_SpecialWedgeRefinement_hpp
#define adapt_TEA_SpecialWedgeRefinement_hpp


#include <adapt/TransitionElementAdapter.hpp>

namespace percept {

  extern bool s_allow_special_wedge_refine;  // for wedge boundary-layer refine

  template<class RefinePredicate>
  class TEA_SpecialWedgeRefinement : public TransitionElementAdapter<RefinePredicate>
  {
    // how to choose to treat parts specially
    stk::mesh::Selector *m_wedge_selector;
    bool m_enable_special_patterns; // if true, try to not refine wedge in vertical direction (don't split quad faces)
    bool m_allow_unrefine;  // this algorithm doesn't yet work well with unrefine, so set this to false if enable_special_patterns is true
    bool m_debug_local;
    bool m_base_debug;
  public:
    typedef TransitionElementAdapter<RefinePredicate> Base;
    typedef std::map<stk::mesh::Entity, int> ChangeRefineFieldRequestType;

    /// A special algorithm for wedge refinement that marks all wedges in a column if the wedge column is adjacent to a marked element
    /// @param wedge_selector :: pass in a selector for which parts to apply this special algorithm to
    /// @param enable_special_patterns :: if true, try to not refine wedge in vertical direction (don't split quad faces)
    /// @param_allow_unrefine :: this algorithm doesn't yet work well with unrefine, so set this to false if enable_special_patterns is true
    ///
    TEA_SpecialWedgeRefinement(RefinePredicate& predicate_refine,
                               percept::PerceptMesh& eMesh, UniformRefinerPatternBase & bp, stk::mesh::FieldBase *proc_rank_field=0, stk::mesh::Selector *wedge_selector=0,
                               bool enable_special_patterns = false, bool allow_unrefine = true, bool debug = false) :
      Base(predicate_refine, eMesh, bp, proc_rank_field, debug), m_wedge_selector(wedge_selector),
      m_enable_special_patterns(enable_special_patterns), m_allow_unrefine(allow_unrefine), m_debug_local(false), m_base_debug(false)
    {
      m_base_debug = Base::m_debug;
    }
    virtual ~TEA_SpecialWedgeRefinement()
    {
    }

    virtual void modify_marks(stk::mesh::Entity element, stk::mesh::EntityRank needed_entity_rank, std::vector<bool>& markInfoVec, bool elem_is_marked)
    {
      if (!s_allow_special_wedge_refine)
        return;
      elem_is_marked = true;
      if (elem_is_marked && m_wedge_selector && (*m_wedge_selector)(Base::m_eMesh.bucket(element)))
        {
          if (needed_entity_rank == Base::m_eMesh.edge_rank())
            {
              // unmark the vertical edges
              for (unsigned iedge = 6; iedge < 9; ++iedge)
                {
                  markInfoVec[iedge] = false;
                }
            }
          else if (needed_entity_rank == Base::m_eMesh.side_rank())
            {
              markInfoVec.assign(markInfoVec.size(), false);
            }
        }
      if (elem_is_marked && Base::m_eMesh.bucket(element).topology() != stk::topology::WEDGE_6
          && ( (   needed_entity_rank == Base::m_eMesh.edge_rank())
               || (needed_entity_rank == Base::m_eMesh.side_rank()) ) )
        {
          stk::mesh::Entity non_wedge = element;
          std::set<stk::mesh::Entity> neighbors;
          Base::m_eMesh.get_node_neighbors(non_wedge, neighbors);

          std::set<stk::mesh::Entity>::iterator it;
          for (it = neighbors.begin(); it != neighbors.end(); ++it)
            {
              stk::mesh::Entity wedge = *it;
              if (m_wedge_selector && (*m_wedge_selector)(Base::m_eMesh.bucket(wedge)))
                {
                  int face_0=0, face_1=0;
                  bool ifn = Base::m_eMesh.is_face_neighbor(non_wedge, wedge, &face_0, &face_1);
                  int edge_0=0, edge_1=0;
                  bool ien = Base::m_eMesh.is_edge_neighbor(non_wedge, wedge, &edge_0, &edge_1);
                  if ((ifn && face_1 < 3) || (ien && edge_1 >= 6))
                    {
                      markInfoVec.assign(markInfoVec.size(), false);
                    }
                }
            }
        }
    }

    void clearWedgeRefineField()
    {
      PerceptMesh& eMesh = Base::m_eMesh;

      std::vector<stk::mesh::Entity> vec;
      VERIFY_OP_ON(m_wedge_selector, !=, 0, "must supply wedge selector to use TEA_SpecialWedgeRefinement");
      stk::mesh::get_selected_entities(*m_wedge_selector , eMesh.get_bulk_data()->buckets(eMesh.element_rank()), vec);

      for (unsigned iElement = 0; iElement < vec.size(); iElement++)
        {
          stk::mesh::Entity element = vec[iElement];

          int *refine_field_elem = stk::mesh::field_data( *Base::m_refine_field , element );

          if (refine_field_elem[0] != 0)
            refine_field_elem[0] = 0;
        }
    }

    virtual void special_processing(const std::string& step, void */*data*/)
    {
      if (step == "tea_pre_refine")
        {
          s_allow_special_wedge_refine = m_enable_special_patterns;
        }
      else if (step == "tea_post_refine")
        {
          s_allow_special_wedge_refine = false;
        }
      else if (step == "tea_pre_unrefine")
        {
          s_allow_special_wedge_refine = m_enable_special_patterns;
          clearWedgeRefineField();

          if (m_allow_unrefine)
            {
              // stage 1
              int max_iter=100;
              int iter=0;
              int64_t did_change=0;
              while ((iter++ < max_iter) && (did_change = this->enforce_boundary_layer_refine_pattern(iter, -1)) )
                {
                  if ((1||m_base_debug) && Base::m_eMesh.get_rank()==0)
                    std::cout << "P[" << Base::m_eMesh.get_rank() << "] tmp srk TEA_SpecialWedgeRefinement::enforce_boundary_layer_refine_pattern unref iter= " << iter << " did_change= " << did_change
                              << std::endl;
                  if (iter > max_iter-5)
                    {
                      throw std::runtime_error("TEAS_1: too many iterations");
                    }
                }
            }
        }
      else if (step == "tea_post_unrefine")
        {
          s_allow_special_wedge_refine = false;
        }
      else if (step == "tea_before_enforce_te_consistency")
        {

          clearWedgeRefineField();

          // stage 1
          {
            int max_iter=100;
            int iter=0;
            int64_t did_change=0;
            while ((iter++ < max_iter) && (did_change = this->enforce_boundary_layer_refine_pattern(iter, 1)) )
              {
                if ((1||m_base_debug) && Base::m_eMesh.get_rank()==0)
                  std::cout << "P[" << Base::m_eMesh.get_rank() << "] tmp srk TEA_SpecialWedgeRefinement::enforce_boundary_layer_refine_pattern iter= " << iter << " did_change= " << did_change
                            << std::endl;
                if (iter > max_iter-5)
                  {
                    throw std::runtime_error("TEAS_2: too many iterations");
                  }
              }
          }
        }
      else if (step == "tea_after_single_pass_enforce_te")
        {
        }
      else if (step == "tea_after_enforce_te_consistency")
        {
        }
      // needed for RunAdaptRun
      else if (step == "refiner_pre_initializeDB")
        {
          s_allow_special_wedge_refine = m_enable_special_patterns;
        }
      else if (step == "refiner_post_initializeDB")
        {
          s_allow_special_wedge_refine = false;
        }
    }

    // this one doesn't iterate neighbors, just checks for any marked neighbors and marks element 
    int64_t mark_neighbor_tri_face_first_pass(stk::mesh::Entity element, int refine_or_unrefine)
    {
      PerceptMesh& eMesh = Base::m_eMesh;

      if (!eMesh.owned(element))
        return 0;

      RefineFieldType *refine_field = Base::m_refine_field;
      bool isParent = (eMesh.hasFamilyTree(element) && eMesh.isParentElement(element));
      if (isParent)
        return 0;

      bool useEdgeNeighbors = true;
      if (eMesh.getProperty("TEA_SpecialWedgeRefinement.useEdgeNeighbors") == "false")
        useEdgeNeighbors = false;

      LocalSetType neighbors;
      int64_t did_change = 0;
      int *refine_field_elem = stk::mesh::field_data( *refine_field , element );

      if (1)
        {
          if (refine_field_elem[0] != 0)
            return 0;
        }

      neighbors.clear();
      Base::get_node_neighbors(element, neighbors);

      for (LocalSetType::iterator it = neighbors.begin(); it != neighbors.end(); ++it)
        {
          stk::mesh::Entity neigh = *it;
          if (neigh == element)
            continue;

          bool isNParent = (eMesh.hasFamilyTree(neigh) && eMesh.isParentElement(neigh));
          if (isNParent)
            continue;
          int *refine_field_neigh = stk::mesh::field_data( *refine_field , neigh );

          bool isFaceN = false;
          int face_0 = -1, face_1 = -1;
          isFaceN = eMesh.is_face_neighbor(element, neigh, &face_0, &face_1);

          bool isEdgeN = false;
          int edge_0=-1, edge_1=-1;
          isEdgeN = eMesh.is_edge_neighbor(element, neigh, &edge_0, &edge_1, 0);

          int markVal = 0;
          if (isFaceN && face_0 >= 3)
            {
              if (refine_or_unrefine == 1 && refine_field_neigh[0] >= 1)
                {
                  VERIFY_OP_ON(refine_field_elem[0], != , refine_or_unrefine, "bad refine_field_elem");
                  markVal = 1;
                }
              else if (refine_or_unrefine == -1 && refine_field_neigh[0] == -1)
                {
                  VERIFY_OP_ON(refine_field_elem[0], != , refine_or_unrefine, "bad refine_field_elem");
                  markVal = -1;
                }
            }

          if (useEdgeNeighbors && isEdgeN && edge_0 < 6 && eMesh.topology(neigh) != stk::topology::WEDGE_6)
            {
              if (refine_or_unrefine == 1 && refine_field_neigh[0] >= 1)
                {
                  VERIFY_OP_ON(refine_field_elem[0], != , refine_or_unrefine, "bad refine_field_elem");
                  markVal = 1;
                }
              else if (refine_or_unrefine == -1 && refine_field_neigh[0] == -1)
                {
                  VERIFY_OP_ON(refine_field_elem[0], != , refine_or_unrefine, "bad refine_field_elem");
                  markVal = -1;
                }
            }

          if (markVal != 0)
            {
              refine_field_elem[0] = markVal;
              ++did_change;
              break;
            }

        }
      return did_change;

    }

    int64_t mark_neighbor_tri_face(stk::mesh::Entity element, int refine_or_unrefine)
    {
      PerceptMesh& eMesh = Base::m_eMesh;

      if (!eMesh.owned(element))
        return 0;

      RefineFieldType *refine_field = Base::m_refine_field;
      bool isParent = (eMesh.hasFamilyTree(element) && eMesh.isParentElement(element));
      if (isParent)
        return 0;

      LocalSetType neighbors;
      int64_t did_change = 0;
      int *refine_field_elem = stk::mesh::field_data( *refine_field , element );

      if (1)
        {
          if (refine_field_elem[0] != 0)
            return 0;
        }

      static std::vector<stk::mesh::Entity> tri_face_neighbors;
      int iter = 0, iter_max=1000;
      while (true)
        {
          ++iter;
          VERIFY_OP_ON(iter, <=, iter_max, "mark_neighbor_tri_face too many iterations");

          refine_field_elem = stk::mesh::field_data( *refine_field , element );
          VERIFY_OP_ON(refine_field_elem[0], ==, 0, "bad refine_field_elem");
          tri_face_neighbors.resize(0);

          neighbors.clear();
          Base::get_node_neighbors(element, neighbors);

          bool element_changed = false;
          for (LocalSetType::iterator it = neighbors.begin(); it != neighbors.end(); ++it)
            {
              stk::mesh::Entity neigh = *it;
              if (neigh == element)
                continue;

              bool isNParent = (eMesh.hasFamilyTree(neigh) && eMesh.isParentElement(neigh));
              if (isNParent)
                continue;
              int *refine_field_neigh = stk::mesh::field_data( *refine_field , neigh );

              bool isFaceN = false;
              int face_0 = -1, face_1 = -1;
              isFaceN = eMesh.is_face_neighbor(element, neigh, &face_0, &face_1);

              int markVal = 0;
              if (isFaceN && face_0 >= 3)
                {
                  if (refine_or_unrefine == 1 && refine_field_neigh[0] >= 1)
                    {
                      markVal = 1;
                    }
                  else if (refine_or_unrefine == -1 && refine_field_neigh[0] == -1)
                    {
                      markVal = -1;
                    }
                  else
                    {
                      tri_face_neighbors.push_back(neigh);
                    }
                }
              if (markVal != 0 && refine_field_elem[0] != markVal)
                {
                  refine_field_elem[0] = markVal;
                  element_changed = true;
                  ++did_change;
                }

            }
          VERIFY_OP_ON(tri_face_neighbors.size(), <=, 2, "bad tri_face_neighbors");
          if (!element_changed)
            break;

          VERIFY_OP_ON(tri_face_neighbors.size(), <=, 1, "bad tri_face_neighbors should be 1");

          bool local_did_change = false;
          // FIXME no need for loop
          for (unsigned ii=0; ii < tri_face_neighbors.size(); ++ii)
            {
              stk::mesh::Entity neigh = tri_face_neighbors[ii];
              if (!eMesh.owned(neigh))
                continue;
              if (eMesh.topology(neigh) != stk::topology::WEDGE_6)
                continue;
              int *refine_field_neigh = stk::mesh::field_data( *refine_field , neigh );
              if (refine_field_neigh[0] == 0)
                {
                  element = neigh;
                  local_did_change = true;
                  break;
                }
              else
                {
                  VERIFY_MSG("bad tri_face_neighbors 3");
                }
            }
          if (!local_did_change)
            break;
        }
      return did_change;
    }

    int64_t enforce_boundary_layer_refine_pattern(int /*iter*/, int refine_or_unrefine)
    {
      PerceptMesh& eMesh = Base::m_eMesh;
      stk::ParallelMachine pm = eMesh.get_bulk_data()->parallel();
      int64_t did_change = 0;
      RefineFieldType *refine_field = Base::m_refine_field;
      if (!refine_field)
        {
          throw std::logic_error("must have refine_field field for transition element refinement");
        }

      TransitionElementType *transition_element_field = Base::m_transition_element_field;
      {
        std::vector< const stk::mesh::FieldBase *> fields;
        fields.push_back(refine_field);
        fields.push_back(transition_element_field);
        stk::mesh::communicate_field_data(eMesh.get_bulk_data()->aura_ghosting(), fields);
      }

      std::vector<stk::mesh::Entity> vec;
      VERIFY_OP_ON(m_wedge_selector, !=, 0, "must supply wedge selector to use TEA_SpecialWedgeRefinement");
      stk::mesh::get_selected_entities(*m_wedge_selector , eMesh.get_bulk_data()->buckets(eMesh.element_rank()), vec);

      did_change = 0;

      for (unsigned iElement = 0; iElement < vec.size(); iElement++)
        {
          stk::mesh::Entity element = vec[iElement];
          if (!eMesh.owned(element))
            continue;

          int *transition_element = stk::mesh::field_data( *transition_element_field , element );
          int *refine_field_elem = stk::mesh::field_data( *refine_field , element );

          bool isParent = (eMesh.hasFamilyTree(element) && eMesh.isParentElement(element));
          if (isParent)
            continue;
          if (transition_element[0] || refine_field_elem[0] != 0)
            continue;

          int64_t local_did_change = mark_neighbor_tri_face_first_pass(element, refine_or_unrefine);
          if (local_did_change)
            did_change += local_did_change;
        }
      int64_t global_did_change = did_change;
      stk::all_reduce( pm, stk::ReduceSum<1>( &global_did_change ) );

      if ((1||m_base_debug) && Base::m_eMesh.get_rank()==0)
        std::cout << "P[" << Base::m_eMesh.get_rank() << "] TEA_SpecialWedgeRefinement::enforce_boundary_layer_refine_pattern nelem= " << vec.size() << " # marked first pass= " << global_did_change << std::endl;

      if (!global_did_change) return 0;

      {
        std::vector< const stk::mesh::FieldBase *> fields;
        fields.push_back(refine_field);
        fields.push_back(transition_element_field);

        stk::mesh::communicate_field_data(eMesh.get_bulk_data()->aura_ghosting(), fields);
      }

      for (unsigned iElement = 0; iElement < vec.size(); iElement++)
        {
          stk::mesh::Entity element = vec[iElement];
          if (!eMesh.owned(element))
            continue;

          int *transition_element = stk::mesh::field_data( *transition_element_field , element );
          int *refine_field_elem = stk::mesh::field_data( *refine_field , element );

          bool isParent = (eMesh.hasFamilyTree(element) && eMesh.isParentElement(element));
          if (m_debug_local)
            {
              std::cout << "P[" << eMesh.get_rank() << "] tmp TEA isParent= " << isParent << " te= " << transition_element[0]
                        << " rf= " << refine_field_elem[0] << std::endl;
            }
          if (isParent)
            continue;
          if (transition_element[0] || refine_field_elem[0] == refine_or_unrefine)
            continue;

          int64_t local_did_change = mark_neighbor_tri_face(element, refine_or_unrefine);
          if (local_did_change)
            did_change += local_did_change;
        }

      {
        std::vector< const stk::mesh::FieldBase *> fields;
        fields.push_back(refine_field);
        fields.push_back(transition_element_field);

        stk::mesh::communicate_field_data(eMesh.get_bulk_data()->aura_ghosting(), fields);
      }

      stk::all_reduce( pm, stk::ReduceSum<1>( &did_change ) );

      return did_change;
    }


  };


}

#endif
