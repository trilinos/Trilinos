// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef adapt_PredicateBasedEdgeAdapter_hpp
#define adapt_PredicateBasedEdgeAdapter_hpp

#include <functional>

#include <adapt/IEdgeAdapter.hpp>

  namespace percept {

    //========================================================================================================================
    //========================================================================================================================
    //========================================================================================================================
    /**
     *  Predicate-based marker
     *
     *  The functor @class RefinePredicate should supply an operator() that returns an entry from AdaptInstruction,
     *    either to do nothing, refine, unrefine, or both refine & unrefine (useful for unit testing, etc.)
     */

    template<class RefinePredicate>
    class PredicateBasedEdgeAdapter : public IEdgeAdapter
    {
      RefinePredicate& m_predicate_refine;

    public:

      PredicateBasedEdgeAdapter(RefinePredicate& predicate_refine,
                                percept::PerceptMesh& eMesh, UniformRefinerPatternBase & bp, stk::mesh::FieldBase *proc_rank_field=0) :
        IEdgeAdapter(eMesh, bp, proc_rank_field), m_predicate_refine(predicate_refine)
      {
      }

      RefinePredicate& getRefinePredicate() { return m_predicate_refine; }

      ///    DO_NOTHING (nothing), DO_REFINE (refine), DO_UNREFINE
      virtual int markEdge(const stk::mesh::Entity element, unsigned which_edge, stk::mesh::Entity node0, stk::mesh::Entity node1,
                           double *coord0, double *coord1, std::vector<int>* existing_edge_marks)
      {
        int mark = m_predicate_refine(element, which_edge, node0, node1, coord0, coord1, existing_edge_marks);
        return mark;
      }


    };


  }

#endif
