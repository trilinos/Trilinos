// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef adapt_IElementBasedAdapterPredicate_hpp
#define adapt_IElementBasedAdapterPredicate_hpp

#include <adapt/IAdapter.hpp>
#include <percept/PerceptMesh.hpp>

#include <functional>

  namespace percept {

    /** Signatures for predicate objects that can be used to select entities (elements, edges, faces,...) for
     *   refinement or unrefinement.  The class must supply an operator() that takes an entity and decides
     *   on whether it should be refined, or unrefined, or ignored.
     *
     * The Selector pattern as shown below is useful for selecting only entities that belong to particular
     *   mesh Parts, for example, or any other definition of Selector.
     *
     * We follow the unary_function pattern to enable the structs to be used in STL algorithms that know about
     *   unary_functions.
     *
     * The following are a couple of examples of what refine and unrefine predicates might look like.
     *
     * Following these examples we show the prototype for the operations that are performed on these predicates.
     */

    // Example
    struct IElementBasedAdapterPredicate : public std::function<int(const stk::mesh::Entity)> {
      PerceptMesh& m_eMesh;
      stk::mesh::Selector * m_eb_selector;
      stk::mesh::FieldBase *m_field;
      double m_tolerance;

      int m_ref_field_stage;


    protected:

    protected:
      IElementBasedAdapterPredicate(PerceptMesh& eMesh, stk::mesh::Selector * selector=0, stk::mesh::FieldBase *field=0, double tolerance=0.0) :
        m_eMesh(eMesh), m_eb_selector(selector), m_field(field), m_tolerance(tolerance),
        m_ref_field_stage(-1), m_mark_centroid_always(false)
      {
      }

    public:
      // This is for making the NodeRegistry consistent when doing multi-refine (hanging-node + transition element approach).
      // It forces a return mark of nothing to allow NodeRegistry looping without new marks.
      void setMarkNone(bool val) { m_eMesh.m_markNone = val; }
      bool getMarkNone() { return m_eMesh.m_markNone; }

      // allow some finer control for refine field values
      void setRefineStage(int x) { m_ref_field_stage = x; }
      int getRefineStage() { return m_ref_field_stage; }

      bool m_mark_centroid_always;

    };


    //get_selected_entities(Selector, rank, EntityVec vec);


    //filter_iterator( Predicate&, vec.begin());
  }

#endif
