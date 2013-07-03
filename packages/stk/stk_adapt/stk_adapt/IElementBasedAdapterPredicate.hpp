#ifndef stk_adapt_IElementBasedAdapterPredicate_hpp
#define stk_adapt_IElementBasedAdapterPredicate_hpp

#include <stk_adapt/IAdapter.hpp>
#include <stk_percept/PerceptMesh.hpp>

#include <functional>

namespace stk {
  namespace adapt {

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
    struct IElementBasedAdapterPredicate : public std::unary_function<const stk::mesh::Entity , int> {
      PerceptMesh& m_eMesh;
      stk::mesh::Selector * m_eb_selector;
      stk::mesh::FieldBase *m_field;
      double m_tolerance;
    protected:
      IElementBasedAdapterPredicate(PerceptMesh& eMesh, stk::mesh::Selector * selector=0, stk::mesh::FieldBase *field=0, double tolerance=0.0) :
        m_eMesh(eMesh), m_eb_selector(selector), m_field(field), m_tolerance(tolerance) {}

    };


    //get_selected_entities(Selector, rank, EntityVec vec);


    //filter_iterator( Predicate&, vec.begin());
  }
}
#endif
