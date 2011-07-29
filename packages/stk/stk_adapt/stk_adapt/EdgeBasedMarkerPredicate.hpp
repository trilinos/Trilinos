#ifndef stk_adapt_EdgeBasedMarkerPredicate_hpp
#define stk_adapt_EdgeBasedMarkerPredicate_hpp

#include <stk_adapt/Refiner.hpp>
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
    struct EdgeBasedMarkerPredicate : public std::unary_function<const stk::mesh::Entity& , bool> {
      stk::mesh::Selector& m_selector;
      stk::mesh::FieldBase *m_field;
      double m_tolerance;
      EdgeBasedMarkerPredicate(stk::mesh::Selector& selector, stk::mesh::FieldBase *field, double tolerance) :
        m_selector(selector), m_field(field), m_tolerance(tolerance) {}
    };

    struct EdgeBasedRefinePredicate : public EdgeBasedMarkerPredicate {
      
      EdgeBasedRefinePredicate(stk::mesh::Selector& selector, stk::mesh::FieldBase *field, double tolerance) :
        EdgeBasedMarkerPredicate(selector, field, tolerance) {}

      /// Return true for refine, false for ignore
      bool operator()(const stk::mesh::Entity& element, unsigned which_edge, stk::mesh::Entity & node0, stk::mesh::Entity & node1,
                      double *coord0, double *coord1, std::vector<int>& existing_edge_marks);

      //double *fdata = stk::mesh::field_data( *static_cast<const ScalarFieldType *>(m_field) , entity );
      //return m_selector(entity) && fdata[0] > 0;
    };

    struct EdgeBasedUnrefinePredicate : public EdgeBasedMarkerPredicate {
      
      EdgeBasedUnrefinePredicate(stk::mesh::Selector& selector, stk::mesh::FieldBase *field, double tolerance) :
        EdgeBasedMarkerPredicate(selector, field, tolerance) {}

      /// Return true for unrefine, false for ignore
      bool operator()(const stk::mesh::Entity& element, unsigned which_edge, stk::mesh::Entity & node0, stk::mesh::Entity & node1,
                      double *coord0, double *coord1);

    };

  }
}
#endif
