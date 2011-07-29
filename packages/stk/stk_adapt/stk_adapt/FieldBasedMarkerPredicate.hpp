#ifndef stk_adapt_FieldBasedMarkerPredicate_hpp
#define stk_adapt_FieldBasedMarkerPredicate_hpp

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
    struct FieldBasedMarkerPredicate : public std::unary_function<const stk::mesh::Entity& , bool> {
      stk::mesh::Selector& m_selector;
      stk::mesh::FieldBase *m_field;
      double m_tolerance;
      FieldBasedMarkerPredicate(stk::mesh::Selector& selector, stk::mesh::FieldBase *field, double tolerance) :
        m_selector(selector), m_field(field), m_tolerance(tolerance) {}
    };

    struct ElementFieldBasedRefinePredicate : public FieldBasedMarkerPredicate {
      
      ElementFieldBasedRefinePredicate(stk::mesh::Selector& selector, stk::mesh::FieldBase *field, double tolerance) :
        FieldBasedMarkerPredicate(selector, field, tolerance) {}

      /// Return true for refine, false for ignore
      bool operator()(const stk::mesh::Entity& entity) {
        double *fdata = stk::mesh::field_data( *static_cast<const ScalarFieldType *>(m_field) , entity );
        return m_selector(entity) && fdata[0] > 0;
      }
    };

    struct ElementFieldBasedUnrefinePredicate : public FieldBasedMarkerPredicate {
      
      ElementFieldBasedUnrefinePredicate(stk::mesh::Selector& selector, stk::mesh::FieldBase *field, double tolerance) :
        FieldBasedMarkerPredicate(selector, field, tolerance) {}

      /// Return true for unrefine, false for ignore
      bool operator()(const stk::mesh::Entity& entity) {
        double *fdata = stk::mesh::field_data( *static_cast<const ScalarFieldType *>(m_field) , entity );
        return m_selector(entity) && fdata[0] < 0;
      }
    };


    /** Example of how the above is used (pseudo-code):
     *
     * class PredicateMarker : public Refiner {
     *   void visit_for_refine_and_unrefin(std::unary_function<stk::mesh::Entity& , bool>& user_predicate_refine,
     *                                     std::unary_function<stk::mesh::Entity& , bool>& user_predicate_unrefine)
     *     {
     *        foreach (Entity entity in list)
     *          {
     *             if (user_predicate_refine(entity)) { do_refine_operation(entity); }
     *          }
     *        foreach (Entity entity in list)
     *          {
     *             if (user_predicate_unrefine(entity)) { do_unrefine_operation(entity); }
     *          }
     *     }
     * };
     */

    //get_selected_entities(Selector, rank, EntityVec vec);


    //filter_iterator( Predicate&, vec.begin());
  }
}
#endif
