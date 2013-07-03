#ifndef stk_adapt_ElementRefinePredicate_hpp
#define stk_adapt_ElementRefinePredicate_hpp

#include <stk_adapt/IElementBasedAdapterPredicate.hpp>

namespace stk {
  namespace adapt {

    // Can be instantiated by the user, or used to define your own
    struct ElementRefinePredicate : public IElementBasedAdapterPredicate {
      //PerceptMesh& m_eMesh;

      ElementRefinePredicate(PerceptMesh& eMesh, stk::mesh::Selector* selector=0, stk::mesh::FieldBase *field=0, double tolerance=0.0) :
        IElementBasedAdapterPredicate(eMesh, selector, field, tolerance) {}

      /// Return DO_REFINE, DO_UNREFINE, DO_NOTHING
      int operator()(const stk::mesh::Entity entity);

      /// wrappers around refine and unrefine
      /// @param enforce_2_to_1_what: {node-, edge-, face-} neighbors
      void refine(IAdapter& breaker, bool enforce_2_to_1_what[3]);
      void unrefine(IAdapter& breaker, bool enforce_2_to_1_what[3]);

      // low-level interface:

      /** return true if mesh is 2-to-1 rule compliant
       * check_what: {node-, edge-, face-} neighbors
       */
      bool check_two_to_one(bool check_what[3]);

      /** after marking using field "refine_field", revisit all elements and upgrade
       * any that need it to enforce two to one
       * @return true if there was a change to refine_field; useful in a while(did_change) loop
       * @param enforce_what: {node-, edge-, face-} neighbors
       */
      bool enforce_two_to_one_refine(bool enforce_what[3]);

      /** after marking using field "refine_field", revisit all elements and downgrade
       * any that need it to enforce two to one during unrefinement
       * @return true if there was a change to refine_field; useful in a while(did_change) loop
       * @param enforce_what: {node-, edge-, face-} neighbors
       */
      bool enforce_two_to_one_unrefine(bool enforce_what[3]);

    private:
      bool is_face_neighbor(stk::mesh::Entity element, int element_level, stk::mesh::Entity neighbor, int neighbor_level);
      bool is_edge_neighbor(stk::mesh::Entity element, int element_level, stk::mesh::Entity neighbor, int neighbor_level);
      bool is_node_neighbor(stk::mesh::Entity element, int element_level, stk::mesh::Entity neighbor, int neighbor_level);
      bool min_max_neighbors_level(stk::mesh::Entity element, int min_max[2], ScalarIntFieldType *refine_level, bool check_what[3] );
      void get_neighbors(stk::mesh::Entity element, ScalarIntFieldType *refine_level, bool get_what[3],
                         std::set<stk::mesh::Entity>& selected_neighbors);
    };

    /** Example of how the above is used (pseudo-code):
     *
     * class PredicateAdapter : public Refiner {
     *   void visit_for_refine_and_unrefin(std::unary_function<stk::mesh::Entity , bool>& user_predicate_refine)
     *     {
     *        foreach (Entity entity in list)
     *          {
     *             if (user_predicate_refine(entity) & DO_REFINE) { do_refine_operation(entity); }
     *          }
     *        foreach (Entity entity in list)
     *          {
     *             if (user_predicate_unrefine(entity) & DO_UNREFINE) { do_unrefine_operation(entity); }
     *          }
     *     }
     * };
     */

  }
}
#endif
