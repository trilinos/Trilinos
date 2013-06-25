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
      int operator()(const stk::mesh::Entity entity) {
        double *fdata = 0;
        if (m_field)
          fdata = m_eMesh.field_data( *static_cast<const ScalarFieldType *>(m_field) , entity );
        bool selected = (m_eb_selector==0 || (*m_eb_selector)(m_eMesh.bucket(entity)));
        bool ref_field_criterion = (fdata  && fdata[0] > 0);
        bool unref_field_criterion = (fdata && fdata[0] < 0);
        int mark = 0;
        if (selected && ref_field_criterion) mark |= DO_REFINE;
        if (selected && unref_field_criterion) mark |= DO_UNREFINE;
        return mark;
      }

      void check_two_to_one(PerceptMesh& eMesh);
      void enforce_two_to_one_refine(PerceptMesh& eMesh);
      void ok_to_unrefine(PerceptMesh& eMesh);



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
