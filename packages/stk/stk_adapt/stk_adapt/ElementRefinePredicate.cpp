
#include <stk_adapt/ElementRefinePredicate.hpp>

namespace stk {
  namespace adapt {

    /// Return DO_REFINE, DO_UNREFINE, DO_NOTHING
    int ElementRefinePredicate::operator()(const stk::mesh::Entity entity) 
    {
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

      // void check_two_to_one(PerceptMesh& eMesh);
      // void enforce_two_to_one_refine(PerceptMesh& eMesh);
      // void ok_to_unrefine(PerceptMesh& eMesh);

    

    bool check_two_to_one(PerceptMesh& eMesh)
    {
      
    }



  }
}

