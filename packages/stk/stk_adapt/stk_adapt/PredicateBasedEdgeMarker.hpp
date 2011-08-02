#ifndef stk_adapt_PredicateBasedEdgeMarker_hpp
#define stk_adapt_PredicateBasedEdgeMarker_hpp

#include <functional>

#include <stk_adapt/EdgeMarker.hpp>

namespace stk {
  namespace adapt {

    //========================================================================================================================
    //========================================================================================================================
    //========================================================================================================================
    /**
     *  Predicate-based marker
     */

    template<class RefinePredicate, class UnrefinePredicate>
    class PredicateBasedEdgeMarker : public EdgeMarker
    {
      RefinePredicate m_predicate_refine;
      UnrefinePredicate m_predicate_unrefine;

    public:

      PredicateBasedEdgeMarker(RefinePredicate predicate_refine, UnrefinePredicate predicate_unrefine,
                               percept::PerceptMesh& eMesh, UniformRefinerPatternBase & bp, stk::mesh::FieldBase *proc_rank_field=0) :
        EdgeMarker(eMesh, bp, proc_rank_field), m_predicate_refine(predicate_refine), m_predicate_unrefine(predicate_unrefine)
      {
      }

      ///    1 (refine), 0 (nothing)
      virtual int mark(const stk::mesh::Entity& element, unsigned which_edge, stk::mesh::Entity & node0, stk::mesh::Entity & node1,
                           double *coord0, double *coord1, std::vector<int>& existing_edge_marks) 
      {
        bool do_ref = m_predicate_refine(element, which_edge, node0, node1, coord0, coord1, existing_edge_marks);
        return do_ref ? 1 : 0;
      }

      ///    -1 (unrefine), 0 (nothing)
      virtual int markUnrefine(const stk::mesh::Entity& element, unsigned which_edge, stk::mesh::Entity & node0, stk::mesh::Entity & node1,
                                double *coord0, double *coord1) 
      {
        bool do_unref = m_predicate_unrefine(element, which_edge, node0, node1, coord0, coord1);
        return do_unref ? -1 : 0;
      }
      
    };


  }
}
#endif
