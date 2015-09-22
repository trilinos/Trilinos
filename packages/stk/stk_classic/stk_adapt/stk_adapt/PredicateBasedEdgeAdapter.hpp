#ifndef stk_adapt_PredicateBasedEdgeAdapter_hpp
#define stk_adapt_PredicateBasedEdgeAdapter_hpp

#include <functional>

#include <stk_adapt/IEdgeAdapter.hpp>

namespace stk_classic {
  namespace adapt {

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
                                percept::PerceptMesh& eMesh, UniformRefinerPatternBase & bp, stk_classic::mesh::FieldBase *proc_rank_field=0) :
        IEdgeAdapter(eMesh, bp, proc_rank_field), m_predicate_refine(predicate_refine)
      {
      }

      RefinePredicate& getRefinePredicate() { return m_predicate_refine; }

      ///    DO_NOTHING (nothing), DO_REFINE (refine), DO_UNREFINE
      virtual int mark(const stk_classic::mesh::Entity& element, unsigned which_edge, stk_classic::mesh::Entity & node0, stk_classic::mesh::Entity & node1,
                           double *coord0, double *coord1, std::vector<int>* existing_edge_marks) 
      {
        int mark = m_predicate_refine(element, which_edge, node0, node1, coord0, coord1, existing_edge_marks);
        return mark;
      }

      
    };


  }
}
#endif
