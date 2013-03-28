#ifndef stk_adapt_IElementAdapter_hpp
#define stk_adapt_IElementAdapter_hpp

#include <stk_adapt/IAdapter.hpp>

namespace stk {
  namespace adapt {

    //========================================================================================================================
    //========================================================================================================================
    //========================================================================================================================
    /**
     * An IElementAdapter is an abstract base class for derived classes that are required to overload the mark method,
     *   which supplies the derived class with the element to be marked for refine, unrefine, or both (@see IAdapter::AdaptInstruction)
     */
    class IElementAdapter : public IAdapter
    {
    public:
      IElementAdapter(percept::PerceptMesh& eMesh, UniformRefinerPatternBase & bp, stk::mesh::FieldBase *proc_rank_field=0)
      : IAdapter(eMesh, bp, proc_rank_field) {}

      virtual ElementUnrefineCollection  buildUnrefineList() ;

    protected:

      /// Client supplies this method - given an element return instruction on what to do to the element:
      ///    DO_NOTHING (nothing), DO_REFINE (refine), DO_UNREFINE
      virtual int mark(const stk::mesh::Entity element) = 0;

      virtual void
      refineMethodApply(NodeRegistry::ElementFunctionPrototype function, const stk::mesh::Entity element,
                                              vector<NeededEntityType>& needed_entity_ranks);

    };
  }
}
#endif
