#ifndef stk_adapt_TestLocalRefiner_hpp
#define stk_adapt_TestLocalRefiner_hpp

#include <stk_adapt/Refiner.hpp>

namespace stk_classic {
  namespace adapt {

    //========================================================================================================================
    //========================================================================================================================
    //========================================================================================================================
    /**
     * A test implementation that does uniform refinement but uses non-uniform methods
     */
    class TestLocalRefiner : public Refiner
    {
    public:
      TestLocalRefiner(percept::PerceptMesh& eMesh, UniformRefinerPatternBase & bp, stk_classic::mesh::FieldBase *proc_rank_field=0);

    protected:
      virtual unsigned
      doForAllElements(stk_classic::mesh::EntityRank rank, NodeRegistry::ElementFunctionPrototype function, 
                       vector< ColorerSetType >& elementColors, unsigned elementType,
                       vector<NeededEntityType>& needed_entity_ranks,
                       bool only_count=false, bool doAllElements=true);


      virtual void 
      refineMethodApply(NodeRegistry::ElementFunctionPrototype function, const stk_classic::mesh::Entity& element, 
            vector<NeededEntityType>& needed_entity_ranks);


    };



  }
}
#endif
