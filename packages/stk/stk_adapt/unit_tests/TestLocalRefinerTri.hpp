#ifndef stk_adapt_TestLocalRefinerTri_hpp
#define stk_adapt_TestLocalRefinerTri_hpp

#include <stk_adapt/Refiner.hpp>

namespace stk {
  namespace adapt {

    //========================================================================================================================
    //========================================================================================================================
    //========================================================================================================================
    /**
     * A test implementation that does uniform refinement but uses non-uniform methods
     */
    class TestLocalRefinerTri : public Refiner
    {
    public:
      TestLocalRefinerTri(percept::PerceptMesh& eMesh, UniformRefinerPatternBase & bp, stk::mesh::FieldBase *proc_rank_field=0);

    protected:

      // not needed
#if 0
      virtual unsigned
      doForAllElements(stk::mesh::EntityRank rank, NodeRegistry::ElementFunctionPrototype function, 
                       vector< ColorerSetType >& elementColors, unsigned elementType,
                       vector<NeededEntityType>& needed_entity_ranks,
                       bool only_count=false, bool doAllElements=true);
#endif


      virtual void 
      applyNodeRegistryFunctionForSubEntities(NodeRegistry::ElementFunctionPrototype function, const stk::mesh::Entity& element, 
                                              vector<NeededEntityType>& needed_entity_ranks);


    };



  }
}
#endif
