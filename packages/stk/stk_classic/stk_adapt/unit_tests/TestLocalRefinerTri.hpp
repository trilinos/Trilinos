#ifndef stk_adapt_TestLocalRefinerTri_hpp
#define stk_adapt_TestLocalRefinerTri_hpp

#include <stk_adapt/Refiner.hpp>

namespace stk_classic {
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
      TestLocalRefinerTri(percept::PerceptMesh& eMesh, UniformRefinerPatternBase & bp, stk_classic::mesh::FieldBase *proc_rank_field=0);

    protected:


      virtual void 
      refineMethodApply(NodeRegistry::ElementFunctionPrototype function, const stk_classic::mesh::Entity& element, 
            vector<NeededEntityType>& needed_entity_ranks);


    };



  }
}
#endif
