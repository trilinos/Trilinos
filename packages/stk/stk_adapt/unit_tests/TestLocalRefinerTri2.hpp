#ifndef stk_adapt_TestLocalRefinerTri2_hpp
#define stk_adapt_TestLocalRefinerTri2_hpp

#include <stk_adapt/Refiner.hpp>

namespace stk {
  namespace adapt {

    //========================================================================================================================
    //========================================================================================================================
    //========================================================================================================================
    /**
     * A test implementation that does uniform refinement but uses non-uniform methods
     */
    class TestLocalRefinerTri2 : public Refiner
    {
      bool m_diagonals;
    public:
      TestLocalRefinerTri2(percept::PerceptMesh& eMesh, UniformRefinerPatternBase & bp, stk::mesh::FieldBase *proc_rank_field=0, bool diagonals=true);

    protected:

      virtual void 
      apply(NodeRegistry::ElementFunctionPrototype function, const stk::mesh::Entity& element, 
            vector<NeededEntityType>& needed_entity_ranks);


    };



  }
}
#endif
