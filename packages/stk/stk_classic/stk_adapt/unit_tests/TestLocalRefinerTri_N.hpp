#ifndef stk_adapt_TestLocalRefinerTri_N_hpp
#define stk_adapt_TestLocalRefinerTri_N_hpp

#include <stk_adapt/Refiner.hpp>

namespace stk_classic {
  namespace adapt {

    //========================================================================================================================
    //========================================================================================================================
    //========================================================================================================================
    /**
     * A test implementation that marks some edges randomly to test RefinerPattern_Tri3_Tri3_N
     */
    class TestLocalRefinerTri_N : public Refiner
    {
    public:
      TestLocalRefinerTri_N(percept::PerceptMesh& eMesh, UniformRefinerPatternBase & bp, stk_classic::mesh::FieldBase *proc_rank_field=0);

      ElementUnrefineCollection  buildTestUnrefineList();

    protected:

      virtual void 
      refineMethodApply(NodeRegistry::ElementFunctionPrototype function, const stk_classic::mesh::Entity& element, 
            vector<NeededEntityType>& needed_entity_ranks);


    };



  }
}
#endif
