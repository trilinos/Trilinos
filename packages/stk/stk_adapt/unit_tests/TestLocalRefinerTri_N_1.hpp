#ifndef stk_adapt_TestLocalRefinerTri_N_1_hpp
#define stk_adapt_TestLocalRefinerTri_N_1_hpp

#include <stk_adapt/Refiner.hpp>

namespace stk {
  namespace adapt {

    //========================================================================================================================
    //========================================================================================================================
    //========================================================================================================================
    /**
     * A test implementation that marks some edges randomly to test RefinerPattern_Tri3_Tri3_N_1
     */
    class TestLocalRefinerTri_N_1 : public Refiner
    {
    public:
      TestLocalRefinerTri_N_1(percept::PerceptMesh& eMesh, UniformRefinerPatternBase & bp, stk::mesh::FieldBase *proc_rank_field=0);

      ElementUnrefineCollection  buildTestUnrefineList();

    protected:

      // not needed

      virtual void 
      apply(NodeRegistry::ElementFunctionPrototype function, const stk::mesh::Entity& element, 
                                              vector<NeededEntityType>& needed_entity_ranks);


    };



  }
}
#endif
