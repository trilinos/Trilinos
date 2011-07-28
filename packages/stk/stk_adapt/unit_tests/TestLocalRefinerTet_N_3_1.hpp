#ifndef stk_adapt_TestLocalRefinerTet_N_3_1_hpp
#define stk_adapt_TestLocalRefinerTet_N_3_1_hpp

#include <stk_adapt/Refiner.hpp>

namespace stk {
  namespace adapt {

    //========================================================================================================================
    //========================================================================================================================
    //========================================================================================================================
    /**
     * A test implementation that marks some edges randomly to test RefinerPattern_Tri3_Tri3_N_3_1
     */
    class TestLocalRefinerTet_N_3_1 : public Refiner
    {
    public:
      TestLocalRefinerTet_N_3_1(percept::PerceptMesh& eMesh, UniformRefinerPatternBase & bp, stk::mesh::FieldBase *proc_rank_field=0, unsigned edge_mark_bitcode=1);

      // ElementUnrefineCollection  buildTestUnrefineList();

    protected:


      virtual void 
      apply(NodeRegistry::ElementFunctionPrototype function, const stk::mesh::Entity& element, 
            vector<NeededEntityType>& needed_entity_ranks);


      unsigned m_edge_mark_bitcode;
    };

  }
}
#endif
