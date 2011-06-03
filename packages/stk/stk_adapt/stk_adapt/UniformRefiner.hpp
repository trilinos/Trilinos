#ifndef stk_adapt_UniformRefiner_hpp
#define stk_adapt_UniformRefiner_hpp

#include <stk_adapt/Refiner.hpp>

namespace stk {
  namespace adapt {

    //========================================================================================================================
    //========================================================================================================================
    //========================================================================================================================
    //template<class UniformRefinerPattern>
    class UniformRefiner : public Refiner
    {
    public:
      UniformRefiner(percept::PerceptMesh& eMesh, UniformRefinerPatternBase & bp, stk::mesh::FieldBase *proc_rank_field=0);

    protected:
    };



  }
}
#endif
