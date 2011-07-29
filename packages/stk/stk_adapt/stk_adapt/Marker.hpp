#ifndef stk_adapt_Marker_hpp
#define stk_adapt_Marker_hpp

#include <stk_adapt/Refiner.hpp>

namespace stk {
  namespace adapt {

    //========================================================================================================================
    //========================================================================================================================
    //========================================================================================================================
    /**
     * A Marker is an abstract base class that provides all functionality of a Refiner but requires sub-classes to 
     *   overload the apply() and buildUnrefineList() methods.
     */
    class Marker : public Refiner
    {
    public:
      Marker(percept::PerceptMesh& eMesh, UniformRefinerPatternBase & bp, stk::mesh::FieldBase *proc_rank_field=0);

      virtual ElementUnrefineCollection  buildUnrefineList() = 0;

    protected:

      virtual void 
      apply(NodeRegistry::ElementFunctionPrototype function, const stk::mesh::Entity& element, 
                                              vector<NeededEntityType>& needed_entity_ranks) = 0;


    };

    Marker::Marker(percept::PerceptMesh& eMesh, UniformRefinerPatternBase &  bp, stk::mesh::FieldBase *proc_rank_field) : 
      Refiner(eMesh, bp, proc_rank_field)
    {
    }


  }
}
#endif
