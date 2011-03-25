#ifndef stk_adapt_UniformRefinerPattern_Beam3_Beam3_2_sierra_hpp
#define stk_adapt_UniformRefinerPattern_Beam3_Beam3_2_sierra_hpp


//#include "UniformRefinerPattern.hpp"
#include <stk_adapt/sierra_element/RefinementTopology.hpp>
#include <stk_adapt/sierra_element/StdMeshObjTopologies.hpp>

namespace stk {
  namespace adapt {

    template <>
    class UniformRefinerPattern<shards::Beam<3>, shards::Beam<3>, 2, SierraPort > : public URP<shards::Beam<3>, shards::Beam<3>  >
    {
    public:

      UniformRefinerPattern(percept::PerceptMesh& eMesh, BlockNamesType block_names = BlockNamesType()) :  URP<shards::Beam<3>, shards::Beam<3>  >(eMesh)
      {
        //         m_primaryEntityRank = mesh::Edge;
        //         if (m_eMesh.getSpatialDim() == 1)
        m_primaryEntityRank = eMesh.element_rank();

        setNeededParts(eMesh, block_names, true);
      }

      virtual void doBreak() {}
      void fillNeededEntities(std::vector<NeededEntityType>& needed_entities)
      {
        needed_entities.resize(1);
        needed_entities[0].first = m_eMesh.edge_rank(); 
        needed_entities[0].second = 3u;
      }

      virtual unsigned getNumNewElemPerElem() { return 2; }

      void 
      createNewElements(percept::PerceptMesh& eMesh, NodeRegistry& nodeRegistry, 
                        Entity& element,  NewSubEntityNodesType& new_sub_entity_nodes, vector<Entity *>::iterator& element_pool,
                        FieldBase *proc_rank_field=0)
      {
        genericRefine_createNewElements(eMesh, nodeRegistry,
                                        element, new_sub_entity_nodes, element_pool,
                                        proc_rank_field);
      }
      
    };

  }
}
#endif
