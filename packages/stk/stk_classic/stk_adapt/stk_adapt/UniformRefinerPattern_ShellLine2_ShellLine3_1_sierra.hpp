#ifndef stk_adapt_UniformRefinerPattern_ShellLine2_ShellLine3_1_sierra_hpp
#define stk_adapt_UniformRefinerPattern_ShellLine2_ShellLine3_1_sierra_hpp


//#include "UniformRefinerPattern.hpp"
#include <stk_adapt/sierra_element/RefinementTopology.hpp>
#include <stk_adapt/sierra_element/StdMeshObjTopologies.hpp>

namespace stk {
  namespace adapt {

    template <>
    class UniformRefinerPattern<shards::ShellLine<2>, shards::ShellLine<3>, 1, SierraPort > : public URP<shards::ShellLine<2>, shards::ShellLine<3>  >
    {
    public:

      UniformRefinerPattern(percept::PerceptMesh& eMesh, BlockNamesType block_names = BlockNamesType()) :  URP<shards::ShellLine<2>, shards::ShellLine<3>  >(eMesh)
      {
        m_primaryEntityRank = m_eMesh.edge_rank();
        if (m_eMesh.get_spatial_dim() == 1)
          m_primaryEntityRank = eMesh.element_rank();

        setNeededParts(eMesh, block_names, false);  // different topologies
        Elem::StdMeshObjTopologies::bootstrap();
      }

      virtual void doBreak() {}

      void setSubPatterns( std::vector<UniformRefinerPatternBase *>& bp, percept::PerceptMesh& eMesh )
      {
        EXCEPTWATCH;
        bp = std::vector<UniformRefinerPatternBase *>(1u, 0);

        if (eMesh.get_spatial_dim() == 1)
          {
            bp[0] = this;
          }
        else 
          {
          }

      }

      void fillNeededEntities(std::vector<NeededEntityType>& needed_entities)
      {
        needed_entities.resize(1);
        //needed_entities[0].first = m_eMesh.edge_rank(); 
        needed_entities[0].first = (m_eMesh.get_spatial_dim() == 1 ? m_eMesh.element_rank() : m_eMesh.edge_rank());
        setToOne(needed_entities);
      }

      virtual unsigned getNumNewElemPerElem() { return 1; }

      void 
      createNewElements(percept::PerceptMesh& eMesh, NodeRegistry& nodeRegistry, 
                        stk::mesh::Entity& element,  NewSubEntityNodesType& new_sub_entity_nodes, vector<stk::mesh::Entity *>::iterator& element_pool,
                        stk::mesh::FieldBase *proc_rank_field=0)
      {
        genericEnrich_createNewElements(eMesh, nodeRegistry,
                                        element, new_sub_entity_nodes, element_pool,
                                        proc_rank_field);
        
      }
      
    };

  }
}
#endif
