#ifndef stk_adapt_RefinerPattern_Hex8_Hex8_N_HangingNode_hpp
#define stk_adapt_RefinerPattern_Hex8_Hex8_N_HangingNode_hpp

#include <stk_adapt/sierra_element/RefinementTopology.hpp>
#include <stk_adapt/sierra_element/StdMeshObjTopologies.hpp>

#include "RefinerPattern_Quad4_Quad4_HangingNode.hpp"

#include <stk_percept/PerceptBoostArray.hpp>

namespace stk {
  namespace adapt {

    template <>
    class RefinerPattern<shards::Hexahedron<8>, shards::Hexahedron<8>, -1 > : public URP<shards::Hexahedron<8>,shards::Hexahedron<8>  >
    {

      RefinerPattern<shards::Quadrilateral<4>, shards::Quadrilateral<4>, -1 > * m_face_breaker;

    public:

      RefinerPattern(percept::PerceptMesh& eMesh, BlockNamesType block_names = BlockNamesType()) :  URP<shards::Hexahedron<8>, shards::Hexahedron<8>  >(eMesh)
      {
        m_primaryEntityRank = stk::mesh::MetaData::ELEMENT_RANK;

        setNeededParts(eMesh, block_names, true);
        Elem::StdMeshObjTopologies::bootstrap();

        m_face_breaker =  new RefinerPattern<shards::Quadrilateral<4>, shards::Quadrilateral<4>, -1 > (eMesh, block_names) ;


      }

      ~RefinerPattern()
      {
        if (m_face_breaker) delete m_face_breaker;
      }

      void setSubPatterns( std::vector<UniformRefinerPatternBase *>& bp, percept::PerceptMesh& eMesh )
      {
        EXCEPTWATCH;
        bp = std::vector<UniformRefinerPatternBase *>(2u, 0);

        bp[0] = this;
        bp[1] = m_face_breaker;
      }

      virtual void doBreak() {}
      void fillNeededEntities(std::vector<NeededEntityType>& needed_entities)
      {
        needed_entities.resize(3);
        needed_entities[0].first = m_eMesh.edge_rank();  
        needed_entities[1].first = m_eMesh.face_rank();
        needed_entities[2].first = stk::mesh::MetaData::ELEMENT_RANK;
        setToOne(needed_entities);

      }

      virtual unsigned getNumNewElemPerElem() { return 8; }

      void 
      createNewElements(percept::PerceptMesh& eMesh, NodeRegistry& nodeRegistry, 
                        stk::mesh::Entity element,  NewSubEntityNodesType& new_sub_entity_nodes, vector<stk::mesh::Entity>::iterator& element_pool,
                        stk::mesh::FieldBase *proc_rank_field=0)
      {
        unsigned num_edges_marked=0;
        for (int iedge = 0; iedge < 8; iedge++)
          {
            unsigned num_nodes_on_edge = new_sub_entity_nodes[m_eMesh.edge_rank()][iedge].size();
            if (num_nodes_on_edge)
              {
                ++num_edges_marked;
              }
          }
        if (num_edges_marked != 8)
          {
            return;
          }

        genericRefine_createNewElements(eMesh, nodeRegistry,
                                        element, new_sub_entity_nodes, element_pool,
                                        proc_rank_field);
      }
      
    };

  }
}
#endif
