#ifndef stk_adapt_UniformRefinerPattern_ShellQuad4_ShellQuad9_1_sierra_hpp
#define stk_adapt_UniformRefinerPattern_ShellQuad4_ShellQuad9_1_sierra_hpp


#include <stk_adapt/sierra_element/RefinementTopology.hpp>
#include <stk_adapt/sierra_element/StdMeshObjTopologies.hpp>

#define EDGE_QU4_QU9_1_BREAKER 1
#if EDGE_QU4_QU9_1_BREAKER
#include "UniformRefinerPattern_ShellLine2_ShellLine3_1_sierra.hpp"
#include "UniformRefinerPattern_Quad4_Quad9_1_sierra.hpp"
#endif

namespace stk {
  namespace adapt {

    template <>
    class UniformRefinerPattern<shards::ShellQuadrilateral<4>, shards::ShellQuadrilateral<9>, 1, SierraPort > : 
      public URP<shards::ShellQuadrilateral<4>,shards::ShellQuadrilateral<9>  >
    {


#if EDGE_QU4_QU9_1_BREAKER
      UniformRefinerPattern<shards::ShellLine<2>, shards::ShellLine<3>, 1, SierraPort > * m_edge_breaker;
      UniformRefinerPattern<shards::Quadrilateral<4>, shards::Quadrilateral<9>, 1, SierraPort > * m_face_breaker;
#endif

    public:

      UniformRefinerPattern(percept::PerceptMesh& eMesh, BlockNamesType block_names = BlockNamesType()) :  
        URP<shards::ShellQuadrilateral<4>, shards::ShellQuadrilateral<9>  >(eMesh)
      {

        if (m_eMesh.getSpatialDim() != 3)
          {
            throw std::runtime_error("can't refine shell elements in 2D");
          }
        m_primaryEntityRank = eMesh.element_rank();

        setNeededParts(eMesh, block_names, false);  // different topologies
        Elem::StdMeshObjTopologies::bootstrap();

#if EDGE_QU4_QU9_1_BREAKER

        //m_edge_breaker = Teuchos::rcp( new UniformRefinerPattern<shards::Line<2>, shards::Line<2>, 2, SierraPort > (eMesh, block_names) );
        if (m_eMesh.getSpatialDim() == 3)
          {
            m_edge_breaker = new UniformRefinerPattern<shards::ShellLine<2>, shards::ShellLine<3>, 1, SierraPort > (eMesh, block_names) ;
            m_face_breaker = new UniformRefinerPattern<shards::Quadrilateral<4>, shards::Quadrilateral<9>, 1, SierraPort > (eMesh, block_names) ;
          }
#endif

      }

      void setSubPatterns( std::vector<UniformRefinerPatternBase *>& bp, percept::PerceptMesh& eMesh )
      {
        EXCEPTWATCH;
#if EDGE_QU4_QU9_1_BREAKER
        bp = std::vector<UniformRefinerPatternBase *>(3u, 0);
#else
        bp = std::vector<UniformRefinerPatternBase *>(1u, 0);
#endif

        if (eMesh.getSpatialDim() == 3)
          {
            bp[0] = this;
#if EDGE_QU4_QU9_1_BREAKER
            bp[1] = m_face_breaker;
            bp[2] = m_edge_breaker;
#endif
          }
        else if (eMesh.getSpatialDim() != 3)
          {
            // FIXME
             std::cout << "ERROR" ;
             throw std::runtime_error("ERROR in shell quad class");
          }

      }

      virtual void doBreak() {}
      void fillNeededEntities(std::vector<NeededEntityType>& needed_entities)
      {
        if (m_eMesh.getSpatialDim() == 2)
          {
            throw std::runtime_error("ERROR in shell quad class fillNeededEntities");
          }

        needed_entities.resize(2);
        needed_entities[0].first  = m_eMesh.edge_rank();    
        needed_entities[0].second = 1u;
        // FIXME - this is a problem with the definition of what a needed entity rank is, and what an stk::mesh::Entity rank is
        needed_entities[1].first  = m_eMesh.face_rank(); // FIXME m_eMesh.element_rank();
        //needed_entities[1].first  = m_eMesh.element_rank();
        needed_entities[1].second = 1u;

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
