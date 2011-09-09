#ifndef stk_adapt_UniformRefinerPattern_ShellQuad8_ShellQuad8_4_sierra_hpp
#define stk_adapt_UniformRefinerPattern_ShellQuad8_ShellQuad8_4_sierra_hpp


#include <stk_adapt/sierra_element/RefinementTopology.hpp>
#include <stk_adapt/sierra_element/StdMeshObjTopologies.hpp>

#define EDGE_QU8_QU8_4_BREAKER 1
#if EDGE_QU8_QU8_4_BREAKER

#include "UniformRefinerPattern_ShellLine3_ShellLine3_2_sierra.hpp"
#include "UniformRefinerPattern_Quad8_Quad8_4_sierra.hpp"
#endif

namespace stk {
  namespace adapt {

    template <>
    class UniformRefinerPattern<shards::ShellQuadrilateral<8>, shards::ShellQuadrilateral<8>, 4, SierraPort > : 
      public URP<shards::ShellQuadrilateral<8>,shards::ShellQuadrilateral<8>  >
    {


#if EDGE_QU8_QU8_4_BREAKER
      UniformRefinerPattern<shards::ShellLine<3>, shards::ShellLine<3>, 2, SierraPort > * m_edge_breaker;
      UniformRefinerPattern<shards::Quadrilateral<8>, shards::Quadrilateral<8>, 4, SierraPort > * m_face_breaker;
#endif

    public:

      UniformRefinerPattern(percept::PerceptMesh& eMesh, BlockNamesType block_names = BlockNamesType()) :  
        URP<shards::ShellQuadrilateral<8>, shards::ShellQuadrilateral<8>  >(eMesh)
      {

        if (m_eMesh.getSpatialDim() != 3)
          {
            throw std::runtime_error("can't refine shell elements in 2D");
          }
        m_primaryEntityRank = eMesh.element_rank();

        setNeededParts(eMesh, block_names, true);
        Elem::StdMeshObjTopologies::bootstrap();

#if EDGE_QU8_QU8_4_BREAKER

        //m_edge_breaker = Teuchos::rcp( new UniformRefinerPattern<shards::Line<2>, shards::Line<2>, 2, SierraPort > (eMesh, block_names) );
        if (m_eMesh.getSpatialDim() == 3)
          {
            m_edge_breaker = new UniformRefinerPattern<shards::ShellLine<3>, shards::ShellLine<3>, 2, SierraPort > (eMesh, block_names) ;
            m_face_breaker = new UniformRefinerPattern<shards::Quadrilateral<8>, shards::Quadrilateral<8>, 4, SierraPort > (eMesh, block_names) ;
          }
#endif

      }

      ~UniformRefinerPattern()
      {
        if (m_edge_breaker) delete m_edge_breaker;
        if (m_face_breaker) delete m_face_breaker;
      }

      void setSubPatterns( std::vector<UniformRefinerPatternBase *>& bp, percept::PerceptMesh& eMesh )
      {
        EXCEPTWATCH;
        bp.resize(0);

        if (eMesh.getSpatialDim() == 3)
          {
            bp.push_back(this);
#if EDGE_QU8_QU8_4_BREAKER
            bp.push_back( m_face_breaker);
            bp.push_back( m_edge_breaker);  
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
        needed_entities[0] = NeededEntityType(m_eMesh.edge_rank(), 3u);    
        needed_entities[1] = NeededEntityType(m_eMesh.face_rank() , 9u);
      }

      virtual unsigned getNumNewElemPerElem() { return 4; }

      void 
      createNewElements(percept::PerceptMesh& eMesh, NodeRegistry& nodeRegistry, 
                        stk::mesh::Entity& element,  NewSubEntityNodesType& new_sub_entity_nodes, vector<stk::mesh::Entity *>::iterator& element_pool,
                        stk::mesh::FieldBase *proc_rank_field=0)
      {
        genericRefine_createNewElements(eMesh, nodeRegistry,
                                        element, new_sub_entity_nodes, element_pool,
                                        proc_rank_field);
      }
      
    };

  }
}
#endif
