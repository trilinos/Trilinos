#ifndef stk_adapt_UniformRefinerPattern_ShellTri3_ShellTri3_4_sierra_hpp
#define stk_adapt_UniformRefinerPattern_ShellTri3_ShellTri3_4_sierra_hpp


//#include "UniformRefinerPattern.hpp"
#include <stk_adapt/sierra_element/RefinementTopology.hpp>
#include <stk_adapt/sierra_element/StdMeshObjTopologies.hpp>

#define EDGE_ST3_ST3_4_BREAKER 1
#if EDGE_ST3_ST3_4_BREAKER
#include "UniformRefinerPattern_ShellLine2_ShellLine2_2_sierra.hpp"
#include "UniformRefinerPattern_Tri3_Tri3_4_sierra.hpp"
#endif

namespace stk {
  namespace adapt {

    template <>
    class UniformRefinerPattern<shards::ShellTriangle<3>, shards::ShellTriangle<3>, 4, SierraPort > : public URP<shards::ShellTriangle<3>,shards::ShellTriangle<3>  >
    {


#if EDGE_ST3_ST3_4_BREAKER
      UniformRefinerPattern<shards::ShellLine<2>, shards::ShellLine<2>, 2, SierraPort > * m_edge_breaker;
      UniformRefinerPattern<shards::Triangle<3>, shards::Triangle<3>, 4, SierraPort > * m_face_breaker;
#endif

    public:

      UniformRefinerPattern(percept::PerceptMesh& eMesh, BlockNamesType block_names = BlockNamesType()) :  URP<shards::ShellTriangle<3>, shards::ShellTriangle<3>  >(eMesh)
      {

        if (m_eMesh.getSpatialDim() != 3)
          {
            throw std::runtime_error("can't refine shell elements in 2D");
          }
        m_primaryEntityRank = eMesh.element_rank();

        setNeededParts(eMesh, block_names, true);
        Elem::StdMeshObjTopologies::bootstrap();

#if EDGE_ST3_ST3_4_BREAKER

        //m_edge_breaker = Teuchos::rcp( new UniformRefinerPattern<shards::Line<2>, shards::Line<2>, 2, SierraPort > (eMesh, block_names) );
        if (m_eMesh.getSpatialDim() == 3)
          {
            m_edge_breaker = new UniformRefinerPattern<shards::ShellLine<2>, shards::ShellLine<2>, 2, SierraPort > (eMesh, block_names) ;
            m_face_breaker = new UniformRefinerPattern<shards::Triangle<3>, shards::Triangle<3>, 4, SierraPort > (eMesh, block_names) ;
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
#if EDGE_ST3_ST3_4_BREAKER
        bp = std::vector<UniformRefinerPatternBase *>(3u, 0);
#else
        bp = std::vector<UniformRefinerPatternBase *>(1u, 0);
#endif

        if (eMesh.getSpatialDim() == 3)
          {
            bp[0] = this;
#if EDGE_ST3_ST3_4_BREAKER
            bp[1] = m_face_breaker;
            bp[2] = m_edge_breaker;
#endif
          }
        else if (eMesh.getSpatialDim() != 3)
          {
            // FIXME
             std::cout << "ERROR" ;
             throw std::runtime_error("ERROR in shell tri class");
          }

      }

      virtual void doBreak() {}
      void fillNeededEntities(std::vector<NeededEntityType>& needed_entities)
      {
        needed_entities.resize(1);
        needed_entities[0].first = m_eMesh.edge_rank();    
        setToOne(needed_entities);
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
