#ifndef stk_adapt_UniformRefinerPattern_Tri3_Tri3_4_sierra_hpp
#define stk_adapt_UniformRefinerPattern_Tri3_Tri3_4_sierra_hpp


//#include "UniformRefinerPattern.hpp"
#include <stk_adapt/sierra_element/RefinementTopology.hpp>
#include <stk_adapt/sierra_element/StdMeshObjTopologies.hpp>

#include "UniformRefinerPattern_Line2_Line2_2_sierra.hpp"

namespace stk {
  namespace adapt {

    template <>
    class UniformRefinerPattern<shards::Triangle<3>, shards::Triangle<3>, 4, SierraPort > : public URP<shards::Triangle<3>,shards::Triangle<3>  >
    {

#define EDGE_BREAKER_T3_T3_4_S 1
#if EDGE_BREAKER_T3_T3_4_S
      UniformRefinerPattern<shards::Line<2>, shards::Line<2>, 2, SierraPort > * m_edge_breaker;
#endif

    public:

      UniformRefinerPattern(percept::PerceptMesh& eMesh, BlockNamesType block_names = BlockNamesType()) :  URP<shards::Triangle<3>, shards::Triangle<3>  >(eMesh),
                                                                                                           m_edge_breaker (0)
      {
        m_primaryEntityRank = m_eMesh.face_rank();
        if (m_eMesh.getSpatialDim() == 2)
          m_primaryEntityRank = eMesh.element_rank();

        setNeededParts(eMesh, block_names, true);
        Elem::StdMeshObjTopologies::bootstrap();

#if EDGE_BREAKER_T3_T3_4_S
        //m_edge_breaker = Teuchos::rcp( new UniformRefinerPattern<shards::Line<2>, shards::Line<2>, 2, SierraPort > (eMesh, block_names) );
        if (m_eMesh.getSpatialDim() == 2)
          {
            m_edge_breaker =  new UniformRefinerPattern<shards::Line<2>, shards::Line<2>, 2, SierraPort > (eMesh, block_names) ;
          }
#endif

      }

      ~UniformRefinerPattern()
      {
        if (m_edge_breaker) delete m_edge_breaker;
      }

      void setSubPatterns( std::vector<UniformRefinerPatternBase *>& bp, percept::PerceptMesh& eMesh )
      {
        EXCEPTWATCH;
        bp = std::vector<UniformRefinerPatternBase *>(2u, 0);

        if (eMesh.getSpatialDim() == 2)
          {
            bp[0] = this;
#if EDGE_BREAKER_T3_T3_4_S
            if (m_eMesh.getSpatialDim() == 2)
              {
                bp[1] = m_edge_breaker;
              }
#endif
          }
        else if (eMesh.getSpatialDim() == 3)
          {
            // FIXME
//             std::cout << "ERROR" ;
//             exit(1);
          }

      }

      virtual void doBreak() {}
      void fillNeededEntities(std::vector<NeededEntityType>& needed_entities)
      {
        needed_entities.resize(1);
        needed_entities[0].first = m_eMesh.edge_rank();    
        //needed_entities[1] = (m_eMesh.getSpatialDim() == 2 ? m_eMesh.element_rank() : m_eMesh.face_rank());
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
