#ifndef stk_adapt_UniformRefinerPattern_Quad9_Quad9_4_sierra_hpp
#define stk_adapt_UniformRefinerPattern_Quad9_Quad9_4_sierra_hpp


//#include "UniformRefinerPattern.hpp"
#include <stk_adapt/sierra_element/RefinementTopology.hpp>
#include <stk_adapt/sierra_element/StdMeshObjTopologies.hpp>

#define EDGE_BREAKER_Q9_Q9 1
#if EDGE_BREAKER_Q9_Q9
#include "UniformRefinerPattern_Line3_Line3_2_sierra.hpp"
#endif

namespace stk {
  namespace adapt {

    template <>
    class UniformRefinerPattern<shards::Quadrilateral<9>, shards::Quadrilateral<9>, 4, SierraPort > : public URP<shards::Quadrilateral<9>,shards::Quadrilateral<9>  >
    {

      //Teuchos::RCP< UniformRefinerPattern<shards::Line<2>, shards::Line<2>, 2, SierraPort >  > m_edge_breaker;
      //UniformRefinerPatternBase * m_edge_breaker;
#if EDGE_BREAKER_Q9_Q9
      UniformRefinerPattern<shards::Line<3>, shards::Line<3>, 2, SierraPort > * m_edge_breaker;
#endif

    public:

      UniformRefinerPattern(percept::PerceptMesh& eMesh, BlockNamesType block_names = BlockNamesType()) :  URP<shards::Quadrilateral<9>, shards::Quadrilateral<9>  >(eMesh)
      {
        m_primaryEntityRank = m_eMesh.face_rank();
        if (m_eMesh.getSpatialDim() == 2)
          m_primaryEntityRank = eMesh.element_rank();

        setNeededParts(eMesh, block_names, true);
        Elem::StdMeshObjTopologies::bootstrap();

#if EDGE_BREAKER_Q9_Q9
        if (m_eMesh.getSpatialDim() == 2)
          m_edge_breaker =  new UniformRefinerPattern<shards::Line<3>, shards::Line<3>, 2, SierraPort > (eMesh, block_names) ;
        else
          m_edge_breaker = 0;
#endif

      }
      ~UniformRefinerPattern()
      {
        if (m_edge_breaker) delete m_edge_breaker;
      }

      void setSubPatterns( std::vector<UniformRefinerPatternBase *>& bp, percept::PerceptMesh& eMesh )
      {
        EXCEPTWATCH;
#if EDGE_BREAKER_Q9_Q9
        bp = std::vector<UniformRefinerPatternBase *>(2u, 0);
#else
        bp = std::vector<UniformRefinerPatternBase *>(1u, 0);
#endif

        if (eMesh.getSpatialDim() == 2)
          {
            bp[0] = this;
#if EDGE_BREAKER_Q9_Q9
            bp[1] = m_edge_breaker;
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
        needed_entities.resize(2);
        needed_entities[0] = NeededEntityType(m_eMesh.edge_rank(), 3u);
        needed_entities[1] = NeededEntityType( (m_eMesh.getSpatialDim() == 2 ? m_eMesh.element_rank() : m_eMesh.face_rank()), 9u);
        //setToOne(needed_entities);
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
