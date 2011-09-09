#ifndef stk_adapt_UniformRefinerPattern_Tet4_Tet4_8_sierra_hpp
#define stk_adapt_UniformRefinerPattern_Tet4_Tet4_8_sierra_hpp


//#include "UniformRefinerPattern.hpp"
#include <stk_adapt/sierra_element/RefinementTopology.hpp>
#include <stk_adapt/sierra_element/StdMeshObjTopologies.hpp>

#define USE_FACE_BREAKER_T4_T4_8 1
#if USE_FACE_BREAKER_T4_T4_8
#include "UniformRefinerPattern_Tri3_Tri3_4_sierra.hpp"
#endif

namespace stk {
  namespace adapt {

    template <>
    class UniformRefinerPattern<shards::Tetrahedron<4>, shards::Tetrahedron<4>, 8, SierraPort > : public URP<shards::Tetrahedron<4>, shards::Tetrahedron<4>  >
    {

#if USE_FACE_BREAKER_T4_T4_8
      // FIXME
      UniformRefinerPattern<shards::Triangle<3>, shards::Triangle<3>, 4, SierraPort > * m_face_breaker;
#endif

    public:

      UniformRefinerPattern(percept::PerceptMesh& eMesh, BlockNamesType block_names = BlockNamesType()) :  URP<shards::Tetrahedron<4>, shards::Tetrahedron<4>  >(eMesh)
      {
        m_primaryEntityRank = eMesh.element_rank();

        setNeededParts(eMesh, block_names, true);
        Elem::StdMeshObjTopologies::bootstrap();

#if USE_FACE_BREAKER_T4_T4_8

        //m_face_breaker = Teuchos::rcp( new UniformRefinerPattern<shards::Line<2>, shards::Line<2>, 2, SierraPort > (eMesh, block_names) );
        m_face_breaker =  new UniformRefinerPattern<shards::Triangle<3>, shards::Triangle<3>, 4, SierraPort > (eMesh, block_names) ;
#endif

      }

      ~UniformRefinerPattern()
      {
        if (m_face_breaker) delete m_face_breaker;
      }

      void setSubPatterns( std::vector<UniformRefinerPatternBase *>& bp, percept::PerceptMesh& eMesh )
      {
        EXCEPTWATCH;
        bp = std::vector<UniformRefinerPatternBase *>(2u, 0);

        if (eMesh.getSpatialDim() == 3)
          {
            bp[0] = this;
#if USE_FACE_BREAKER_T4_T4_8
            bp[1] = m_face_breaker;
#endif
          }
        else if (eMesh.getSpatialDim() == 2)
          {
            // FIXME
            std::cout << "ERROR" ;
            exit(1);
          }

      }

      virtual void doBreak() {}
      void fillNeededEntities(std::vector<NeededEntityType>& needed_entities)
      {
        needed_entities.resize(1);
        needed_entities[0].first = m_eMesh.edge_rank();   
        setToOne(needed_entities);
      }

      virtual unsigned getNumNewElemPerElem() { return 8; }

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
