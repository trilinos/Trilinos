#ifndef stk_adapt_UniformRefinerPattern_Tet10_Tet10_8_sierra_hpp
#define stk_adapt_UniformRefinerPattern_Tet10_Tet10_8_sierra_hpp


//#include "UniformRefinerPattern.hpp"
#include <stk_adapt/sierra_element/RefinementTopology.hpp>
#include <stk_adapt/sierra_element/StdMeshObjTopologies.hpp>

#define FACE_BREAKER_T10_T10 1
#if FACE_BREAKER_T10_T10
#include "UniformRefinerPattern_Tri6_Tri6_4_sierra.hpp"
#endif

namespace stk {
  namespace adapt {

    template <>
    class UniformRefinerPattern<shards::Tetrahedron<10>, shards::Tetrahedron<10>, 8, SierraPort > : public URP<shards::Tetrahedron<10>, shards::Tetrahedron<10>  >
    {

      //Teuchos::RCP< UniformRefinerPattern<shards::Line<2>, shards::Line<2>, 2, SierraPort >  > m_face_breaker;
      //UniformRefinerPatternBase * m_face_breaker;
#if FACE_BREAKER_T10_T10
      UniformRefinerPattern<shards::Triangle<6>, shards::Triangle<6>, 4, SierraPort > * m_face_breaker;
#endif

    public:

      UniformRefinerPattern(percept::PerceptMesh& eMesh, BlockNamesType block_names = BlockNamesType()) :  URP<shards::Tetrahedron<10>, shards::Tetrahedron<10>  >(eMesh)
      {
        m_primaryEntityRank = eMesh.element_rank();

        setNeededParts(eMesh, block_names, true);
        Elem::StdMeshObjTopologies::bootstrap();

#if FACE_BREAKER_T10_T10
        
        m_face_breaker =  new UniformRefinerPattern<shards::Triangle<6>, shards::Triangle<6>, 4, SierraPort > (eMesh, block_names);

#endif

      }
      ~UniformRefinerPattern()
      {
        if (m_face_breaker) delete m_face_breaker;
      }

      void setSubPatterns( std::vector<UniformRefinerPatternBase *>& bp, percept::PerceptMesh& eMesh )
      {
        EXCEPTWATCH;
#if FACE_BREAKER_T10_T10
        bp = std::vector<UniformRefinerPatternBase *>(2u, 0);
#else
        bp = std::vector<UniformRefinerPatternBase *>(1u, 0);
#endif

        bp[0] = this;
#if FACE_BREAKER_T10_T10
        bp[1] = m_face_breaker;
#endif
      }

      virtual void doBreak() {}
      void fillNeededEntities(std::vector<NeededEntityType>& needed_entities)
      {
        needed_entities.resize(3);
        // 4 vertices
        needed_entities[0] = NeededEntityType(m_eMesh.edge_rank(), 3u); // 18
        needed_entities[1] = NeededEntityType(m_eMesh.face_rank(), 3u); // 12
        needed_entities[2] = NeededEntityType(m_eMesh.element_rank(), 1u); // 1
        //setToOne(needed_entities);
      }

      virtual unsigned getNumNewElemPerElem() { return 8; }

      //typedef boost::array<unsigned, ToTopology::node_count > refined_element_type;


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
