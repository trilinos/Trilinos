#ifndef stk_adapt_UniformRefinerPattern_Hex20_Hex20_8_sierra_hpp
#define stk_adapt_UniformRefinerPattern_Hex20_Hex20_8_sierra_hpp


//#include "UniformRefinerPattern.hpp"
#include <stk_adapt/sierra_element/RefinementTopology.hpp>
#include <stk_adapt/sierra_element/StdMeshObjTopologies.hpp>

#define FACE_BREAKER_H20_H20 1
#if FACE_BREAKER_H20_H20
#include "UniformRefinerPattern_Quad8_Quad8_4_sierra.hpp"
#endif

namespace stk {
  namespace adapt {

    template <>
    class UniformRefinerPattern<shards::Hexahedron<20>, shards::Hexahedron<20>, 8, SierraPort > : public URP<shards::Hexahedron<20>, shards::Hexahedron<20>  >
    {

      //Teuchos::RCP< UniformRefinerPattern<shards::Line<2>, shards::Line<2>, 2, SierraPort >  > m_face_breaker;
      //UniformRefinerPatternBase * m_face_breaker;
#if FACE_BREAKER_H20_H20
      UniformRefinerPattern<shards::Quadrilateral<8>, shards::Quadrilateral<8>, 4, SierraPort > * m_face_breaker;
#endif

    public:

      UniformRefinerPattern(percept::PerceptMesh& eMesh, BlockNamesType block_names = BlockNamesType()) :  URP<shards::Hexahedron<20>, shards::Hexahedron<20>  >(eMesh)
      {
        m_primaryEntityRank = eMesh.element_rank();

        setNeededParts(eMesh, block_names, true);
        Elem::StdMeshObjTopologies::bootstrap();

#if FACE_BREAKER_H20_H20
        
        m_face_breaker =  new UniformRefinerPattern<shards::Quadrilateral<8>, shards::Quadrilateral<8>, 4, SierraPort > (eMesh, block_names);

#endif

      }
      ~UniformRefinerPattern()
      {
        if (m_face_breaker) delete m_face_breaker;
      }

      void setSubPatterns( std::vector<UniformRefinerPatternBase *>& bp, percept::PerceptMesh& eMesh )
      {
        EXCEPTWATCH;
#if FACE_BREAKER_H20_H20
        bp = std::vector<UniformRefinerPatternBase *>(2u, 0);
#else
        bp = std::vector<UniformRefinerPatternBase *>(1u, 0);
#endif

        bp[0] = this;
#if FACE_BREAKER_H20_H20
        bp[1] = m_face_breaker;
#endif
      }

      virtual void doBreak() {}
      void fillNeededEntities(std::vector<NeededEntityType>& needed_entities)
      {
        needed_entities.resize(3);
        needed_entities[0] = NeededEntityType(m_eMesh.edge_rank(), 3u);
//         needed_entities[1] = NeededEntityType(m_eMesh.face_rank(), 5u);
//         needed_entities[2] = NeededEntityType(m_eMesh.element_rank(), 7u);
        needed_entities[1] = NeededEntityType(m_eMesh.face_rank(), 9u);
        needed_entities[2] = NeededEntityType(m_eMesh.element_rank(), 27u);
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
