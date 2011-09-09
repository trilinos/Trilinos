#ifndef stk_adapt_UniformRefinerPattern_Hex8_Hex8_8_sierra_hpp
#define stk_adapt_UniformRefinerPattern_Hex8_Hex8_8_sierra_hpp

#include <stk_adapt/sierra_element/RefinementTopology.hpp>
#include <stk_adapt/sierra_element/StdMeshObjTopologies.hpp>

#include "UniformRefinerPattern_Quad4_Quad4_4_sierra.hpp"

#include <boost/array.hpp>

namespace stk {
  namespace adapt {

    template <>
    class UniformRefinerPattern<shards::Hexahedron<8>, shards::Hexahedron<8>, 8, SierraPort > : public URP<shards::Hexahedron<8>,shards::Hexahedron<8>  >
    {

#define FACE_BREAKER_H8_H8_8 1
#if FACE_BREAKER_H8_H8_8
      UniformRefinerPattern<shards::Quadrilateral<4>, shards::Quadrilateral<4>, 4, SierraPort > * m_face_breaker;
#endif

    public:

      UniformRefinerPattern(percept::PerceptMesh& eMesh, BlockNamesType block_names = BlockNamesType()) :  URP<shards::Hexahedron<8>, shards::Hexahedron<8>  >(eMesh)
      {
        m_primaryEntityRank = eMesh.element_rank();

        setNeededParts(eMesh, block_names, true);
        Elem::StdMeshObjTopologies::bootstrap();

#if FACE_BREAKER_H8_H8_8

        m_face_breaker =  new UniformRefinerPattern<shards::Quadrilateral<4>, shards::Quadrilateral<4>, 4, SierraPort > (eMesh, block_names) ;
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

        bp[0] = this;
#if FACE_BREAKER_H8_H8_8
        bp[1] = m_face_breaker;
#endif
      }

      virtual void doBreak() {}
      void fillNeededEntities(std::vector<NeededEntityType>& needed_entities)
      {
        needed_entities.resize(3);
        needed_entities[0].first = m_eMesh.edge_rank();  
        needed_entities[1].first = m_eMesh.face_rank();
        needed_entities[2].first = m_eMesh.element_rank();
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
