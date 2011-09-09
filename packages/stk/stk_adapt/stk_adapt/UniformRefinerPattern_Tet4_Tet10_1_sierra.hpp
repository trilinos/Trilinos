#ifndef stk_adapt_UniformRefinerPattern_Tet4_Tet10_1_sierra_hpp
#define stk_adapt_UniformRefinerPattern_Tet4_Tet10_1_sierra_hpp


#include <stk_adapt/sierra_element/RefinementTopology.hpp>
#include <stk_adapt/sierra_element/StdMeshObjTopologies.hpp>

#include <boost/array.hpp>

#define FACE_BREAKER_TET4_TET10_1 1
#if FACE_BREAKER_TET4_TET10_1
#include "UniformRefinerPattern_Tri3_Tri6_1_sierra.hpp"
#endif

namespace stk {
  namespace adapt {

    template <>
    class UniformRefinerPattern< shards::Tetrahedron<4>, shards::Tetrahedron<10>, 1, SierraPort > : public URP<shards::Tetrahedron<4> , shards::Tetrahedron<10> >
    {

#if FACE_BREAKER_TET4_TET10_1
      UniformRefinerPattern<shards::Triangle<3>, shards::Triangle<6>, 1, SierraPort > * m_subDim_breaker;
#endif

    public:

      UniformRefinerPattern(percept::PerceptMesh& eMesh, BlockNamesType block_names = BlockNamesType()) : URP<shards::Tetrahedron<4> , shards::Tetrahedron<10> >(eMesh)
      {
        m_primaryEntityRank = eMesh.element_rank();

        setNeededParts(eMesh, block_names, false);
        Elem::StdMeshObjTopologies::bootstrap();
#if FACE_BREAKER_TET4_TET10_1

        m_subDim_breaker =  new UniformRefinerPattern<shards::Triangle<3>, shards::Triangle<6>, 1, SierraPort > (eMesh, block_names) ;
#endif

      }

      ~UniformRefinerPattern()
      {
        if (m_subDim_breaker) delete m_subDim_breaker;
      }


      virtual void doBreak() {}

      void setSubPatterns( std::vector<UniformRefinerPatternBase *>& bp, percept::PerceptMesh& eMesh )
      {
        EXCEPTWATCH;
        bp = std::vector<UniformRefinerPatternBase *>(2u, 0);

        bp[0] = this;
#if FACE_BREAKER_TET4_TET10_1
        bp[1] = m_subDim_breaker;
#endif

      }

      void fillNeededEntities(std::vector<NeededEntityType>& needed_entities)
      {
        needed_entities.resize(1);
        needed_entities[0].first = m_eMesh.edge_rank();   
        setToOne(needed_entities);
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

  } // namespace adapt
} // namespace stk
#endif
