#ifndef stk_adapt_UniformRefinerPattern_Wedge6_Wedge18_1_sierra_hpp
#define stk_adapt_UniformRefinerPattern_Wedge6_Wedge18_1_sierra_hpp


#include <stk_adapt/sierra_element/RefinementTopology.hpp>
#include <stk_adapt/sierra_element/StdMeshObjTopologies.hpp>

#include <boost/array.hpp>

#define FACE_BREAKER_W6_W18_1 1
#if FACE_BREAKER_W6_W18_1
#include "UniformRefinerPattern_Tri3_Tri6_1_sierra.hpp"
#include "UniformRefinerPattern_Quad4_Quad9_1_sierra.hpp"
#endif

namespace stk {
  namespace adapt {

    template <>
    class UniformRefinerPattern< shards::Wedge<6>, shards::Wedge<18>, 1, SierraPort > : public URP<shards::Wedge<6> , shards::Wedge<18> >
    {

#if FACE_BREAKER_W6_W18_1
      UniformRefinerPattern<shards::Triangle<3>, shards::Triangle<6>, 1, SierraPort > * m_subDim_breaker;
      UniformRefinerPattern<shards::Quadrilateral<4>, shards::Quadrilateral<9>, 1, SierraPort > * m_subDim_breaker_quad;
#endif

    public:

      UniformRefinerPattern(percept::PerceptMesh& eMesh, BlockNamesType block_names = BlockNamesType()) : URP<shards::Wedge<6> , shards::Wedge<18> >(eMesh)
      {
        m_primaryEntityRank = eMesh.element_rank();

        setNeededParts(eMesh, block_names, false);
        Elem::StdMeshObjTopologies::bootstrap();
#if FACE_BREAKER_W6_W18_1

        m_subDim_breaker =  new UniformRefinerPattern<shards::Triangle<3>, shards::Triangle<6>, 1, SierraPort > (eMesh, block_names) ;
        m_subDim_breaker_quad = new UniformRefinerPattern<shards::Quadrilateral<4>, shards::Quadrilateral<9>, 1, SierraPort > (eMesh, block_names);
#endif

      }
      ~UniformRefinerPattern()
      {
        if (m_subDim_breaker) delete m_subDim_breaker;
        if (m_subDim_breaker_quad) delete m_subDim_breaker_quad;
      }


      virtual void doBreak() {}

      void setSubPatterns( std::vector<UniformRefinerPatternBase *>& bp, percept::PerceptMesh& eMesh )
      {
        EXCEPTWATCH;
        bp = std::vector<UniformRefinerPatternBase *>(3u, 0);

        bp[0] = this;
#if FACE_BREAKER_W6_W18_1
        bp[1] = m_subDim_breaker;
        bp[2] = m_subDim_breaker_quad;
#endif

      }

      void fillNeededEntities(std::vector<NeededEntityType>& needed_entities)
      {
        needed_entities.resize(2);
        needed_entities[0] = NeededEntityType(m_eMesh.edge_rank(), 1u);
        needed_entities[1] = NeededEntityType(m_eMesh.face_rank(), 1u); 
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
