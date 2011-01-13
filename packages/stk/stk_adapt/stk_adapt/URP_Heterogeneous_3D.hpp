#ifndef stk_adapt_URP_Heterogeneous_3D_hpp
#define stk_adapt_URP_Heterogeneous_3D_hpp


#include <stk_adapt/sierra_element/RefinementTopology.hpp>
#include <stk_adapt/sierra_element/StdMeshObjTopologies.hpp>

#define FACE_BREAKER_HETERO_3D 0
#if FACE_BREAKER_HETERO_3D
#include "UniformRefinerPattern_Quad9_Quad9_4_sierra.hpp"
#endif

namespace stk {
  namespace adapt {

    class URP_Heterogeneous_3D : public UniformRefinerPatternBase
    {

      std::vector<UniformRefinerPatternBase *> m_bp;

#if FACE_BREAKER_HETERO_3D
      UniformRefinerPattern<shards::Quadrilateral<9>, shards::Quadrilateral<9>, 4, SierraPort > * m_face_breaker;
#endif

    protected:

      percept::PerceptMesh& m_eMesh;

    public:

      URP_Heterogeneous_3D(percept::PerceptMesh& eMesh, BlockNamesType block_names = BlockNamesType()) :  m_eMesh(eMesh) 
      {
        m_primaryEntityRank = mesh::Element;

        //!setNeededParts(eMesh, block_names, true);
        Elem::StdMeshObjTopologies::bootstrap();

        m_bp.resize(3);
        m_bp[0] = new UniformRefinerPattern<shards::Hexahedron<8>,    shards::Hexahedron<8>,    8, SierraPort > (eMesh, block_names);
        m_bp[1] = new UniformRefinerPattern<shards::Wedge<6>,         shards::Wedge<6>,         8, SierraPort > (eMesh, block_names);
        m_bp[2] = new UniformRefinerPattern<shards::Tetrahedron<4>,   shards::Tetrahedron<4>,   8, SierraPort > (eMesh, block_names);

#if FACE_BREAKER_HETERO_3D
        
        m_face_breaker =  new UniformRefinerPattern<shards::Quadrilateral<9>, shards::Quadrilateral<9>, 4, SierraPort > (eMesh, block_names);

#endif

      }

      void setSubPatterns( std::vector<UniformRefinerPatternBase *>& bp, percept::PerceptMesh& eMesh )
      {
        EXCEPTWATCH;
#if FACE_BREAKER_HETERO_3D
        bp = std::vector<UniformRefinerPatternBase *>(2u, 0);
#else
        bp = std::vector<UniformRefinerPatternBase *>(3u, 0);   // FIXME
#endif

        bp[0] = m_bp[0];
        bp[1] = m_bp[1];
        bp[2] = m_bp[2];

#if FACE_BREAKER_HETERO_3D
        bp[1] = m_face_breaker;
#endif
      }

      virtual void doBreak() 
      {
        throw std::runtime_error("shouldn't call URP_Heterogeneous_3D::doBreak()");

      }
      virtual unsigned getFromTypeKey()
      {
        throw std::runtime_error("shouldn't call URP_Heterogeneous_3D::getFromTypeKey()");

      }

      void fillNeededEntities(std::vector<NeededEntityType>& needed_entities)
      {
        throw std::runtime_error("shouldn't call URP_Heterogeneous_3D::fillNeededEntities()");
      }

      virtual unsigned getNumNewElemPerElem() 
      { 
        throw std::runtime_error("shouldn't call URP_Heterogeneous_3D::getNumNewElemPerElem()");
        return 8; 
      }

      void 
      createNewElements(percept::PerceptMesh& eMesh, NodeRegistry& nodeRegistry, 
                        Entity& element,  NewSubEntityNodesType& new_sub_entity_nodes, vector<Entity *>::iterator& element_pool,
                        FieldBase *proc_rank_field=0)
      {
        throw std::runtime_error("shouldn't call URP_Heterogeneous_3D::createNewElements()");
      }      
    };

  }
}
#endif
