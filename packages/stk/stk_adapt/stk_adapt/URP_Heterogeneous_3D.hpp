#ifndef stk_adapt_URP_Heterogeneous_3D_hpp
#define stk_adapt_URP_Heterogeneous_3D_hpp


#include <stk_adapt/sierra_element/RefinementTopology.hpp>
#include <stk_adapt/sierra_element/StdMeshObjTopologies.hpp>

#include <Shards_CellTopologyData.h>

#define FACE_BREAKER_HETERO_3D 1
#if FACE_BREAKER_HETERO_3D
#include "UniformRefinerPattern_Quad4_Quad4_4_sierra.hpp"
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

        // list all types of known break patterns to be used here
        m_bp.resize(0);
        m_bp.push_back(  new UniformRefinerPattern<shards::Hexahedron<8>,    shards::Hexahedron<8>,    8, SierraPort > (eMesh, block_names) );
        m_bp.push_back(  new UniformRefinerPattern<shards::Wedge<6>,         shards::Wedge<6>,         8, SierraPort > (eMesh, block_names) );
        m_bp.push_back(  new UniformRefinerPattern<shards::Tetrahedron<4>,   shards::Tetrahedron<4>,   8, SierraPort > (eMesh, block_names) );

        m_bp.push_back(  new UniformRefinerPattern<shards::ShellQuadrilateral<4>,   shards::ShellQuadrilateral<4>,   4, SierraPort > (eMesh, block_names) );

#if FACE_BREAKER_HETERO_3D
        
        m_bp.push_back(  new UniformRefinerPattern<shards::Quadrilateral<4>, shards::Quadrilateral<4>, 4, SierraPort > (eMesh, block_names) );

#endif

      }

      void setSubPatterns( std::vector<UniformRefinerPatternBase *>& bp, percept::PerceptMesh& eMesh )
      {
        EXCEPTWATCH;

        //bp.resize( m_bp.size() );
        bp = m_bp;

#if 0
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

      virtual const CellTopologyData * const getFromTopology()
      {
        throw std::runtime_error("shouldn't call URP_Heterogeneous_3D::getFromTopology()");
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
