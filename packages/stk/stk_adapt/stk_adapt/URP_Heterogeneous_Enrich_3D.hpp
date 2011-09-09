#ifndef stk_adapt_URP_Heterogeneous_Enrich_3D_hpp
#define stk_adapt_URP_Heterogeneous_Enrich_3D_hpp


#include <stk_adapt/sierra_element/RefinementTopology.hpp>
#include <stk_adapt/sierra_element/StdMeshObjTopologies.hpp>

#include <Shards_CellTopologyData.h>


#define FACE_BREAKER_HETERO_ENRICH_3D 1

#if FACE_BREAKER_HETERO_ENRICH_3D
#include "UniformRefinerPattern_Quad4_Quad4_4_sierra.hpp"
#include "UniformRefinerPattern_Tri3_Tri3_4_sierra.hpp"
#endif

namespace stk {
  namespace adapt {

    class URP_Heterogeneous_Enrich_3D : public UniformRefinerPatternBase
    {

      std::vector<UniformRefinerPatternBase *> m_bp;

#if FACE_BREAKER_HETERO_ENRICH_3D
      UniformRefinerPattern<shards::Quadrilateral<4>, shards::Quadrilateral<4>, 4, SierraPort > * m_face_breaker;
#endif

    protected:

      percept::PerceptMesh& m_eMesh;

    public:

      URP_Heterogeneous_Enrich_3D(percept::PerceptMesh& eMesh, BlockNamesType block_names = BlockNamesType()) :  m_eMesh(eMesh) 
      {
        m_primaryEntityRank = eMesh.element_rank();

        Elem::StdMeshObjTopologies::bootstrap();

        // list all types of known break patterns to be used here
        m_bp.resize(0);

        int spatialDim = eMesh.getSpatialDim();

        // refine

        // put them in reverse topological rank order

        if (spatialDim != 3)
          {
            throw std::runtime_error("URP_Heterogeneous_Enrich_3D is only for 3D meshes");
          }

//         // refine
//         m_bp.push_back(  new UniformRefinerPattern<shards::Hexahedron<8>,    shards::Hexahedron<8>,    8, SierraPort > (eMesh, block_names) );
//         m_bp.push_back(  new UniformRefinerPattern<shards::Wedge<6>,         shards::Wedge<6>,         8, SierraPort > (eMesh, block_names) );
//         m_bp.push_back(  new UniformRefinerPattern<shards::Tetrahedron<4>,   shards::Tetrahedron<4>,   8, SierraPort > (eMesh, block_names) );


        // enrich
        m_bp.push_back ( new UniformRefinerPattern< shards::Wedge<6>, shards::Wedge<15>, 1, SierraPort >                 (eMesh, block_names) );
        m_bp.push_back ( new UniformRefinerPattern<shards::Tetrahedron<4>,   shards::Tetrahedron<10>,  1, SierraPort >   (eMesh, block_names) );
        m_bp.push_back ( new UniformRefinerPattern<shards::Hexahedron<8>,    shards::Hexahedron<27>,   1, SierraPort >   (eMesh, block_names) );

        //m_bp.push_back(  new UniformRefinerPattern<shards::ShellQuadrilateral<4>,  shards::ShellQuadrilateral<9>,   1, SierraPort > (eMesh, block_names) );
        //m_bp.push_back(  new UniformRefinerPattern<shards::ShellTriangle<3>,       shards::ShellTriangle<3>,   4, SierraPort > (eMesh, block_names) );


#if FACE_BREAKER_HETERO_ENRICH_3D
        
        m_bp.push_back(  new UniformRefinerPattern<shards::Quadrilateral<4>, shards::Quadrilateral<4>, 4, SierraPort > (eMesh, block_names) );
        m_bp.push_back(  new UniformRefinerPattern<shards::Triangle<3>, shards::Triangle<3>, 4, SierraPort > (eMesh, block_names) );

#endif

      }

      ~URP_Heterogeneous_Enrich_3D()
      {
        for (unsigned ibp=0; ibp < m_bp.size(); ibp++)
          {
            if (m_bp[ibp]) delete m_bp[ibp];
          }
      }

      void setSubPatterns( std::vector<UniformRefinerPatternBase *>& bp, percept::PerceptMesh& eMesh )
      {
        EXCEPTWATCH;

        bp = m_bp;

      }

      virtual void doBreak() 
      {
        throw std::runtime_error("shouldn't call URP_Heterogeneous_Enrich_3D::doBreak()");

      }
      virtual unsigned getFromTypeKey()
      {
        throw std::runtime_error("shouldn't call URP_Heterogeneous_Enrich_3D::getFromTypeKey()");

      }

      virtual const CellTopologyData * const getFromTopology()
      {
        throw std::runtime_error("shouldn't call URP_Heterogeneous_Enrich_3D::getFromTopology()");
      }
      virtual const CellTopologyData * const getToTopology()
      {
        throw std::runtime_error("shouldn't call URP_Heterogeneous_3D::getToTopology()");
      }


      void fillNeededEntities(std::vector<NeededEntityType>& needed_entities)
      {
        throw std::runtime_error("shouldn't call URP_Heterogeneous_Enrich_3D::fillNeededEntities()");
      }

      virtual unsigned getNumNewElemPerElem() 
      { 
        throw std::runtime_error("shouldn't call URP_Heterogeneous_Enrich_3D::getNumNewElemPerElem()");
        return 8; 
      }

      void 
      createNewElements(percept::PerceptMesh& eMesh, NodeRegistry& nodeRegistry, 
                        stk::mesh::Entity& element,  NewSubEntityNodesType& new_sub_entity_nodes, vector<stk::mesh::Entity *>::iterator& element_pool,
                        stk::mesh::FieldBase *proc_rank_field=0)
      {
        throw std::runtime_error("shouldn't call URP_Heterogeneous_Enrich_3D::createNewElements()");
      }      
    };

  }
}
#endif
