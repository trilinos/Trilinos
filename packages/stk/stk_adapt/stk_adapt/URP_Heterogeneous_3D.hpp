#ifndef stk_adapt_URP_Heterogeneous_3D_hpp
#define stk_adapt_URP_Heterogeneous_3D_hpp


#include <stk_adapt/sierra_element/RefinementTopology.hpp>
#include <stk_adapt/sierra_element/StdMeshObjTopologies.hpp>

#include <Shards_CellTopologyData.h>

namespace stk {
  namespace adapt {

    class URP_Heterogeneous_3D : public UniformRefinerPatternBase
    {

      std::vector<UniformRefinerPatternBase *> m_bp;

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

        // refine

        // put them in reverse topological rank order

        m_bp.push_back(new  UniformRefinerPattern<shards::Tetrahedron<4>,   shards::Tetrahedron<4>,   8, SierraPort > (eMesh, block_names));
        m_bp.push_back(new  UniformRefinerPattern<shards::Hexahedron<8>,    shards::Hexahedron<8>,    8, SierraPort > (eMesh, block_names));
        m_bp.push_back(new  UniformRefinerPattern<shards::Wedge<6>,         shards::Wedge<6>,         8, SierraPort > (eMesh, block_names));
        m_bp.push_back(new  UniformRefinerPattern<shards::Hexahedron<27>,   shards::Hexahedron<27>,   8, SierraPort > (eMesh, block_names));
        m_bp.push_back(new  UniformRefinerPattern<shards::Hexahedron<20>,   shards::Hexahedron<20>,   8, SierraPort > (eMesh, block_names));
        m_bp.push_back(new  UniformRefinerPattern<shards::Tetrahedron<10>,  shards::Tetrahedron<10>,  8, SierraPort > (eMesh, block_names));
        m_bp.push_back(new  UniformRefinerPattern<shards::Wedge<15>,        shards::Wedge<15>,        8, SierraPort > (eMesh, block_names));


        m_bp.push_back(new  UniformRefinerPattern<shards::ShellQuadrilateral<4>, shards::ShellQuadrilateral<4>, 4, SierraPort > (eMesh, block_names) );
        m_bp.push_back(new  UniformRefinerPattern<shards::ShellTriangle<3>,      shards::ShellTriangle<3>,      4, SierraPort > (eMesh, block_names) );
        //m_bp.push_back(new  UniformRefinerPattern<shards::ShellTriangle<6>,      shards::ShellTriangle<6>,      4, SierraPort > (eMesh, block_names) );
        m_bp.push_back(new  UniformRefinerPattern<shards::ShellQuadrilateral<8>, shards::ShellQuadrilateral<8>, 4, SierraPort > (eMesh, block_names));

        m_bp.push_back(new  UniformRefinerPattern<shards::Quadrilateral<4>, shards::Quadrilateral<4>, 4, SierraPort > (eMesh, block_names) );
        m_bp.push_back(new  UniformRefinerPattern<shards::Triangle<3>,      shards::Triangle<3>,      4, SierraPort > (eMesh, block_names) );
        //m_bp.push_back(new  UniformRefinerPattern<shards::Triangle<6>,      shards::Triangle<6>,      4, SierraPort > (eMesh, block_names));
        m_bp.push_back(new  UniformRefinerPattern<shards::Quadrilateral<9>, shards::Quadrilateral<9>, 4, SierraPort > (eMesh, block_names));
        m_bp.push_back(new  UniformRefinerPattern<shards::Quadrilateral<8>, shards::Quadrilateral<8>, 4, SierraPort > (eMesh, block_names));

#if 0
        m_bp.push_back(new  UniformRefinerPattern<shards::ShellLine<2>,     shards::ShellLine<2>,     2, SierraPort > (eMesh, block_names));

        m_bp.push_back(new  UniformRefinerPattern<shards::Line<2>,          shards::Line<2>,          2, SierraPort > (eMesh, block_names));
        m_bp.push_back(new  UniformRefinerPattern<shards::Line<3>,          shards::Line<3>,          2, SierraPort > (eMesh, block_names));
#endif

      }

      void setSubPatterns( std::vector<UniformRefinerPatternBase *>& bp, percept::PerceptMesh& eMesh )
      {
        EXCEPTWATCH;

        bp = m_bp;
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
