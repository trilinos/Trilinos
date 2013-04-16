#ifndef stk_adapt_RefinerPattern_Quad4_Quad4_HangingNode_hpp
#define stk_adapt_RefinerPattern_Quad4_Quad4_HangingNode_hpp


//#include "UniformRefinerPattern.hpp"
#include <stk_adapt/sierra_element/RefinementTopology.hpp>
#include <stk_adapt/sierra_element/StdMeshObjTopologies.hpp>

#include "RefinerPattern_Line2_Line2_N.hpp"

namespace stk {
  namespace adapt {

    template <>
    class RefinerPattern<shards::Quadrilateral<4>, shards::Quadrilateral<4>, -1 > : public URP<shards::Quadrilateral<4>,shards::Quadrilateral<4>  >
    {

      RefinerPattern<shards::Line<2>, shards::Line<2>, -1 > * m_edge_breaker;

    public:

      RefinerPattern(percept::PerceptMesh& eMesh, BlockNamesType block_names = BlockNamesType()) :  URP<shards::Quadrilateral<4>, shards::Quadrilateral<4>  >(eMesh)
      {
        m_primaryEntityRank = eMesh.face_rank();
        if (m_eMesh.get_spatial_dim() == 2)
          m_primaryEntityRank = stk::mesh::MetaData::ELEMENT_RANK;

        setNeededParts(eMesh, block_names, true);
        Elem::StdMeshObjTopologies::bootstrap();

        if (m_eMesh.get_spatial_dim() == 2)
          m_edge_breaker =  new RefinerPattern<shards::Line<2>, shards::Line<2>, -1 > (eMesh, block_names) ;
        else
          m_edge_breaker = 0;

      }
      ~RefinerPattern()
      {
        if (m_edge_breaker) delete m_edge_breaker;
      }

      void setSubPatterns( std::vector<UniformRefinerPatternBase *>& bp, percept::PerceptMesh& eMesh )
      {
        EXCEPTWATCH;
        bp = std::vector<UniformRefinerPatternBase *>(2u, 0);

        if (eMesh.get_spatial_dim() == 2)
          {
            bp[0] = this;
            bp[1] = m_edge_breaker;
          }
        else if (eMesh.get_spatial_dim() == 3)
          {
            // FIXME
            //             std::cout << "ERROR" ;
            //             exit(1);
          }

      }

      virtual void doBreak() {}
      void fillNeededEntities(std::vector<NeededEntityType>& needed_entities)
      {
        needed_entities.resize(2);
        needed_entities[0].first = m_eMesh.edge_rank();
        needed_entities[1].first = (m_eMesh.get_spatial_dim() == 2 ? stk::mesh::MetaData::ELEMENT_RANK : m_eMesh.face_rank());
        setToOne(needed_entities);
      }

      virtual unsigned getNumNewElemPerElem() { return 4; }

      void
      createNewElements(percept::PerceptMesh& eMesh, NodeRegistry& nodeRegistry,
                        stk::mesh::Entity element,  NewSubEntityNodesType& new_sub_entity_nodes, vector<stk::mesh::Entity>::iterator& element_pool,
                        stk::mesh::FieldBase *proc_rank_field=0)
      {
        unsigned num_edges_marked=0;
        for (int iedge = 0; iedge < 4; iedge++)
          {
            unsigned num_nodes_on_edge = new_sub_entity_nodes[m_eMesh.edge_rank()][iedge].size();
            if (num_nodes_on_edge)
              {
                ++num_edges_marked;
              }
          }
        if (num_edges_marked != 4)
          {
            return;
          }

        genericRefine_createNewElements(eMesh, nodeRegistry,
                                        element, new_sub_entity_nodes, element_pool,
                                        proc_rank_field);
      }
    };

  }
}
#endif
