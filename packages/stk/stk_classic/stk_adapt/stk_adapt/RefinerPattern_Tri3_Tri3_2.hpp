#ifndef stk_adapt_RefinerPattern_Tri3_Tri3_2_sierra_hpp
#define stk_adapt_RefinerPattern_Tri3_Tri3_2_sierra_hpp


//#include "UniformRefinerPattern.hpp"
#include <stk_adapt/sierra_element/RefinementTopology.hpp>
#include <stk_adapt/sierra_element/StdMeshObjTopologies.hpp>

#include "UniformRefinerPattern_Line2_Line2_2_sierra.hpp"

namespace stk_classic {
  namespace adapt {

    /// this is for testing only - or could be unsed in future for a bisection-based refinement scheme

    template <>
    class RefinerPattern<shards::Triangle<3>, shards::Triangle<3>, 2 > : public URP<shards::Triangle<3>,shards::Triangle<3>  >
    {

      UniformRefinerPattern<shards::Line<2>, shards::Line<2>, 2, SierraPort > * m_edge_breaker;

    public:

      RefinerPattern(percept::PerceptMesh& eMesh, BlockNamesType block_names = BlockNamesType()) :  URP<shards::Triangle<3>, shards::Triangle<3>  >(eMesh),
                                                                                                    m_edge_breaker(0)
      {
        m_primaryEntityRank = m_eMesh.face_rank();
        if (m_eMesh.get_spatial_dim() == 2)
          m_primaryEntityRank = eMesh.element_rank();

        setNeededParts(eMesh, block_names, true);
        Elem::StdMeshObjTopologies::bootstrap();

        if (m_eMesh.get_spatial_dim() == 2)
          {
            m_edge_breaker =  new UniformRefinerPattern<shards::Line<2>, shards::Line<2>, 2, SierraPort > (eMesh, block_names) ;
          }

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
            if (m_eMesh.get_spatial_dim() == 2)
              {
                bp[1] = m_edge_breaker;
              }
          }

      }

      virtual void doBreak() {}
      void fillNeededEntities(std::vector<NeededEntityType>& needed_entities)
      {
        needed_entities.resize(1);
        needed_entities[0].first = m_eMesh.edge_rank();
        needed_entities[0].second = 1u;
      }

      virtual unsigned getNumNewElemPerElem() { return 2; }

      void
      createNewElements(percept::PerceptMesh& eMesh, NodeRegistry& nodeRegistry,
                        stk_classic::mesh::Entity& element,  NewSubEntityNodesType& new_sub_entity_nodes, vector<stk_classic::mesh::Entity *>::iterator& element_pool,
                        stk_classic::mesh::FieldBase *proc_rank_field=0)
      {
        const CellTopologyData * const cell_topo_data = stk_classic::percept::PerceptMesh::get_cell_topology(element);
        typedef boost::tuple<stk_classic::mesh::EntityId, stk_classic::mesh::EntityId, stk_classic::mesh::EntityId> tri_tuple_type;
        static vector<tri_tuple_type> elems(2);

        shards::CellTopology cell_topo(cell_topo_data);
        const stk_classic::mesh::PairIterRelation elem_nodes = element.relations(stk_classic::mesh::fem::FEMMetaData::NODE_RANK);

        std::vector<stk_classic::mesh::Part*> add_parts;
        std::vector<stk_classic::mesh::Part*> remove_parts;
        add_parts = m_toParts;


        unsigned num_edges_marked=0;
        for (int iedge = 0; iedge < 3; iedge++)
          {
            unsigned num_nodes_on_edge = new_sub_entity_nodes[m_eMesh.edge_rank()][iedge].size();
            if (num_nodes_on_edge)
              {
                ++num_edges_marked;
              }
          }
        if (num_edges_marked > 1)
          {
            throw std::runtime_error("RefinerPattern_Tri3_Tri3_2 can only refine element with one marked edge");
          }
        if (num_edges_marked == 0)
          return;

        double tmp_x[3];
        for (int iedge = 0; iedge < 3; iedge++)
          {
            double * mp = midPoint(EDGE_COORD(iedge,0), EDGE_COORD(iedge,1), eMesh.get_spatial_dim(), tmp_x);

            unsigned num_nodes_on_edge = new_sub_entity_nodes[m_eMesh.edge_rank()][iedge].size();
            if (num_nodes_on_edge)
              {
                eMesh.createOrGetNode(EDGE_N(iedge), mp);
                elems[0] = tri_tuple_type(VERT_N(iedge), EDGE_N(iedge), VERT_N((iedge+2)%3) );
                elems[1] = tri_tuple_type(EDGE_N(iedge), VERT_N((iedge+1)%3), VERT_N((iedge+2)%3) );
              }
          }

        //nodeRegistry.makeCentroidCoords(*const_cast<stk_classic::mesh::Entity *>(&element), m_eMesh.element_rank(), 0u);


        for (unsigned ielem=0; ielem < elems.size(); ielem++)
          {
            stk_classic::mesh::Entity& newElement = *(*element_pool);
            //std::cout << "tmp newElement id = " << newElement.identifier() << std::endl;

            if (proc_rank_field)
              {
                double *fdata = stk_classic::mesh::field_data( *static_cast<const ScalarFieldType *>(proc_rank_field) , newElement );
                //fdata[0] = double(m_eMesh.get_rank());
                fdata[0] = double(newElement.owner_rank());
              }

            eMesh.get_bulk_data()->change_entity_parts( newElement, add_parts, remove_parts );

            {
              if (!elems[ielem].get<0>())
                {
                  std::cout << "P[" << eMesh.get_rank() << "] nid = 0 << " << std::endl;
                  //exit(1);
                }
            }

            // 3 nodes of the new tris
            eMesh.get_bulk_data()->declare_relation(newElement, eMesh.createOrGetNode(elems[ielem].get<0>()), 0);
            eMesh.get_bulk_data()->declare_relation(newElement, eMesh.createOrGetNode(elems[ielem].get<1>()), 1);
            eMesh.get_bulk_data()->declare_relation(newElement, eMesh.createOrGetNode(elems[ielem].get<2>()), 2);

            set_parent_child_relations(eMesh, element, newElement, ielem);

            element_pool++;

          }


      }

    };

  }
}
#endif
