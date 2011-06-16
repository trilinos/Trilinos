#ifndef stk_adapt_RefinerPattern_Tri3_Tri3_N_sierra_hpp
#define stk_adapt_RefinerPattern_Tri3_Tri3_N_sierra_hpp


//#include "UniformRefinerPattern.hpp"
#include <stk_adapt/sierra_element/RefinementTopology.hpp>
#include <stk_adapt/sierra_element/StdMeshObjTopologies.hpp>

#include "UniformRefinerPattern_Line2_Line2_2_sierra.hpp"

namespace stk {
  namespace adapt {

    /// general refinement pattern
    
    // the "-1" here signifies the number of elements created is not fixed, depends on the marking pattern
    template <>
    class RefinerPattern<shards::Triangle<3>, shards::Triangle<3>, -1 > : public URP<shards::Triangle<3>,shards::Triangle<3>  >
    {

      UniformRefinerPattern<shards::Line<2>, shards::Line<2>, 2, SierraPort > * m_edge_breaker;

    public:

      RefinerPattern(percept::PerceptMesh& eMesh, BlockNamesType block_names = BlockNamesType()) :  URP<shards::Triangle<3>, shards::Triangle<3>  >(eMesh)
      {
        m_primaryEntityRank = m_eMesh.face_rank();
        if (m_eMesh.getSpatialDim() == 2)
          m_primaryEntityRank = eMesh.element_rank();

        setNeededParts(eMesh, block_names, true);
        Elem::StdMeshObjTopologies::bootstrap();

        if (m_eMesh.getSpatialDim() == 2)
          {
            m_edge_breaker =  new UniformRefinerPattern<shards::Line<2>, shards::Line<2>, 2, SierraPort > (eMesh, block_names) ;
          }

      }

      void setSubPatterns( std::vector<UniformRefinerPatternBase *>& bp, percept::PerceptMesh& eMesh )
      {
        EXCEPTWATCH;
        bp = std::vector<UniformRefinerPatternBase *>(2u, 0);

        if (eMesh.getSpatialDim() == 2)
          {
            bp[0] = this;
            if (m_eMesh.getSpatialDim() == 2)
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

      // FIXME - for now, create more than we need (to fix this right we need a change to the Refiner.cpp interface)
      virtual unsigned getNumNewElemPerElem() { return 4; }

      void 
      createNewElements(percept::PerceptMesh& eMesh, NodeRegistry& nodeRegistry, 
                        stk::mesh::Entity& element,  NewSubEntityNodesType& new_sub_entity_nodes, vector<stk::mesh::Entity *>::iterator& element_pool,
                        stk::mesh::FieldBase *proc_rank_field=0)
      {
        const CellTopologyData * const cell_topo_data = stk::percept::PerceptMesh::get_cell_topology(element);
        typedef boost::tuple<stk::mesh::EntityId, stk::mesh::EntityId, stk::mesh::EntityId> tri_tuple_type;
        typedef boost::tuple<int, int, int> tri_tuple_type_int;
        static vector<tri_tuple_type> elems(4);
        //int num_new_elems=0;

        CellTopology cell_topo(cell_topo_data);
        const stk::mesh::PairIterRelation elem_nodes = element.relations(stk::mesh::fem::FEMMetaData::NODE_RANK);
        VectorFieldType* coordField = eMesh.getCoordinatesField();

        std::vector<stk::mesh::Part*> add_parts;
        std::vector<stk::mesh::Part*> remove_parts;
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

        //         if (num_edges_marked > 1)
        //           {
        //             throw std::runtime_error("RefinerPattern_Tri3_Tri3_2 can only refine element with one marked edge");
        //           }

//         if (num_edges_marked == 0)
//           return;
        
        // if (num_edges_marked == 2)
        //  return;

        if (num_edges_marked == 3)
          {
            elems.resize(4);

            elems[0] = tri_tuple_type( VERT_N(0),    EDGE_N(0), EDGE_N(2) );
            elems[1] = tri_tuple_type( VERT_N(1),    EDGE_N(1), EDGE_N(0) );
            elems[2] = tri_tuple_type( VERT_N(2),    EDGE_N(2), EDGE_N(1) );
            elems[3] = tri_tuple_type( EDGE_N(0),    EDGE_N(1), EDGE_N(2) );
          }
        else if (num_edges_marked == 2)
          {
            // Note: this will form the basis of triangulating faces in 3D, so it should be generalized to a 
            //   generic method: triangulate_face(node_coordinates, node_ids, edge_coords, edge_ids, returned_triangles)
            // Note: code below doesn't orient the face except for a rotation - we need a polarity flip check as
            //   well for the general, 3D face case
            //
            /**
             *
             *   case 1: jedge == max length edge
             *
             *                i2  
             *                o
             *               /|\
             *              / | \
             *             /  |  \
             *            /   |   * jedgep
             *           /    |  / \
             *          /     | /   \
             *         /      |/     \ 
             *        o-------*-------o
             *       i0      jedge     i1
             *
             *
             *   case 2: jedge+1 == max length edge
             *
             *                i2 
             *                o
             *               / \ 
             *              /   \
             *             /     \
             *            /     _.* jedgep
             *           /   _.* / \
             *          / _.*   /   \
             *         /.*     /     \
             *        o-------*-------o
             *       i0      jedge     i1
             *
             */

            elems.resize(3);

            // find first of two marked edges in sequence (get in "standard" orientation), and longest marked edge
            int jedge = -1;
            int jedge_max_edge = -1;
            double max_edge_length = -1.0;
            unsigned id_diff_0 = 0u;
            unsigned id_diff_1 = 0u;
            for (int iedge = 0; iedge < 3; iedge++)
              {

                unsigned num_nodes_on_edge = new_sub_entity_nodes[m_eMesh.edge_rank()][iedge].size();
                unsigned num_nodes_on_edge_p = new_sub_entity_nodes[m_eMesh.edge_rank()][(iedge+1)%3].size();
                if (num_nodes_on_edge && num_nodes_on_edge_p)
                  {
                    jedge = iedge;
                  }

                if (num_nodes_on_edge)
                  {
                    stk::mesh::Entity * node_0 = elem_nodes[cell_topo_data->edge[iedge].node[0]].entity();
                    stk::mesh::Entity * node_1 = elem_nodes[cell_topo_data->edge[iedge].node[1]].entity();

                    bool reverse = false;
                    // ensure edge_len is computed identically, independent of edge orientation
                    if (node_0->identifier() > node_1->identifier())
                      {
                        reverse = true;
                        stk::mesh::Entity *node_temp = node_0;
                        node_0 = node_1;
                        node_1 = node_temp;
                      }

                    double * const coord_0 = stk::mesh::field_data( *coordField , *node_0 );
                    double * const coord_1 = stk::mesh::field_data( *coordField , *node_1 );
                    double edge_len_squared = 0.0;

                    edge_len_squared = 
                      (coord_0[0] - coord_1[0])*(coord_0[0] - coord_1[0])+
                      (coord_0[1] - coord_1[1])*(coord_0[1] - coord_1[1])+
                      (m_eMesh.getSpatialDim() == 2 ? 0 : 
                       (coord_0[2] - coord_1[2])*(coord_0[2] - coord_1[2]) );

                    if (edge_len_squared > max_edge_length)
                      {
                        id_diff_0 = node_0->identifier();
                        id_diff_1 = node_1->identifier();
                        max_edge_length = edge_len_squared;
                        jedge_max_edge = iedge;
                      }
                    // intentional floating-point comparison (tie-break)
                    else if (edge_len_squared == max_edge_length)
                      {
                        unsigned loc_id_diff_0 = node_0->identifier();
                        unsigned loc_id_diff_1 = node_1->identifier();
                        bool lexical_less = false;
                        if (loc_id_diff_0 < id_diff_0)
                          {
                            lexical_less = true;
                          }
                        else if (loc_id_diff_0 == id_diff_0 && loc_id_diff_1 < id_diff_1)
                          {
                            lexical_less = true;
                          }
                        if (!lexical_less)
                          {
                            max_edge_length = edge_len_squared;
                            jedge_max_edge = iedge;
                          }
                      }
                  }
              }

            if (jedge < 0 || jedge_max_edge < 0) 
              {
                std::cout << "jedge = " << jedge << " jedge_max_edge = " << jedge_max_edge << std::endl;
                throw std::runtime_error("RefinerPattern_Tri3_Tri3_N jedge < 0");
              }

            //stk::mesh::Entity & node0 = *elem_nodes[iii].entity();
            int i0 = cell_topo_data->edge[jedge].node[0];
            if (i0 != jedge)
              {
                std::cout << "i0 = " << i0 << " jedge= " << jedge << std::endl;
                throw std::runtime_error("RefinerPattern_Tri3_Tri3_N i0 != jedge");
              }

            int i1 = (i0+1)%3;
            int i2 = (i0+2)%3;
            int jedgep = (jedge+1)%3;
            if (jedge_max_edge == jedge)
              {
                elems[0] = tri_tuple_type( VERT_N(i0),    EDGE_N(jedge),         VERT_N(i2)       );
                elems[1] = tri_tuple_type( EDGE_N(jedge), EDGE_N( jedgep ),      VERT_N(i2)       );
                elems[2] = tri_tuple_type( EDGE_N(jedge), VERT_N(i1),            EDGE_N( jedgep ) );
              }
            else
              {
                elems[0] = tri_tuple_type( VERT_N(i0),    EDGE_N(jedge),         EDGE_N( jedgep ) );
                elems[1] = tri_tuple_type( VERT_N(i0),    EDGE_N( jedgep ),      VERT_N(i2)       );
                elems[2] = tri_tuple_type( EDGE_N(jedge), VERT_N(i1),            EDGE_N( jedgep ) );
              }
          }
        else if (num_edges_marked == 1)
          {
            elems.resize(2);
            for (int iedge = 0; iedge < 3; iedge++)
              {
                unsigned num_nodes_on_edge = new_sub_entity_nodes[m_eMesh.edge_rank()][iedge].size();
                if (num_nodes_on_edge)
                  {
                    elems[0] = tri_tuple_type(VERT_N(iedge), EDGE_N(iedge), VERT_N((iedge+2)%3) );
                    elems[1] = tri_tuple_type(EDGE_N(iedge), VERT_N((iedge+1)%3), VERT_N((iedge+2)%3) );
                    break;
                  }
              }
          }
        else if (num_edges_marked == 0)
          {
            elems.resize(1);
            elems[0] = tri_tuple_type(VERT_N(0), VERT_N(1), VERT_N(2) );
          }

        //nodeRegistry.makeCentroidCoords(*const_cast<stk::mesh::Entity *>(&element), m_eMesh.element_rank(), 0u);
        
        for (unsigned ielem=0; ielem < elems.size(); ielem++)
          {
            stk::mesh::Entity& newElement = *(*element_pool);

            if (proc_rank_field)
              {
                double *fdata = stk::mesh::field_data( *static_cast<const ScalarFieldType *>(proc_rank_field) , newElement );
                //fdata[0] = double(m_eMesh.getRank());
                fdata[0] = double(newElement.owner_rank());
              }

            eMesh.getBulkData()->change_entity_parts( newElement, add_parts, remove_parts );

            set_parent_child_relations(eMesh, element, newElement, ielem);

            {
              if (!elems[ielem].get<0>())
                {
                  std::cout << "P[" << eMesh.getRank() << "] nid = 0 << " << std::endl;
                  //exit(1);
                }
            }

            // 3 nodes of the new tris
            eMesh.getBulkData()->declare_relation(newElement, eMesh.createOrGetNode(elems[ielem].get<0>()), 0);
            eMesh.getBulkData()->declare_relation(newElement, eMesh.createOrGetNode(elems[ielem].get<1>()), 1);
            eMesh.getBulkData()->declare_relation(newElement, eMesh.createOrGetNode(elems[ielem].get<2>()), 2);

            element_pool++;

          }

      
      }
      
    };

  }
}
#endif
