#ifndef stk_adapt_UniformRefinerPattern_Quad4_Tri3_6_hpp
#define stk_adapt_UniformRefinerPattern_Quad4_Tri3_6_hpp


//#include "UniformRefinerPattern.hpp"

namespace stk {
  namespace adapt {

    template <>
    class UniformRefinerPattern<shards::Quadrilateral<4>, shards::Triangle<3>, 6 > : public URP<shards::Quadrilateral<4> , shards::Triangle<3> >
    {
    public:

      UniformRefinerPattern(percept::PerceptMesh& eMesh, BlockNamesType block_names = BlockNamesType()) : URP<shards::Quadrilateral<4> , shards::Triangle<3> >(eMesh)
       {
         m_primaryEntityRank = mesh::Face;
         if (m_eMesh.getSpatialDim() == 2)
           m_primaryEntityRank = mesh::Element;

         setNeededParts(eMesh, block_names, false);
       }

      virtual void doBreak() {}
      void fillNeededEntities(std::vector<NeededEntityType>& needed_entities)
      {
        needed_entities.resize(1);
        needed_entities[0].first = stk::mesh::Edge; // edges have 2 nodes
        setToOne(needed_entities);
      }

      virtual unsigned getNumNewElemPerElem() { return 6; }

      void 
      createNewElements(percept::PerceptMesh& eMesh, NodeRegistry& nodeRegistry, 
                        Entity& element,  NewSubEntityNodesType& new_sub_entity_nodes, vector<Entity *>::iterator& element_pool,
                        FieldBase *proc_rank_field=0)
      {
        const CellTopologyData * const cell_topo_data = stk::percept::PerceptMesh::get_cell_topology(element);
        typedef boost::tuple<EntityId, EntityId, EntityId> tri_tuple_type;
        static vector<tri_tuple_type> elems(6);

        CellTopology cell_topo(cell_topo_data);
        const mesh::PairIterRelation elem_nodes = element.relations(Node);

        //stk::mesh::Part & active = mesh->ActivePart();
        //stk::mesh::Part & quad4  = mesh->QuadPart();

        std::vector<stk::mesh::Part*> add_parts;
        std::vector<stk::mesh::Part*> remove_parts;

        //add_parts.push_back( &active );
        //FIXME 
        //add_parts.push_back( const_cast<mesh::Part*>( eMesh.getPart(m_toTopoPartName) ));
        add_parts = m_toParts;
        
        /**
           \node[above] at (p4.side 1){2};
           \node[left] at (p4.side 2){3};
           \node[below] at (p4.side 3){0};
           \node[right] at (p4.side 4){1};
        */



        double tmp_x[3];
        for (int iedge = 0; iedge < 4; iedge++)
          {
            double * mp = midPoint(EDGE_COORD(iedge,0), EDGE_COORD(iedge,1), eMesh.getSpatialDim(), tmp_x);

            if (!EDGE_N(iedge))
              {
                std::cout << "P[" << eMesh.getRank() << " nid ## = 0 << " << std::endl;
              }
            eMesh.createOrGetNode(EDGE_N(iedge), mp);

          }


        elems[0] = tri_tuple_type(VERT_N(0), EDGE_N(0), EDGE_N(3));
        elems[1] = tri_tuple_type(VERT_N(1), EDGE_N(1), EDGE_N(0));
        elems[2] = tri_tuple_type(EDGE_N(0), EDGE_N(1), EDGE_N(3));

        elems[3] = tri_tuple_type(VERT_N(2), EDGE_N(2), EDGE_N(1));
        elems[4] = tri_tuple_type(VERT_N(3), EDGE_N(3), EDGE_N(2));
        elems[5] = tri_tuple_type(EDGE_N(2), EDGE_N(3), EDGE_N(1));

        // write a diagram of the refinement pattern as a vtk file, or a latex/tikz/pgf file
#define WRITE_DIAGRAM 0
#if WRITE_DIAGRAM

#endif
        
        for (unsigned ielem=0; ielem < elems.size(); ielem++)
          {
            //Entity& newElement = eMesh.getBulkData()->declare_entity(Element, *element_id_pool, eMesh.getPart(interface_table::shards_Triangle_3) );
            //Entity& newElement = eMesh.getBulkData()->declare_entity(Element, *element_id_pool, eMesh.getPart(interface_table::shards_Triangle_3) );

            Entity& newElement = *(*element_pool);

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
                  std::cout << "P[" << eMesh.getRank() << " nid = 0 << " << std::endl;
                  exit(1);
                }

            }
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
