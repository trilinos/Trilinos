#ifndef stk_adapt_UniformRefinerPattern_Quad4_Quad4_4_hpp
#define stk_adapt_UniformRefinerPattern_Quad4_Quad4_4_hpp


//#include "UniformRefinerPattern.hpp"

namespace stk_classic {
  namespace adapt {

    template <>
    class UniformRefinerPattern<shards::Quadrilateral<4>, shards::Quadrilateral<4>, 4 > : public URP<shards::Quadrilateral<4>,shards::Quadrilateral<4>  >
    {
    public:

//       UniformRefinerPattern(percept::PerceptMesh& eMesh, std::string fromTopoPartName="block_1", std::string toTopoPartName="block_quad_4")
//       {
//         setNeededParts(eMesh, fromTopoPartName, toTopoPartName);
//       }
      UniformRefinerPattern(percept::PerceptMesh& eMesh, BlockNamesType block_names = BlockNamesType()) : URP<shards::Quadrilateral<4>, shards::Quadrilateral<4>  >(eMesh)
      {
        m_primaryEntityRank = eMesh.face_rank(); 
        if (m_eMesh.get_spatial_dim() == 2)
          m_primaryEntityRank = eMesh.element_rank();

        setNeededParts(eMesh, block_names, true);
      }

      virtual void doBreak() {}
      void fillNeededEntities(std::vector<NeededEntityType>& needed_entities)
      {
        needed_entities.resize(2);
        needed_entities[0].first = m_eMesh.edge_rank();    // edges have 2 nodes
        needed_entities[1].first = m_eMesh.element_rank(); 
        setToOne(needed_entities);
      }

      virtual unsigned getNumNewElemPerElem() { return 4; }

      void 
      createNewElements(percept::PerceptMesh& eMesh, NodeRegistry& nodeRegistry, 
                        stk_classic::mesh::Entity& element,  NewSubEntityNodesType& new_sub_entity_nodes, vector<stk_classic::mesh::Entity *>::iterator& element_pool,
                        stk_classic::mesh::FieldBase *proc_rank_field=0)
      {
        const CellTopologyData * const cell_topo_data = stk_classic::percept::PerceptMesh::get_cell_topology(element);
        typedef boost::tuple<stk_classic::mesh::EntityId, stk_classic::mesh::EntityId, stk_classic::mesh::EntityId, stk_classic::mesh::EntityId> quad_tuple_type;
        static vector<quad_tuple_type> elems(4);

        CellTopology cell_topo(cell_topo_data);
        const stk_classic::mesh::PairIterRelation elem_nodes = element.relations(stk_classic::mesh::fem::FEMMetaData::NODE_RANK);

        //stk_classic::mesh::Part & active = mesh->ActivePart();
        //stk_classic::mesh::Part & quad4  = mesh->QuadPart();

        std::vector<stk_classic::mesh::Part*> add_parts;
        std::vector<stk_classic::mesh::Part*> remove_parts;

        //add_parts.push_back( &active );
        //FIXME 
        //add_parts.push_back( const_cast<mesh::Part*>( eMesh.getPart(m_toTopoPartName) ));
        add_parts = m_toParts;
        
        double tmp_x[3];
        for (int iedge = 0; iedge < 4; iedge++)
          {
            double * mp = midPoint(EDGE_COORD(iedge,0), EDGE_COORD(iedge,1), eMesh.get_spatial_dim(), tmp_x);

            if (!EDGE_N(iedge))
              {
                std::cout << "P[" << eMesh.get_rank() << " nid ## = 0 << " << std::endl;
              }
            eMesh.createOrGetNode(EDGE_N(iedge), mp);

          }

        nodeRegistry.makeCentroidCoords(*const_cast<stk_classic::mesh::Entity *>(&element), m_eMesh.element_rank(), 0u);


// new_sub_entity_nodes[i][j]
#define CENTROID_N NN(m_primaryEntityRank,0)  

        elems[0] = quad_tuple_type(VERT_N(0), EDGE_N(0), CENTROID_N, EDGE_N(3));
        elems[1] = quad_tuple_type(VERT_N(1), EDGE_N(1), CENTROID_N, EDGE_N(0));
        elems[2] = quad_tuple_type(VERT_N(2), EDGE_N(2), CENTROID_N, EDGE_N(1));
        elems[3] = quad_tuple_type(VERT_N(3), EDGE_N(3), CENTROID_N, EDGE_N(2));

#undef CENTROID_N

        // write a diagram of the refinement pattern as a vtk file, or a latex/tikz/pgf file
#define WRITE_DIAGRAM 0
#if WRITE_DIAGRAM

#endif
        
        for (unsigned ielem=0; ielem < elems.size(); ielem++)
          {
            stk_classic::mesh::Entity& newElement = *(*element_pool);

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
                  std::cout << "P[" << eMesh.get_rank() << " nid = 0 << " << std::endl;
                  exit(1);
                }

            }
            eMesh.get_bulk_data()->declare_relation(newElement, eMesh.createOrGetNode(elems[ielem].get<0>()), 0);
            eMesh.get_bulk_data()->declare_relation(newElement, eMesh.createOrGetNode(elems[ielem].get<1>()), 1);
            eMesh.get_bulk_data()->declare_relation(newElement, eMesh.createOrGetNode(elems[ielem].get<2>()), 2);
            eMesh.get_bulk_data()->declare_relation(newElement, eMesh.createOrGetNode(elems[ielem].get<3>()), 3);

            set_parent_child_relations(eMesh, element, newElement, ielem);


            element_pool++;

          }

      }
      
    };

  }
}
#endif
