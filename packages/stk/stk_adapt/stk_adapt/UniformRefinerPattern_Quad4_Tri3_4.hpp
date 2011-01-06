#ifndef stk_adapt_UniformRefinerPattern_Quad4_Tri3_4_hpp
#define stk_adapt_UniformRefinerPattern_Quad4_Tri3_4_hpp


//#include "UniformRefinerPattern.hpp"

namespace stk {
  namespace adapt {

    struct Specialization {};

    template <>
    class UniformRefinerPattern<shards::Quadrilateral<4>, shards::Triangle<3>, 4, Specialization > : public URP<shards::Quadrilateral<4> , shards::Triangle<3> >
    {
    public:

      UniformRefinerPattern(percept::PerceptMesh& eMesh, BlockNamesType block_names = BlockNamesType()) : URP<shards::Quadrilateral<4> , shards::Triangle<3> >(eMesh)
       {
         EXCEPTWATCH;
         m_primaryEntityRank = mesh::Face;
         if (m_eMesh.getSpatialDim() == 2)
           m_primaryEntityRank = mesh::Element;

         setNeededParts(eMesh, block_names, false);

       }

      virtual void doBreak() {}

      void fillNeededEntities(std::vector<NeededEntityType>& needed_entities)
      {
        needed_entities.resize(1);
        needed_entities[0].first = (m_eMesh.getSpatialDim() == 2 ? stk::mesh::Element :  stk::mesh::Face);
        setToOne(needed_entities);
      }

      virtual unsigned getNumNewElemPerElem() { return 4; }

      virtual StringStringMap fixSurfaceAndEdgeSetNamesMap()
      {
        StringStringMap str_map;
        str_map["hex8"] = "tet4";
        str_map["quad4"] = "tri3";
        return str_map;
      }

      void 
      createNewElements(percept::PerceptMesh& eMesh, NodeRegistry& nodeRegistry, 
                        Entity& element,  NewSubEntityNodesType& new_sub_entity_nodes, vector<Entity *>::iterator& element_pool,
                        FieldBase *proc_rank_field=0)
      {
        const CellTopologyData * const cell_topo_data = get_cell_topology(element);
        typedef boost::tuple<EntityId, EntityId, EntityId> tri_tuple_type;
        static vector<tri_tuple_type> elems(4);

        CellTopology cell_topo(cell_topo_data);
        const mesh::PairIterRelation elem_nodes = element.relations(Node);

        //stk::mesh::Part & active = mesh->ActivePart();
        //stk::mesh::Part & quad4  = mesh->QuadPart();

        std::vector<stk::mesh::Part*> add_parts;
        std::vector<stk::mesh::Part*> remove_parts;

        add_parts = m_toParts;
        
        //std::cout << "P["<< m_eMesh.getRank() << "] add_parts = " << add_parts << std::endl;

        EntityRank my_rank = m_primaryEntityRank;

        //         EntityRank my_rank = stk::mesh::Face;
        //         if (m_eMesh.getSpatialDim() == 2)
        //           my_rank = stk::mesh::Element;

        nodeRegistry.makeCentroid(*const_cast<Entity *>(&element), my_rank, 0u);
        nodeRegistry.addToExistingParts(*const_cast<Entity *>(&element), my_rank, 0u);
        nodeRegistry.interpolateFields(*const_cast<Entity *>(&element), my_rank, 0u);
        
#define CENTROID_N NN(m_primaryEntityRank, 0)  

        elems[0] = tri_tuple_type(VERT_N(0), VERT_N(1), CENTROID_N);
        elems[1] = tri_tuple_type(VERT_N(1), VERT_N(2), CENTROID_N);
        elems[2] = tri_tuple_type(VERT_N(2), VERT_N(3), CENTROID_N);
        elems[3] = tri_tuple_type(VERT_N(3), VERT_N(0), CENTROID_N);

#undef CENTROID_N

        // write a diagram of the refinement pattern as a vtk file, or a latex/tikz/pgf file
#define WRITE_DIAGRAM 0
#if WRITE_DIAGRAM

        /**
           \node[above] at (p4.side 1){2};
           \node[left] at (p4.side 2){3};
           \node[below] at (p4.side 3){0};
           \node[right] at (p4.side 4){1};
        */

#endif
        
        for (unsigned ielem=0; ielem < elems.size(); ielem++)
          {
            Entity& newElement = *(*element_pool);

            //std::cout << "P["<< m_eMesh.getRank() << "] urp tmp 3 "  << proc_rank_field << std::endl;
            if (proc_rank_field && element.entity_rank() == mesh::Element)
              {
                double *fdata = stk::mesh::field_data( *static_cast<const ScalarFieldType *>(proc_rank_field) , newElement );
                fdata[0] = double(newElement.owner_rank());
              }

            //std::cout << "P["<< m_eMesh.getRank() << "] urp tmp 4 "  << std::endl;
            change_entity_parts(eMesh, element, newElement);

            //std::cout << "P["<< m_eMesh.getRank() << "] urp tmp 5 "  << std::endl;

            {
              if (!elems[ielem].get<0>())
                {
                  std::cout << "P[" << eMesh.getRank() << " nid = 0 << " << std::endl;
                  exit(1);
                }

            }
            //std::cout << "P["<< m_eMesh.getRank() << "] urp tmp 6 "  << std::endl;

            eMesh.getBulkData()->declare_relation(newElement, eMesh.createOrGetNode(elems[ielem].get<0>()), 0);
            eMesh.getBulkData()->declare_relation(newElement, eMesh.createOrGetNode(elems[ielem].get<1>()), 1);
            eMesh.getBulkData()->declare_relation(newElement, eMesh.createOrGetNode(elems[ielem].get<2>()), 2);

            //std::cout << "P["<< m_eMesh.getRank() << "] urp tmp 7 "  << std::endl;

            element_pool++;

          }

      }
      
    };

  }
}
#endif
