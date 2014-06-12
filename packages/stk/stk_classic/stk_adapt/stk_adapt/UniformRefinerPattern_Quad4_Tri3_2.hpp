#ifndef stk_adapt_UniformRefinerPattern_Quad4_Tri3_2_hpp
#define stk_adapt_UniformRefinerPattern_Quad4_Tri3_2_hpp


//#include "UniformRefinerPattern.hpp"

namespace stk_classic {
  namespace adapt {

    //struct Specialization {};

    template <>
    class UniformRefinerPattern<shards::Quadrilateral<4>, shards::Triangle<3>, 2  > : public URP<shards::Quadrilateral<4> , shards::Triangle<3> >
    {
    public:

      UniformRefinerPattern(percept::PerceptMesh& eMesh, BlockNamesType block_names = BlockNamesType()) : URP<shards::Quadrilateral<4> , shards::Triangle<3> >(eMesh)
       {
         EXCEPTWATCH;
         m_primaryEntityRank = eMesh.face_rank();
         if (m_eMesh.get_spatial_dim() == 2)
           m_primaryEntityRank = eMesh.element_rank();

         setNeededParts(eMesh, block_names, false);

       }

      virtual void doBreak() {}

      void fillNeededEntities(std::vector<NeededEntityType>& needed_entities)
      {
        needed_entities.resize(0);
        //needed_entities[0] = (m_eMesh.get_spatial_dim() == 2 ? m_eMesh.element_rank() :  m_eMesh.face_rank());
        setToOne(needed_entities);
      }

      virtual unsigned getNumNewElemPerElem() { return 2; }

      virtual StringStringMap fixSurfaceAndEdgeSetNamesMap()
      {
        StringStringMap str_map;
        str_map["hex8"] = "tet4";
        str_map["quad4"] = "tri3";
        return str_map;
      }

      void 
      createNewElements(percept::PerceptMesh& eMesh, NodeRegistry& nodeRegistry, 
                        stk_classic::mesh::Entity& element,  NewSubEntityNodesType& new_sub_entity_nodes, vector<stk_classic::mesh::Entity *>::iterator& element_pool,
                        stk_classic::mesh::FieldBase *proc_rank_field=0)
      {
        const CellTopologyData * const cell_topo_data = stk_classic::percept::PerceptMesh::get_cell_topology(element);
        typedef boost::tuple<stk_classic::mesh::EntityId, stk_classic::mesh::EntityId, stk_classic::mesh::EntityId> tri_tuple_type;
        static vector<tri_tuple_type> elems(2);

        CellTopology cell_topo(cell_topo_data);
        const stk_classic::mesh::PairIterRelation elem_nodes = element.relations(stk_classic::mesh::fem::FEMMetaData::NODE_RANK);

        //std::cout << "tmp quad elem= " << element << std::endl;

        //stk_classic::mesh::Part & active = mesh->ActivePart();
        //stk_classic::mesh::Part & quad4  = mesh->QuadPart();

        std::vector<stk_classic::mesh::Part*> add_parts;
        std::vector<stk_classic::mesh::Part*> remove_parts;

        add_parts = m_toParts;
        
        //std::cout << "P["<< m_eMesh.get_rank() << "] add_parts = " << add_parts << std::endl;

        //stk_classic::mesh::EntityRank my_rank = m_primaryEntityRank;

        //nodeRegistry.makeCentroid(*const_cast<stk_classic::mesh::Entity *>(&element), my_rank, 0u);
        //nodeRegistry.addToExistingParts(*const_cast<stk_classic::mesh::Entity *>(&element), my_rank, 0u);
        //nodeRegistry.interpolateFields(*const_cast<stk_classic::mesh::Entity *>(&element), my_rank, 0u);
        
        {
          unsigned globalIqf  = VERT_N(0);
          unsigned minVal     = globalIqf;
          unsigned indxMinVal = 0;
          for (unsigned iFaceNodeOrd=1; iFaceNodeOrd < 4; iFaceNodeOrd++)
            {
              globalIqf = VERT_N(iFaceNodeOrd);
              if (globalIqf < minVal)
                {
                  minVal = globalIqf;
                  indxMinVal = iFaceNodeOrd;
                }
            }

          unsigned istart = indxMinVal;
          elems[0] = tri_tuple_type(VERT_N((0 + istart) % 4), VERT_N((1 + istart) % 4), VERT_N((2 + istart) % 4));
          elems[1] = tri_tuple_type(VERT_N((0 + istart) % 4), VERT_N((2 + istart) % 4), VERT_N((3 + istart) % 4));
        }

        for (unsigned ielem=0; ielem < elems.size(); ielem++)
          {
            stk_classic::mesh::Entity& newElement = *(*element_pool);

            if (proc_rank_field && element.entity_rank() == m_eMesh.element_rank())
              {
                double *fdata = stk_classic::mesh::field_data( *static_cast<const ScalarFieldType *>(proc_rank_field) , newElement );
                fdata[0] = double(newElement.owner_rank());
              }

            change_entity_parts(eMesh, element, newElement);

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

            set_parent_child_relations(eMesh, element, newElement, ielem);

            element_pool++;

          }

      }
      
    };

  }
}
#endif
