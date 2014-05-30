#ifndef stk_adapt_UniformRefinerPattern_ShellLine2_ShellLine2_2_sierra_hpp
#define stk_adapt_UniformRefinerPattern_ShellLine2_ShellLine2_2_sierra_hpp


//#include "UniformRefinerPattern.hpp"
#include <stk_adapt/sierra_element/RefinementTopology.hpp>
#include <stk_adapt/sierra_element/StdMeshObjTopologies.hpp>

namespace stk_classic {
  namespace adapt {

    template <>
    class UniformRefinerPattern<shards::ShellLine<2>, shards::ShellLine<2>, 2, SierraPort > : public URP<shards::ShellLine<2>, shards::ShellLine<2>  >
    {
    public:

      UniformRefinerPattern(percept::PerceptMesh& eMesh, BlockNamesType block_names = BlockNamesType()) :  URP<shards::ShellLine<2>, shards::ShellLine<2>  >(eMesh)
      {
        m_primaryEntityRank = m_eMesh.edge_rank();
        if (m_eMesh.get_spatial_dim() == 1)
          m_primaryEntityRank = eMesh.element_rank();

        setNeededParts(eMesh, block_names, true);
        Elem::StdMeshObjTopologies::bootstrap();

      }

      virtual void doBreak() {}
      void fillNeededEntities(std::vector<NeededEntityType>& needed_entities)
      {
        needed_entities.resize(1);
        needed_entities[0].first = m_eMesh.edge_rank(); 
        setToOne(needed_entities);
      }

      virtual unsigned getNumNewElemPerElem() { return 2; }

      void 
      createNewElements(percept::PerceptMesh& eMesh, NodeRegistry& nodeRegistry, 
                        stk_classic::mesh::Entity& element,  NewSubEntityNodesType& new_sub_entity_nodes, vector<stk_classic::mesh::Entity *>::iterator& element_pool,
                        stk_classic::mesh::FieldBase *proc_rank_field=0)
      {
        const CellTopologyData * const cell_topo_data = stk_classic::percept::PerceptMesh::get_cell_topology(element);
        typedef boost::tuple<stk_classic::mesh::EntityId, stk_classic::mesh::EntityId> line_tuple_type;
        static vector<line_tuple_type> elems(2);

        CellTopology cell_topo(cell_topo_data);
        const stk_classic::mesh::PairIterRelation elem_nodes = element.relations(stk_classic::mesh::fem::FEMMetaData::NODE_RANK);

        std::vector<stk_classic::mesh::Part*> add_parts;
        std::vector<stk_classic::mesh::Part*> remove_parts;

        add_parts = m_toParts;
        
        double coord_x[3];
        for (int iedge = 0; iedge < 1; iedge++)
          {
            //double * mp = midPoint(EDGE_COORD(iedge,0), EDGE_COORD(iedge,1), eMesh.get_spatial_dim(), coord_x);
            //double * mp = midPoint(FACE_COORD(iedge,0), FACE_COORD(iedge,1), eMesh.get_spatial_dim(), coord_x);
            double * mp = midPoint(VERT_COORD(0), VERT_COORD(1), eMesh.get_spatial_dim(), coord_x);

            if (!EDGE_N(iedge))
              {
                std::cout << "P[" << eMesh.get_rank() << " nid ## = 0  " << std::endl;
              }

            eMesh.createOrGetNode(EDGE_N(iedge), mp);
          }

        // FIXME
        nodeRegistry.makeCentroidCoords(*const_cast<stk_classic::mesh::Entity *>(&element), m_primaryEntityRank, 0u);
        nodeRegistry.addToExistingParts(*const_cast<stk_classic::mesh::Entity *>(&element), m_primaryEntityRank, 0u);

        nodeRegistry.interpolateFields(*const_cast<stk_classic::mesh::Entity *>(&element), m_primaryEntityRank, 0u);

        Elem::CellTopology elem_celltopo = Elem::getCellTopology< FromTopology >();
        const Elem::RefinementTopology* ref_topo_p = Elem::getRefinementTopology(elem_celltopo);
        const Elem::RefinementTopology& ref_topo = *ref_topo_p;

#ifndef NDEBUG
        unsigned num_child = ref_topo.num_child();
        VERIFY_OP(num_child, == , 2, "createNewElements num_child problem");
        bool homogeneous_child = ref_topo.homogeneous_child();
        VERIFY_OP(homogeneous_child, ==, true, "createNewElements homogeneous_child");
#endif

        // new_sub_entity_nodes[i][j]
        //const UInt * const * child_nodes() const {
        //const UInt * child_node_0 = ref_topo.child_node(0);

        typedef Elem::StdMeshObjTopologies::RefTopoX RefTopoX;
        RefTopoX& l2 = Elem::StdMeshObjTopologies::RefinementTopologyExtra< FromTopology > ::refinement_topology;

#define CENTROID_N NN(m_primaryEntityRank,0)  

        for (unsigned iChild = 0; iChild < 2; iChild++)
          {
            unsigned EN[2];
            for (unsigned jNode = 0; jNode < 2; jNode++)
              {
                unsigned childNodeIdx = ref_topo.child_node(iChild)[jNode];
#ifndef NDEBUG
                unsigned childNodeIdxCheck = l2[childNodeIdx].ordinal_of_node;
                VERIFY_OP(childNodeIdx, ==, childNodeIdxCheck, "childNodeIdxCheck");
#endif
                unsigned inode=0;

                if (l2[childNodeIdx].rank_of_subcell == 0)
                  inode = VERT_N(l2[childNodeIdx].ordinal_of_subcell);
                else if (l2[childNodeIdx].rank_of_subcell == 1)
                  inode = EDGE_N(l2[childNodeIdx].ordinal_of_subcell);

                //                 else if (l2[childNodeIdx].rank_of_subcell == 2)
                //                   inode = CENTROID_N;

                EN[jNode] = inode;
              }
            elems[iChild] = line_tuple_type(EN[0], EN[1]);
          }

#undef CENTROID_N

        // write a diagram of the refinement pattern as a vtk file, or a latex/tikz/pgf file
#define WRITE_DIAGRAM 0
#if WRITE_DIAGRAM

#endif

        for (unsigned ielem=0; ielem < elems.size(); ielem++)
          {
            stk_classic::mesh::Entity& newElement = *(*element_pool);

#if 0
            if (proc_rank_field && proc_rank_field->rank() == m_eMesh.edge_rank()) //&& m_eMesh.get_spatial_dim()==1)
              {
                double *fdata = stk_classic::mesh::field_data( *static_cast<const ScalarFieldType *>(proc_rank_field) , newElement );
                //fdata[0] = double(m_eMesh.get_rank());
                fdata[0] = double(newElement.owner_rank());
              }
#endif
            stk_classic::mesh::FieldBase * proc_rank_field_edge = m_eMesh.get_field("proc_rank_edge");
            if (proc_rank_field_edge)
              {
                double *fdata = stk_classic::mesh::field_data( *static_cast<const ScalarFieldType *>(proc_rank_field_edge) , newElement );
                fdata[0] = double(newElement.owner_rank());
                //fdata[0] = 1234.56;
                if (0)
                std::cout << "P[" << m_eMesh.get_rank() << "] tmp set proc_rank_field_edge to value = " << newElement.owner_rank() 
                          << " for side element = " << newElement.identifier()
                          << std::endl;
              }

            //eMesh.get_bulk_data()->change_entity_parts( newElement, add_parts, remove_parts );

            change_entity_parts(eMesh, element, newElement);

            {
              if (!elems[ielem].get<0>())
                {
                  std::cout << "P[" << eMesh.get_rank() << " nid = 0  " << std::endl;
                  exit(1);
                }

            }

            eMesh.get_bulk_data()->declare_relation(newElement, eMesh.createOrGetNode(elems[ielem].get<0>()), 0);
            eMesh.get_bulk_data()->declare_relation(newElement, eMesh.createOrGetNode(elems[ielem].get<1>()), 1);

            set_parent_child_relations(eMesh, element, newElement, ielem);


            element_pool++;

          }

      }
      
    };

  }
}
#endif
