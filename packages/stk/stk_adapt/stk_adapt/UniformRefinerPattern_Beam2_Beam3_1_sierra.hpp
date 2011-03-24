#ifndef stk_adapt_UniformRefinerPattern_Beam2_Beam3_1_sierra_hpp
#define stk_adapt_UniformRefinerPattern_Beam2_Beam3_1_sierra_hpp


//#include "UniformRefinerPattern.hpp"
#include <stk_adapt/sierra_element/RefinementTopology.hpp>
#include <stk_adapt/sierra_element/StdMeshObjTopologies.hpp>

#include <boost/array.hpp>


namespace stk {
  namespace adapt {

    template <>
    class UniformRefinerPattern< shards::Beam<2>, shards::Beam<3>, 1, SierraPort > : public URP<shards::Beam<2> , shards::Beam<3> >
    {

    public:

      UniformRefinerPattern(percept::PerceptMesh& eMesh, BlockNamesType block_names = BlockNamesType()) : URP<shards::Beam<2> , shards::Beam<3> >(eMesh)
      {
        //m_primaryEntityRank = mesh::Edge;
        //if (m_eMesh.getSpatialDim() == 1)
        m_primaryEntityRank = mesh::Element;

        setNeededParts(eMesh, block_names, false);
        Elem::StdMeshObjTopologies::bootstrap();

      }


      virtual void doBreak() {}

      void setSubPatterns( std::vector<UniformRefinerPatternBase *>& bp, percept::PerceptMesh& eMesh )
      {
        EXCEPTWATCH;
        bp = std::vector<UniformRefinerPatternBase *>(1u, 0);
        bp[0] = this;
      }

      void fillNeededEntities(std::vector<NeededEntityType>& needed_entities)
      {
        needed_entities.resize(1);
        needed_entities[0].first = stk::mesh::Element;
        setToOne(needed_entities);
      }

      virtual unsigned getNumNewElemPerElem() { return 1; }

      void 
      createNewElements(percept::PerceptMesh& eMesh, NodeRegistry& nodeRegistry, 
                        Entity& element,  NewSubEntityNodesType& new_sub_entity_nodes, vector<Entity *>::iterator& element_pool,
                        FieldBase *proc_rank_field=0)
      {
        const CellTopologyData * const cell_topo_data = stk::percept::PerceptMesh::get_cell_topology(element);
        typedef boost::array<unsigned,3> quadratic_type;
        static vector<quadratic_type> elems(1);

        CellTopology cell_topo(cell_topo_data);
        const mesh::PairIterRelation elem_nodes = element.relations(Node);

        std::vector<stk::mesh::Part*> add_parts;
        std::vector<stk::mesh::Part*> remove_parts;

        add_parts = m_toParts;
        nodeRegistry.addToExistingParts(*const_cast<Entity *>(&element), m_primaryEntityRank, 0u);

#define CENTROID_N NN(m_primaryEntityRank,0)  

        {
          quadratic_type& EN = elems[0];

          for (unsigned ind = 0; ind < 2; ind++)
            {
              unsigned inode = VERT_N(ind);
              EN[ind] = inode;
            }

          EN[2] = CENTROID_N;
        }

#undef CENTROID_N

        for (unsigned ielem=0; ielem < elems.size(); ielem++)
          {
            Entity& newElement = *(*element_pool);

            // FIXME
            if (0 && proc_rank_field)
              {
                double *fdata = stk::mesh::field_data( *static_cast<const ScalarFieldType *>(proc_rank_field) , newElement );
                fdata[0] = double(newElement.owner_rank());
              }

            change_entity_parts(eMesh, element, newElement);

            set_parent_child_relations(eMesh, element, newElement, ielem);

            {
              if (!elems[ielem][0])
                {
                  std::cout << "P[" << eMesh.getRank() << " nid = 0  " << std::endl;
                  exit(1);
                }

            }
            for (int inode=0; inode < 3; inode++)
              {
                mesh::EntityId eid = elems[ielem][inode];
                mesh::Entity& node = eMesh.createOrGetNode(eid);
                eMesh.get_bulkData()->declare_relation(newElement, node, inode);
              }

            element_pool++;

          }

      }
      
    };

  }
}
#endif
