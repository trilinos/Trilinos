// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef adapt_UniformRefinerPattern_Quad4_Tri3_4_hpp
#define adapt_UniformRefinerPattern_Quad4_Tri3_4_hpp


  namespace percept {

    struct Specialization {};

    template <>
    class UniformRefinerPattern<shards::Quadrilateral<4>, shards::Triangle<3>, 4, Specialization > : public URP<shards::Quadrilateral<4> , shards::Triangle<3> >
    {
    public:

      virtual bool edgeMarkIsEnough() override { return false; }
      UniformRefinerPattern(percept::PerceptMesh& eMesh, BlockNamesType block_names = BlockNamesType()) : URP<shards::Quadrilateral<4> , shards::Triangle<3> >(eMesh)
       {
         EXCEPTWATCH;
         m_primaryEntityRank = eMesh.face_rank();
         if (m_eMesh.get_spatial_dim() == 2)
           m_primaryEntityRank = stk::topology::ELEMENT_RANK;

         setNeededParts(eMesh, block_names, false);

       }

      virtual void doBreak() override {}

      void fillNeededEntities(std::vector<NeededEntityType>& needed_entities) override
      {
        needed_entities.resize(1);
        needed_entities[0].first = (m_eMesh.get_spatial_dim() == 2 ? stk::topology::ELEMENT_RANK :  m_eMesh.face_rank());
        setToOne(needed_entities);
      }

      virtual unsigned getNumNewElemPerElem() override { return 4; }

      virtual StringStringMap fixSurfaceAndEdgeSetNamesMap() override
      {
        StringStringMap str_map;
        str_map["hex8"] = "tet4";
        str_map["quad4"] = "tri3";
        return str_map;
      }

      void
      createNewElements(percept::PerceptMesh& eMesh, NodeRegistry& nodeRegistry,
                        stk::mesh::Entity element,  NewSubEntityNodesType& new_sub_entity_nodes, vector<stk::mesh::Entity>::iterator& element_pool,
                        vector<stk::mesh::Entity>::iterator& ft_element_pool,
                        stk::mesh::FieldBase *proc_rank_field=0) override
      {
        const CellTopologyData * const cell_topo_data = m_eMesh.get_cell_topology(element);
        typedef std::array<stk::mesh::EntityId, 3> tri_tuple_type;
        static vector<tri_tuple_type> elems(4);

        CellTopology cell_topo(cell_topo_data);
        const percept::MyPairIterRelation elem_nodes (m_eMesh, element, stk::topology::NODE_RANK);

        //stk::mesh::Part & active = mesh->ActivePart();
        //stk::mesh::Part & quad4  = mesh->QuadPart();

        std::vector<stk::mesh::Part*> add_parts;
        std::vector<stk::mesh::Part*> remove_parts;

        add_parts = m_toParts;

        //std::cout << "P["<< m_eMesh.get_rank() << "] add_parts = " << add_parts << std::endl;

        stk::mesh::EntityRank my_rank = m_primaryEntityRank;

        nodeRegistry.prolongateCoords(element, my_rank, 0u);
        nodeRegistry.addToExistingParts(element, my_rank, 0u);
        nodeRegistry.prolongateFields(element, my_rank, 0u);

#define CENTROID_N NN(m_primaryEntityRank, 0)

        elems[0] = {VERT_N(0), VERT_N(1), CENTROID_N};
        elems[1] = {VERT_N(1), VERT_N(2), CENTROID_N};
        elems[2] = {VERT_N(2), VERT_N(3), CENTROID_N};
        elems[3] = {VERT_N(3), VERT_N(0), CENTROID_N};

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

        bool use_declare_element_side = UniformRefinerPatternBase::USE_DECLARE_ELEMENT_SIDE &&  m_primaryEntityRank == eMesh.side_rank();

        for (unsigned ielem=0; ielem < elems.size(); ielem++)
          {
            stk::mesh::Entity newElement = stk::mesh::Entity();
            if (!use_declare_element_side)
              newElement = *element_pool;

            stk::mesh::Entity nodes[3] = {
              eMesh.createOrGetNode(elems[ielem][0]),
              eMesh.createOrGetNode(elems[ielem][1]),
              eMesh.createOrGetNode(elems[ielem][2])};

            create_side_element(eMesh, use_declare_element_side, nodes, 3, newElement);

            //std::cout << "P["<< m_eMesh.get_rank() << "] urp tmp 3 "  << proc_rank_field << std::endl;
            if (proc_rank_field && m_eMesh.entity_rank(element) == stk::topology::ELEMENT_RANK)
              {
                double *fdata = stk::mesh::field_data( *static_cast<const ScalarFieldType *>(proc_rank_field) , newElement );
                fdata[0] = double(eMesh.owner_rank(newElement));
              }

            //std::cout << "P["<< m_eMesh.get_rank() << "] urp tmp 4 "  << std::endl;
            change_entity_parts(eMesh, element, newElement);

            //std::cout << "P["<< m_eMesh.get_rank() << "] urp tmp 7 "  << std::endl;
            set_parent_child_relations(eMesh, element, newElement, *ft_element_pool, ielem);

            ft_element_pool++;
            if (!use_declare_element_side)
              element_pool++;

          }

      }

    };

  }

#endif
