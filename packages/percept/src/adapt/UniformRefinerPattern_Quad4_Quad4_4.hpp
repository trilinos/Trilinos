// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef adapt_UniformRefinerPattern_Quad4_Quad4_4_hpp
#define adapt_UniformRefinerPattern_Quad4_Quad4_4_hpp


//#include "UniformRefinerPattern.hpp"

  namespace percept {

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
          m_primaryEntityRank = stk::topology::ELEMENT_RANK;

        setNeededParts(eMesh, block_names, true);
      }

      virtual void doBreak() override {}
      void fillNeededEntities(std::vector<NeededEntityType>& needed_entities) override
      {
        needed_entities.resize(2);
        needed_entities[0].first = m_eMesh.edge_rank();    // edges have 2 nodes
        needed_entities[1].first = stk::topology::ELEMENT_RANK;
        setToOne(needed_entities);
      }

      virtual unsigned getNumNewElemPerElem() override { return 4; }

      void
      createNewElements(percept::PerceptMesh& eMesh, NodeRegistry& nodeRegistry,
                        stk::mesh::Entity element,  NewSubEntityNodesType& new_sub_entity_nodes, vector<stk::mesh::Entity>::iterator& element_pool,
                        vector<stk::mesh::Entity>::iterator& ft_element_pool,
                        stk::mesh::FieldBase *proc_rank_field=0) override
      {
        const CellTopologyData * const cell_topo_data = m_eMesh.get_cell_topology(element);
        typedef std::array<stk::mesh::EntityId, 4> quad_tuple_type;
        static vector<quad_tuple_type> elems(4);

        CellTopology cell_topo(cell_topo_data);
        const percept::MyPairIterRelation elem_nodes (m_eMesh, element,stk::topology::NODE_RANK);

        //stk::mesh::Part & active = mesh->ActivePart();
        //stk::mesh::Part & quad4  = mesh->QuadPart();

        std::vector<stk::mesh::Part*> add_parts;
        std::vector<stk::mesh::Part*> remove_parts;

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

        nodeRegistry.prolongateCoords(element, stk::topology::ELEMENT_RANK, 0u);


// new_sub_entity_nodes[i][j]
#define CENTROID_N NN(m_primaryEntityRank,0)

        elems[0] = {VERT_N(0), EDGE_N(0), CENTROID_N, EDGE_N(3)};
        elems[1] = {VERT_N(1), EDGE_N(1), CENTROID_N, EDGE_N(0)};
        elems[2] = {VERT_N(2), EDGE_N(2), CENTROID_N, EDGE_N(1)};
        elems[3] = {VERT_N(3), EDGE_N(3), CENTROID_N, EDGE_N(2)};

#undef CENTROID_N

        // write a diagram of the refinement pattern as a vtk file, or a latex/tikz/pgf file
#define WRITE_DIAGRAM 0
#if WRITE_DIAGRAM

#endif
        bool use_declare_element_side = UniformRefinerPatternBase::USE_DECLARE_ELEMENT_SIDE &&  m_primaryEntityRank == eMesh.side_rank();

        for (unsigned ielem=0; ielem < elems.size(); ielem++)
          {
            stk::mesh::Entity newElement = stk::mesh::Entity();
            if (!use_declare_element_side)
              newElement = *element_pool;

            stk::mesh::Entity nodes[4] = {
              eMesh.createOrGetNode(elems[ielem][0]),
              eMesh.createOrGetNode(elems[ielem][1]),
              eMesh.createOrGetNode(elems[ielem][2]),
              eMesh.createOrGetNode(elems[ielem][3]) };

            create_side_element(eMesh, use_declare_element_side, nodes, 4, newElement);

            if (proc_rank_field)
              {
                double *fdata = stk::mesh::field_data( *static_cast<const ScalarFieldType *>(proc_rank_field) , newElement );
                //fdata[0] = double(m_eMesh.get_rank());
                fdata[0] = double(eMesh.owner_rank(newElement));
              }

            eMesh.get_bulk_data()->change_entity_parts( newElement, add_parts, remove_parts );

            set_parent_child_relations(eMesh, element, newElement, *ft_element_pool, ielem);

            ft_element_pool++;
            if (!use_declare_element_side)
              element_pool++;

          }

      }

    };

  }

#endif
