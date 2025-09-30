// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef adapt_UniformRefinerPattern_Line2_Line2_2_sierra_hpp
#define adapt_UniformRefinerPattern_Line2_Line2_2_sierra_hpp


//#include "UniformRefinerPattern.hpp"
#include <adapt/sierra_element/RefinementTopology.hpp>
#include <adapt/sierra_element/StdMeshObjTopologies.hpp>

  namespace percept {

    template <>
    class UniformRefinerPattern<shards::Line<2>, shards::Line<2>, 2, SierraPort > : public URP<shards::Line<2>, shards::Line<2>  >
    {
    public:

      UniformRefinerPattern(percept::PerceptMesh& eMesh, BlockNamesType block_names = BlockNamesType()) :  URP<shards::Line<2>, shards::Line<2>  >(eMesh)
      {
        m_primaryEntityRank = m_eMesh.edge_rank();
        if (m_eMesh.get_spatial_dim() == 1)
          m_primaryEntityRank = stk::topology::ELEMENT_RANK;

        setNeededParts(eMesh, block_names, true);
        Elem::StdMeshObjTopologies::bootstrap();

      }

      virtual void doBreak() override {}
      void fillNeededEntities(std::vector<NeededEntityType>& needed_entities) override
      {
        needed_entities.resize(1);
        needed_entities[0].first = m_eMesh.edge_rank();
        setToOne(needed_entities);
      }

      virtual unsigned getNumNewElemPerElem() override { return 2; }

      void
      createNewElements(percept::PerceptMesh& eMesh, NodeRegistry& nodeRegistry,
                        stk::mesh::Entity element,  NewSubEntityNodesType& new_sub_entity_nodes, vector<stk::mesh::Entity>::iterator& element_pool,
                        vector<stk::mesh::Entity>::iterator& ft_element_pool,
                        stk::mesh::FieldBase */*proc_rank_field*/=0) override
      {
        const CellTopologyData * const cell_topo_data = m_eMesh.get_cell_topology(element);
        typedef std::array<stk::mesh::EntityId, 2> line_tuple_type;
        static vector<line_tuple_type> elems(2);

        CellTopology cell_topo(cell_topo_data);
        const percept::MyPairIterRelation elem_nodes (m_eMesh, element, stk::topology::NODE_RANK);

        std::vector<stk::mesh::Part*> add_parts;
        std::vector<stk::mesh::Part*> remove_parts;

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
        nodeRegistry.prolongateCoords(element, m_primaryEntityRank, 0u);
        nodeRegistry.addToExistingParts(element, m_primaryEntityRank, 0u);

        nodeRegistry.prolongateFields(element, m_primaryEntityRank, 0u);

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

        RefTopoX& l2 = Elem::StdMeshObjTopologies::RefinementTopologyExtra< FromTopology > ::refinement_topology;

#define CENTROID_N NN(m_primaryEntityRank,0)

        for (unsigned iChild = 0; iChild < 2; iChild++)
          {
            stk::mesh::EntityId EN[2];
            for (unsigned jNode = 0; jNode < 2; jNode++)
              {
                unsigned childNodeIdx = ref_topo.child_node(iChild)[jNode];
#ifndef NDEBUG
                unsigned childNodeIdxCheck = l2[childNodeIdx].ordinal_of_node;
                VERIFY_OP(childNodeIdx, ==, childNodeIdxCheck, "childNodeIdxCheck");
#endif
                stk::mesh::EntityId inode=0;

                if (l2[childNodeIdx].rank_of_subcell == 0)
                  inode = VERT_N(l2[childNodeIdx].ordinal_of_subcell);
                else if (l2[childNodeIdx].rank_of_subcell == 1)
                  inode = EDGE_N(l2[childNodeIdx].ordinal_of_subcell);

                //                 else if (l2[childNodeIdx].rank_of_subcell == 2)
                //                   inode = CENTROID_N;

                EN[jNode] = inode;
              }
            elems[iChild] = {EN[0], EN[1]};
          }

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

            stk::mesh::Entity nodes[2] = {eMesh.createOrGetNode(elems[ielem][0]), eMesh.createOrGetNode(elems[ielem][1])};
            create_side_element(eMesh, use_declare_element_side, nodes, 2, newElement);

            change_entity_parts(eMesh, element, newElement);

            set_parent_child_relations(eMesh, element, newElement, *ft_element_pool, ielem);

            ft_element_pool++;
            if (!use_declare_element_side)
              element_pool++;

          }

      }

    };

  }

#endif
