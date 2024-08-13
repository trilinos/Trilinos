// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef adapt_UniformRefinerPattern_Tet4_Tet4_8_sierra_hpp
#define adapt_UniformRefinerPattern_Tet4_Tet4_8_sierra_hpp


//#include "UniformRefinerPattern.hpp"
#include <adapt/sierra_element/RefinementTopology.hpp>
#include <adapt/sierra_element/StdMeshObjTopologies.hpp>
#include <adapt/CompareCoordinates.hpp>
#include <adapt/Percept_MOAB_SimplexTemplateRefiner.hpp>

#define USE_FACE_BREAKER_T4_T4_8 1
#if USE_FACE_BREAKER_T4_T4_8
#include "UniformRefinerPattern_Tri3_Tri3_4_sierra.hpp"
#endif

#define USE_PERCEPT_MOAB_TET_REFINE 1

  namespace percept {

    template <>
    class UniformRefinerPattern<shards::Tetrahedron<4>, shards::Tetrahedron<4>, 8, SierraPort > : public URP<shards::Tetrahedron<4>, shards::Tetrahedron<4>  >
    {

#if USE_FACE_BREAKER_T4_T4_8
      // FIXME
      UniformRefinerPattern<shards::Triangle<3>, shards::Triangle<3>, 4, SierraPort > * m_face_breaker;
#endif

    public:

      UniformRefinerPattern(percept::PerceptMesh& eMesh, BlockNamesType block_names = BlockNamesType()) :  URP<shards::Tetrahedron<4>, shards::Tetrahedron<4>  >(eMesh)
      {
        m_primaryEntityRank = stk::topology::ELEMENT_RANK;

        setNeededParts(eMesh, block_names, true);
        Elem::StdMeshObjTopologies::bootstrap();

#if USE_FACE_BREAKER_T4_T4_8

        //m_face_breaker = Teuchos::rcp( new UniformRefinerPattern<shards::Line<2>, shards::Line<2>, 2, SierraPort > (eMesh, block_names) );
        m_face_breaker =  new UniformRefinerPattern<shards::Triangle<3>, shards::Triangle<3>, 4, SierraPort > (eMesh, block_names) ;
#endif

      }

      ~UniformRefinerPattern()
      {
        if (m_face_breaker) delete m_face_breaker;
      }

      void setSubPatterns( std::vector<UniformRefinerPatternBase *>& bp, percept::PerceptMesh& eMesh )
      {
        EXCEPTWATCH;
        bp.resize(2);

        if (eMesh.get_spatial_dim() == 3)
          {
            bp[0] = this;
#if USE_FACE_BREAKER_T4_T4_8
            bp[1] = m_face_breaker;
#endif
          }
        else if (eMesh.get_spatial_dim() == 2)
          {
            // FIXME
            std::cout << "ERROR" ;
            exit(1);
          }

      }

      virtual void doBreak() {}
      void fillNeededEntities(std::vector<NeededEntityType>& needed_entities)
      {
        needed_entities.resize(1);
        needed_entities[0].first = m_eMesh.edge_rank();
        setToOne(needed_entities);
      }

      virtual unsigned getNumNewElemPerElem() { return 8; }

#if USE_PERCEPT_MOAB_TET_REFINE
      typedef std::array<unsigned, 4> TetTupleTypeLocal;
      typedef std::array<stk::mesh::EntityId, 4> TetTupleType;

#define TET_VERT_N(i) (i)
#define TET_EDGE_N(i) ((i)+4)
#define TET_TRI_CV_EV(iface,i) ( i < 3 ? tbl_tet_face_nodes[iface][i] : tbl_tet_face_edge_map[iface][i-3] )
#define TET_CENTROID_NODE (4+6)
#define TET_CV_EV(i) ( i < 4 ? VERT_N(i) : EDGE_N(i-4) )

      static void triangulate_tet(PerceptMesh& eMesh, stk::mesh::Entity tet_elem_nodes[4], unsigned edge_marks[6],
                                  std::vector<TetTupleTypeLocal>& tets)
      {
//         static const unsigned tbl_tet_face_edge_map[4][3]  = { {0, 4, 3}, {1, 5, 4}, {3, 5, 2}, {2, 1, 0} };
//         static const unsigned tbl_tet_face_nodes[4][3]     = { {0, 1, 3}, {1, 2, 3}, {0, 3, 2}, {0, 2, 1} };
//         static const unsigned tbl_tet_edge_nodes[6][2]     = { {0, 1}, {1, 2}, {2, 0}, {0, 3}, {1, 3}, {2, 3} };

        const CellTopologyData * const cell_topo_data = shards::getCellTopologyData< shards::Tetrahedron<4> >();

        shards::CellTopology cell_topo(cell_topo_data);
        //CoordinatesFieldType* coordField = eMesh.get_coordinates_field();

        //unsigned num_edges_marked=6;

        // uniform refinement
        tets.resize(8);

        std::array<double *,4> node_coords;
        for (int inode=0; inode < 4; inode++)
          {
            node_coords[inode] = static_cast<double*>(stk::mesh::field_data( *eMesh.get_coordinates_field() , tet_elem_nodes[inode] ));
            if (0) std::cout << "tmp RP node_coords= "
                             << node_coords[inode][0] << " "
                             << node_coords[inode][1] << " "
                             << node_coords[inode][2] << std::endl;
          }

        const std::array<int,4> node_rank = get_rank_of_nodes_based_on_coordinates(node_coords);
        std::vector<moab::TetTupleInt> new_tets;
        bool choose_best_tets = true;
        moab::SimplexTemplateRefiner str(choose_best_tets);
        str.refine_3_simplex(new_tets, edge_marks, 1,
                             node_coords[0], 0, node_rank[0],
                             node_coords[1], 0, node_rank[1],
                             node_coords[2], 0, node_rank[2],
                             node_coords[3], 0, node_rank[3] );

        if (0)
          {
            for (int inode=0; inode < 4; inode++)
              {
                node_coords[inode] = static_cast<double*>(stk::mesh::field_data( *eMesh.get_coordinates_field() , tet_elem_nodes[inode] ));
                std::cout << "tmp RefPatt::createNewElements node_coords after= "
                          << node_coords[inode][0] << " "
                          << node_coords[inode][1] << " "
                          << node_coords[inode][2] << std::endl;
              }
          }

        tets.resize(new_tets.size());
        for (unsigned i = 0; i < new_tets.size(); i++)
          {
            tets[i] = {(unsigned)new_tets[i][0],
                       (unsigned)new_tets[i][1],
                       (unsigned)new_tets[i][2],
                       (unsigned)new_tets[i][3]};
            if (0)
              std::cout << "tmp RefPatt::createNewElements new tet= " << tets[i] << std::endl;

          }
      }

      void
      createNewElements(percept::PerceptMesh& eMesh, NodeRegistry& nodeRegistry,
                        stk::mesh::Entity element,  NewSubEntityNodesType& new_sub_entity_nodes, std::vector<stk::mesh::Entity>::iterator& element_pool,
                                        vector<stk::mesh::Entity>::iterator& ft_element_pool,
                        stk::mesh::FieldBase *proc_rank_field=0)
      {
        const CellTopologyData * const cell_topo_data = eMesh.get_cell_topology(element);
        static std::vector<TetTupleType> elems(8);
        static std::vector<TetTupleTypeLocal> elems_local(8);
        unsigned num_new_elems=0;

        shards::CellTopology cell_topo(cell_topo_data);
        const percept::MyPairIterRelation elem_nodes (eMesh, element, stk::topology::NODE_RANK);

        std::vector<stk::mesh::Part*> add_parts;
        std::vector<stk::mesh::Part*> remove_parts;
        add_parts = m_toParts;

        unsigned edge_marks[6] = {0,0,0,0,0,0};
        for (int iedge = 0; iedge < 6; iedge++)
          {
            edge_marks[iedge] = 1;
          }

        stk::mesh::Entity elem_nodes_local[4] = {stk::mesh::Entity()};
        for (int inode=0; inode < 4; inode++)
          {
            elem_nodes_local[inode] = elem_nodes[inode].entity();
          }
        triangulate_tet(eMesh, elem_nodes_local, edge_marks, elems_local);

        num_new_elems = elems_local.size();
        elems.resize(num_new_elems);
        for (unsigned ielem=0; ielem < num_new_elems; ielem++)
          {
            elems[ielem] = { TET_CV_EV(elems_local[ielem][0] ),
                             TET_CV_EV(elems_local[ielem][1] ),
                             TET_CV_EV(elems_local[ielem][2] ),
                             TET_CV_EV(elems_local[ielem][3] ) };
            if (0)
              std::cout << "tmp RefPatt::createNewElements new tet= " << elems[ielem] << std::endl;

          }

        //std::cout << "tmp RefinerPattern_Tet4_Tet4_N::num_edges_marked= " << num_edges_marked << std::endl;

        //nodeRegistry.prolongateCoords(*const_cast<stk::mesh::Entity>(&element), stk::topology::ELEMENT_RANK, 0u);

        for (unsigned ielem=0; ielem < elems.size(); ielem++)
          {
            stk::mesh::Entity newElement = *element_pool;

            if (proc_rank_field)
              {
                double *fdata = stk::mesh::field_data( *static_cast<const ScalarFieldType *>(proc_rank_field) , newElement );
                //fdata[0] = double(eMesh.get_rank());
                fdata[0] = double(eMesh.owner_rank(newElement));
              }

            //eMesh.get_bulk_data()->change_entity_parts( newElement, add_parts, remove_parts );
            change_entity_parts(eMesh, element, newElement);

            // 4 nodes of the new tets
            eMesh.get_bulk_data()->declare_relation(newElement, {eMesh.createOrGetNode(elems[ielem][0]),
                                                                 eMesh.createOrGetNode(elems[ielem][1]),
                                                                 eMesh.createOrGetNode(elems[ielem][2]),
                                                                 eMesh.createOrGetNode(elems[ielem][3])});

            set_parent_child_relations(eMesh, element, newElement, *ft_element_pool, ielem);

            std::vector<stk::mesh::Entity> elements(1,element);
            eMesh.prolongateElementFields( elements, newElement);

            if (0)
              {
                std::cout << "tmp RefPatt::createNewElements eMesh.identifier(element)= " << eMesh.identifier(element)
                          << " newElement= " << eMesh.identifier(newElement) << std::endl;

              }

            ft_element_pool++;
            element_pool++;

          }


      }

#undef TET_VERT_N
#undef TET_EDGE_N
#undef TET_TRI_CV_EV
#undef TET_CENTROID_NODE
#undef TET_CV_EV

#else
      void
      createNewElements(percept::PerceptMesh& eMesh, NodeRegistry& nodeRegistry,
                        stk::mesh::Entity element,  NewSubEntityNodesType& new_sub_entity_nodes, vector<stk::mesh::Entity>::iterator& element_pool,
                        vector<stk::mesh::Entity>::iterator& ft_element_pool,
                        stk::mesh::FieldBase *proc_rank_field=0)
      {
        genericRefine_createNewElements(eMesh, nodeRegistry,
                                        element, new_sub_entity_nodes, element_pool, ft_element_pool,
                                        proc_rank_field);
      }
#endif

    };

  }

#endif
