// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef adapt_RefinerPattern_Tet4_Tet4_N_sierra_hpp
#define adapt_RefinerPattern_Tet4_Tet4_N_sierra_hpp


//#include "UniformRefinerPattern.hpp"
#include <adapt/sierra_element/RefinementTopology.hpp>
#include <adapt/sierra_element/StdMeshObjTopologies.hpp>
#include <adapt/CompareCoordinates.hpp>
#include <adapt/RefinerPattern_Tri3_Tri3_N.hpp>

#include <adapt/Percept_MOAB_SimplexTemplateRefiner.hpp>

//#include "UniformRefinerPattern_Line2_Line2_2_sierra.hpp"
#include "RefinerPattern_Line2_Line2_N.hpp"

/** NOTE: A lot of the following code is unfinished (greedy triangulation scheme).
 *  The call to triangulate_tet_generic has been replaced with a call to the
 *  MOAB SimplexTemplateRefiner as modified for Percept/Adapt.
 *
 *  We are leaving this unfinished code in place until the MOAB testing is verified.
 */

#define PERCEPT_USE_MOAB_REFINER 1

  namespace percept {

    /*---------------------------------------------------------------------*/
    /** From Shards_BasicTopologies.hpp
     *
     * typedef MakeTypeList< IndexList< 0 , 1 , 4 > ,
     *                       IndexList< 1 , 2 , 5 > ,
     *                       IndexList< 2 , 0 , 6 > ,
     *                       IndexList< 0 , 3 , 7 > ,
     *                       IndexList< 1 , 3 , 8 > ,
     *                       IndexList< 2 , 3 , 9 > >::type
     *   TetrahedronEdgeNodeMap ;
     *
     * typedef MakeTypeList< IndexList< 0 , 1 , 3 ,   4 , 8 , 7 > ,
     *                       IndexList< 1 , 2 , 3 ,   5 , 9 , 8 > ,
     *                       IndexList< 0 , 3 , 2 ,   7 , 9 , 6 > ,
     *                       IndexList< 0 , 2 , 1 ,   6 , 5 , 4 > >::type
     *   TetrahedronSideNodeMap ;
     *
     *
     */

    /*---------------------------------------------------------------------*/
    /**
     *                         PARENT 4-Node Tetrahedron Object Nodes
     *              3
     *              o
     *             /|\
     *            / | \       (PARENT) 4-Node Tetrahedron Object
     *           /  |  \               Edge Node Map:
     *          /   |   \         0      1       2       3       4       5
     *         /    |    \    { {0, 1}, {1, 2}, {2, 0}, {0, 3}, {1, 3}, {2, 3} };
     *      0 o-----|-----o 2
     *         \    |    /
     *          \   |   /         Face Node Map:
     *           \  |  /
     *            \ | /       { {0, 1, 3}, {1, 2, 3}, {0, 3, 2}, {0, 2, 1}  }
     *             \|/
     *              o
     *              1
     *
     *              3
     *              o
     *             /|\
     *            / | \
     *         7 *  |  * 9
     *          /   |   \
     *         /   6|    \
     *      0 o----*|-----o 2
     *         \    *8   /
     *          \   |   /
     *         4 *  |  * 5
     *            \ | /
     *             \|/
     *              o
     *              1
     */
    /*---------------------------------------------------------------------*/
    /** Edge (j) of face (i) maps to edge (k) of tet:  ==> tet_face_edge_map[4][3] == tet_face_edge_map[i][j]
     *
     *  face 0 {0, 1, 3}:  {0, 4, 3}
     *  face 1 {1, 2, 3}:  {1, 5, 4}
     *  face 2 {0, 3, 2}:  {3, 5, 2}
     *  face 3 {0, 2, 1}:  {2, 1, 0}
     */

    //static unsigned tbl_tet_face_edge_map[4][3]  = { {0, 4, 3}, {1, 5, 4}, {3, 5, 2}, {2, 1, 0} };
    //static unsigned tbl_tet_face_nodes[4][3]     = { {0, 1, 3}, {1, 2, 3}, {0, 3, 2}, {0, 2, 1} };
    //static unsigned tbl_tet_edge_nodes[6][2]     = { {0, 1}, {1, 2}, {2, 0}, {0, 3}, {1, 3}, {2, 3} };

    typedef std::array<unsigned, 4> TetTupleTypeLocal;
    typedef std::array<stk::mesh::EntityId, 4> TetTupleType;

    /// general refinement pattern

    // the "-1" here signifies the number of elements created is not fixed, depends on the marking pattern
    template <>
    class RefinerPattern<shards::Tetrahedron<4>, shards::Tetrahedron<4>, -1 > : public URP<shards::Tetrahedron<4>,shards::Tetrahedron<4>  >
    {

      RefinerPattern<shards::Triangle<3>, shards::Triangle<3>, -1 > * m_face_breaker;
      RefinerPattern<shards::Line<2>, shards::Line<2>, -1 > * m_edge_breaker;
      TransitionElementType *m_transition_element_field;
      bool m_transition_element_field_set;

    public:

      RefinerPattern(percept::PerceptMesh& eMesh, BlockNamesType block_names = BlockNamesType()) :  URP<shards::Tetrahedron<4>, shards::Tetrahedron<4>  >(eMesh),
                                                                                                    m_face_breaker(0), m_edge_breaker(0),
                                                                                                    m_transition_element_field(0), m_transition_element_field_set(false)
      {
        m_primaryEntityRank = stk::topology::ELEMENT_RANK;

        setNeededParts(eMesh, block_names, true);
        Elem::StdMeshObjTopologies::bootstrap();

        m_face_breaker =  new RefinerPattern<shards::Triangle<3>, shards::Triangle<3>, -1 > (eMesh, block_names) ;
        m_edge_breaker =  new RefinerPattern<shards::Line<2>, shards::Line<2>, -1 > (eMesh, block_names) ;

      }

      ~RefinerPattern()
      {
        if (m_face_breaker) delete m_face_breaker;
        if (m_edge_breaker) delete m_edge_breaker;
      }

      void setSubPatterns( std::vector<UniformRefinerPatternBase *>& bp, percept::PerceptMesh& /*eMesh*/ ) override
      {
        EXCEPTWATCH;
        bp.resize(3);
        bp[0] = this;
        bp[1] = m_face_breaker;
        bp[2] = m_edge_breaker;
      }

      virtual void doBreak() override {}
      void fillNeededEntities(std::vector<NeededEntityType>& needed_entities) override
      {
        if (!m_mark_centroid_always)
          {
            needed_entities.resize(1);
            needed_entities[0].first = m_eMesh.edge_rank();
            //needed_entities[1].first = m_eMesh.face_rank();
            setToOne(needed_entities);
          }
        else
          {
            needed_entities.resize(2);
            needed_entities[0].first = m_eMesh.edge_rank();
            needed_entities[1].first = m_eMesh.face_rank();
            setToOne(needed_entities);
          }
      }

      // FIXME - for now, create more than we need (to fix this right we need a change to the Refiner.cpp interface)
      virtual unsigned getNumNewElemPerElem() override { return 8; }

      /**
       *
       *   Convention: input is the element's nodes and the marks on the 6 edges.  Output is an array
       *     of "elements" defined as local id's of nodes forming those elements, where {0,1,2,3} represent
       *     the original vertices and {4,..,9} are the edges:
       *
       *
       *              3
       *              o
       *             /|\
       *            / | \
       *         7 *  |  * 9
       *          /   |   \
       *         /   6|    \
       *      0 o----*|-----o 2
       *         \    *8   /
       *          \   |   /
       *         4 *  |  * 5
       *            \ | /
       *             \|/
       *              o
       *              1
       */



#define TET_VERT_N(i) (i)
#define TET_EDGE_N(i) ((i)+4)

      static void triangulate_tet(PerceptMesh& eMesh, stk::mesh::Entity tet_elem_nodes[4], unsigned edge_marks[6],
                                  std::vector<TetTupleTypeLocal>& tets)
      {
        tets.resize(0);
        const CellTopologyData * const cell_topo_data = shards::getCellTopologyData< shards::Tetrahedron<4> >();

        shards::CellTopology cell_topo(cell_topo_data);
        //CoordinatesFieldType* coordField = eMesh.get_coordinates_field();

        unsigned num_edges_marked=0;
        for (int iedge = 0; iedge < 6; iedge++)
          {
            unsigned num_nodes_on_edge = edge_marks[iedge];
            if (num_nodes_on_edge)
              {
                ++num_edges_marked;
              }
          }

        if (0)
          std::cout << "tmp RefinerPattern_Tet4_Tet4_N::num_edges_marked= " << num_edges_marked << std::endl;

        if (num_edges_marked == 0)
          {
            return;
          }
        // general case (now includes uniform refinement (all edges marked))
        else
          {
            std::array<double *, 4> node_coords;
            for (int inode=0; inode < 4; inode++)
              {
                node_coords[inode] = static_cast<double*>(stk::mesh::field_data( *eMesh.get_coordinates_field() , tet_elem_nodes[inode] ));
                if (0) std::cout << "tmp RP node_coords= "
                                 << node_coords[inode][0] << " "
                                 << node_coords[inode][1] << " "
                                 << node_coords[inode][2] << std::endl;
              }

            const std::array<int, 4> node_rank = get_rank_of_nodes_based_on_coordinates(node_coords);
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
              }
          }
      }

      void
      createNewElements(percept::PerceptMesh& eMesh, NodeRegistry& /*nodeRegistry*/,
                        stk::mesh::Entity element,  NewSubEntityNodesType& new_sub_entity_nodes, std::vector<stk::mesh::Entity>::iterator& element_pool,
                        vector<stk::mesh::Entity>::iterator& ft_element_pool,
                        stk::mesh::FieldBase *proc_rank_field=0) override
      {
        const CellTopologyData * const cell_topo_data = eMesh.get_cell_topology(element);
        static std::vector<TetTupleType> elems(8);
        static std::vector<TetTupleTypeLocal> elems_local(8);
        unsigned num_new_elems=0;

        shards::CellTopology cell_topo(cell_topo_data);
        const percept::MyPairIterRelation elem_nodes (eMesh, element,stk::topology::NODE_RANK);
        //CoordinatesFieldType* coordField = eMesh.get_coordinates_field();

        std::vector<stk::mesh::Part*> add_parts;
        std::vector<stk::mesh::Part*> remove_parts;
        add_parts = m_toParts;

        unsigned edge_marks[6] = {0,0,0,0,0,0};
        unsigned num_edges_marked=0;
        for (int iedge = 0; iedge < 6; iedge++)
          {
            unsigned num_nodes_on_edge = new_sub_entity_nodes[eMesh.edge_rank()][iedge].size();
            if (num_nodes_on_edge)
              {
                edge_marks[iedge] = 1;
                ++num_edges_marked;
              }
          }
        if (!m_transition_element_field_set)
          {
            m_transition_element_field_set = true;
            m_transition_element_field = eMesh.get_transition_element_field();
            if(m_transition_element_field && m_transition_element_field->name() != "transition_element_3") m_transition_element_field = nullptr;
          }

        if (num_edges_marked == 0)
          return;

        stk::mesh::Entity elem_nodes_local[4] = {stk::mesh::Entity()};
        for (int inode=0; inode < 4; inode++)
          {
            elem_nodes_local[inode] = elem_nodes[inode].entity();
          }
        triangulate_tet(eMesh, elem_nodes_local, edge_marks, elems_local);

#define TET_CV_EV(i) ( i < 4 ? VERT_N(i) : EDGE_N(i-4) )

        num_new_elems = elems_local.size();
        elems.resize(num_new_elems);
        for (unsigned ielem=0; ielem < num_new_elems; ielem++)
          {
            elems[ielem] = { TET_CV_EV(elems_local[ielem][0] ),
                             TET_CV_EV(elems_local[ielem][1] ),
                             TET_CV_EV(elems_local[ielem][2] ),
                             TET_CV_EV(elems_local[ielem][3] ) };
          }

        for (unsigned ielem=0; ielem < elems.size(); ielem++)
          {
            stk::mesh::Entity newElement = *element_pool;

            if (proc_rank_field)
              {
                double *fdata = static_cast<double*>(stk::mesh::field_data( *proc_rank_field, newElement ));
                //fdata[0] = double(eMesh.get_rank());
                fdata[0] = double(eMesh.owner_rank(newElement));
              }

            if (m_transition_element_field)
              {
                int *transition_element = 0;
                transition_element = stk::mesh::field_data( *m_transition_element_field , newElement );
                if (num_edges_marked == 6)
                  {
                    transition_element[0] = 0;
                  }
                else
                  {
                    transition_element[0] = 1;
                  }
              }

            change_entity_parts(eMesh, element, newElement);

            // 4 nodes of the new tets
            eMesh.get_bulk_data()->declare_relation(newElement, eMesh.createOrGetNode(elems[ielem][0]), 0);
            eMesh.get_bulk_data()->declare_relation(newElement, eMesh.createOrGetNode(elems[ielem][1]), 1);
            eMesh.get_bulk_data()->declare_relation(newElement, eMesh.createOrGetNode(elems[ielem][2]), 2);
            eMesh.get_bulk_data()->declare_relation(newElement, eMesh.createOrGetNode(elems[ielem][3]), 3);
            
            set_parent_child_relations(eMesh, element, newElement, *ft_element_pool, ielem);

            std::vector<stk::mesh::Entity> elements(1,element);
            eMesh.prolongateElementFields( elements, newElement);

            ft_element_pool++;
            element_pool++;

          }


      }

    };

  }

#endif
