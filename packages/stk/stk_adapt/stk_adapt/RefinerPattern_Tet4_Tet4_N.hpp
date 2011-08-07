#ifndef stk_adapt_RefinerPattern_Tet4_Tet4_N_sierra_hpp
#define stk_adapt_RefinerPattern_Tet4_Tet4_N_sierra_hpp


//#include "UniformRefinerPattern.hpp"
#include <stk_adapt/sierra_element/RefinementTopology.hpp>
#include <stk_adapt/sierra_element/StdMeshObjTopologies.hpp>

#include <stk_adapt/RefinerPattern_Tri3_Tri3_N.hpp>

#include <stk_adapt/Percept_MOAB_SimplexTemplateRefiner.hpp>

#include "UniformRefinerPattern_Line2_Line2_2_sierra.hpp"

/** NOTE: A lot of the following code is unfinished (greedy triangulation scheme).  
 *  The call to triangulate_tet_generic has been replaced with a call to the 
 *  MOAB SimplexTemplateRefiner as modified for Percept/Adapt.
 *
 *  We are leaving this unfinished code in place until the MOAB testing is verified.  
 */

#define PERCEPT_USE_MOAB_REFINER 1

namespace stk {
  namespace adapt {

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

    static unsigned tbl_tet_face_edge_map[4][3]  = { {0, 4, 3}, {1, 5, 4}, {3, 5, 2}, {2, 1, 0} };
    static unsigned tbl_tet_face_nodes[4][3]     = { {0, 1, 3}, {1, 2, 3}, {0, 3, 2}, {0, 2, 1} };
    static unsigned tbl_tet_edge_nodes[6][2]     = { {0, 1}, {1, 2}, {2, 0}, {0, 3}, {1, 3}, {2, 3} };

    typedef boost::tuple<unsigned, unsigned, unsigned, unsigned> TetTupleTypeLocal;
    typedef boost::tuple<stk::mesh::EntityId, stk::mesh::EntityId, stk::mesh::EntityId, stk::mesh::EntityId> TetTupleType;

    /// general refinement pattern
    
    // the "-1" here signifies the number of elements created is not fixed, depends on the marking pattern
    template <>
    class RefinerPattern<shards::Tetrahedron<4>, shards::Tetrahedron<4>, -1 > : public URP<shards::Tetrahedron<4>,shards::Tetrahedron<4>  >
    {

      UniformRefinerPattern<shards::Line<2>, shards::Line<2>, 2, SierraPort > * m_edge_breaker;
      RefinerPattern<shards::Triangle<3>, shards::Triangle<3>, -1 > * m_face_breaker;

    public:

      RefinerPattern(percept::PerceptMesh& eMesh, BlockNamesType block_names = BlockNamesType()) :  URP<shards::Tetrahedron<4>, shards::Tetrahedron<4>  >(eMesh)
      {
        m_primaryEntityRank = eMesh.element_rank();

        setNeededParts(eMesh, block_names, true);
        Elem::StdMeshObjTopologies::bootstrap();

        m_edge_breaker =  new UniformRefinerPattern<shards::Line<2>, shards::Line<2>, 2, SierraPort > (eMesh, block_names) ;
        m_face_breaker =  new RefinerPattern<shards::Triangle<3>, shards::Triangle<3>, -1 > (eMesh, block_names) ;

      }

      void setSubPatterns( std::vector<UniformRefinerPatternBase *>& bp, percept::PerceptMesh& eMesh )
      {
        EXCEPTWATCH;
        bp = std::vector<UniformRefinerPatternBase *>(3u, 0);
        bp[0] = this;
        bp[1] = m_face_breaker;
        bp[2] = m_edge_breaker;

      }

      virtual void doBreak() {}
      void fillNeededEntities(std::vector<NeededEntityType>& needed_entities)
      {
#if 0
        // FIXME - tmp, for now to create a centroid node which we delete later
        needed_entities.resize(2);
        needed_entities[0].first = m_eMesh.edge_rank();    
        needed_entities[0].second = 1u;
        needed_entities[1].first = m_eMesh.element_rank();    
        needed_entities[1].second = 1u;
#else
        needed_entities.resize(1);
        needed_entities[0].first = m_eMesh.edge_rank();    
        needed_entities[0].second = 1u;
#endif

      }

      // FIXME - for now, create more than we need (to fix this right we need a change to the Refiner.cpp interface)
      virtual unsigned getNumNewElemPerElem() { return 8; }

      static stk::mesh::Entity* get_new_node(percept::PerceptMesh& eMesh)
      {
        //stk::mesh::PartVector empty ;
        std::vector<stk::mesh::Entity *> requested_entities;
        eMesh.createEntities(stk::mesh::fem::FEMMetaData::NODE_RANK, 1, requested_entities);
        return requested_entities[0];
      }

      // rotate in place so that the minimum node index is first
      static void normalize(tri_tuple_type_local& face)
      {
        //tri_tuple_type_local new_face;
        int nodes[3] = {face.get<0>(), face.get<1>(), face.get<2>() };
        int min_node_val = 1000;
        int min_node_index = -1;
        for (int i = 0; i < 3; i++)
          {
            if (nodes[i] < min_node_val)
              {
                min_node_val = nodes[i];
                min_node_index = i;
              }
          }
        face.get<0>() = nodes[min_node_index];
        face.get<1>() = nodes[(min_node_index+1) % 3];
        face.get<2>() = nodes[(min_node_index+2) % 3];
      }


      static bool is_valid_tet(std::vector<tri_tuple_type_local>& three_faces, TetTupleTypeLocal& valid_tet)
      {
        bool is_valid = false;
        std::set<unsigned> node_set;
        for (unsigned i=0; i < 3; i++)
          {
            node_set.insert(three_faces[i].get<0>());
            node_set.insert(three_faces[i].get<1>());
            node_set.insert(three_faces[i].get<2>());
          }
        if (node_set.size() == 4)
          {
            is_valid = true;
            valid_tet.get<0>() = three_faces[0].get<0>();
            valid_tet.get<1>() = three_faces[0].get<1>();
            valid_tet.get<2>() = three_faces[0].get<2>();
            node_set.erase(three_faces[0].get<0>());
            node_set.erase(three_faces[0].get<1>());
            node_set.erase(three_faces[0].get<2>());
            valid_tet.get<3>() = *node_set.begin();

            // check_volume(valid_tet,...
          }
        return is_valid;
      }


      static bool check_tet(std::vector<tri_tuple_type_local>& three_faces, TetTupleTypeLocal& valid_tet)
      {
        // FIXME
        return false;
      }

      static bool add_face(std::set<tri_tuple_type_local>& tet_faces_local_set, std::vector<tri_tuple_type_local>& tet_faces_local, std::vector<TetTupleTypeLocal>& tets)
      {
        unsigned nface=tet_faces_local.size();
        if (nface < 3)
          {
            throw std::logic_error("RefinerPattern_Tet4_Tet4_N::add_face: found nface < 3");
          }
        std::vector<tri_tuple_type_local> three_faces(3);
        for (unsigned i = 0; i < nface-2; i++)
          {
            //
            // double best_quality = 1e30; ... FIXME
            //
            for (unsigned j = i+1; j < nface-1; j++)
              {
                for (unsigned k = j+1; k < nface; k++)
                  {
                    three_faces[0] = tet_faces_local[i];
                    three_faces[1] = tet_faces_local[j];
                    three_faces[2] = tet_faces_local[k];
                    TetTupleTypeLocal valid_tet;
                    bool is_valid_tet = check_tet(three_faces, valid_tet);
                    if (is_valid_tet)
                      {
                        tets.push_back(valid_tet);
                        unsigned nodes[4] = {valid_tet.get<0>(), valid_tet.get<1>(), valid_tet.get<2>(), valid_tet.get<3>() };
                        for (int iface = 0; iface < 4; iface++)
                          {
                            tri_tuple_type_local tet_face(nodes[tbl_tet_face_nodes[iface][0]], 
                                                          nodes[tbl_tet_face_nodes[iface][1]],
                                                          nodes[tbl_tet_face_nodes[iface][2]] );
                            normalize(tet_face);
                            tet_faces_local_set.insert(tet_face);
                          }
                        for (int iface = 0; iface < 3; iface++)
                          {
                            tet_faces_local_set.erase(three_faces[iface]);
                          }
                        tet_faces_local = std::vector<tri_tuple_type_local>(tet_faces_local_set.begin(), tet_faces_local_set.end());
                        return true;
                      }
                  }
              }
          }
        return false;
      }

      static void triangulate_tet_generic(std::vector<tri_tuple_type_local>& tet_faces, double * centroid_coord, double *node_coords[10], std::vector<TetTupleTypeLocal>& tets)
      {
#if PERCEPT_USE_MOAB_REFINER

#else
        // algorithm is a greedy triangulation scheme:
        //   1. form triples of external faces (initialize tet_faces_local from @param tet_faces)
        //   2. if a triple forms a valid tet (i.e. it has positive volume and shares exactly 4 nodes):
        //         a. output it into @param tets
        //         b. remove its faces from tet_faces_local, add any face not originally tet_faces_local to tet_faces_local to
        //                maintain tet_faces_local as the list of external faces
        //   
        // Note: forming triples is brute-force below, but could easily be optimized by using a node-neighbors map
        tets.resize(0);
        TetTupleTypeLocal valid_tet(0,0,0,0);

        std::set<tri_tuple_type_local> tet_faces_local_set(tet_faces.begin(), tet_faces.end());
        std::vector<tri_tuple_type_local> tet_faces_local = tet_faces;

        int iloop=0;
        while (tet_faces_local_set.size())
          {
            add_face(tet_faces_local_set, tet_faces_local, tets);
            ++iloop;
            if (iloop > 1000)
              throw std::logic_error("RefinerPattern_Tet4_Tet4_N::add_face: iloop");
          }
#endif
      }
      

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

      static void triangulate_tet(PerceptMesh& eMesh, stk::mesh::Entity *tet_elem_nodes[4], unsigned edge_marks[6], 
                                  std::vector<TetTupleTypeLocal>& tets)
      {

        static stk::mesh::Entity *centroid_node = 0;
        static stk::mesh::Entity *temp_edge_nodes[6] = {0,0,0,0,0,0};
        if (!centroid_node)
          {
            centroid_node = get_new_node(eMesh);
            if (!centroid_node)
              throw std::logic_error("RefinerPattern_Tet4_Tet4_N::triangulate_tet: centroid_node is null");
            for (int i = 0; i < 6; i++)
              temp_edge_nodes[i] = get_new_node(eMesh);
          }

        const CellTopologyData * const cell_topo_data = shards::getCellTopologyData< shards::Tetrahedron<4> >();

        CellTopology cell_topo(cell_topo_data);
        //VectorFieldType* coordField = eMesh.getCoordinatesField();

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

        // uniform refinement
        if (num_edges_marked == 6)
          {
            tets.resize(8);

            // FIXME - use MOAB

            // FIXME - use StdMeshObjTopologies directly
            //stk::adapt::Elem::MeshObjTopology mesh_obj_topo;
            //const RefinementTopology *ref_topo = mesh_obj_topo.getRefinementTopology(Elem::CellTopology(cell_topo));

            /**
             * 
             *   CHILD 4-Node Tetrahedron 3D Object Node Maps:
             * |
             * | static const UInt child_0[] = { 0, 4, 6, 7 };  // srkenno 091410 fixed (used to be {0, 4, 8, 7} )
             * | static const UInt child_1[] = { 4, 1, 5, 8 };
             * | static const UInt child_2[] = { 6, 5, 2, 9 };
             * | static const UInt child_3[] = { 7, 8, 9, 3 };
             * | static const UInt child_4[] = { 8, 7, 6, 4 };
             * | static const UInt child_5[] = { 6, 9, 8, 5 };
             * | static const UInt child_6[] = { 9, 8, 7, 6 };
             * | static const UInt child_7[] = { 5, 6, 4, 8 };
             * |
             */

            tets[0] = TetTupleTypeLocal( 0, 4, 6, 7 );
            tets[1] = TetTupleTypeLocal( 4, 1, 5, 8 );
            tets[2] = TetTupleTypeLocal( 6, 5, 2, 9 );
            tets[3] = TetTupleTypeLocal( 7, 8, 9, 3 );
            tets[4] = TetTupleTypeLocal( 8, 7, 6, 4 );
            tets[5] = TetTupleTypeLocal( 6, 9, 8, 5 );
            tets[6] = TetTupleTypeLocal( 9, 8, 7, 6 );
            tets[7] = TetTupleTypeLocal( 5, 6, 4, 8 );
          }
        else if (num_edges_marked == 0)
          {
            return;
          }
        // general case
        else
          {
            if (PERCEPT_USE_MOAB_REFINER)
              {
                double * node_coords[4];
                for (int inode=0; inode < 4; inode++)
                  {
                    node_coords[inode] = stk::mesh::field_data( *eMesh.getCoordinatesField() , *tet_elem_nodes[inode] );
                    if (0) std::cout << "tmp RP node_coords= " 
                                     << node_coords[inode][0] << " "
                                     << node_coords[inode][1] << " "
                                     << node_coords[inode][2] << std::endl;
                  }

                std::vector<moab::TetTupleInt> new_tets;
                moab::SimplexTemplateRefiner str;
                str.refine_3_simplex(new_tets, edge_marks, 1,
                                     node_coords[0], 0, tet_elem_nodes[0]->identifier(),
                                     node_coords[1], 0, tet_elem_nodes[1]->identifier(),
                                     node_coords[2], 0, tet_elem_nodes[2]->identifier(),
                                     node_coords[3], 0, tet_elem_nodes[3]->identifier() );

                if (0)
                  {
                    for (int inode=0; inode < 4; inode++)
                      {
                        node_coords[inode] = stk::mesh::field_data( *eMesh.getCoordinatesField() , *tet_elem_nodes[inode] );
                        std::cout << "tmp RefPatt::createNewElements node_coords after= " 
                                  << node_coords[inode][0] << " "
                                  << node_coords[inode][1] << " "
                                  << node_coords[inode][2] << std::endl;
                      }
                  }

                tets.resize(new_tets.size());
                for (unsigned i = 0; i < new_tets.size(); i++)
                  {
                    tets[i] = TetTupleTypeLocal((unsigned)new_tets[i].get<0>(),
                                                (unsigned)new_tets[i].get<1>(),
                                                (unsigned)new_tets[i].get<2>(),
                                                (unsigned)new_tets[i].get<3>() );
                    if (0)
                      std::cout << "tmp RefPatt::createNewElements new tet= " << tets[i] << std::endl;

                  }
                return;
              }


            tets.resize(0);
            std::vector<tri_tuple_type_local> tet_faces;

            // algorithm: triangulate each face, create tets from centroid, collapse shortest edge

            for (int iface = 0; iface < 4; iface++)
              {
                unsigned tri_face_edge_marks[3] = {0,0,0};
                for (int iedge = 0; iedge < 3; iedge++)
                  {
                    tri_face_edge_marks[iedge] = edge_marks[tbl_tet_face_edge_map[iface][iedge]];
                  }

                stk::mesh::Entity *tri_elem_nodes_local[3] = {0,0,0};
                for (int inode=0; inode < 3; inode++)
                  {
                    tri_elem_nodes_local[inode] = tet_elem_nodes[tbl_tet_face_nodes[iface][inode]];
                  }

                std::vector<tri_tuple_type_local> tri_face_elems_local;
                RefinerPattern<shards::Triangle<3>,      shards::Triangle<3>,     -1  >::triangulate_face(eMesh, tri_elem_nodes_local, tri_face_edge_marks, tri_face_elems_local);


                // triangulate_face returns tri_face_elems_local, a vector of 3-tuples with local face/edge numbers: 
                //     {0,1,2} for the vertices, and {3,4,5} for the edges
                // These need to be converted to the parent tet's vertices {0-3}, and edges {4-9}, done by TET_TRI_CV_EV 

#define TET_TRI_CV_EV(iface,i) ( i < 3 ? tbl_tet_face_nodes[iface][i] : tbl_tet_face_edge_map[iface][i-3] )

                for (unsigned itri = 0; itri < tri_face_elems_local.size(); itri++)
                  {
                    normalize(tri_face_elems_local[itri]);

                    // FIXME - do we need to reverse two nodes for positive volumes since faces have right-hand-rule pointing outward?
                    // The current scheme gives all outward-pointing faces, so when tets are created, there should be a reversal, or
                    // we could do it here.
                    tri_tuple_type_local tet_face(TET_TRI_CV_EV(iface, tri_face_elems_local[itri].get<0>()), 
                                                  TET_TRI_CV_EV(iface, tri_face_elems_local[itri].get<1>()), 
                                                  TET_TRI_CV_EV(iface, tri_face_elems_local[itri].get<2>()) );
                    normalize(tet_face);
                    tet_faces.push_back(tet_face);
                  }
              }

#define TET_CENTROID_NODE (4+6)

            // compute centroid and temp_edge_nodes
            double * node_coords[10];  // all node coords in order, vertices {0-3}, edges {4-9}
            double * centroid_coord = stk::mesh::field_data( *eMesh.getCoordinatesField() , *centroid_node);
            for (int dim=0; dim < 3; dim++) centroid_coord[dim] = 0.0;
            for (int inode=0; inode < 4; inode++)
              {
                double * coord = stk::mesh::field_data( *eMesh.getCoordinatesField() , *tet_elem_nodes[inode] );
                for (int dim=0; dim < 3; dim++) centroid_coord[dim] += coord[dim]/4.;
                node_coords[inode] = coord;
              }
            int jnode = 4;
            for (int iedge=0; iedge < 6; iedge++)
              {
                double * coord0 = stk::mesh::field_data( *eMesh.getCoordinatesField() , *tet_elem_nodes[tbl_tet_edge_nodes[iedge][0]] );
                double * coord1 = stk::mesh::field_data( *eMesh.getCoordinatesField() , *tet_elem_nodes[tbl_tet_edge_nodes[iedge][1]] );
                double * temp_node_coord = stk::mesh::field_data( *eMesh.getCoordinatesField() , *temp_edge_nodes[iedge] );
                node_coords[jnode++] = temp_node_coord;
                for (int dim=0; dim < 3; dim++) temp_node_coord[dim] = 0.5*(coord0[dim] + coord1[dim]);
              }
                
            triangulate_tet_generic(tet_faces, centroid_coord, node_coords, tets);
          }
      }


      void 
      createNewElements(percept::PerceptMesh& eMesh, NodeRegistry& nodeRegistry, 
                        stk::mesh::Entity& element,  NewSubEntityNodesType& new_sub_entity_nodes, std::vector<stk::mesh::Entity *>::iterator& element_pool,
                        stk::mesh::FieldBase *proc_rank_field=0)
      {
        const CellTopologyData * const cell_topo_data = stk::percept::PerceptMesh::get_cell_topology(element);
        static std::vector<TetTupleType> elems(8);
        static std::vector<TetTupleTypeLocal> elems_local(8);
        unsigned num_new_elems=0;

        CellTopology cell_topo(cell_topo_data);
        const stk::mesh::PairIterRelation elem_nodes = element.relations(stk::mesh::fem::FEMMetaData::NODE_RANK);
        //VectorFieldType* coordField = eMesh.getCoordinatesField();

        std::vector<stk::mesh::Part*> add_parts;
        std::vector<stk::mesh::Part*> remove_parts;
        add_parts = m_toParts;
        
        unsigned edge_marks[6] = {0,0,0,0,0,0};
        unsigned num_edges_marked=0;
        for (int iedge = 0; iedge < 6; iedge++)
          {
            unsigned num_nodes_on_edge = new_sub_entity_nodes[m_eMesh.edge_rank()][iedge].size();
            if (num_nodes_on_edge)
              {
                edge_marks[iedge] = 1;
                ++num_edges_marked;
              }
          }
        if (num_edges_marked == 0)
          return;

        stk::mesh::Entity *elem_nodes_local[4] = {0,0,0,0};
        for (int inode=0; inode < 4; inode++)
          {
            elem_nodes_local[inode] = elem_nodes[inode].entity();
          }
        triangulate_tet(eMesh, elem_nodes_local, edge_marks, elems_local);
        
        //#define TET_CV_EV(i) ( i < 4 ? VERT_N(i) : (i < TET_CENTROID_NODE ? EDGE_N(i-4) : -1) )
#define TET_CV_EV(i) ( i < 4 ? VERT_N(i) : EDGE_N(i-4) )

        num_new_elems = elems_local.size();
        elems.resize(num_new_elems);
        for (unsigned ielem=0; ielem < num_new_elems; ielem++)
          {
            elems[ielem] = TetTupleType( TET_CV_EV(elems_local[ielem].get<0>() ), 
                                         TET_CV_EV(elems_local[ielem].get<1>() ), 
                                         TET_CV_EV(elems_local[ielem].get<2>() ), 
                                         TET_CV_EV(elems_local[ielem].get<3>() ) );
            if (0)
              std::cout << "tmp RefPatt::createNewElements new tet= " << elems[ielem] << std::endl;

          }

        //std::cout << "tmp RefinerPattern_Tet4_Tet4_N::num_edges_marked= " << num_edges_marked << std::endl;

        //nodeRegistry.makeCentroidCoords(*const_cast<stk::mesh::Entity *>(&element), m_eMesh.element_rank(), 0u);
        
        for (unsigned ielem=0; ielem < elems.size(); ielem++)
          {
            stk::mesh::Entity& newElement = *(*element_pool);

            if (proc_rank_field)
              {
                double *fdata = stk::mesh::field_data( *static_cast<const ScalarFieldType *>(proc_rank_field) , newElement );
                //fdata[0] = double(m_eMesh.getRank());
                fdata[0] = double(newElement.owner_rank());
              }

            eMesh.getBulkData()->change_entity_parts( newElement, add_parts, remove_parts );

            set_parent_child_relations(eMesh, element, newElement, ielem);

            interpolateElementFields(eMesh, element, newElement);

            {
              if (!elems[ielem].get<0>())
                {
                  std::cout << "P[" << eMesh.getRank() << "] nid = 0 << " << std::endl;
                  //exit(1);
                }
            }

            // 4 nodes of the new tets
            eMesh.getBulkData()->declare_relation(newElement, eMesh.createOrGetNode(elems[ielem].get<0>()), 0);
            eMesh.getBulkData()->declare_relation(newElement, eMesh.createOrGetNode(elems[ielem].get<1>()), 1);
            eMesh.getBulkData()->declare_relation(newElement, eMesh.createOrGetNode(elems[ielem].get<2>()), 2);
            eMesh.getBulkData()->declare_relation(newElement, eMesh.createOrGetNode(elems[ielem].get<3>()), 3);

            if (0)
              {
                std::cout << "tmp RefPatt::createNewElements element.identifier()= " << element.identifier() 
                          << " newElement= " << newElement.identifier() << std::endl;
                
              }

            element_pool++;

          }

      
      }
      
    };

  }
}
#endif
