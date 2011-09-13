#ifndef stk_adapt_UniformRefinerPattern_Hex8_Tet4_24_hpp
#define stk_adapt_UniformRefinerPattern_Hex8_Tet4_24_hpp

#include "UniformRefinerPattern.hpp"
#include "UniformRefinerPattern_Quad4_Tri3_4.hpp"

namespace stk {
  namespace adapt {

    /** From Shards_BasicTopologies.hpp
     * 
     * \brief Topological traits: Dimension = 3, Sides = 6, Edges = 12,              
     *         Vertices = 8, and Nodes = 8, 20, or 27.                               
     *                                                                               
     *  <PRE>                                                                        
     *   Linear 8-Node Hexahedron node locations.                                    
     *                                                                               
     *          7                    6                                               
     *           o------------------o                                                
     *          /|                 /|                                                
     *         / |                / |                                                
     *        /  |               /  |                                                
     *       /   |              /   |                                                
     *      /    |             /    |                                                
     *     /     |            /     |                                                
     *  4 /      |         5 /      |                                                
     *   o------------------o       |                                                
     *   |       |          |       |                                                
     *   |     3 o----------|-------o 2                                              
     *   |      /           |      /                                                 
     *   |     /            |     /                                                  
     *   |    /             |    /                                                   
     *   |   /              |   /                                                    
     *   |  /               |  /                                                     
     *   | /                | /                                                      
     *   |/                 |/                                                       
     *   o------------------o                                                        
     *  0                    1                                                       
     *                                                                               
     *                                                                               
     *   face numbering for symmetric hex to tet break pattern                      |   typedef                                                                  
     *                                                                              |     MakeTypeList< IndexList< 0 , 1 ,   8 > ,                               
     *           7                                                                  |                   IndexList< 1 , 2 ,   9 > ,                               
     *            o------------------o 6                                            |                   IndexList< 2 , 3 ,  10 > ,                               
     *           /|                 /|                                              |                   IndexList< 3 , 0 ,  11 > ,                               
     *          / |                / |                                              |                   IndexList< 4 , 5 ,  16 > ,                               
     *         /  |   13          /  |                                              |                   IndexList< 5 , 6 ,  17 > ,                               
     *        /   |    o         /   |                                              |                   IndexList< 6 , 7 ,  18 > ,                               
     *       /    |       o10   /    |     Node #14 is at centroid of element       |                   IndexList< 7 , 4 ,  19 > ,                               
     *      /     |            /     |                                              |                   IndexList< 0 , 4 ,  12 > ,                               
     *   4 /      |         5 /      |     "2D surface" containing nodes            |                   IndexList< 1 , 5 ,  13 > ,                               
     *    o------------------o    9  |      0,1,5,4 has node 25 at center....       |                   IndexList< 2 , 6 ,  14 > ,                               
     *    | 11o   | 3        |   o   |                                              |                   IndexList< 3 , 7 ,  15 > >::type                         
     *    |       o----------|-------o 2                                            |     HexahedronEdgeNodeMap ;                                                
     *    |      /           |      /                                               |                                                                            
     *    |     /   8        |     /                                                |   typedef                                                                  
     *    |    /    o        |    /                                                 |     MakeTypeList< IndexList< 0, 1, 5, 4,   8, 13, 16, 12,   25 > ,         
     *    |   /        o12   |   /                                                  |                   IndexList< 1, 2, 6, 5,   9, 14, 17, 13,   24 > ,         
     *    |  /               |  /                                                   |                   IndexList< 2, 3, 7, 6,  10, 15, 18, 14,   26 > ,         
     *    | /                | /                                                    |                   IndexList< 0, 4, 7, 3,  12, 19, 15, 11,   23 > ,         
     *    |/                 |/                                                     |                   IndexList< 0, 3, 2, 1,  11, 10,  9,  8,   21 > ,         
     *    o------------------o                                                      |                   IndexList< 4, 5, 6, 7,  16, 17, 18, 19,   22 > >::type   
     *   0                    1                                                     |     HexahedronFaceNodeMap ;                                                
     *                                                                              |
     * </PRE>                                                                       |
     *
     *
     */


    template <>
    class UniformRefinerPattern<shards::Hexahedron<8>, shards::Tetrahedron<4>, 24 > : public URP<shards::Hexahedron<8>, shards::Tetrahedron<4> >
    {
    private:
      UniformRefinerPattern<shards::Quadrilateral<4>, shards::Triangle<3>, 4, Specialization > *m_face_breaker;

    public:

      UniformRefinerPattern(percept::PerceptMesh& eMesh, BlockNamesType block_names = BlockNamesType()) : URP<shards::Hexahedron<8>, shards::Tetrahedron<4> >(eMesh)
       {
         EXCEPTWATCH;

         m_primaryEntityRank = eMesh.element_rank();

         setNeededParts(eMesh, block_names, false);

#define USE_FACE_BREAKER 1
#if USE_FACE_BREAKER
         m_face_breaker = new UniformRefinerPattern<shards::Quadrilateral<4>, shards::Triangle<3>, 4, Specialization > (eMesh, block_names);
#endif
       }
      ~UniformRefinerPattern()
      {
        if (m_face_breaker) delete m_face_breaker;
      }


      virtual void doBreak() {}

      void fillNeededEntities(std::vector<NeededEntityType>& needed_entities)
      {
        needed_entities.resize(2);
        needed_entities[0].first = m_eMesh.face_rank();
        needed_entities[1].first = m_eMesh.element_rank();  
        setToOne(needed_entities);

      }

      void setSubPatterns( std::vector<UniformRefinerPatternBase *>& bp, percept::PerceptMesh& eMesh )
      {
        EXCEPTWATCH;
        bp = std::vector<UniformRefinerPatternBase *>(1u, 0);
        bp[0] = this;

#if USE_FACE_BREAKER
        bp.push_back( m_face_breaker);
#endif
      }

      virtual unsigned getNumNewElemPerElem() { return 24u; }

      void 
      createNewElements(percept::PerceptMesh& eMesh, NodeRegistry& nodeRegistry,
                        stk::mesh::Entity& element,  NewSubEntityNodesType& new_sub_entity_nodes, vector<stk::mesh::Entity *>::iterator& element_pool,
                        stk::mesh::FieldBase *proc_rank_field=0)
      {
        EXCEPTWATCH;
        const CellTopologyData * const cell_topo_data = stk::percept::PerceptMesh::get_cell_topology(element);
        typedef boost::tuple<stk::mesh::EntityId, stk::mesh::EntityId, stk::mesh::EntityId, stk::mesh::EntityId> tet_tuple_type;
        static vector<tet_tuple_type> elems(24);

        shards::CellTopology cell_topo(cell_topo_data);
        const stk::mesh::PairIterRelation elem_nodes = element.relations(stk::mesh::fem::FEMMetaData::NODE_RANK);


        // FIXME - maybe the computation of node coorinates should go in the calling code?
        double tmp_x[3];
        for (unsigned i_face = 0; i_face < 6; i_face++)
          {
            double * pts[4] = {FACE_COORD(i_face, 0), FACE_COORD(i_face, 1), FACE_COORD(i_face, 2), FACE_COORD(i_face, 3)};

            double * mp = getCentroid(pts, 4, eMesh.getSpatialDim(), tmp_x);

#if 0
            std::cout << "pts = \n" 

                      << pts[0][0] << " " << pts[0][1] << " " << pts[0][2] << " \n" 
                      << pts[1][0] << " " << pts[1][1] << " " << pts[1][2] << " \n" 
                      << pts[2][0] << " " << pts[2][1] << " " << pts[2][2] << " \n" 
                      << pts[3][0] << " " << pts[3][1] << " " << pts[3][2] << " \n" 
                      << " centroid = " 
                      << mp[0] << " " << mp[1] << " " << mp[2] <<  std::endl;
#endif

            //stk::mesh::Entity& new_node =eMesh.createOrGetNode(FACE_N(i_face), mp);

            eMesh.createOrGetNode(FACE_N(i_face), mp);
            nodeRegistry.addToExistingParts(*const_cast<stk::mesh::Entity *>(&element), m_eMesh.face_rank(), i_face);
            nodeRegistry.interpolateFields(*const_cast<stk::mesh::Entity *>(&element), m_eMesh.face_rank(), i_face);

          }

        nodeRegistry.makeCentroidCoords(*const_cast<stk::mesh::Entity *>(&element), m_eMesh.element_rank(), 0u);
        nodeRegistry.addToExistingParts(*const_cast<stk::mesh::Entity *>(&element), m_eMesh.element_rank(), 0u);
        nodeRegistry.interpolateFields(*const_cast<stk::mesh::Entity *>(&element), m_eMesh.element_rank(), 0u);
        

        //#define C 14

// new_sub_entity_nodes[i][j]
#define CENTROID_N NN(m_eMesh.element_rank(), 0)  

        unsigned iele = 0;
        for (unsigned i_face = 0; i_face < 6; i_face++)
          {
            for (unsigned i_tri_on_face = 0; i_tri_on_face < 4; i_tri_on_face++)
              {
                unsigned itf  = cell_topo_data->side[i_face].node[ i_tri_on_face ];
                unsigned itfp = cell_topo_data->side[i_face].node[ (i_tri_on_face + 1) % 4];
                if ( itf > 7)
                  {
                    throw std::logic_error("UniformRefinerPattern_Hex8_Tet4_20 logic err");
                  }

                 elems[iele++] = tet_tuple_type(CENTROID_N, VERT_N(itf), VERT_N(itfp), FACE_N(i_face));
              }
          }

#undef CENTROID_N

        // write a diagram of the refinement pattern as a vtk file, or a latex/tikz/pgf file
#define WRITE_DIAGRAM 0
#if WRITE_DIAGRAM

#endif
        
        for (unsigned ielem=0; ielem < elems.size(); ielem++)
          {
            //stk::mesh::Entity& newElement = eMesh.getBulkData()->declare_entity(Element, *element_id_pool, eMesh.getPart(interface_table::shards_Triangle_3) );
            //stk::mesh::Entity& newElement = eMesh.getBulkData()->declare_entity(Element, *element_id_pool, eMesh.getPart(interface_table::shards_Triangle_3) );

            stk::mesh::Entity& newElement = *(*element_pool);

            if (proc_rank_field)
              {
                double *fdata = stk::mesh::field_data( *static_cast<const ScalarFieldType *>(proc_rank_field) , newElement );
                fdata[0] = double(newElement.owner_rank());
              }

            change_entity_parts(eMesh, element, newElement);
            set_parent_child_relations(eMesh, element, newElement, ielem);

            {
              if (!elems[ielem].get<0>())
                {
                  std::cout << "P[" << eMesh.getRank() << " nid = 0 << " << std::endl;
                  exit(1);
                }
            }
            eMesh.getBulkData()->declare_relation(newElement, eMesh.createOrGetNode(elems[ielem].get<0>()), 0);
            eMesh.getBulkData()->declare_relation(newElement, eMesh.createOrGetNode(elems[ielem].get<1>()), 1);
            eMesh.getBulkData()->declare_relation(newElement, eMesh.createOrGetNode(elems[ielem].get<2>()), 2);
            eMesh.getBulkData()->declare_relation(newElement, eMesh.createOrGetNode(elems[ielem].get<3>()), 3);


            element_pool++;

          }

      }
      
    };

  }
}
#endif
