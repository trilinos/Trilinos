#ifndef stk_adapt_UniformRefinerPattern_Hex8_Tet4_6_12_hpp
#define stk_adapt_UniformRefinerPattern_Hex8_Tet4_6_12_hpp

#include "UniformRefinerPattern.hpp"
#include "UniformRefinerPattern_Quad4_Tri3_2.hpp"

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
    class UniformRefinerPattern<shards::Hexahedron<8>, shards::Tetrahedron<4>, 6 > : public URP<shards::Hexahedron<8>, shards::Tetrahedron<4> >
    {
    private:
#define USE_FACE_BREAKER_H8_T4_6 1
#if USE_FACE_BREAKER_H8_T4_6
      UniformRefinerPattern<shards::Quadrilateral<4>, shards::Triangle<3>, 2 > *m_face_breaker;
#endif

    public:

      UniformRefinerPattern(percept::PerceptMesh& eMesh, BlockNamesType block_names = BlockNamesType()) : URP<shards::Hexahedron<8>, shards::Tetrahedron<4> >(eMesh)
      {
        EXCEPTWATCH;

        m_primaryEntityRank = eMesh.element_rank();

        setNeededParts(eMesh, block_names, false);

#if USE_FACE_BREAKER_H8_T4_6
        m_face_breaker = new UniformRefinerPattern<shards::Quadrilateral<4>, shards::Triangle<3>, 2 > (eMesh, block_names);
#endif
      }
      ~UniformRefinerPattern()
      {
        if (m_face_breaker) delete m_face_breaker;
      }


      virtual void doBreak() {}

      void fillNeededEntities(std::vector<NeededEntityType>& needed_entities)
      {
        needed_entities.resize(1);
        needed_entities[0].first = m_eMesh.element_rank();  
        setToOne(needed_entities);
      }

      void setSubPatterns( std::vector<UniformRefinerPatternBase *>& bp, percept::PerceptMesh& eMesh )
      {
        EXCEPTWATCH;
        bp = std::vector<UniformRefinerPatternBase *>(2u, 0);
        bp[0] = this;

#if USE_FACE_BREAKER_H8_T4_6
        if (1)
          {
            //             UniformRefinerPattern<shards::Quadrilateral<4>, shards::Triangle<3>, 4, Specialization > *face_breaker = 
            //               new UniformRefinerPattern<shards::Quadrilateral<4>, shards::Triangle<3>, 4, Specialization > (eMesh);
            bp[1] = m_face_breaker;
          }
#endif
      }

      /// NOTE: we create additional un-used elements if the Hex8 can be broken into 6 tets
      virtual unsigned getNumNewElemPerElem() { return 12u; }

      void 
      createNewElements(percept::PerceptMesh& eMesh, NodeRegistry& nodeRegistry,
                        stk::mesh::Entity& element,  NewSubEntityNodesType& new_sub_entity_nodes, vector<stk::mesh::Entity *>::iterator& element_pool,
                        stk::mesh::FieldBase *proc_rank_field=0)
      {
        EXCEPTWATCH;
        const CellTopologyData * const cell_topo_data = stk::percept::PerceptMesh::get_cell_topology(element);
        typedef boost::tuple<stk::mesh::EntityId, stk::mesh::EntityId, stk::mesh::EntityId, stk::mesh::EntityId> tet_tuple_type;
        vector<tet_tuple_type> new_elements(6);

        shards::CellTopology cell_topo(cell_topo_data);
        const stk::mesh::PairIterRelation elem_nodes = element.relations(stk::mesh::fem::FEMMetaData::NODE_RANK);

        // for cases that have a single center node, we just compute the new node's quantities here instead of globally
        //stk::mesh::Entity * node = getBulkData()->get_entity( Node, node_id );

#define CENTROID_N NN(m_eMesh.element_rank(), 0)  

#if STK_ADAPT_URP_LOCAL_NODE_COMPS
        nodeRegistry.makeCentroidCoords(*const_cast<stk::mesh::Entity *>(&element), Element, 0u);
        nodeRegistry.addToExistingParts(*const_cast<stk::mesh::Entity *>(&element), Element, 0u);
        nodeRegistry.interpolateFields(*const_cast<stk::mesh::Entity *>(&element), Element, 0u);
#endif

        // following code is from SweepMesher::breakElement, modified here for stk_mesh
        // from here------------------------------------------->>>>>>
        {
          static unsigned loc_trifaces[6][2][3];  // iQuadFace, kTriOnQuadFace, jTriNodeIndex
          static unsigned loc_qfaces[6][4]; // iHexFaceOrd, iFaceNodeOrd

          static unsigned element_globalIds[8] = {0,0,0,0, 0,0,0,0};
          //static std::vector<unsigned> element_globalIds(8);
          const stk::mesh::PairIterRelation elem_nodes = element.relations( stk::mesh::fem::FEMMetaData::NODE_RANK );

          //std::cout << "tmp hex elem= " << element << std::endl;
          for (int inode=0; inode < 8; inode++)
            {
              stk::mesh::Entity & node = * elem_nodes[inode].entity();
              element_globalIds[inode] = node.identifier();
            }

          // build quad face break patterns - use simple min node index as first node
          const CellTopologyData *hex_topo = shards::getCellTopologyData<shards::Hexahedron<8> >();

          // node valences
          unsigned valences[8]={0,0,0,0,0,0,0,0};
          for (unsigned iHexFaceOrd = 0; iHexFaceOrd < 6; iHexFaceOrd++)
            {
              const CellTopologyData_Subcell& face = hex_topo->side[iHexFaceOrd];
              unsigned iqf        = face.node[0];
              unsigned globalIqf  = element_globalIds[iqf];
              unsigned minVal     = globalIqf;
              unsigned indxMinVal = 0;
              for (unsigned iFaceNodeOrd=1; iFaceNodeOrd < 4; iFaceNodeOrd++)
                {
                  iqf = face.node[iFaceNodeOrd];
                  globalIqf = element_globalIds[iqf];
                  if (globalIqf < minVal)
                    {
                      minVal = globalIqf;
                      indxMinVal = iFaceNodeOrd;
                    }
                }
              // permute to make min node index come first
              for (unsigned iFaceNodeOrd=0; iFaceNodeOrd < 4; iFaceNodeOrd++)
                {
                  unsigned jFaceNodeOrd = (iFaceNodeOrd + indxMinVal) % 4;
                  //qfaces[iHexFaceOrd][iFaceNodeOrd] = element_globalIds[face.node[jFaceNodeOrd]];
                  loc_qfaces[iHexFaceOrd][iFaceNodeOrd] = face.node[jFaceNodeOrd];
                }
              if (0)
                std::cout << "tmp hex face[" << iHexFaceOrd << "] = " 
                          << element_globalIds[loc_qfaces[iHexFaceOrd][0]] << " "
                          << element_globalIds[loc_qfaces[iHexFaceOrd][1]] << " "
                          << element_globalIds[loc_qfaces[iHexFaceOrd][2]] << " "
                          << element_globalIds[loc_qfaces[iHexFaceOrd][3]] << std::endl;

              // each quad face is now broken into tri faces as {0,1,2}, {0,2,3}
              loc_trifaces[iHexFaceOrd][0][0] = loc_qfaces[iHexFaceOrd][0];
              loc_trifaces[iHexFaceOrd][0][1] = loc_qfaces[iHexFaceOrd][1];
              loc_trifaces[iHexFaceOrd][0][2] = loc_qfaces[iHexFaceOrd][2];
              loc_trifaces[iHexFaceOrd][1][0] = loc_qfaces[iHexFaceOrd][0];
              loc_trifaces[iHexFaceOrd][1][1] = loc_qfaces[iHexFaceOrd][2];
              loc_trifaces[iHexFaceOrd][1][2] = loc_qfaces[iHexFaceOrd][3];

              valences[loc_trifaces[iHexFaceOrd][0][0]]++;
              valences[loc_trifaces[iHexFaceOrd][0][1]]++;
              valences[loc_trifaces[iHexFaceOrd][0][2]]++;
              valences[loc_trifaces[iHexFaceOrd][1][0]]++;
              valences[loc_trifaces[iHexFaceOrd][1][1]]++;
              valences[loc_trifaces[iHexFaceOrd][1][2]]++;
            }

          // find max valence 
          unsigned vmaxIndx = 0;
          unsigned vmax = valences[0];
          for (unsigned iv = 1; iv < 8; iv++)
            {
              if (valences[iv] > vmax)
                {
                  vmax = valences[iv];
                  vmaxIndx = iv;
                }
            }

          if (vmax != 6)
            {
              /// Rare case - create tets by joining centroid to each face - for now, just throw an exception to see how often this
              /// case occurs - FIXME - take this exception out later
              if (1)
                {
                  std::ostringstream msg;
                  msg << "shouldn't ever get this exception - please inform srkenno@sandia.gov\n";
                  throw std::runtime_error( msg.str() );
                }
              new_elements.resize(12);

              //               boost::array<double,3> centroid = {{0,0,0}};
              //               for (unsigned iv = 0; iv < 8; iv++)
              //                 {
              //                   centroid[0] += m_node_coords[element_globalIds[iv]][0]/8.0;
              //                   centroid[1] += m_node_coords[element_globalIds[iv]][1]/8.0;
              //                   centroid[2] += m_node_coords[element_globalIds[iv]][2]/8.0;
              //                 }
              //               m_node_coords.push_back(centroid);

              //               unsigned newNodeIndex = m_node_coords.size() - 1;

              int iele=0;
              for (unsigned iQuadFace = 0; iQuadFace < 6; iQuadFace++)
                {
                  for (unsigned kTriOnQuadFace = 0; kTriOnQuadFace < 2; kTriOnQuadFace++)
                    {
                      //                       m_elems[shards_Tetrahedron_4].push_back(newNodeIndex);
                      //                       for (unsigned jTriNodeIndex = 0; jTriNodeIndex < 3; jTriNodeIndex++)
                      //                         {
                      //                           m_elems[shards_Tetrahedron_4].push_back(element_globalIds[loc_trifaces[iQuadFace][kTriOnQuadFace][jTriNodeIndex]]);
                      //                         }

                      new_elements[iele++] = tet_tuple_type(CENTROID_N, 
                                                            element_globalIds[loc_trifaces[iQuadFace][kTriOnQuadFace][0]],
                                                            element_globalIds[loc_trifaces[iQuadFace][kTriOnQuadFace][1]],
                                                            element_globalIds[loc_trifaces[iQuadFace][kTriOnQuadFace][2]]);

                    }
                }
            }
          else
            {
              /// normal case - connect max valence node to other faces without that node (should be 6 tets always)
              /// Note: there is a 5-tet configuration that exists for some face diagonal configurations - FIXME - could add this case later
              ///    The 5-tet case consists of an interior tet with no boundary faces, and 4 corner tets; the boundary faces have
              ///      to each have alternating diagonals along the 3 axis directions for this configuration to exist
              unsigned count=0;
              for (unsigned iQuadFace = 0; iQuadFace < 6; iQuadFace++)
                {
                  for (unsigned kTriOnQuadFace = 0; kTriOnQuadFace < 2; kTriOnQuadFace++)
                    {
                      bool isEqual = false;
                      for (unsigned jTriNodeIndex = 0; jTriNodeIndex < 3; jTriNodeIndex++)
                        {
                          if (vmaxIndx == loc_trifaces[iQuadFace][kTriOnQuadFace][jTriNodeIndex])
                            {
                              isEqual = true;
                              break;
                            }
                        }
                      if (not isEqual)
                        {
                          //  m_elems[shards_Tetrahedron_4].push_back(element_globalIds[vmaxIndx]);
                          //  for (unsigned jTriNodeIndex = 0; jTriNodeIndex < 3; jTriNodeIndex++)
                          //    {
                          //       m_elems[shards_Tetrahedron_4].push_back(element_globalIds[loc_trifaces[iQuadFace][kTriOnQuadFace][jTriNodeIndex]]);
                          //    }
                          new_elements[count] = tet_tuple_type(element_globalIds[vmaxIndx],
                                                               element_globalIds[loc_trifaces[iQuadFace][kTriOnQuadFace][0]],
                                                               element_globalIds[loc_trifaces[iQuadFace][kTriOnQuadFace][1]],
                                                               element_globalIds[loc_trifaces[iQuadFace][kTriOnQuadFace][2]]);

                          ++count;

                        }
                    }
                }
              assert(count == 6);

            }
        }
        // to here <<<<<< -----------------------

        for (unsigned ielem=0; ielem < new_elements.size(); ielem++)
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

            unsigned nchild = new_elements.size();
            set_parent_child_relations(eMesh, element, newElement, ielem, &nchild);

            {
              if (!new_elements[ielem].get<0>())
                {
                  std::cout << "P[" << eMesh.getRank() << " nid = 0 << " << std::endl;
                  exit(123);
                }
            }
            eMesh.getBulkData()->declare_relation(newElement, eMesh.createOrGetNode(new_elements[ielem].get<0>()), 0);
            eMesh.getBulkData()->declare_relation(newElement, eMesh.createOrGetNode(new_elements[ielem].get<1>()), 1);
            eMesh.getBulkData()->declare_relation(newElement, eMesh.createOrGetNode(new_elements[ielem].get<2>()), 2);
            eMesh.getBulkData()->declare_relation(newElement, eMesh.createOrGetNode(new_elements[ielem].get<3>()), 3);

            element_pool++;
          }

        if (1 && new_elements.size() == 6)
          {
            for (unsigned ielem=0; ielem < 6; ielem++)
              {
                // destroy un-needed elems
                // elems_to_destroy.push_back(*element_pool);  ++element_pool;
                eMesh.getBulkData()->destroy_entity(*element_pool);
                ++element_pool;
              }
            //nodes_to_destroy.push_back(CENTROID_N)
            stk::mesh::Entity * node = eMesh.getBulkData()->get_entity( stk::mesh::fem::FEMMetaData::NODE_RANK, CENTROID_N);
            if (!node)
              {
                throw std::logic_error("UniformRefinerPattern_Hex8_Tet4_6_12:: node is null");
              }
            else
              {
                eMesh.getBulkData()->destroy_entity(node);
              }
          }
      }
      
    };
#undef CENTROID_N

  }
}
#endif
