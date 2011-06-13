#include <stdexcept>
#include <iostream>
#include <utility>
#include <set>

#undef NDEBUG
//#define NDEBUG
#include <cassert>

#include "SweepMesher.hpp"
#include "TransformPath.hpp"

#include <Shards_BasicTopologies.hpp>
#include <Shards_CellTopologyData.h>
#include <stk_util/parallel/Parallel.hpp>

#include <stk_mesh/fem/FEMMetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/GetEntities.hpp>

#include <stk_mesh/fem/FEMHelpers.hpp>
#include <stk_mesh/fem/Stencils.hpp>

#include <stk_io/MeshReadWriteUtils.hpp>
#include <stk_io/IossBridge.hpp>
#include <init/Ionit_Initializer.h>


// template<typename Element>
//   void clone(const VectorOfInt& oldElems, VectorOfInt& newElems, VectorOfCoord& nodePool, VectorOfCoord& newNodes, Transform& xform)
// {
  
// }

SHARDS_ARRAY_DIM_TAG_SIMPLE_IMPLEMENTATION( Tag1 )
  SHARDS_ARRAY_DIM_TAG_SIMPLE_IMPLEMENTATION( Tag2 )
  SHARDS_ARRAY_DIM_TAG_SIMPLE_IMPLEMENTATION( Tag3 )
  SHARDS_ARRAY_DIM_TAG_SIMPLE_IMPLEMENTATION( Tag4 )


namespace stk
{
  namespace percept 
  {
    using namespace util;

    enum { SpatialDim = 3 };
    enum { node_count = 21 };
    enum { number_hex = 3 };
    enum { number_wedge = 3 };
    enum { number_tetra = 3 };
    enum { number_pyramid = 2 };
    enum { number_shell_quad = 3 };
    enum { number_shell_tri = 3 };



    static const double node_coord_data[ node_count ][ SpatialDim ] = {
      { 0 , 0 , 0 } , { 1 , 0 , 0 } , { 2 , 0 , 0 } , { 3 , 0 , 0 } ,
      { 0 , 1 , 0 } , { 1 , 1 , 0 } , { 2 , 1 , 0 } , { 3 , 1 , 0 } ,
      { 0 , 2 , 0 } , { 1 , 2 , 0 } ,
      { 0 , 0 , -1 } , { 1 , 0 , -1 } , { 2 , 0 , -1 } , { 3 , 0 , -1 } ,
      { 0 , 1 , -1 } , { 1 , 1 , -1 } , { 2 , 1 , -1 } , { 3 , 1 , -1 } ,
      { 0 , 2 , -1 } , { 1 , 2 , -1 } ,
      { 1 , 1 , -2 } };

    static const stk::mesh::EntityId hex_node_ids[3][ shards::Hexahedron<> ::node_count ] = {
      { 1 , 2 , 12 , 11 , 5 , 6 , 16 , 15 } ,
      { 2 , 3 , 13 , 12 , 6 , 7 , 17 , 16 } ,
      { 3 , 4 , 14 , 13 , 7 , 8 , 18 , 17 } };

    static const stk::mesh::EntityId wedge_node_ids[3][ shards::Wedge<> ::node_count ] = {
      { 15 , 16 , 19 ,  5 ,  6 ,  9 } ,
      { 10 ,  9 ,  6 , 20 , 19 , 16 } ,
      { 16 , 17 , 20 ,  6 ,  7 , 10 } };

    static const stk::mesh::EntityId tetra_node_ids[3][ shards::Tetrahedron<> ::node_count ] = {
      { 15 , 19 , 16 , 21 } ,
      { 19 , 20 , 16 , 21 } ,
      { 16 , 20 , 17 , 21 } };

    static const stk::mesh::EntityId pyramid_node_ids[2][ shards::Pyramid<> ::node_count ] = {
      { 11 , 15 , 16 , 12 , 21 } ,
      { 12 , 16 , 17 , 13 , 21 } };

    static const stk::mesh::EntityId shell_quad_node_ids[3][ shards::ShellQuadrilateral<> ::node_count ]={
      { 9 , 6 , 16 , 19 } ,
      { 6 , 7 , 17 , 16 } ,
      { 7 , 8 , 18 , 17 } };

    static const stk::mesh::EntityId shell_tri_node_ids[3][ shards::ShellTriangle<> ::node_count ] ={
      { 19 , 16 , 21 } ,
      { 16 , 17 , 21 } ,
      { 17 , 13 , 21 } };


      


    template<> void SweepMesher::breakElement<shards_Quadrilateral_4, shards_Triangle_3>(unsigned elemIndex)
    {
      unsigned* elem = &m_elems[shards_Quadrilateral_4][elemIndex*m_elemInfo[shards_Quadrilateral_4].vertex_count];
      VectorOfInt newElem = VectorOfInt(3);
      newElem[0] = elem[0];
      newElem[1] = elem[1];
      newElem[2] = elem[2];
      push_back(m_elems[shards_Triangle_3], newElem);
      newElem[0] = elem[0];
      newElem[1] = elem[2];
      newElem[2] = elem[3];
      push_back(m_elems[shards_Triangle_3], newElem);
      if (m_deleteAfterBreak)
        {
          VectorOfInt::iterator pos = m_elems[shards_Quadrilateral_4].begin()+(elemIndex*m_elemInfo[shards_Quadrilateral_4].vertex_count);
          m_elems[shards_Quadrilateral_4].erase(pos, pos + m_elemInfo[shards_Quadrilateral_4].vertex_count);
        }
    }


#if 0
    template<> void SweepMesher::breakAllElements<SweepMesher::shards_Quadrilateral_4, SweepMesher::shards_Triangle_3>()
    {
      unsigned numElems = m_elems[shards_Quadrilateral_4].size();
      for (unsigned elemIndex = 0; elemIndex < numElems; elemIndex++)
        {
          breakElement<SweepMesher::shards_Quadrilateral_4, SweepMesher::shards_Triangle_3 > ( elemIndex);
        }
    }
#endif


    /// Break a wedge element into 3 (common case) or 8 (rare case) tets using a constructive algorithm.
    template<> void SweepMesher::breakElement<shards_Wedge_6, shards_Tetrahedron_4>(unsigned elemIndex)
    {
      static unsigned loc_trifaces[4][2][3];  // iPseudoQuadFace, kTriOnQuadFace, jTriNodeIndex
      static unsigned loc_qfaces[3][4]; // iWedgeFaceOrd, iFaceNodeOrd
      unsigned* elem = &m_elems[shards_Wedge_6][elemIndex * m_elemInfo[shards_Wedge_6].vertex_count];

      // build quad face break patterns - use simple smallest node index as first node
      //shards::CellTopology wedge_topo(shards::getCellTopologyData<shards::Wedge<6> >() );
      const CellTopologyData *wedge_topo = shards::getCellTopologyData<shards::Wedge<6> >();
      unsigned numFaces = wedge_topo->side_count; // should be 5
      if (numFaces != 5)
        {
          assert(5 == numFaces);
        }

      // quad faces are first in the face list
      unsigned valences[6]={0,0,0,0,0,0};
      for (unsigned iWedgeFaceOrd = 0; iWedgeFaceOrd < 3; iWedgeFaceOrd++)
        {
          const CellTopologyData_Subcell& face = wedge_topo->side[iWedgeFaceOrd];
          unsigned iqf = face.node[0];
          unsigned globalIqf = elem[iqf];
          unsigned minVal = globalIqf;
          unsigned indxMinVal = 0;
          for (unsigned iFaceNodeOrd=1; iFaceNodeOrd < 4; iFaceNodeOrd++)
            {
              //qfaces[iWedgeFaceOrd][iFaceNodeOrd] = elem[face.node[iFaceNodeOrd]];
              iqf = face.node[iFaceNodeOrd];
              globalIqf = elem[iqf];
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
              //qfaces[iWedgeFaceOrd][iFaceNodeOrd] = elem[face.node[jFaceNodeOrd]];
              loc_qfaces[iWedgeFaceOrd][iFaceNodeOrd] = face.node[jFaceNodeOrd];
            }
          // each quad face is now broken into tri faces as {0,1,2}, {0,2,3}
          loc_trifaces[iWedgeFaceOrd][0][0] = loc_qfaces[iWedgeFaceOrd][0];
          loc_trifaces[iWedgeFaceOrd][0][1] = loc_qfaces[iWedgeFaceOrd][1];
          loc_trifaces[iWedgeFaceOrd][0][2] = loc_qfaces[iWedgeFaceOrd][2];
          loc_trifaces[iWedgeFaceOrd][1][0] = loc_qfaces[iWedgeFaceOrd][0];
          loc_trifaces[iWedgeFaceOrd][1][1] = loc_qfaces[iWedgeFaceOrd][2];
          loc_trifaces[iWedgeFaceOrd][1][2] = loc_qfaces[iWedgeFaceOrd][3];

          valences[loc_trifaces[iWedgeFaceOrd][0][0]]++;
          valences[loc_trifaces[iWedgeFaceOrd][0][1]]++;
          valences[loc_trifaces[iWedgeFaceOrd][0][2]]++;
          valences[loc_trifaces[iWedgeFaceOrd][1][0]]++;
          valences[loc_trifaces[iWedgeFaceOrd][1][1]]++;
          valences[loc_trifaces[iWedgeFaceOrd][1][2]]++;
        }
      // add in the top and bottom face as a new pseudo quad face
      loc_trifaces[3][0][0] = wedge_topo->side[3].node[0];
      loc_trifaces[3][0][1] = wedge_topo->side[3].node[1];
      loc_trifaces[3][0][2] = wedge_topo->side[3].node[2];
      loc_trifaces[3][1][0] = wedge_topo->side[4].node[0];
      loc_trifaces[3][1][1] = wedge_topo->side[4].node[1];
      loc_trifaces[3][1][2] = wedge_topo->side[4].node[2];

      // find max valence 
      unsigned vmaxIndx = 0;
      unsigned vmax = valences[0];
      bool all3 = true;
      for (unsigned iv = 1; iv < 6; iv++)
        {
          if (valences[iv] > vmax)
            {
              if (valences[iv] != 3)
                all3 = false;
              vmax = valences[iv];
              vmaxIndx = iv;
            }
        }

      if (vmax == 3)
        {
          /// Rare case where all valences are 3 (each face is broken in same twisting direction - this is the classic 
          /// "un-tetrahedralizable" configuration, the Schonhardt prism) - in this case, we have to add a Steiner point.
          /// The point can be added along one face's diagonal midpoint, but this introduces a need to investigate neighbors,
          /// so, we simply choose to create more tets by using the centroid as the Steiner point.
          /// (cf. http://www.ams.org/journals/spmj/2005-16-04/S1061-0022-05-00872-1/S1061-0022-05-00872-1.pdf
          /// St. Petersburg Math. J. Tom. 16 (2004), vyp. 4	Vol. 16 (2005), No. 4, Pages 673–690 S 1061-0022(05)00872-1 
          /// Article electronically published on June 24, 2005
          /// REGULAR TRIANGULATIONS AND STEINER POINTS
          /// M. YU. ZVAGELSKI I, A. V. PROSKURNIKOV, AND YU. R. ROMANOVSKI I
          /// )

          //assert(0);
          assert(all3);
          if (1)
            {
              std::ostringstream msg;
              msg << "shouldn't get here, but if we do, let's exit for now " << 1 << " " << 1.e-10 << "\n";
              throw std::runtime_error( msg.str() );
            }

          //exit(1);
          boost::array<double,3> centroid = {{0,0,0}};
          for (unsigned iv = 0; iv < 6; iv++)
            {
              centroid[0] += m_node_coords[elem[iv]][0]/6.0;
              centroid[1] += m_node_coords[elem[iv]][0]/6.0;
              centroid[2] += m_node_coords[elem[iv]][0]/6.0;
            }
          m_node_coords.push_back(centroid);
          unsigned newNodeIndex = m_node_coords.size() - 1;

          for (unsigned iPseudoQuadFace = 0; iPseudoQuadFace < 4; iPseudoQuadFace++)
            {
              for (unsigned kTriOnQuadFace = 0; kTriOnQuadFace < 2; kTriOnQuadFace++)
                {
                  m_elems[shards_Tetrahedron_4].push_back(newNodeIndex);
                  for (unsigned jTriNodeIndex = 0; jTriNodeIndex < 3; jTriNodeIndex++)
                    {
                      m_elems[shards_Tetrahedron_4].push_back(elem[loc_trifaces[iPseudoQuadFace][kTriOnQuadFace][jTriNodeIndex]]);
                    }
                }
            }
        }
      else
        {
          /// normal case - connect max valence node to other faces without that node (should be 3 tets always)
          unsigned count=0;
          for (unsigned iPseudoQuadFace = 0; iPseudoQuadFace < 4; iPseudoQuadFace++)
            {
              for (unsigned kTriOnQuadFace = 0; kTriOnQuadFace < 2; kTriOnQuadFace++)
                {
                  bool isEqual = false;
                  for (unsigned jTriNodeIndex = 0; jTriNodeIndex < 3; jTriNodeIndex++)
                    {
                      if (vmaxIndx == loc_trifaces[iPseudoQuadFace][kTriOnQuadFace][jTriNodeIndex])
                        {
                          isEqual = true;
                          break;
                        }
                    }
                  if (not isEqual)
                    {
                      ++count;
                      m_elems[shards_Tetrahedron_4].push_back(elem[vmaxIndx]);
                      for (unsigned jTriNodeIndex = 0; jTriNodeIndex < 3; jTriNodeIndex++)
                        {
                          m_elems[shards_Tetrahedron_4].push_back(elem[loc_trifaces[iPseudoQuadFace][kTriOnQuadFace][jTriNodeIndex]]);
                        }
                    }
                }
            }
          assert(count == 3);

        }
      if (m_deleteAfterBreak)
        {
          VectorOfInt::iterator pos = m_elems[shards_Wedge_6].begin()+(elemIndex * m_elemInfo[shards_Wedge_6].vertex_count);
          m_elems[shards_Wedge_6].erase(pos, pos + m_elemInfo[shards_Wedge_6].vertex_count);
        }

    }

    /// Break a hex element into 6 (common case) or 12 (rare case) tets using a constructive algorithm.
    /// Note: the 5-tet case cannot be constructed with this algorithm - a table lookup scheme would be more efficient
    ///   and allow for the 5-tet case
    template<> void SweepMesher::breakElement<shards_Hexahedron_8, shards_Tetrahedron_4>(unsigned elemIndex)
    {
      static unsigned loc_trifaces[6][2][3];  // iQuadFace, kTriOnQuadFace, jTriNodeIndex
      static unsigned loc_qfaces[6][4]; // iHexFaceOrd, iFaceNodeOrd
      unsigned* elem = &m_elems[shards_Hexahedron_8][elemIndex * m_elemInfo[shards_Hexahedron_8].vertex_count];

      // build quad face break patterns - use simple min node index as first node
      const CellTopologyData *hex_topo = shards::getCellTopologyData<shards::Hexahedron<8> >();
      unsigned numFaces = hex_topo->side_count; // should be 6
      assert(6 == numFaces);

      // node valences
      unsigned valences[8]={0,0,0,0,0,0,0,0};
      for (unsigned iHexFaceOrd = 0; iHexFaceOrd < 6; iHexFaceOrd++)
        {
          const CellTopologyData_Subcell& face = hex_topo->side[iHexFaceOrd];
          unsigned iqf = face.node[0];
          unsigned globalIqf = elem[iqf];
          unsigned minVal = globalIqf;
          unsigned indxMinVal = 0;
          for (unsigned iFaceNodeOrd=1; iFaceNodeOrd < 4; iFaceNodeOrd++)
            {
              iqf = face.node[iFaceNodeOrd];
              globalIqf = elem[iqf];
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
              //qfaces[iHexFaceOrd][iFaceNodeOrd] = elem[face.node[jFaceNodeOrd]];
              loc_qfaces[iHexFaceOrd][iFaceNodeOrd] = face.node[jFaceNodeOrd];
            }
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
              msg << "shouldn't get here, but if we do, let's exit for now \n";
              throw std::runtime_error( msg.str() );
            }

          //exit(1);
          boost::array<double,3> centroid = {{0,0,0}};
          for (unsigned iv = 0; iv < 8; iv++)
            {
              centroid[0] += m_node_coords[elem[iv]][0]/8.0;
              centroid[1] += m_node_coords[elem[iv]][0]/8.0;
              centroid[2] += m_node_coords[elem[iv]][0]/8.0;
            }
          m_node_coords.push_back(centroid);
          unsigned newNodeIndex = m_node_coords.size() - 1;

          for (unsigned iQuadFace = 0; iQuadFace < 6; iQuadFace++)
            {
              for (unsigned kTriOnQuadFace = 0; kTriOnQuadFace < 2; kTriOnQuadFace++)
                {
                  m_elems[shards_Tetrahedron_4].push_back(newNodeIndex);
                  for (unsigned jTriNodeIndex = 0; jTriNodeIndex < 3; jTriNodeIndex++)
                    {
                      m_elems[shards_Tetrahedron_4].push_back(elem[loc_trifaces[iQuadFace][kTriOnQuadFace][jTriNodeIndex]]);
                    }
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
                      ++count;
                      m_elems[shards_Tetrahedron_4].push_back(elem[vmaxIndx]);
                      for (unsigned jTriNodeIndex = 0; jTriNodeIndex < 3; jTriNodeIndex++)
                        {
                          m_elems[shards_Tetrahedron_4].push_back(elem[loc_trifaces[iQuadFace][kTriOnQuadFace][jTriNodeIndex]]);
                        }
                    }
                }
            }
          assert(count == 6);

        }
      if (m_deleteAfterBreak)
        {
          VectorOfInt::iterator pos = m_elems[shards_Hexahedron_8].begin()+(elemIndex * m_elemInfo[shards_Hexahedron_8].vertex_count);
          m_elems[shards_Hexahedron_8].erase(pos, pos + m_elemInfo[shards_Hexahedron_8].vertex_count);
        }
    }

    /// based on UseCase_3 in stk_mesh/use_cases - creates nodes and elements in stk::mesh database
    void SweepMesher::stkMeshCreate(stk::ParallelMachine& comm)
    {
      stkMeshCreateMetaNoCommit(comm);
      m_metaData->commit();
      stkMeshCreateBulkAfterMetaCommit(comm);
    }

    static std::vector<std::string> get_entity_rank_names(unsigned dim)
    {
      std::vector<std::string> names = stk::mesh::fem::entity_rank_names(dim);
#if PERCEPT_USE_FAMILY_TREE
      names.push_back("FAMILY_TREE");
#endif
      return names;
    }

    void SweepMesher::stkMeshCreateMetaNoCommit(stk::ParallelMachine& comm)
    {
      //m_metaData = new stk::mesh::fem::FEMMetaData(3); //  stk::mesh::fem::fem_entity_rank_names() );  // FAMILY_TREE search
      m_metaData = new stk::mesh::fem::FEMMetaData(3, get_entity_rank_names(3u) ); //  stk::mesh::fem::fem_entity_rank_names() );
      //m_metaData = & stk::mesh::fem::FEMMetaData::get_meta_data(*m_metaData);
      m_bulkData = new stk::mesh::BulkData( stk::mesh::fem::FEMMetaData::get_meta_data(*m_metaData) , comm );
      m_parts.resize(NUM_ELEM_TYPES);

      for (unsigned ieletype = 0; ieletype < NUM_ELEM_TYPES; ieletype++)
        {
          if (m_elems[ieletype].size() > 0)
            {
              m_parts[ieletype] = &(m_metaData->declare_part( std::string("block_").append(std::string(m_elemInfo[ieletype].name)) , m_metaData->element_rank() ));
            }
        }
      m_coordinates_field = &m_metaData->declare_field< VectorFieldType >( "coordinates" );

      m_element_node_coordinates_field = &m_metaData->declare_field< ElementNodePointerFieldType >( "elem_node_coord" );

      // set cell topology
      for (unsigned ieletype = 0; ieletype < NUM_ELEM_TYPES; ieletype++)
        {
          if (m_elems[ieletype].size() > 0)
            {
              m_elemInfo[ieletype].setCellTopoFptr(*m_parts[ieletype]);
              stk::io::put_io_part_attribute(*m_parts[ieletype]);
            }
        }

      // Field restrictions:
      stk::mesh::Part & universal = m_metaData->universal_part();

      put_field( *m_coordinates_field , stk::mesh::fem::FEMMetaData::NODE_RANK , universal );
  
      m_metaData->declare_field_relation(
                                         *m_element_node_coordinates_field ,
                                         stk::mesh::fem::get_element_node_stencil(3) ,
                                         *m_coordinates_field 
                                         );

      for (unsigned ieletype = 0; ieletype < NUM_ELEM_TYPES; ieletype++)
        {
          if (m_elems[ieletype].size() > 0)
            {
#if 0
              std::cout << "shards::Hexahedron<> ::node_count = " << shards::Hexahedron<> ::node_count 
                        << " " << m_elemInfo[ieletype].node_count << std::endl;
              std::cout << "shards::Wedge<> ::node_count = " << shards::Wedge<> ::node_count 
                        << " " << m_elemInfo[ieletype].node_count << std::endl;
#endif
              put_field( *m_element_node_coordinates_field, m_metaData->element_rank(), *m_parts[ieletype], m_elemInfo[ieletype].node_count);
            }
        }
    }

    void SweepMesher::stkMeshCreateBulkAfterMetaCommit(stk::ParallelMachine& comm)
    {

      //stk::mesh::BulkData & bulkData = modifiableBulkData();
      stk::mesh::BulkData & bulkData = *m_bulkData;
      bulkData.modification_begin();

      stk::mesh::EntityId elem_id = 1 ;

      for (unsigned ieletype = 0; ieletype < NUM_ELEM_TYPES; ieletype++)
        {
          if (m_elems[ieletype].size() > 0)
            {
              static stk::mesh::EntityId node_ids[27];  // FIXME - do we have more than 27 nodes?
              stk::mesh::Part & part =  *m_parts[ieletype];
              unsigned nodes_per_elem = m_elemInfo[ieletype].vertex_count;
              unsigned numElems = m_elems[ieletype].size()/nodes_per_elem;
                
              //std::cout << " elems[" << m_elemInfo[ieletype].name << "] = " << numElems << std::endl;

              for (unsigned jelem = 0; jelem < numElems; jelem++, ++elem_id)
                {
                  for (unsigned inode = 0; inode < m_elemInfo[ieletype].vertex_count; inode++)
                    {
                      node_ids[inode] = 1+m_elems[ieletype][jelem*m_elemInfo[ieletype].vertex_count + inode];
                    }

                  //std::cout << "elem_id = " << elem_id << std::endl;
                  stk::mesh::fem::declare_element( bulkData , part , elem_id , node_ids );
                }
            }
        }

      unsigned node_count_1 = m_node_coords.size();
      for ( unsigned i = 0 ; i < node_count_1 ; ++i ) {
        stk::mesh::Entity * const node = m_bulkData->get_entity( stk::mesh::fem::FEMMetaData::NODE_RANK , i + 1 );
        double * const coord = field_data( *m_coordinates_field , *node );

        coord[0] = m_node_coords[i][0];
        coord[1] = m_node_coords[i][1];
        coord[2] = m_node_coords[i][2];
      }

      bulkData.modification_end();
    }


    void SweepMesher::dumpSTK()
    {
      const stk::mesh::BulkData & bulkData = *m_bulkData;

      for (unsigned ieletype = 0; ieletype < NUM_ELEM_TYPES; ieletype++)
        {
          unsigned nodes_per_elem = m_elemInfo[ieletype].vertex_count;
          if (m_elems[ieletype].size() > 0 && nodes_per_elem > 0)
            {
              //std::cout << nodes_per_elem << std::endl;
              unsigned numElems = m_elems[ieletype].size()/nodes_per_elem;

              //std::cout << " elems[" << m_elemInfo[ieletype].name << "] = " << numElems << std::endl;

              stk::mesh::Part & part = *m_parts[ieletype];

              const unsigned expected_num_nodes = m_node_coords.size();
              const unsigned expected_num_elems = numElems;
              //const unsigned expected_num_edges = 0;
              // const unsigned expected_num_faces = 0;
  
              bool result = true;
              std::vector<unsigned> entity_counts;
              stk::mesh::Selector selector(part);
              stk::mesh::count_entities( selector, bulkData , entity_counts );
              if (0) std::cout << "num_nodes = " << entity_counts[0] << " " << expected_num_nodes << std::endl;
              if (0) std::cout << "num_elems = " << entity_counts[m_metaData->element_rank()] << " " << expected_num_elems << std::endl;
                

              if (
                  //(entity_counts[stk::mesh::fem::FEMMetaData::NODE_RANK] != expected_num_nodes) ||  
                  //(entity_counts[Edge] != expected_num_edges) ||
                  //(entity_counts[Face] != expected_num_faces) ||
                  (entity_counts[m_metaData->element_rank()] != expected_num_elems)
                  ) {
                std::cerr<< "Error, the  entity counts are incorrect!" << std::endl;
                result = false;

              }
            }
        }

    }

    void SweepMesher::writeSTKMesh(const char *filename)
    {
      //const std::string out_filename("tp2.e");
      const std::string out_filename(filename);

      std::vector< stk::mesh::Part * > parts;
      for (unsigned ielemType = 0; ielemType < NUM_ELEM_TYPES; ielemType++)
        {
          if (m_elems[ielemType].size())
            {
              parts.push_back(m_parts[ielemType]);
            }
        }

      const stk::ParallelMachine& comm = m_bulkData->parallel();

      Ioss::Init::Initializer init_db;
      stk::io::MeshData mesh;
      stk::io::create_output_mesh(out_filename, comm, *m_bulkData, mesh);
    }

    void SweepMesher::sweep(const VectorOfCoord& path, const VectorOfCoord& dir)
    {
      unsigned npoints = path.size();
      std::vector<Transform *> xforms(npoints-1);
      for (unsigned i = 0; i < npoints-1; i++)
        {
          xforms[i] = new TransformPath(path[i], dir[i], path[i+1], dir[i+1]);
        }
      sweep(xforms);
      for (unsigned i = 0; i < npoints-1; i++)
        {
          delete xforms[i];
        }
    }


  }//namespace percept
}//namespace stk

