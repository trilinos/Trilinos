// @HEADER
// *****************************************************************************
//                Shards : Shared Discretization Tools
//
// Copyright 2008-2011 NTESS and the Shards contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
 

/** \file
    \brief  Examples of cell topology access using shards::CellTopology class.
    \author Created by P. Bochev, H. Carter Edwards and D. Ridzal
*/
#include <iostream>
#include "Shards_CellTopology.hpp"
#include "Shards_CellTopologyManagedData.hpp"


using namespace std;
using namespace shards;


int main(int argc, char *argv[]) {
  std::cout \
  << "===============================================================================\n" \
  << "|                                                                             |\n" \
  << "|     Example use of the Shards package:                                      |\n" \
  << "|                                                                             |\n" \
  << "|    1) Definition of CellTopology objects                                    |\n" \
  << "|    2) Definition of custom cell topologies                                  |\n" \
  << "|                                                                             |\n" \
  << "|  Questions? Contact  Pavel Bochev (pbboche@sandia.gov)                      |\n" \
  << "|                      H. Carter Edwards (hcedwar@sandia.gov)                 |\n" \
  << "|                      Denis Ridzal (dridzal@sandia.gov).                     |\n" \
  << "|                                                                             |\n" \
  << "|  shards's website:   http://trilinos.sandia.gov/packages/shards             |\n" \
  << "|  Trilinos website:   http://trilinos.sandia.gov                             |\n" \
  << "|                                                                             |\n" \
  << "===============================================================================\n\n";
  
    
  std::cout \
    << "===============================================================================\n"      
    << "| EXAMPLE 1: Creating CellTopology objects from CellTopologyData structs      |\n"      
    << "===============================================================================\n\n" ; 
  
  // Creates invalid cell with NULL topology and NULL base topology
  shards::CellTopology emptyCell( NULL );

  
  // Creates CellTopology from raw CellTopologyData struct to provide safe access to the latter 
  shards::CellTopology my_triangle(getCellTopologyData<Triangle<6> > () );
  std::cout << my_triangle << "\n";
  
  shards::CellTopology my_shell_quad(getCellTopologyData<ShellQuadrilateral<8> >() );
  std::cout << my_shell_quad << "\n";
  
  
  // Shards has some predefined polygon topologies:
  shards::CellTopology pentagon(getCellTopologyData<Pentagon<5> >() );
  std::cout << pentagon << "\n";
  
  shards::CellTopology hexagon(getCellTopologyData<Hexagon<6> >() );
  std::cout << hexagon << "\n";
  
  
  
  std::cout \
    << "===============================================================================\n"      
    << "| EXAMPLE 2: Creating custom 1D and 2D CellTopology objects                   |\n"      
    << "===============================================================================\n\n" ; 
  
  
  /*************************************************************************************************
    *  Line with 5 nodes                                                                           *
    ***********************************************************************************************/
  
  shards::CellTopologyManagedData customLineManaged("customLine", 4);
  shards::CellTopology customLine( & customLineManaged );
  std::cout << customLine << "\n";
  
  
  /*************************************************************************************************
    *  Pentagon_10:                                                                                *
    *  extended topology version of the basic Pentagon<5> topology in which each edge has 3 nodes  *
    ***********************************************************************************************/
  
  // Flat array with 5 identical edge topologies
  std::vector<const CellTopologyData * > homog_edges;
  homog_edges.resize(5);
  for(int i = 0; i < 5; i++){
    homog_edges[i] = getCellTopologyData<Line<3> >();
  }
  
  // Flat array with edge to node map: vertices numbered first, mid-edge nodes second
  unsigned nodeList[] = {0, 5, 1,   1, 6, 2,   2, 7, 3,   3, 8, 4,   4, 9, 0}; 
  std::vector<unsigned>  homog_edges_node_map;
  homog_edges_node_map.assign(nodeList, nodeList + 15);

  // Custom 2D cell constructor: provide Pentagon<5> as base topology
  shards::CellTopologyManagedData Pentagon_10_managed("Pentagon_10",
                                   5,
                                   10,
                                   homog_edges,
                                   homog_edges_node_map,
                                   getCellTopologyData<Pentagon<5> >() );
  shards::CellTopology Pentagon_10( & Pentagon_10_managed);
  
  std::cout << Pentagon_10 
    << "Edge homogeneity = " << Pentagon_10.getSubcellHomogeneity(1) << "\n\n";

  
  
  /*************************************************************************************************
    *  Pentagon_7:                                                                                 *
    *  extended topology version of the basic Pentagon<5> topology in which edges 0 and 2 have     *
    *  2 nodes and all other edges have 2 nodes                                                    *
    ***********************************************************************************************/
  
  // Flat array with 5 inhomogeneous edges
  std::vector<const CellTopologyData * > inhomog_edges;
  inhomog_edges.resize(5);
  inhomog_edges[0] = getCellTopologyData<Line<3> >();
  inhomog_edges[1] = getCellTopologyData<Line<2> >();
  inhomog_edges[2] = getCellTopologyData<Line<3> >();
  inhomog_edges[3] = getCellTopologyData<Line<2> >();
  inhomog_edges[4] = getCellTopologyData<Line<2> >();
  
  // Flat array with edge to node map: vertices numbered first, mid-edge nodes second
  unsigned inhomogNodeList[] = {0, 5, 1,   1, 2,   2, 6, 3,   3, 4,   4, 0}; 
  std::vector<unsigned>  inhomog_edges_node_map;
  inhomog_edges_node_map.assign(inhomogNodeList, inhomogNodeList + 12);
  
  shards::CellTopologyManagedData Pentagon_7_managed("Pentagon_7",
                                  5,
                                  7,
                                  inhomog_edges,
                                  inhomog_edges_node_map,
                                  getCellTopologyData<Pentagon<5> >() );
  shards::CellTopology Pentagon_7( & Pentagon_7_managed );
  
  std::cout << Pentagon_7 
    << "Edge homogeneity = " << Pentagon_7.getSubcellHomogeneity(1) << "\n\n";
  

  
  std::cout \
    << "===============================================================================\n"      
    << "| EXAMPLE 3: Creating custom 3D CellTopology objects                          |\n"      
    << "===============================================================================\n\n" ; 
  
  
  
  /*************************************************************************************************
    *  Beehive_12                                                                                  *
    *  prismatic cell with base Hexagon<6>, sides Quad<4> and self-referential base topology       *
    ***********************************************************************************************/
  
  // Array for edge topologies
  std::vector<const CellTopologyData *> beehive_edges;
  beehive_edges.resize(18);
  for(unsigned i = 0; i < 18; i++){
    beehive_edges[i] = getCellTopologyData<Line<2> >(); 
  }
  
  // edges are numbered on the bottom face, top face and on the sides
  unsigned beehive_edge_node_list[] = 
    {
      0, 1,   1, 2,   2, 3,   3,  4,    4,  5,   5, 0,          // bottom face edges
      6, 7,   7, 8,   8, 9,   9, 10,   10, 11,  11, 6,          // top face edges
      0, 6,   1, 7,   2, 8,   3,  9,    4, 10,   5,11           // side edges
    };
  std::vector<unsigned> beehive_edge_node_map;
  beehive_edge_node_map.resize(36);
  beehive_edge_node_map.assign(beehive_edge_node_list, beehive_edge_node_list + 36);
  
  // Array for face topologies
  std::vector<const CellTopologyData *> beehive_faces;
  beehive_faces.resize(8);
  
  // Top and bottom faces are Hexagon<6>
  beehive_faces[0] = getCellTopologyData<Hexagon<6> >();
  beehive_faces[7] = beehive_faces[0];
  
  // All other faces are Quadrilateral<4>
  for(unsigned i = 1; i < 7; i++) {
    beehive_faces[i] = getCellTopologyData<Quadrilateral<4> >();
  }
  
  unsigned beehive_face_node_list[] = 
    {
     0, 1, 2, 3, 4, 5,
     0, 1, 6, 7,   1, 2, 8, 7,   2, 3, 9, 8,  3, 4, 10, 9,  4, 5, 11, 10,  5, 0, 6, 11,  
     6, 7, 8, 9, 10, 11 
    };
  std::vector<unsigned> beehive_face_node_map;
  beehive_face_node_map.resize(36);
  beehive_face_node_map.assign(beehive_face_node_list, beehive_face_node_list + 36);
  
  
  shards::CellTopologyManagedData Beehive_12_managed( "Beehive_12",
                                   12,
                                   12,
                                   beehive_edges ,
                                   beehive_edge_node_map ,
                                   beehive_faces ,
                                   beehive_face_node_map, 0);
  shards::CellTopology Beehive_12( & Beehive_12_managed );
  
  std::cout << Beehive_12 << "\n"
    << "Edge homogeneity = " << Beehive_12.getSubcellHomogeneity(1) << "\n\n";
  
  return 0;
}


