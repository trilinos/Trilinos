// @HEADER
// ************************************************************************
//
//                           shards Package
//                 Copyright (2007) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Pavel Bochev      (pbboche@sandia.gov)
//                    H. Carter Edwards (hcedwar@sandia.gov)
//                    Denis Ridzal      (dridzal@sandia.gov).
//
// ************************************************************************
// @HEADER
 

/** \file
\brief  Examples of cell topology access using shards::CellTopology class.
\author Created by P. Bochev, H. Carter Edwards and D. Ridzal
*/
#include "Teuchos_Time.hpp"

#include "Shards_CellTopology.hpp"
#include <Shards_CellTopologyTraits.hpp>
#include <Shards_BasicTopologies.hpp>

using namespace std;
using namespace shards;


int main(int argc, char *argv[]) {
  std::cout \
  << "===============================================================================\n" \
  << "|                                                                             |\n" \
  << "|     Example use of the CellTopologyData struct and CellTopology class       |\n" \
  << "|                                                                             |\n" \
  << "|    1) Timing tests of different access methods                              |\n" \
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
  
  // Confuse compiler which cell we will use: 
   const CellTopologyData & myCellData = 
    (argc < 100) ? *getCellTopologyData<Hexahedron<27> >() : *getCellTopologyData<Triangle<3> >();
  
  // Create an instance of CellTopology using topology data
  CellTopology myCell(&myCellData);
  
  // Same for an index to access vertices in one of the tests that measures vertex_count access
  int vertIndex = (argc < 100) ? 0 : -1;
    
 // Default arguments
  int sbcDim = 2;
  int sbcOrd = 0;
  int numCycles1 = 100000;
  int numCycles2 = 100000;
  
  // command line arguments: sbcDim, sbcOrd., numCycles1, numCycles2
  if(argc == 2) {
    sbcDim = atoi(argv[1]);
  }
  if(argc == 3) {
    sbcDim = atoi(argv[1]);
    sbcOrd = atoi(argv[2]);
  }
  if(argc == 4) {
    sbcDim = atoi(argv[1]);
    sbcOrd = atoi(argv[2]);
    numCycles1 = atoi(argv[3]);
  }
  if(argc == 5) {
    sbcDim = atoi(argv[1]);
    sbcOrd = atoi(argv[2]);
    numCycles1 = atoi(argv[3]);
    numCycles2 = atoi(argv[4]);
  }
  
  std::cout \
    << "===============================================================================\n"      
    << "| EXAMPLE 1: Timing of different access methods to cell topology              |\n"      
    << "===============================================================================\n"  
    << "\t  argc = " << argc <<"\n"
    << "\t  subcell dimension  = " << sbcDim <<"\n"
    << "\t  subcell ordinal    = " << sbcOrd <<"\n"
    << "\t  num cycles 1       = " << numCycles1 <<"\n"
    << "\t  num cycles 2       = " << numCycles2 <<"\n\n";
  
  // Initialize timer array
  double timerResults[4][4];
  for(int i = 0; i < 4; i++){
    for(int j = 0; j < 4; j++){
      timerResults[i][j] = 0.0;
    }
  }
  
  //////////////////////////////////////////////////////////////////////////////////////////////////
  //                                                                                              //
  //                          Timing of access to subcell key                                     //
  //                                                                                              //
  //////////////////////////////////////////////////////////////////////////////////////////////////
  
  // Access using CellTopologyData struct
  {
    Teuchos::Time timer("Timer - using cell topology data struct");
    timer.start();
    for(int i = 0; i < numCycles1; i++) {
      int cellDim = myCellData.subcell[sbcDim][sbcOrd].topology -> key;    
    }
    timer.stop();
    timerResults[0][0] = timer.totalElapsedTime();
  }
  
  // Access using getCellTopologyData<Traits>() ->
  {
    Teuchos::Time timer("Timer - using getCellTopologyData<Traits>() ->");
    timer.start();
    for(int i = 0; i < numCycles1; i++) {
      unsigned cellDim = getCellTopologyData< Hexahedron<8> >() -> subcell[sbcDim][sbcOrd].topology -> key;
    }
    timer.stop();
    timerResults[0][1] = timer.totalElapsedTime();
  }
  
  // Access using CellTopology object
  {
    Teuchos::Time timer("Timer - using CellTopology object");
    timer.start();
    for(int i = 0; i < numCycles1; i++) {
      unsigned cellDim = myCell.getKey(sbcDim, sbcOrd);
    }
    timer.stop();
    timerResults[0][2] = timer.totalElapsedTime();
  }
  
  // Access using static methods of CellTopology
  {
    Teuchos::Time timer("Timer - using CellTopology static method");
    /*
    timer.start();
    for(int i = 0; i < numCycles1; i++) {
      int cellDim = CellTopology::getKey(&myCellData, sbcDim, sbcOrd);
    }
    timer.stop();
    timerResults[0][3] = timer.totalElapsedTime();
     */
    timerResults[0][3] = -1.0;
  }
  
  
  //////////////////////////////////////////////////////////////////////////////////////////////////
  //                                                                                              //
  //                       Timing of access to number of subcells                                 //
  //                                                                                              //
  //////////////////////////////////////////////////////////////////////////////////////////////////
  
  {
    Teuchos::Time timer("Timer - using cell topology data struct");
    int cellDim = myCellData.dimension;
    timer.start();
    for(int i = 0; i < numCycles1; i++) {
      for(int dim = 0; dim < cellDim; dim++){
        unsigned subcCount = myCellData.subcell_count[dim];   
      }
    }
    timer.stop();
    timerResults[1][0] = timer.totalElapsedTime();
  }
  
  {
    Teuchos::Time timer("Timer -  using getCellTopologyData<Traits>() ->");
    int cellDim = getCellTopologyData<Hexahedron<8> >() -> dimension;
    timer.start();
    for(int i = 0; i < numCycles1; i++) {
      for(int dim = 0; dim < cellDim; dim++){
        unsigned subcCount = getCellTopologyData<Hexahedron<8> >() -> subcell_count[dim];   
      }
    }
    timer.stop();
    timerResults[1][1] = timer.totalElapsedTime();
  }
  
  {
    Teuchos::Time timer("Timer - using CellTopology object");
    int cellDim = myCell.getDimension();
    timer.start();
    for(int i = 0; i < numCycles1; i++) {
      for(int dim = 0; dim < cellDim; dim++){
        unsigned subcCount = myCell.getSubcellCount(dim);   
      }
    }
    timer.stop();
    timerResults[1][2] = timer.totalElapsedTime();
  }
  
  {
    Teuchos::Time timer("Timer - using CellTopology static method");
    /*
    int cellDim = myCellData.dimension;
    timer.start();
    for(int i = 0; i < numCycles1; i++) {
      for(int dim = 0; dim < cellDim; dim++){
        int subcCount = CellTopology::getSubcellCount(&myCellData, dim);   
      }
    }
    timer.stop();
    timerResults[1][3] = timer.totalElapsedTime();
     */
    timerResults[1][3] = -1.0;
  }
  
  
  //////////////////////////////////////////////////////////////////////////////////////////////////
  //                                                                                              //
  //                         Timing of access to all subcell nodes                                //
  //                                                                                              //
  //////////////////////////////////////////////////////////////////////////////////////////////////
  
  {
    Teuchos::Time timer("Timer - using cell topology data struct");
    timer.start();
    unsigned numNodes = myCellData.subcell[sbcDim][sbcOrd].topology -> node_count;
    for(int i = 0; i < numCycles1; i++) {
      for(unsigned n = 0; n < numNodes; n++){
        unsigned subcNode = myCellData.subcell[sbcDim][sbcOrd].node[n];   
      }
    }
    timer.stop();
    timerResults[2][0] = timer.totalElapsedTime();
  }
  
  {
    Teuchos::Time timer("Timer -  using getCellTopologyData<Traits>() ->");
    timer.start();
    unsigned numNodes = getCellTopologyData<Hexahedron<8> >() -> subcell[sbcDim][sbcOrd].topology -> node_count;
    for(int i = 0; i < numCycles1; i++) {
      for(unsigned n = 0; n < numNodes; n++){
        unsigned subcNode = getCellTopologyData<Hexahedron<8> >() -> subcell[sbcDim][sbcOrd].node[n];   
      }
    }
    timer.stop();
    timerResults[2][1] = timer.totalElapsedTime();
  }
  
  {
    Teuchos::Time timer("Timer - using CellTopology object");
    timer.start();
    unsigned numNodes = myCell.getNodeCount(sbcDim, sbcOrd);
    for(int i = 0; i < numCycles1; i++) {
      for(unsigned n = 0; n < numNodes; n++){
        unsigned subcNode = myCell.getNodeMap(sbcDim, sbcOrd, n);   
      }
    }
    timer.stop();
    timerResults[2][2] = timer.totalElapsedTime();
  }

  
  {
    Teuchos::Time timer("Timer - using a static method of CellTopology");
    /*
    timer.start();
    int numNodes = CellTopology::getNodeCount(&myCellData, sbcDim, sbcOrd);
    for(int i = 0; i < numCycles1; i++) {
      for(int n = 0; n < numNodes; n++){
        int subcNode = CellTopology::getNodeMap(&myCellData, sbcDim, sbcOrd, n);   
      }
    }
    timer.stop();
    timerResults[2][3] = timer.totalElapsedTime();
     */
    timerResults[2][3] = -1.0;
  }
  
  
  //////////////////////////////////////////////////////////////////////////////////////////////////
  //                                                                                              //
  //                       Timing of access to vertex ordinals of all subcells                    //
  //                                                                                              //
  //////////////////////////////////////////////////////////////////////////////////////////////////
  
  {
    Teuchos::Time timer("Timer - using cell topology data struct");
    timer.start();
    for(int i = 0; i < numCycles2; i++) {
      
      // subcell dimension
      for(unsigned dim = 0; dim < myCellData.dimension; dim++){
        
        // number of subcells of specified dimension
        for(unsigned ord = 0; ord < myCellData.subcell_count[dim]; ord++) {
          
          // number of vertices in the subcell
          int numVert = myCellData.subcell[dim][ord].topology -> vertex_count;
          for(int i = 0; i < numVert; i++){
            unsigned vertexOrd = myCellData.subcell[dim][ord].node[i]; 
          }
        }
      }
    }
    timer.stop();
    timerResults[3][0] = timer.totalElapsedTime();
  }
  
  {
    Teuchos::Time timer("Timer -  using getCellTopologyData<Traits>() ->");
    timer.start();
    for(int i = 0; i < numCycles2; i++) {
      
      // subcell dimension
      int cellDim = getCellTopologyData<Hexahedron<8> >() -> dimension;
      for(int dim = 0; dim < cellDim; dim++){
        
        // number of subcells of specified dimension
        int subcCount = getCellTopologyData<Hexahedron<8> >() -> subcell_count[dim];
        for(int ord = 0; ord < subcCount; ord++) {
          
          // number of vertices in the subcell
          int numVert = getCellTopologyData<Hexahedron<8> >() -> subcell[dim][ord].topology -> vertex_count;
          for(int i = 0; i < numVert; i++){
            unsigned vertexOrd = getCellTopologyData<Hexahedron<8> >() -> subcell[dim][ord].node[i]; 
          }
        }
      }
    }
    timer.stop();
    timerResults[3][1] = timer.totalElapsedTime();
  }
  
  {
    Teuchos::Time timer("Timer - using CellTopology object");
    timer.start();
    for(int i = 0; i < numCycles2; i++) {

      // subcell dimension
      for(unsigned dim = 0; dim < myCell.getDimension(); dim++){
        
        // number of subcells of specified dimension
        for(unsigned ord = 0; ord < myCell.getSubcellCount(dim); ord++) {
          
          // number of vertices in the subcell
          int numVert = myCell.getVertexCount(dim,ord);
          for(int i = 0; i < numVert; i++){
            unsigned vertexOrd = myCell.getNodeMap(dim, ord, i); 
          }
        }
      }
    }
    timer.stop();
    timerResults[3][2] = timer.totalElapsedTime();
  }

  {
    /*
    Teuchos::Time timer("Timer - using a static method of CellTopology");
    timer.start();
    for(int i = 0; i < numCycles2; i++) {
      
      // subcell dimension
      for(unsigned dim = 0; dim < CellTopology::getDimension(&myCellData); dim++){
        
        // number of subcells of specified dimension
        for(unsigned ord = 0; ord < CellTopology::getSubcellCount(&myCellData, dim); ord++) {
          
          // number of vertices in the subcell
          int numVert = CellTopology::getVertexCount(&myCellData, dim, ord);
          for(int i = 0; i < numVert; i++){
            unsigned vertexOrd = CellTopology::getNodeMap(&myCellData, dim, ord, i); 
          }
        }
      }
    }
    timer.stop();
    timerResults[3][3] = timer.totalElapsedTime();
     */
    timerResults[3][3] = -1.0;
  }

  std::cout \
    << "===========================================================================================\n" 
    << "   TEST/ACCESS    |   CellTopoData  | getCellTopoData |  CellTopo object | CellTopo static|\n"
    << "===========================================================================================\n" 
    << " cell dimension   | " 
    << std::setw(14) << timerResults[0][0] << "  | " 
    << std::setw(14) << timerResults[0][1] << "  | "  
    << std::setw(14) << timerResults[0][2] << "  | "  
    << std::setw(14) << timerResults[0][3] << "  |\n"  
    << " num. subcells    | "
    << std::setw(14) << timerResults[1][0] << "  | " 
    << std::setw(14) << timerResults[1][1] << "  | "  
    << std::setw(14) << timerResults[1][2] << "  | "  
    << std::setw(14) << timerResults[1][3] << "  |\n"  
    << " all subc. nodes  | "
    << std::setw(14) << timerResults[2][0] << "  | " 
    << std::setw(14) << timerResults[2][1] << "  | "  
    << std::setw(14) << timerResults[2][2] << "  | "  
    << std::setw(14) << timerResults[2][3] << "  |\n"
    << " subc. vert. ord. | "
    << std::setw(14) << timerResults[3][0] << "  | " 
    << std::setw(14) << timerResults[3][1] << "  | "  
    << std::setw(14) << timerResults[3][2] << "  | "  
    << std::setw(14) << timerResults[3][3] << "  |\n"
    << "===========================================================================================\n\n";
    
  std::cout \
    << "===============================================================================\n"      
    << "| EXAMPLE 2: Creating CellTopology objects from CellTopologyData structs      |\n"      
    << "===============================================================================\n\n" ; 
  
  // Default ctor creates invalid cell with NULL topology and NULL base topology
  shards::CellTopology emptyCell;
  
  // Creates CellTopology from raw CellTopologyData struct to provide safe access to the latter 
  shards::CellTopology triangle(getCellTopologyData<Triangle<6> > () );
  std::cout << triangle << "\n";
  
  // Available polygon topologies (no extended polygons provided, can skip template parameter)
  shards::CellTopology pentagon(getCellTopologyData<Pentagon<5> >() );
  std::cout << pentagon << "\n";
  
  shards::CellTopology hexagon(getCellTopologyData<Hexagon<6> >() );
  std::cout << hexagon << "\n";
  
  
  
  std::cout \
    << "===============================================================================\n"      
    << "| EXAMPLE 3: Creating custom 1D and 2D CellTopology objects                   |\n"      
    << "===============================================================================\n\n" ; 
  
  
  // Line with 5 nodes
  shards::CellTopology customLine("customLine", 4);
  std::cout << customLine << "\n";
  
  
  // 
  // Pentagon_10: 
  // extended topology version of the basic Pentagon<5> topology in which each edge has 3 nodes
  // 
  
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
  shards::CellTopology Pentagon_10("Pentagon_10",
                                   5,
                                   10,
                                   homog_edges,
                                   homog_edges_node_map,
                                   getCellTopologyData<Pentagon<5> >() );
  
  std::cout << Pentagon_10 
    << "Edge homogeneity = " << Pentagon_10.getSubcellHomogeneity(1) << "\n\n";

  
  
  //
  //  Pentagon_7: 
  //  extended topology version of the basic Pentagon<5> topology in which edges 0 and 2 have
  //  2 nodes and all other edges have 2 nodes
  //
  
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
  
  shards::CellTopology Pentagon_7("Pentagon_7",
                                  5,
                                  7,
                                  inhomog_edges,
                                  inhomog_edges_node_map,
                                  getCellTopologyData<Pentagon<5> >() );
  
  std::cout << Pentagon_7 
    << "Edge homogeneity = " << Pentagon_7.getSubcellHomogeneity(1) << "\n\n";
  

  std::cout \
    << "===============================================================================\n"      
    << "| EXAMPLE 4: Creating custom 3D CellTopology objects                          |\n"      
    << "===============================================================================\n\n" ; 
  
  // Beehive_12
  // prismatic cell with base Hexagon<6>, sides Quad<4> and self-referential base topology
  
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
  
  /*
  for(int i=0; i<36; i++){
    std::cout << beehive_edge_node_map[i] << "  ";
  }
  */
  
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
  
  
  
  shards::CellTopology Beehive_12("Beehive_12",
                                   12,
                                   12,
                                   beehive_edges ,
                                   beehive_edge_node_map ,
                                   beehive_faces ,
                                   beehive_face_node_map);
  
  std::cout << Beehive_12 << "\n"
    << "Edge homogeneity = " << Beehive_12.getSubcellHomogeneity(1) << "\n\n";

  
  
  
  /*
  std::cout 
    << "Triangle: \n" << *shards::getCellTopologyData<Triangle<> >() <<"\n" 
    << "Triangle: \n" << triangle.getTopology() <<"\n" 
    << "Triangle: \n" << *triangle.getTopology() <<"\n" 
    << "Triangle: \n" << triangle <<"\n" 
    << "Custom line : \n " << customLine << "\n"
    << "Pentagon : \n " << *shards::getCellTopologyData<Pentagon<5> >() << "\n";
    //    << "Custom 1-cell: \n" << emptyCell.getBaseTopology() <<"\n" ;
    //<< "Base topology of custom 1-cell: \n" << Line5.getBaseTopology(1,0) -> name << "\n"; 
  */
  
  /**
    // accessing cell topology using compile time traits
   cout << " number 0-subcells = " << Triangle<3>::subcell<0>::count <<"\n";
   cout << " number 1-subcells = " << Triangle<3>::subcell<1>::count <<"\n";
   cout << " number 2-subcells = " << Triangle<3>::subcell<2>::count <<"\n";
   
   // access subcell-1 with ordinal 0 (edge 0)
   cout << " num nodes in edge 0 = " << Triangle<3>::subcell<1,0>::topology::node_count <<"\n";
   cout << " num nodes in edge 0 = " << Triangle<3>::subcell<1,0>::topology::subcell<0>::count <<"\n";
   
   // access node numbers:
   cout << " node 0 of edge  2 = " << Triangle<3>::subcell<1,2,0>::node <<"\n";
   cout << " node 1 of edge  2 = " << Triangle<3>::subcell<1,2,1>::node <<"\n";
   
   
   // runtime access:
   const phdmesh::CellTopology & myTri = *cell_topology<Triangle<3> >();
   const phdmesh::CellTopology & myOtherTri = *cell_topology<Triangle<> >();
   
   
   cout << " topology name   = " <<  myTri.name << "\n";
   cout << " dimension       = " <<  myTri.dimension << "\n";
   cout << " number of nodes = " <<  myTri.subcell_count[0] <<"\n";
   cout << " number of edges = " <<  myTri.subcell_count[1] <<"\n";
   
   cell_topology<Triangle<3> >() -> subcell_count[0];
   
   
   cout << cell_topology<Triangle<> >()[0] << "\n";
   cout << cell_topology<Triangle<3> >()[0] << "\n";
   cout << cell_topology<Triangle<6> >()[0] << "\n\n";
   cout << cell_topology<ShellTriangle<> >()[0] << "\n";
   cout << cell_topology<ShellTriangle<3> >()[0] << "\n\n";
   
   cout << cell_topology<Tetrahedron<> >()[0] << "\n";
   cout << cell_topology<Tetrahedron<4> >()[0] << "\n";
   
   **/
  
  
  return 0;
}




























