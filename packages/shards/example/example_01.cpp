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
#include "Teuchos_Time.hpp"
#include "Teuchos_Assert.hpp"
#include "Shards_CellTopology.hpp"
#include <limits>


using namespace std;
using namespace shards;


int main(int argc, char *argv[]) {
  std::cout \
  << "===============================================================================\n" \
  << "|                                                                             |\n" \
  << "|     Example use of the Shards package:                                      |\n" \
  << "|                                                                             |\n" \
  << "|    1) Timing tests of different access methods                              |\n" \
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
  double timerResults[4][3];
  for(int i = 0; i < 4; i++){
    for(int j = 0; j < 3; j++){
      timerResults[i][j] = 0.0;
    }
  }
  
  // Define variables used in tests
  unsigned vertexOrd = std::numeric_limits<unsigned>::max();
  unsigned subcNode = std::numeric_limits<unsigned>::max();
  unsigned subcCount = std::numeric_limits<unsigned>::max();
  unsigned cellDim = std::numeric_limits<unsigned>::max();
  
  /*************************************************************************************************
    *                                                                                              *
    *                          Timing of access to subcell key                                     *
    *                                                                                              *
    ***********************************************************************************************/
  
  // Access using CellTopologyData struct
  {
    Teuchos::Time timer("Timer - using cell topology data struct");
    timer.start();
    for(int i = 0; i < numCycles1; i++) {
      cellDim = myCellData.subcell[sbcDim][sbcOrd].topology -> key;    
    }
    timer.stop();
    timerResults[0][0] = timer.totalElapsedTime();
  }
  
  // Access using getCellTopologyData<Traits>() ->
  {
    Teuchos::Time timer("Timer - using getCellTopologyData<Traits>() ->");
    timer.start();
    for(int i = 0; i < numCycles1; i++) {
      cellDim = getCellTopologyData< Hexahedron<8> >() -> subcell[sbcDim][sbcOrd].topology -> key;
    }
    timer.stop();
    timerResults[0][1] = timer.totalElapsedTime();
  }
  
  // Access using CellTopology object
  {
    Teuchos::Time timer("Timer - using CellTopology object");
    timer.start();
    for(int i = 0; i < numCycles1; i++) {
      cellDim = myCell.getKey(sbcDim, sbcOrd);
    }
    timer.stop();
    timerResults[0][2] = timer.totalElapsedTime();
  }
  
  
  
  /*************************************************************************************************
    *                                                                                              *
    *                       Timing of access to number of subcells                                 *
    *                                                                                              *
    ***********************************************************************************************/
  
  {
    Teuchos::Time timer("Timer - using cell topology data struct");
    int cellDim = myCellData.dimension;
    timer.start();
    for(int i = 0; i < numCycles1; i++) {
      for(int dim = 0; dim < cellDim; dim++){
        subcCount = myCellData.subcell_count[dim];   
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
        subcCount = getCellTopologyData<Hexahedron<8> >() -> subcell_count[dim];   
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
        subcCount = myCell.getSubcellCount(dim);   
      }
    }
    timer.stop();
    timerResults[1][2] = timer.totalElapsedTime();
  }
  
  
  /*************************************************************************************************
    *                                                                                              *
    *                        Timing of access to all subcell nodes                                 *
    *                                                                                              *
    ***********************************************************************************************/
  
  {
    Teuchos::Time timer("Timer - using cell topology data struct");
    timer.start();
    unsigned numNodes = myCellData.subcell[sbcDim][sbcOrd].topology -> node_count;
    for(int i = 0; i < numCycles1; i++) {
      for(unsigned n = 0; n < numNodes; n++){
        subcNode = myCellData.subcell[sbcDim][sbcOrd].node[n];   
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
        subcNode = getCellTopologyData<Hexahedron<8> >() -> subcell[sbcDim][sbcOrd].node[n];   
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
        subcNode = myCell.getNodeMap(sbcDim, sbcOrd, n);   
      }
    }
    timer.stop();
    timerResults[2][2] = timer.totalElapsedTime();
  }

  
  
  /*************************************************************************************************
    *                                                                                              *
    *                      Timing of access to vertex ordinals of all subcells                     *
    *                                                                                              *
    ***********************************************************************************************/
  
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
            vertexOrd = myCellData.subcell[dim][ord].node[i]; 
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
            vertexOrd = getCellTopologyData<Hexahedron<8> >() -> subcell[dim][ord].node[i]; 
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
            vertexOrd = myCell.getNodeMap(dim, ord, i); 
          }
        }
      }
    }
    timer.stop();
    timerResults[3][2] = timer.totalElapsedTime();
  }

  // fix compiler warnings: set but unused
  TEUCHOS_ASSERT(vertexOrd < std::numeric_limits<unsigned>::max());
  TEUCHOS_ASSERT(subcNode < std::numeric_limits<unsigned>::max());
  TEUCHOS_ASSERT(subcCount < std::numeric_limits<unsigned>::max());
  TEUCHOS_ASSERT(cellDim < std::numeric_limits<unsigned>::max());

  std::cout \
    << "===============================================================================\n"      
    << "   TEST/ACCESS          |  CellTopoData   | getCellTopoData | CellTopo object |\n"
    << "===============================================================================\n"      
    << " cell dimension         | " 
    << std::setw(13) << timerResults[0][0] << "   | " 
    << std::setw(13) << timerResults[0][1] << "   | "  
    << std::setw(13) << timerResults[0][2] << "   |\n"  
    << " number of subcells     | "
    << std::setw(13) << timerResults[1][0] << "   | " 
    << std::setw(13) << timerResults[1][1] << "   | "  
    << std::setw(13) << timerResults[1][2] << "   |\n"  
    << " all subcell nodes      | "
    << std::setw(13) << timerResults[2][0] << "   | " 
    << std::setw(13) << timerResults[2][1] << "   | "  
    << std::setw(13) << timerResults[2][2] << "   |\n"
    << " subcell vertex ordinal | "
    << std::setw(13) << timerResults[3][0] << "   | " 
    << std::setw(13) << timerResults[3][1] << "   | "  
    << std::setw(13) << timerResults[3][2] << "   |\n"
    << "===============================================================================\n\n";      
  
  return 0;
}


