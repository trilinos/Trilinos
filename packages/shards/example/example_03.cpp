// @HEADER
// *****************************************************************************
//                Shards : Shared Discretization Tools
//
// Copyright 2008-2011 NTESS and the Shards contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


/** \file
    \brief  Example of the CellTools class.
    \author Created by P. Bochev, H. Carter Edwards and D. Ridzal
*/
#include <iostream>
#include <iomanip>
#include "Shards_CellTopology.hpp"


using namespace std;
using namespace shards;

/** \brief  Prints the vector with the selected topologies.

    \param  topologies      [in]    - vector containing the selected topologies
    \param  cellType        [in]    - enum for the selected cell type 
    \param  topologyType    [in]    - enum for the selected topology type
  */
void printSelectTopologies(const std::vector<CellTopology>&   topologies,
                     const ECellType                    cellType = ALL_CELLS,
                     const ETopologyType                topologyType = ALL_TOPOLOGIES);




int main(int argc, char *argv[]) {
  std::cout \
  << "===============================================================================\n" \
  << "|                                                                             |\n" \
  << "|     Example use of the Shards package:                                      |\n" \
  << "|                                                                             |\n" \
  << "|    1) Query of the available cell topologies                                |\n" \
  << "|                                                                             |\n" \
  << "|  Questions? Contact  Pavel Bochev (pbboche@sandia.gov)                      |\n" \
  << "|                      H. Carter Edwards (hcedwar@sandia.gov)                 |\n" \
  << "|                      Denis Ridzal (dridzal@sandia.gov).                     |\n" \
  << "|                                                                             |\n" \
  << "|  shards's website:   http://trilinos.sandia.gov/packages/shards             |\n" \
  << "|  Trilinos website:   http://trilinos.sandia.gov                             |\n" \
  << "|                                                                             |\n" \
  << "===============================================================================\n\n";
  
  
  /*************************************************************************************************
    *                                                                                              *
    *  Cell type enums:     ALL_CELLS,       STANDARD_CELL,  NONSTANDARD_CELL                      *
    *  Topology type enums: ALL_TOPOLOGIES,  BASE_TOPOLOGY,  EXTENDED_TOPOLOGY                     *
    *                                                                                              *
    ***********************************************************************************************/
  std::vector<CellTopology> topologies;
  
  std::cout \
    << "===============================================================================\n"\
    << "| EXAMPLE 1: Queries of predefined 0D, 1D, 2D and 3D cell topologies          |\n"\
    << "===============================================================================\n";
  
  for(int cellDim = 0; cellDim < 4; cellDim++){
    
    std::cout << "************ Selected cell dimension = " << cellDim << " ************\n";
    
    // All cells 
    shards::getTopologies(topologies, cellDim);
    printSelectTopologies(topologies);
    
    
    // All standard cells
    shards::getTopologies(topologies, cellDim, STANDARD_CELL);
    printSelectTopologies(topologies,          STANDARD_CELL);
    
    // All standard cells with base topology
    shards::getTopologies(topologies, cellDim, STANDARD_CELL, BASE_TOPOLOGY);
    printSelectTopologies(topologies,          STANDARD_CELL, BASE_TOPOLOGY);
    
    // All standard cells with extended topology
    shards::getTopologies(topologies, cellDim, STANDARD_CELL, EXTENDED_TOPOLOGY);
    printSelectTopologies(topologies,          STANDARD_CELL, EXTENDED_TOPOLOGY);
    
    
    
    // All non-standard cells
    shards::getTopologies(topologies, cellDim, NONSTANDARD_CELL);
    printSelectTopologies(topologies,          NONSTANDARD_CELL);
    
    // All non-standard 0D cells with base topology
    shards::getTopologies(topologies, cellDim, NONSTANDARD_CELL, BASE_TOPOLOGY);
    printSelectTopologies(topologies,          NONSTANDARD_CELL, BASE_TOPOLOGY);
    
    
    // All non-standard cells with extended topology
    shards::getTopologies(topologies, cellDim,  NONSTANDARD_CELL, EXTENDED_TOPOLOGY);
    printSelectTopologies(topologies,           NONSTANDARD_CELL, EXTENDED_TOPOLOGY);
    
    
  }
  
 
  
  std::cout \
    << "===============================================================================\n"\
    << "| EXAMPLE 2: Query of all predefined cell topologies                          |\n"\
    << "===============================================================================\n";
  
  // This query uses default argument values for all input arguments:
  shards::getTopologies(topologies);
  printSelectTopologies(topologies);
  
return 0;
}


/***************************************************************************************************
  *                                                                                                *
  *    Helper function                                                                             *
  *                                                                                                *
  *************************************************************************************************/
void printSelectTopologies(const std::vector<CellTopology>&   topologies,
                     const ECellType                          cellType,
                     const ETopologyType                      topologyType)
{
  
  std::cout << "List of " << shards::ECellTypeToString(cellType) << " ";
  
  // If topologies contains all 33 predefined topologies do not print cell dimension
  if( topologies.size() == 33 ) {
    std::cout << "cells and ";
  }
  else if ( ! topologies.empty() ) {
    std::cout << topologies[0].getDimension() << "D cells and ";
 
  }
  std::cout << shards::ETopologyTypeToString(topologyType) << " topology types  (total of " 
  << topologies.size() << " cells)\n";

  std::cout << "-------------------------------------------------------------------------------\n";
  std::cout << setw(25) << " Cell Topology " << setw(25) << " Base topology" << setw(30) << "|\n";
  std::cout << "-------------------------------------------------------------------------------\n";
  
  for(unsigned i = 0; i < topologies.size(); i++){
    std::cout << setw(25) << topologies[i].getName() << setw(25) << topologies[i].getBaseName() << "\n"; 
  }
  std::cout << "===============================================================================\n\n";
}  

























