// @HEADER
// ************************************************************************
//
//                           Intrepid Package
//                 Copyright (2007) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Pavel Bochev  (pbboche@sandia.gov)
//                    Denis Ridzal  (dridzal@sandia.gov), or
//                    Kara Peterson (kjpeter@sandia.gov)
//
// ************************************************************************
// @HEADER


/** \file
    \brief  Test of the CellTools class.
    \author Created by P. Bochev, D. Ridzal and K. Peterson
*/
#include "Intrepid_CellTools.hpp"
#include "Intrepid_FieldContainer.hpp"
#include "Intrepid_DefaultCubatureFactory.hpp"
#include "Shards_CellTopology.hpp"

#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_ScalarTraits.hpp"

using namespace std;
using namespace Intrepid;
using namespace shards;
  
#define INTREPID_TEST_COMMAND( S , throwCounter, nException )                                                              \
{                                                                                                                          \
  ++nException;                                                                                                            \
    try {                                                                                                                    \
      S ;                                                                                                                    \
    }                                                                                                                        \
    catch (std::logic_error err) {                                                                                           \
      ++throwCounter;                                                                                                      \
        *outStream << "Expected Error " << nException << " -------------------------------------------------------------\n"; \
          *outStream << err.what() << '\n';                                                                                    \
            *outStream << "-------------------------------------------------------------------------------" << "\n\n";           \
    };                                                                                                                       \
}



int main(int argc, char *argv[]) {

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  typedef CellTools<double>       CellTools;
  typedef shards::CellTopology    CellTopology;
  
  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  
  Teuchos::RCP<std::ostream> outStream;
  Teuchos::oblackholestream bhs; // outputs nothing
  
  if (iprint > 0)
    outStream = Teuchos::rcp(&std::cout, false);
  else
    outStream = Teuchos::rcp(&bhs, false);
  
  // Save the format state of the original std::cout.
  Teuchos::oblackholestream oldFormatState;
  oldFormatState.copyfmt(std::cout);
  
  *outStream \
    << "===============================================================================\n" \
    << "|                                                                             |\n" \
    << "|                              Unit Test CellTools                            |\n" \
    << "|                                                                             |\n" \
    << "|     1) Mapping to and from reference cells with base and extended topologies|\n" \
    << "|        using default initial guesses when computing the inverse F^{-1}      |\n" \
    << "|     2) Repeat all tests from 1) using user-defined initial guess for F^{-1} |\n" \
    << "|     3) Exception testing                                                    |\n" \
    << "|                                                                             |\n" \
    << "|  Questions? Contact  Pavel Bochev (pbboche@sandia.gov)                      |\n" \
    << "|                      Denis Ridzal (dridzal@sandia.gov), or                  |\n" \
    << "|                      Kara Peterson (kjpeter@sandia.gov)                     |\n" \
    << "|                                                                             |\n" \
    << "|  Intrepid's website: http://trilinos.sandia.gov/packages/intrepid           |\n" \
    << "|  Trilinos website:   http://trilinos.sandia.gov                             |\n" \
    << "|                                                                             |\n" \
    << "===============================================================================\n";
  
  int errorFlag  = 0;

  // Collect all supported cell topologies
  std::vector<shards::CellTopology> supportedTopologies;
  supportedTopologies.push_back(shards::getCellTopologyData<Triangle<3> >() );
  supportedTopologies.push_back(shards::getCellTopologyData<Triangle<6> >() );
  supportedTopologies.push_back(shards::getCellTopologyData<Quadrilateral<4> >() );
  supportedTopologies.push_back(shards::getCellTopologyData<Quadrilateral<9> >() );
  supportedTopologies.push_back(shards::getCellTopologyData<Tetrahedron<4> >() );
  supportedTopologies.push_back(shards::getCellTopologyData<Tetrahedron<10> >() );
  supportedTopologies.push_back(shards::getCellTopologyData<Hexahedron<8> >() );
  supportedTopologies.push_back(shards::getCellTopologyData<Hexahedron<27> >() );
  supportedTopologies.push_back(shards::getCellTopologyData<Wedge<6> >() );
  supportedTopologies.push_back(shards::getCellTopologyData<Wedge<18> >() );
  
  // Declare iterator to loop over the cell topologies
  std::vector<shards::CellTopology>::iterator topo_iterator;

  // Test 1 scope
  try{

    *outStream \
    << "\n"
    << "===============================================================================\n"\
    << "| Test 1: computing F(x) and F^{-1}(x) using default initial guesses.         |\n"\
    << "===============================================================================\n\n";
    /*
     *  Test summary:
     *
     *    A reference point set is mapped to physical frame and then back to reference frame.
     *    Test passes if the final set of points matches the first set of points. The cell workset
     *    is generated by perturbing randomly the cellWorkset of a reference cell with the specified 
     *    cell topology. 
     *
     */
    // Declare arrays for cell workset and point sets. We will have 10 cells in the wset. and 10 pts per pt. set
    FieldContainer<double> cellWorkset;                 // physical cell workset
    FieldContainer<double> refPoints;                   // reference point set(s) 
    FieldContainer<double> physPoints;                  // physical point set(s)
    FieldContainer<double> controlPoints;               // preimages: physical points mapped back to ref. frame
    
    // We will use cubature factory to get some points on the reference cells. Declare necessary arrays
    DefaultCubatureFactory<double>  cubFactory;   
    FieldContainer<double> cubPoints;
    FieldContainer<double> cubWeights;

    // Initialize number of cells in the cell workset
    int numCells  = 10;
    

    // Loop over cell topologies, make cell workset for each one by perturbing the cellWorkset & test methods
    for(topo_iterator = supportedTopologies.begin(); topo_iterator != supportedTopologies.end(); ++topo_iterator){
      
      // 1.   Define a single reference point set using cubature factory with order 4 cubature
      Teuchos::RCP<Cubature<double> > cellCubature = cubFactory.create( (*topo_iterator), 4); 
      int cubDim = cellCubature -> getDimension();
      int numPts = cellCubature -> getNumPoints();
      cubPoints.resize(numPts, cubDim);
      cubWeights.resize(numPts);
      cellCubature -> getCubature(cubPoints, cubWeights);
             
      // 2.   Define a cell workset by perturbing the cellWorkset of the reference cell with the specified topology
      // 2.1  Resize dimensions of the rank-3 (C,N,D) cell workset array for the current topology
      int numNodes = (*topo_iterator).getNodeCount();
      int cellDim  = (*topo_iterator).getDimension();
      cellWorkset.resize(numCells, numNodes, cellDim);
      
      // 2.2  Copy cellWorkset of the reference cell with the same topology to temp rank-2 (N,D) array
      FieldContainer<double> refCellNodes(numNodes, cellDim );
      CellTools::getReferenceSubcellNodes(refCellNodes, cellDim, 0, (*topo_iterator) );
      
      // 2.3  Create randomly perturbed version of the reference cell and save in the cell workset array
      for(int cellOrd = 0; cellOrd < numCells; cellOrd++){
        
        // Move vertices +/-0.125 along their axes. Gives nondegenerate cells for base and extended topologies 
        for(int nodeOrd = 0; nodeOrd < numNodes; nodeOrd++){
          for(int d = 0; d < cellDim; d++){
            double delta = Teuchos::ScalarTraits<double>::random()/16.0;
            cellWorkset(cellOrd, nodeOrd, d) = refCellNodes(nodeOrd, d) + delta;
          } // d
        }// nodeOrd           
      }// cellOrd
      /* 
       * 3.1 Test 1: single point set to single physical cell: map ref. point set in rank-2 (P,D) array
       *      to a physical point set in rank-2 (P,D) array for a specified cell ordinal. Use the cub.
       *      points array for this test. Resize physPoints and controlPoints to rank-2 (P,D) arrays.
       */
      physPoints.resize(numPts, cubDim);
      controlPoints.resize(numPts, cubDim);
      
      *outStream 
        << " Mapping a set of " << numPts << " points to one cell in a workset of " << numCells << " " 
        << (*topo_iterator).getName() << " cells. \n"; 
      
      for(int cellOrd = 0; cellOrd < numCells; cellOrd++){
        
        // Forward map:: requires cell ordinal
        CellTools::mapToPhysicalFrame(physPoints, cubPoints, cellWorkset, (*topo_iterator), cellOrd);
        // Inverse map: requires cell ordinal
        CellTools::mapToReferenceFrame(controlPoints, physPoints, cellWorkset, (*topo_iterator), cellOrd);

        // Points in controlPoints should match the originals in cubPoints up to a tolerance
        for(int pt = 0; pt < numPts; pt++){
          for(int d = 0; d < cellDim; d++){
            
            if( abs( controlPoints(pt, d) - cubPoints(pt, d) ) > 100.0*INTREPID_TOL ){
              errorFlag++;
              *outStream
                << std::setw(70) << "^^^^----FAILURE!" << "\n"
                << " Mapping a single point set to a single physical cell in a workset failed for: \n"
                << "                    Cell Topology = " << (*topo_iterator).getName() << "\n"
                << " Physical cell ordinal in workset = " << cellOrd << "\n"
                << "          Reference point ordinal = " << setprecision(12) << pt << "\n"
                << "    At reference point coordinate = " << setprecision(12) << d << "\n"
                << "                   Original value = " << cubPoints(pt, d) << "\n"
                << "                     F^{-1}F(P_d) = " << controlPoints(pt, d) <<"\n";
            }
          }// d
        }// pt
      }// cellOrd
      /* 
       * 3.2  Test 2: single point set to multiple physical cells: map ref. point set in rank-2 (P,D) array
       *      to a physical point set in rank-3 (C, P,D) array for all cell ordinals. Use the cub.
       *      points array for this test. Resize physPoints and controlPoints to rank-3 (C,P,D) arrays.
       */
      physPoints.clear(); 
      controlPoints.clear();
      physPoints.resize(numCells, numPts, cubDim);
      controlPoints.resize(numCells, numPts, cubDim);
      
      *outStream 
        << " Mapping a set of " << numPts << " points to all cells in workset of " << numCells << " " 
        << (*topo_iterator).getName() << " cells. \n"; 
      
      // Forward map: do not specify cell ordinal
      CellTools::mapToPhysicalFrame(physPoints, cubPoints, cellWorkset, (*topo_iterator));
      // Inverse map: do not specify cell ordinal
      CellTools::mapToReferenceFrame(controlPoints, physPoints, cellWorkset, (*topo_iterator));
      
      // Check: points in controlPoints should match the originals in cubPoints up to a tolerance
      for(int cellOrd = 0; cellOrd < numCells; cellOrd++){
        for(int pt = 0; pt < numPts; pt++){
          for(int d = 0; d < cellDim; d++){
            
            if( abs( controlPoints(cellOrd, pt, d) - cubPoints(pt, d) ) > 100.0*INTREPID_TOL ){
              errorFlag++;
              *outStream
                << std::setw(70) << "^^^^----FAILURE!" << "\n"
                << " Mapping a single point set to all physical cells in a workset failed for: \n"
                << "                    Cell Topology = " << (*topo_iterator).getName() << "\n"
                << " Physical cell ordinal in workset = " << cellOrd << "\n"
                << "          Reference point ordinal = " << setprecision(12) << pt << "\n"
                << "    At reference point coordinate = " << setprecision(12) << d << "\n"
                << "                   Original value = " << cubPoints(pt, d) << "\n"
                << "                     F^{-1}F(P_d) = " << controlPoints(cellOrd, pt, d) <<"\n";
            }
          }// d
        }// pt
      }// cellOrd      
      /* 
       * 3.3 Test 3: multiple point sets to multiple physical cells: map ref. point sets in rank-3 (C,P,D) array
       *     to physical point sets in rank-3 (C, P,D) array for all cell ordinals. The (C,P,D) array
       *     with reference point sets is obtained by cloning the cubature point array.
       */
      physPoints.clear(); 
      controlPoints.clear();
      refPoints.resize(numCells, numPts, cubDim);
      physPoints.resize(numCells, numPts, cubDim);
      controlPoints.resize(numCells, numPts, cubDim);
      
      // Clone cubature points in refPoints:
      for(int c = 0; c < numCells; c++){
        for(int pt = 0; pt < numPts; pt++){
          for(int d = 0; d < cellDim; d++){
            refPoints(c, pt, d) = cubPoints(pt, d);
          }// d
        }// pt
      }// c
      
      *outStream 
        << " Mapping " << numCells << " sets of " << numPts << " points to corresponding cells in workset of " << numCells << " " 
        << (*topo_iterator).getName() << " cells. \n"; 
      
      // Forward map: do not specify cell ordinal
      CellTools::mapToPhysicalFrame(physPoints, refPoints, cellWorkset, (*topo_iterator));
      // Inverse map: do not specify cell ordinal
      CellTools::mapToReferenceFrame(controlPoints, physPoints, cellWorkset, (*topo_iterator));
      
      // Check: points in controlPoints should match the originals in cubPoints up to a tolerance
      for(int cellOrd = 0; cellOrd < numCells; cellOrd++){
        for(int pt = 0; pt < numPts; pt++){
          for(int d = 0; d < cellDim; d++){
            
            if( abs( controlPoints(cellOrd, pt, d) - cubPoints(pt, d) ) > 100.0*INTREPID_TOL ){
              errorFlag++;
              *outStream
                << std::setw(70) << "^^^^----FAILURE!" << "\n"
                << " Mapping multiple point sets to corresponding physical cells in a workset failed for: \n"
                << "                    Cell Topology = " << (*topo_iterator).getName() << "\n"
                << " Physical cell ordinal in workset = " << cellOrd << "\n"
                << "          Reference point ordinal = " << setprecision(12) << pt << "\n"
                << "    At reference point coordinate = " << setprecision(12) << d << "\n"
                << "                   Original value = " << refPoints(cellOrd, pt, d) << "\n"
                << "                     F^{-1}F(P_d) = " << controlPoints(cellOrd, pt, d) <<"\n";
            }
          }// d
        }// pt
      }// cellOrd
    }// topo_iterator
  }// try using default initial guesses branch
  
  /*************************************************************************************************
    *         Wrap up test: check if the test broke down unexpectedly due to an exception          *
    ************************************************************************************************/

    catch (std::logic_error err) {
    *outStream << err.what() << "\n";
    errorFlag = -1000;
  };
  
  
  // Test 2: repeat all of the above using the original points as user-defined initial guess
  try{
    
    *outStream \
    << "\n"
    << "===============================================================================\n"\
    << "| Test 2: computing F(x) and F^{-1}(x) using user-defined initial guess.      |\n"\
    << "===============================================================================\n\n";
    /*
     *  Test summary:
     *     
     *    Repeats all parts of Test 1 using user-defined initial guesses in the computation of F^{-1}.
     *    The guesses are simply the exact solutions (the original set of points we started with)
     *    and so, this tests runs much faster because Newton converges in a single iteration.
     *
     */
    
    // Declare arrays for cell workset and point sets. We will have 10 cells in the wset. and 10 pts per pt. set
    FieldContainer<double> cellWorkset;                 // physical cell workset
    FieldContainer<double> physPoints;                  // physical point set(s)
    FieldContainer<double> controlPoints;               // preimages: physical points mapped back to ref. frame
    FieldContainer<double> initialGuess;                // User-defined initial guesses for F^{-1}
    
    // We will use cubature factory to get some points on the reference cells. Declare necessary arrays
    DefaultCubatureFactory<double>  cubFactory;   
    FieldContainer<double> cubPoints;
    FieldContainer<double> cubWeights;
    
    // Initialize number of cells in the cell workset
    int numCells  = 10;
    
    
    // Loop over cell topologies, make cell workset for each one by perturbing the cellWorkset & test methods
    for(topo_iterator = supportedTopologies.begin(); topo_iterator != supportedTopologies.end(); ++topo_iterator){
      
      // 1.   Define a single reference point set using cubature factory with order 6 cubature
      Teuchos::RCP<Cubature<double> > cellCubature = cubFactory.create( (*topo_iterator), 4); 
      int cubDim = cellCubature -> getDimension();
      int numPts = cellCubature -> getNumPoints();
      cubPoints.resize(numPts, cubDim);
      cubWeights.resize(numPts);
      cellCubature -> getCubature(cubPoints, cubWeights);
      
      // 2.   Define a cell workset by perturbing the cellWorkset of the reference cell with the specified topology
      // 2.1  Resize dimensions of the rank-3 (C,N,D) cell workset array for the current topology
      int numNodes = (*topo_iterator).getNodeCount();
      int cellDim  = (*topo_iterator).getDimension();
      cellWorkset.resize(numCells, numNodes, cellDim);
      
      // 2.2  Copy cellWorkset of the reference cell with the same topology to temp rank-2 (N,D) array
      FieldContainer<double> refCellNodes(numNodes, cellDim );
      CellTools::getReferenceSubcellNodes(refCellNodes, cellDim, 0, (*topo_iterator) );
      
      // 2.3  Create randomly perturbed version of the reference cell and save in the cell workset array
      for(int cellOrd = 0; cellOrd < numCells; cellOrd++){
        
        // Move vertices +/-0.125 along their axes. Gives nondegenerate cells for base and extended topologies 
        for(int nodeOrd = 0; nodeOrd < numNodes; nodeOrd++){
          for(int d = 0; d < cellDim; d++){
            double delta = Teuchos::ScalarTraits<double>::random()/16.0;
            cellWorkset(cellOrd, nodeOrd, d) = refCellNodes(nodeOrd, d) + delta;
          } // d
        }// nodeOrd           
      }// cellOrd
      /* 
       * 3.1 Test 1: single point set to single physical cell: map ref. point set in rank-2 (P,D) array
       *      to a physical point set in rank-2 (P,D) array for a specified cell ordinal. Use the cub.
       *      points array for this test. Resize physPoints and controlPoints to rank-2 (P,D) arrays.
       */
      physPoints.resize(numPts, cubDim);
      controlPoints.resize(numPts, cubDim);
      
      *outStream 
        << " Mapping a set of " << numPts << " points to one cell in a workset of " << numCells << " " 
        << (*topo_iterator).getName() << " cells. \n"; 
      
      for(int cellOrd = 0; cellOrd < numCells; cellOrd++){
        
        // Forward map:: requires cell ordinal
        CellTools::mapToPhysicalFrame(physPoints, cubPoints, cellWorkset, (*topo_iterator), cellOrd);
        // Inverse map: requires cell ordinal. Use cubPoints as initial guess
        CellTools::mapToReferenceFrameInitGuess(controlPoints, cubPoints,  physPoints, cellWorkset, (*topo_iterator), cellOrd);
        
        // Points in controlPoints should match the originals in cubPoints up to a tolerance
        for(int pt = 0; pt < numPts; pt++){
          for(int d = 0; d < cellDim; d++){
            
            if( abs( controlPoints(pt, d) - cubPoints(pt, d) ) > 100.0*INTREPID_TOL ){
              errorFlag++;
              *outStream
                << std::setw(70) << "^^^^----FAILURE!" << "\n"
                << " Mapping a single point set to a single physical cell in a workset failed for: \n"
                << "                    Cell Topology = " << (*topo_iterator).getName() << "\n"
                << " Physical cell ordinal in workset = " << cellOrd << "\n"
                << "          Reference point ordinal = " << setprecision(12) << pt << "\n"
                << "    At reference point coordinate = " << setprecision(12) << d << "\n"
                << "                   Original value = " << cubPoints(pt, d) << "\n"
                << "                     F^{-1}F(P_d) = " << controlPoints(pt, d) <<"\n";
            }
          }// d
        }// pt
      }// cellOrd
      /* 
       * 3.2  Test 2: single point set to multiple physical cells: map ref. point set in rank-2 (P,D) array
       *      to a physical point set in rank-3 (C, P,D) array for all cell ordinals. Use the cub.
       *      points array for this test. Resize physPoints and controlPoints to rank-3 (C,P,D) arrays.
       */
      physPoints.clear(); 
      controlPoints.clear();
      physPoints.resize(numCells, numPts, cubDim);
      controlPoints.resize(numCells, numPts, cubDim);
    
      // Clone cubature points in initialGuess:
      initialGuess.resize(numCells, numPts, cubDim);
      for(int c = 0; c < numCells; c++){
        for(int pt = 0; pt < numPts; pt++){
          for(int d = 0; d < cellDim; d++){
            initialGuess(c, pt, d) = cubPoints(pt, d);
          }// d
        }// pt
      }// c
      
      *outStream 
        << " Mapping a set of " << numPts << " points to all cells in workset of " << numCells << " " 
        << (*topo_iterator).getName() << " cells. \n"; 
      
      // Forward map: do not specify cell ordinal
      CellTools::mapToPhysicalFrame(physPoints, cubPoints, cellWorkset, (*topo_iterator));
      // Inverse map: do not specify cell ordinal
      CellTools::mapToReferenceFrameInitGuess(controlPoints, initialGuess, physPoints, cellWorkset, (*topo_iterator));
      
      // Check: points in controlPoints should match the originals in cubPoints up to a tolerance
      for(int cellOrd = 0; cellOrd < numCells; cellOrd++){
        for(int pt = 0; pt < numPts; pt++){
          for(int d = 0; d < cellDim; d++){
            
            if( abs( controlPoints(cellOrd, pt, d) - cubPoints(pt, d) ) > 100.0*INTREPID_TOL ){
              errorFlag++;
              *outStream
                << std::setw(70) << "^^^^----FAILURE!" << "\n"
                << " Mapping a single point set to all physical cells in a workset failed for: \n"
                << "                    Cell Topology = " << (*topo_iterator).getName() << "\n"
                << " Physical cell ordinal in workset = " << cellOrd << "\n"
                << "          Reference point ordinal = " << setprecision(12) << pt << "\n"
                << "    At reference point coordinate = " << setprecision(12) << d << "\n"
                << "                   Original value = " << cubPoints(pt, d) << "\n"
                << "                     F^{-1}F(P_d) = " << controlPoints(cellOrd, pt, d) <<"\n";
            }
          }// d
        }// pt
      }// cellOrd
      /* 
       * 3.3 Test 3: multiple point sets to multiple physical cells: map ref. point sets in rank-3 (C,P,D) array
       *     to physical point sets in rank-3 (C, P,D) array for all cell ordinals. The initialGuess
       *     array from last test is used as the required (C,P,D) array of reference points for the
       *     forward map and as the user-defined initial guess array for the inverse map
       */
      physPoints.clear(); 
      controlPoints.clear();
      physPoints.resize(numCells, numPts, cubDim);
      controlPoints.resize(numCells, numPts, cubDim);
            
      *outStream 
        << " Mapping " << numCells << " sets of " << numPts << " points to corresponding cells in workset of " << numCells << " " 
        << (*topo_iterator).getName() << " cells. \n"; 
      
      // Forward map: do not specify cell ordinal
      CellTools::mapToPhysicalFrame(physPoints, initialGuess, cellWorkset, (*topo_iterator));
      // Inverse map: do not specify cell ordinal
      CellTools::mapToReferenceFrameInitGuess(controlPoints, initialGuess, physPoints, cellWorkset, (*topo_iterator));
      
      // Check: points in controlPoints should match the originals in cubPoints up to a tolerance
      for(int cellOrd = 0; cellOrd < numCells; cellOrd++){
        for(int pt = 0; pt < numPts; pt++){
          for(int d = 0; d < cellDim; d++){
            
            if( abs( controlPoints(cellOrd, pt, d) - cubPoints(pt, d) ) > 100.0*INTREPID_TOL ){
              errorFlag++;
              *outStream
                << std::setw(70) << "^^^^----FAILURE!" << "\n"
                << " Mapping multiple point sets to corresponding physical cells in a workset failed for: \n"
                << "                    Cell Topology = " << (*topo_iterator).getName() << "\n"
                << " Physical cell ordinal in workset = " << cellOrd << "\n"
                << "          Reference point ordinal = " << setprecision(12) << pt << "\n"
                << "    At reference point coordinate = " << setprecision(12) << d << "\n"
                << "                   Original value = " << initialGuess(cellOrd, pt, d) << "\n"
                << "                     F^{-1}F(P_d) = " << controlPoints(cellOrd, pt, d) <<"\n";
            }
          }// d
        }// pt
      }// cellOrd
    } //topo-iterator
  }// try user-defined initial guess    
  
  /*************************************************************************************************
    *         Wrap up test: check if the test broke down unexpectedly due to an exception          *
    ************************************************************************************************/
  
  catch (std::logic_error err) {
    *outStream << err.what() << "\n";
    errorFlag = -1000;
  };
 
  *outStream \
    << "\n"
    << "===============================================================================\n"\
    << "| Test 3: Exception testing - only when HAVE_INTREPID_DEBUG is defined.       |\n"\
    << "===============================================================================\n\n";
  /*
   *  Test summary:
   *    Calls methods of CellTools class with incorrectly configured arguments. This test is run only
   *    in debug mode because many of the exceptions are checked only in that mode.
   *
   */
  
  // Initialize throw counter for exception testing
  int nException     = 0;
  int throwCounter   = 0;  
  
  try {
    
#ifdef HAVE_INTREPID_DEBUG
    // Some arbitrary dimensions
    int C = 10;
    int P = 21;
    int N;
    int D;
    int V;
    
    // Array arguments
    FieldContainer<double> jacobian;
    FieldContainer<double> jacobianInv;
    FieldContainer<double> jacobianDet;
    FieldContainer<double> points;
    FieldContainer<double> cellWorkset;
    FieldContainer<double> physPoints;
    FieldContainer<double> refPoints;
    FieldContainer<double> initGuess;
        
    /***********************************************************************************************
      *                          Exception tests for setJacobian method                            *
      **********************************************************************************************/
    
    // Use the second cell topology for these tests (Triangle<6>)
    topo_iterator = supportedTopologies.begin() + 1;
    D = (*topo_iterator).getDimension();
    N = (*topo_iterator).getNodeCount();
    V = (*topo_iterator).getVertexCount();

    // 1. incorrect jacobian rank
    jacobian.resize(C, P, D);
    points.resize(P, D);
    cellWorkset.resize(C, N, D);
    INTREPID_TEST_COMMAND( CellTools::setJacobian(jacobian, points, cellWorkset, (*topo_iterator) ), 
                           throwCounter, nException );
    
    // 2. Incorrect cellWorkset rank
    cellWorkset.resize(C, D);
    INTREPID_TEST_COMMAND( CellTools::setJacobian(jacobian, points, cellWorkset, (*topo_iterator) ), 
                           throwCounter, nException );
    
    // 3. Incorrect points rank
    cellWorkset.resize(C, N, D);
    points.resize(D);
    INTREPID_TEST_COMMAND( CellTools::setJacobian(jacobian, points, cellWorkset, (*topo_iterator) ), 
                           throwCounter, nException );
    
    // 4. points rank incompatible with whichCell = valid cell ordinal
    points.resize(C, P, D);
    INTREPID_TEST_COMMAND( CellTools::setJacobian(jacobian, points, cellWorkset, (*topo_iterator), 0 ), 
                           throwCounter, nException );
    
    // 5. Non-matching dim
    jacobian.resize(C, P, D, D);
    points.resize(C, P, D - 1);
    INTREPID_TEST_COMMAND( CellTools::setJacobian(jacobian, points, cellWorkset, (*topo_iterator) ), 
                           throwCounter, nException );
 
    // 6. Non-matching dim
    jacobian.resize(C, P, D, D);
    points.resize(C, P - 1, D);
    INTREPID_TEST_COMMAND( CellTools::setJacobian(jacobian, points, cellWorkset, (*topo_iterator) ), 
                           throwCounter, nException );
    
    // 7. Non-matching dim
    jacobian.resize(C, P, D, D);
    points.resize(C - 1, P, D);
    INTREPID_TEST_COMMAND( CellTools::setJacobian(jacobian, points, cellWorkset, (*topo_iterator) ), 
                           throwCounter, nException );
 
    // 8. Non-matching dim
    jacobian.resize(C, P, D, D);
    points.resize(C, P, D);
    cellWorkset.resize(C, N, D - 1);
    INTREPID_TEST_COMMAND( CellTools::setJacobian(jacobian, points, cellWorkset, (*topo_iterator) ), 
                           throwCounter, nException );
    
    // 9. Non-matching dim
    jacobian.resize(C, P, D, D);
    points.resize(C, P, D);
    cellWorkset.resize(C - 1, N, D);
    INTREPID_TEST_COMMAND( CellTools::setJacobian(jacobian, points, cellWorkset, (*topo_iterator) ), 
                           throwCounter, nException );
    
    // 10. Incompatible ranks
    jacobian.resize(C, D, D);
    points.resize(C, P, D);
    cellWorkset.resize(C, N, D);
    INTREPID_TEST_COMMAND( CellTools::setJacobian(jacobian, points, cellWorkset, (*topo_iterator) ), 
                           throwCounter, nException );
    
    /***********************************************************************************************
      *                          Exception tests for setJacobianInv method                         *
      **********************************************************************************************/
    
    // 11. incompatible ranks
    jacobian.resize(C, P, D, D);
    jacobianInv.resize(P, D, D);
    INTREPID_TEST_COMMAND( CellTools::setJacobianInv(jacobianInv, jacobian), 
                           throwCounter, nException );
    
    // 12. incorrect ranks
    jacobian.resize(D, D);
    jacobianInv.resize(D, D);
    INTREPID_TEST_COMMAND( CellTools::setJacobianInv(jacobianInv, jacobian), 
                           throwCounter, nException );

    // 13. nonmatching dimensions
    jacobian.resize(C, P, D, D - 1);
    jacobianInv.resize(C, P, D, D);
    INTREPID_TEST_COMMAND( CellTools::setJacobianInv(jacobianInv, jacobian), 
                           throwCounter, nException );
    
    // 14. nonmatching dimensions
    jacobian.resize(C, P, D - 1, D);
    jacobianInv.resize(C, P, D, D);
    INTREPID_TEST_COMMAND( CellTools::setJacobianInv(jacobianInv, jacobian), 
                           throwCounter, nException );
    
    // 15. nonmatching dimensions
    jacobian.resize(C, P - 1, D, D);
    jacobianInv.resize(C, P, D, D);
    INTREPID_TEST_COMMAND( CellTools::setJacobianInv(jacobianInv, jacobian), 
                           throwCounter, nException );
    
    // 16. nonmatching dimensions
    jacobian.resize(C - 1, P, D, D);
    jacobianInv.resize(C, P, D, D);
    INTREPID_TEST_COMMAND( CellTools::setJacobianInv(jacobianInv, jacobian), 
                           throwCounter, nException );
    
    /***********************************************************************************************
      *                          Exception tests for setJacobianDet method                         *
      **********************************************************************************************/
    
    // 17. Incompatible ranks
    jacobian.resize(C, P, D, D);
    jacobianDet.resize(C, P, D);
    INTREPID_TEST_COMMAND( CellTools::setJacobianDet(jacobianDet, jacobian), 
                           throwCounter, nException );

    // 18. Incompatible ranks
    jacobian.resize(P, D, D);
    jacobianDet.resize(C, P);
    INTREPID_TEST_COMMAND( CellTools::setJacobianDet(jacobianDet, jacobian), 
                           throwCounter, nException );
    
    // 19. Incorrect rank
    jacobian.resize(D, D);
    jacobianDet.resize(C, P);
    INTREPID_TEST_COMMAND( CellTools::setJacobianDet(jacobianDet, jacobian), 
                           throwCounter, nException );
    
    // 20. Incorrect rank
    jacobian.resize(C, P, D, D);
    jacobianDet.resize(C);
    INTREPID_TEST_COMMAND( CellTools::setJacobianDet(jacobianDet, jacobian), 
                           throwCounter, nException );
    
    // 21. Incorrect dimension
    jacobian.resize(C, P, D, D);
    jacobianDet.resize(C, P-1);
    INTREPID_TEST_COMMAND( CellTools::setJacobianDet(jacobianDet, jacobian), 
                           throwCounter, nException );

    // 22. Incorrect dimension
    jacobian.resize(C - 1, P, D, D);
    jacobianDet.resize(C, P);
    INTREPID_TEST_COMMAND( CellTools::setJacobianDet(jacobianDet, jacobian), 
                           throwCounter, nException );
    
    /***********************************************************************************************
      *                        Exception tests for mapToPhysicalFrame method                       *
      **********************************************************************************************/
    
    // 23. Incorrect refPoint rank
    refPoints.resize(P);
    physPoints.resize(P, D);
    cellWorkset.resize(C, N, D);
    INTREPID_TEST_COMMAND( CellTools::mapToPhysicalFrame(physPoints, refPoints, cellWorkset, (*topo_iterator) ),
                           throwCounter, nException );
    // 24. Incorrect workset rank
    cellWorkset.resize(P, D);
    refPoints.resize(P, D);
    physPoints.resize(C, P, D);
    INTREPID_TEST_COMMAND( CellTools::mapToPhysicalFrame(physPoints, refPoints, cellWorkset, (*topo_iterator) ),
                           throwCounter, nException );
    
    // 25. Incompatible ranks
    refPoints.resize(C, P, D);
    physPoints.resize(P, D);
    cellWorkset.resize(C, N, D);
    INTREPID_TEST_COMMAND( CellTools::mapToPhysicalFrame(physPoints, refPoints, cellWorkset, (*topo_iterator) ),
                           throwCounter, nException );
    
    // 26. Incompatible dimensions
    refPoints.resize(C, P, D);
    physPoints.resize(C, P, D - 1);
    INTREPID_TEST_COMMAND( CellTools::mapToPhysicalFrame(physPoints, refPoints, cellWorkset, (*topo_iterator) ),
                           throwCounter, nException );

    // 27. Incompatible dimensions
    refPoints.resize(C, P, D);
    physPoints.resize(C, P - 1, D);
    INTREPID_TEST_COMMAND( CellTools::mapToPhysicalFrame(physPoints, refPoints, cellWorkset, (*topo_iterator) ),
                           throwCounter, nException );
    
    // 28. Incompatible dimensions
    refPoints.resize(C, P, D);
    physPoints.resize(C - 1, P, D);
    INTREPID_TEST_COMMAND( CellTools::mapToPhysicalFrame(physPoints, refPoints, cellWorkset, (*topo_iterator) ),
                           throwCounter, nException );
    
    // 29. Incorrect physPoints rank when whichCell is valid cell ordinal
    refPoints.resize(P, D);
    physPoints.resize(C, P, D);
    INTREPID_TEST_COMMAND( CellTools::mapToPhysicalFrame(physPoints, refPoints, cellWorkset, (*topo_iterator), 0 ),
                           throwCounter, nException );
    
    /***********************************************************************************************
      *          Exception tests for mapToReferenceFrame method (with default initial guesses)     *
      **********************************************************************************************/
    
    // 30. incompatible ranks
    refPoints.resize(C, P, D);
    physPoints.resize(P, D);
    cellWorkset.resize(C, N, D);
    INTREPID_TEST_COMMAND( CellTools::mapToReferenceFrame(refPoints, physPoints, cellWorkset, (*topo_iterator) ),
                           throwCounter, nException );
    
    // 31. Incompatible ranks with whichCell = valid cell ordinal
    refPoints.resize(C, P, D);
    physPoints.resize(C, P, D);
    INTREPID_TEST_COMMAND( CellTools::mapToReferenceFrame(refPoints, physPoints, cellWorkset, (*topo_iterator), 0 ),
                           throwCounter, nException );

    // 32. Incompatible ranks with whichCell = -1 (default)
    refPoints.resize(P, D);
    physPoints.resize(P, D);
    INTREPID_TEST_COMMAND( CellTools::mapToReferenceFrame(refPoints, physPoints, cellWorkset, (*topo_iterator)),
                           throwCounter, nException );
    
    // 33. Nonmatching dimensions
    refPoints.resize(C, P, D - 1);
    physPoints.resize(C, P, D);
    INTREPID_TEST_COMMAND( CellTools::mapToReferenceFrame(refPoints, physPoints, cellWorkset, (*topo_iterator)),
                           throwCounter, nException );
    
    // 34. Nonmatching dimensions
    refPoints.resize(C, P - 1, D);
    physPoints.resize(C, P, D);
    INTREPID_TEST_COMMAND( CellTools::mapToReferenceFrame(refPoints, physPoints, cellWorkset, (*topo_iterator)),
                           throwCounter, nException );

    // 35. Nonmatching dimensions
    refPoints.resize(C - 1, P, D);
    physPoints.resize(C, P, D);
    INTREPID_TEST_COMMAND( CellTools::mapToReferenceFrame(refPoints, physPoints, cellWorkset, (*topo_iterator)),
                           throwCounter, nException );
    
    // 36. Incorrect rank for cellWorkset
    refPoints.resize(C, P, D);
    physPoints.resize(C, P, D);
    cellWorkset.resize(C, N);    
    INTREPID_TEST_COMMAND( CellTools::mapToReferenceFrame(refPoints, physPoints, cellWorkset, (*topo_iterator)),
                           throwCounter, nException );
    
    /***********************************************************************************************
      *   Exception tests for mapToReferenceFrameInitGuess method (initial guess is a parameter)   *
      **********************************************************************************************/
    
    // 37. Incompatible ranks
    refPoints.resize(C, P, D);
    physPoints.resize(C, P, D);
    initGuess.resize(P, D);
    cellWorkset.resize(C, N, D);
    INTREPID_TEST_COMMAND( CellTools::mapToReferenceFrameInitGuess(refPoints, initGuess, physPoints, cellWorkset, (*topo_iterator) ),
                           throwCounter, nException );

    // 38. Incompatible ranks when whichCell is valid ordinal
    refPoints.resize(P, D);
    physPoints.resize(P, D);
    initGuess.resize(C, P, D);
    INTREPID_TEST_COMMAND( CellTools::mapToReferenceFrameInitGuess(refPoints, initGuess, physPoints, cellWorkset, (*topo_iterator), 0),
                           throwCounter, nException );

    // 39. Nonmatching dimensions
    refPoints.resize(C, P, D);
    physPoints.resize(C, P, D);
    initGuess.resize(C, P, D - 1);
    INTREPID_TEST_COMMAND( CellTools::mapToReferenceFrameInitGuess(refPoints, initGuess, physPoints, cellWorkset, (*topo_iterator)),
                           throwCounter, nException );
    
    // 40. Nonmatching dimensions
    initGuess.resize(C, P - 1, D);
    INTREPID_TEST_COMMAND( CellTools::mapToReferenceFrameInitGuess(refPoints, initGuess, physPoints, cellWorkset, (*topo_iterator)),
                           throwCounter, nException );
    
    // 41. Nonmatching dimensions
    initGuess.resize(C - 1, P, D);
    INTREPID_TEST_COMMAND( CellTools::mapToReferenceFrameInitGuess(refPoints, initGuess, physPoints, cellWorkset, (*topo_iterator)),
                           throwCounter, nException );
        
    /***********************************************************************************************
      *                        Exception tests for mapToReferenceSubcell method                    *
      **********************************************************************************************/
    
    FieldContainer<double> refSubcellPoints;
    FieldContainer<double> paramPoints;
    int subcellDim = 2;
    int subcellOrd = 0;
    
    // This should set cell topology to Tetrahedron<10> so that we have real edges and faces.
    topo_iterator += 5;
    D = (*topo_iterator).getDimension();
    
    // 42. Incorrect array rank
    refSubcellPoints.resize(P,3);
    paramPoints.resize(P);
    INTREPID_TEST_COMMAND( CellTools::mapToReferenceSubcell(refSubcellPoints, paramPoints, subcellDim, subcellOrd, (*topo_iterator)),
                           throwCounter, nException );
   
    // 43. Incorrect array rank
    refSubcellPoints.resize(P);
    paramPoints.resize(P, 2);
    INTREPID_TEST_COMMAND( CellTools::mapToReferenceSubcell(refSubcellPoints, paramPoints, subcellDim, subcellOrd, (*topo_iterator)),
                           throwCounter, nException );
    
    // 44. Incorrect array dimension for face of 3D cell (should be 3)
    refSubcellPoints.resize(P, 2);
    paramPoints.resize(P, 2);
    INTREPID_TEST_COMMAND( CellTools::mapToReferenceSubcell(refSubcellPoints, paramPoints, subcellDim, subcellOrd, (*topo_iterator)),
                           throwCounter, nException );
    
    // 45. Incorrect array dimension for parametrization domain of a face of 3D cell (should be 2)
    refSubcellPoints.resize(P, 3);
    paramPoints.resize(P, 3);
    INTREPID_TEST_COMMAND( CellTools::mapToReferenceSubcell(refSubcellPoints, paramPoints, subcellDim, subcellOrd, (*topo_iterator)),
                           throwCounter, nException );
    
    /***********************************************************************************************
      *                        Exception tests for getReferenceEdgeTangent method                  *
      **********************************************************************************************/
    
    FieldContainer<double> refEdgeTangent;

    // 46. Incorrect rank
    refEdgeTangent.resize(C,P,D);
    INTREPID_TEST_COMMAND( CellTools::getReferenceEdgeTangent(refEdgeTangent, 0, (*topo_iterator)),
                           throwCounter, nException );
    
    // 47. Incorrect dimension D for Tet<10> cell
    refEdgeTangent.resize(2);
    INTREPID_TEST_COMMAND( CellTools::getReferenceEdgeTangent(refEdgeTangent, 0, (*topo_iterator)),
                           throwCounter, nException );
    
    // 48. Invalid edge ordinal for Tet<10>
    refEdgeTangent.resize(C,P,D);
    INTREPID_TEST_COMMAND( CellTools::getReferenceEdgeTangent(refEdgeTangent, 10, (*topo_iterator)),
                           throwCounter, nException );
    
    /***********************************************************************************************
      *                        Exception tests for getReferenceFaceTangents method                 *
      **********************************************************************************************/
    
    FieldContainer<double> refFaceTanU;
    FieldContainer<double> refFaceTanV;
    
    // 49. Incorrect rank
    refFaceTanU.resize(P, D);
    refFaceTanV.resize(D);
    INTREPID_TEST_COMMAND( CellTools::getReferenceFaceTangents(refFaceTanU, refFaceTanV, 0, (*topo_iterator)),
                           throwCounter, nException );
    
    // 50. Incorrect rank
    refFaceTanU.resize(D);
    refFaceTanV.resize(P, D);
    INTREPID_TEST_COMMAND( CellTools::getReferenceFaceTangents(refFaceTanU, refFaceTanV, 0, (*topo_iterator)),
                           throwCounter, nException );

    // 51. Incorrect dimension for 3D cell
    refFaceTanU.resize(D - 1);
    refFaceTanV.resize(D);
    INTREPID_TEST_COMMAND( CellTools::getReferenceFaceTangents(refFaceTanU, refFaceTanV, 0, (*topo_iterator)),
                           throwCounter, nException );

    // 52. Incorrect dimension for 3D cell
    refFaceTanU.resize(D);
    refFaceTanV.resize(D - 1);
    INTREPID_TEST_COMMAND( CellTools::getReferenceFaceTangents(refFaceTanU, refFaceTanV, 0, (*topo_iterator)),
                           throwCounter, nException );
    
    // 53. Invalid face ordinal
    refFaceTanU.resize(D);
    refFaceTanV.resize(D);
    INTREPID_TEST_COMMAND( CellTools::getReferenceFaceTangents(refFaceTanU, refFaceTanV, 10, (*topo_iterator)),
                           throwCounter, nException );
    
    /***********************************************************************************************
      *                        Exception tests for getReferenceSide/FaceNormal methods             *
      **********************************************************************************************/
    
    FieldContainer<double> refSideNormal;
    
    // 54-55. Incorrect rank
    refSideNormal.resize(C,P);
    INTREPID_TEST_COMMAND( CellTools::getReferenceSideNormal(refSideNormal, 0, (*topo_iterator)),
                           throwCounter, nException );
    INTREPID_TEST_COMMAND( CellTools::getReferenceFaceNormal(refSideNormal, 0, (*topo_iterator)),
                           throwCounter, nException );
    
    // 56-57. Incorrect dimension for 3D cell 
    refSideNormal.resize(D - 1);
    INTREPID_TEST_COMMAND( CellTools::getReferenceSideNormal(refSideNormal, 0, (*topo_iterator)),
                           throwCounter, nException );
    INTREPID_TEST_COMMAND( CellTools::getReferenceFaceNormal(refSideNormal, 0, (*topo_iterator)),
                           throwCounter, nException );
    
    // 58-59. Invalid side ordinal for Tet<10>
    refSideNormal.resize(D);
    INTREPID_TEST_COMMAND( CellTools::getReferenceSideNormal(refSideNormal, 10, (*topo_iterator)),
                           throwCounter, nException );
    INTREPID_TEST_COMMAND( CellTools::getReferenceFaceNormal(refSideNormal, 10, (*topo_iterator)),
                           throwCounter, nException );
    
    // 60. Incorrect dimension for 2D cell: reset topo_iterator to the first cell in supportedTopologies which is Tri<3> 
    topo_iterator = supportedTopologies.begin();
    D = (*topo_iterator).getDimension();
    refSideNormal.resize(D - 1);
    INTREPID_TEST_COMMAND( CellTools::getReferenceSideNormal(refSideNormal, 0, (*topo_iterator)),
                           throwCounter, nException );
    
    // 61. Invalid side ordinal for Tri<3>
    refSideNormal.resize(D);
    INTREPID_TEST_COMMAND( CellTools::getReferenceSideNormal(refSideNormal, 10, (*topo_iterator)),
                           throwCounter, nException );
    
    // 62. Cannot call the "face" method for 2D cells
    refSideNormal.resize(D);
    INTREPID_TEST_COMMAND( CellTools::getReferenceFaceNormal(refSideNormal, 0, (*topo_iterator)),
                           throwCounter, nException );
    
    /***********************************************************************************************
      *          Exception tests for checkPoint/Pointset/PointwiseInclusion methods        *
      **********************************************************************************************/
    points.resize(2,3,3,4);
    FieldContainer<int> inCell;
    
    // 63. Point dimension does not match cell topology
    double * point = 0;
    INTREPID_TEST_COMMAND(CellTools::checkPointInclusion(point, (*topo_iterator).getDimension() + 1, (*topo_iterator) ),
                          throwCounter, nException );
    
    // 64. Invalid cell topology
    CellTopology pentagon_5(shards::getCellTopologyData<shards::Pentagon<> >() );
    INTREPID_TEST_COMMAND(CellTools::checkPointInclusion(point, pentagon_5.getDimension(), pentagon_5 ),
                          throwCounter, nException );
        
    // 65. Incorrect spatial dimension of points
    points.resize(10, 10, (*topo_iterator).getDimension() + 1);
    INTREPID_TEST_COMMAND(CellTools::checkPointsetInclusion(points, (*topo_iterator) ),
                          throwCounter, nException );
    
    // 66. Incorrect rank of input array
    points.resize(10,10,10,3);
    INTREPID_TEST_COMMAND(CellTools::checkPointsetInclusion(points, (*topo_iterator) ),
                          throwCounter, nException );
    
    // 67. Incorrect rank of output array
    points.resize(10,10,(*topo_iterator).getDimension() );
    inCell.resize(10);
    INTREPID_TEST_COMMAND(CellTools::checkPointwiseInclusion(inCell, points, (*topo_iterator) ),
                          throwCounter, nException );
  
    // 68. Incorrect rank of output array
    points.resize(10, (*topo_iterator).getDimension() );
    inCell.resize(10, 10);
    INTREPID_TEST_COMMAND(CellTools::checkPointwiseInclusion(inCell, points, (*topo_iterator) ),
                          throwCounter, nException );
    
    // 69. Incorrect rank of output array
    points.resize((*topo_iterator).getDimension() );
    inCell.resize(10, 10);
    INTREPID_TEST_COMMAND(CellTools::checkPointwiseInclusion(inCell, points, (*topo_iterator) ),
                          throwCounter, nException );
    
    // 70. Incorrect dimension of output array
    points.resize(10, 10, (*topo_iterator).getDimension() );
    inCell.resize(10, 9);
    INTREPID_TEST_COMMAND(CellTools::checkPointwiseInclusion(inCell, points, (*topo_iterator) ),
                          throwCounter, nException );
    
    // 71. Incorrect dimension of output array
    points.resize(10, 10, (*topo_iterator).getDimension() );
    inCell.resize(9, 10);
    INTREPID_TEST_COMMAND(CellTools::checkPointwiseInclusion(inCell, points, (*topo_iterator) ),
                          throwCounter, nException );
    
    // 72. Incorrect dimension of output array
    points.resize(10, (*topo_iterator).getDimension() );
    inCell.resize(9);
    INTREPID_TEST_COMMAND(CellTools::checkPointwiseInclusion(inCell, points, (*topo_iterator) ),
                          throwCounter, nException );
    
    // 73. Incorrect spatial dimension of input array
    points.resize(10, 10, (*topo_iterator).getDimension() + 1);
    inCell.resize(10, 10);
    INTREPID_TEST_COMMAND(CellTools::checkPointwiseInclusion(inCell, points, (*topo_iterator) ),
                          throwCounter, nException );
    
    // 74. Incorrect rank of input array.
    points.resize(10,10,10,3);
    INTREPID_TEST_COMMAND(CellTools::checkPointwiseInclusion(inCell, points, (*topo_iterator) ),
                          throwCounter, nException );
    
        
    physPoints.resize(C, P, D);
    inCell.resize(C, P);
    // 75. Invalid rank of cellWorkset
    cellWorkset.resize(C, N, D, D);
    INTREPID_TEST_COMMAND(CellTools::checkPointwiseInclusion(inCell, physPoints, cellWorkset, (*topo_iterator) ),
                          throwCounter, nException );
    
    // 76. Invalid dimension 1 (node count) of cellWorkset
    cellWorkset.resize(C, N + 1, D);
    INTREPID_TEST_COMMAND(CellTools::checkPointwiseInclusion(inCell, physPoints, cellWorkset, (*topo_iterator) ),
                          throwCounter, nException );
    
    // 77. Invalid dimension 2 (spatial dimension) of cellWorkset
    cellWorkset.resize(C, N, D + 1);
    INTREPID_TEST_COMMAND(CellTools::checkPointwiseInclusion(inCell, physPoints, cellWorkset, (*topo_iterator) ),
                          throwCounter, nException );
    
    // 78. Invalid whichCell value (exceeds cell count in the workset)
    cellWorkset.resize(C, N, D);
    INTREPID_TEST_COMMAND(CellTools::checkPointwiseInclusion(inCell, physPoints, cellWorkset, (*topo_iterator), C + 1 ),
                          throwCounter, nException );
    
    // 79. Invalid whichCell for rank-3 physPoints (must be -1, here it is valid cell ordinal)
    INTREPID_TEST_COMMAND(CellTools::checkPointwiseInclusion(inCell, physPoints, cellWorkset, (*topo_iterator), 0 ),
                          throwCounter, nException );
    
    // 80. Invalid whichCell for rank-2 physPoints (must be a valid cell ordinal, here it is the default -1)
    physPoints.resize(P, D);
    inCell.resize(P);
    INTREPID_TEST_COMMAND(CellTools::checkPointwiseInclusion(inCell, physPoints, cellWorkset, (*topo_iterator) ),
                          throwCounter, nException );
    
    // 81. Incompatible ranks of I/O arrays
    physPoints.resize(C, P, D);
    inCell.resize(P);
    INTREPID_TEST_COMMAND(CellTools::checkPointwiseInclusion(inCell, physPoints, cellWorkset, (*topo_iterator)),
                          throwCounter, nException );
    
    // 82. Incompatible ranks of I/O arrays
    physPoints.resize(P, D);
    inCell.resize(C, P);
    INTREPID_TEST_COMMAND(CellTools::checkPointwiseInclusion(inCell, physPoints, cellWorkset, (*topo_iterator), 0),
                          throwCounter, nException );
    
    // 83. Incompatible dimensions of I/O arrays
    physPoints.resize(C, P, D);
    inCell.resize(C, P + 1);
    INTREPID_TEST_COMMAND(CellTools::checkPointwiseInclusion(inCell, physPoints, cellWorkset, (*topo_iterator)),
                          throwCounter, nException );

    // 84. Incompatible dimensions of I/O arrays: rank-3 Input
    physPoints.resize(C + 1, P, D);
    inCell.resize(C, P);
    INTREPID_TEST_COMMAND(CellTools::checkPointwiseInclusion(inCell, physPoints, cellWorkset, (*topo_iterator)),
                          throwCounter, nException );
    
    // 85. Incompatible dimensions of I/O arrays: rank-2 Input
    physPoints.resize(P, D);
    inCell.resize(P + 1);
    INTREPID_TEST_COMMAND(CellTools::checkPointwiseInclusion(inCell, physPoints, cellWorkset, (*topo_iterator), 0 ),
                          throwCounter, nException );
    
    
    /***********************************************************************************************
      *               Exception tests for getReferenceVertex/vertices/Node/Nodes methods           *
      **********************************************************************************************/
    
    FieldContainer<double> subcellNodes;
    
    // 86-89. Cell does not have reference cell
    INTREPID_TEST_COMMAND(CellTools::getReferenceVertex(pentagon_5, 0), throwCounter, nException);
    INTREPID_TEST_COMMAND(CellTools::getReferenceNode(pentagon_5, 0), throwCounter, nException);    
    INTREPID_TEST_COMMAND(CellTools::getReferenceSubcellVertices(subcellNodes, 0, 0, pentagon_5), throwCounter, nException);
    INTREPID_TEST_COMMAND(CellTools::getReferenceSubcellNodes(subcellNodes, 0, 0, pentagon_5), throwCounter, nException);

    // Use last cell topology (Wedge<18>) for these tests
    topo_iterator = supportedTopologies.end() - 1;
    D = (*topo_iterator).getDimension();
    int subcDim = D - 1;
    int S = (*topo_iterator).getSubcellCount(subcDim);
    V = (*topo_iterator).getVertexCount(subcDim, S - 1);
    subcellNodes.resize(V, D);
    // 90. subcell ordinal out of range
    INTREPID_TEST_COMMAND(CellTools::getReferenceSubcellVertices(subcellNodes, subcDim, S + 1, (*topo_iterator)), 
                          throwCounter, nException);
    
    // 91. subcell dim out of range
    INTREPID_TEST_COMMAND(CellTools::getReferenceSubcellVertices(subcellNodes, D + 1, S, (*topo_iterator)), 
                          throwCounter, nException);
    
    // 92. Incorrect rank for subcellNodes 
    subcellNodes.resize(V, D, D); 
    INTREPID_TEST_COMMAND(CellTools::getReferenceSubcellVertices(subcellNodes, subcDim, S - 1, (*topo_iterator)), 
                          throwCounter, nException);

    // 93. Incorrect dimension for subcellNodes 
    subcellNodes.resize(V - 1, D); 
    INTREPID_TEST_COMMAND(CellTools::getReferenceSubcellVertices(subcellNodes, subcDim, S - 1, (*topo_iterator)), 
                          throwCounter, nException);
    
    // 94. Incorrect dimension for subcellNodes 
    subcellNodes.resize(V, D - 1); 
    INTREPID_TEST_COMMAND(CellTools::getReferenceSubcellVertices(subcellNodes, subcDim, S - 1, (*topo_iterator)), 
                          throwCounter, nException);
    
          
    N = (*topo_iterator).getNodeCount(subcDim, S - 1);
    subcellNodes.resize(N, D);
    // 95. subcell ordinal out of range
    INTREPID_TEST_COMMAND(CellTools::getReferenceSubcellNodes(subcellNodes, subcDim, S + 1, (*topo_iterator)), 
                          throwCounter, nException);
    
    // 96. subcell dim out of range
    INTREPID_TEST_COMMAND(CellTools::getReferenceSubcellNodes(subcellNodes, D + 1, S, (*topo_iterator)), 
                          throwCounter, nException);
    
    // 97. Incorrect rank for subcellNodes 
    subcellNodes.resize(N, D, D); 
    INTREPID_TEST_COMMAND(CellTools::getReferenceSubcellNodes(subcellNodes, subcDim, S - 1, (*topo_iterator)), 
                          throwCounter, nException);
    
    // 98. Incorrect dimension for subcellNodes 
    subcellNodes.resize(N - 1, D); 
    INTREPID_TEST_COMMAND(CellTools::getReferenceSubcellNodes(subcellNodes, subcDim, S - 1, (*topo_iterator)), 
                          throwCounter, nException);
    
    // 99. Incorrect dimension for subcellNodes 
    subcellNodes.resize(N, D - 1); 
    INTREPID_TEST_COMMAND(CellTools::getReferenceSubcellNodes(subcellNodes, subcDim, S - 1, (*topo_iterator)), 
                          throwCounter, nException);
    
#endif    
  } // try exception testing
  
  /*************************************************************************************************
    *         Wrap up test: check if the test broke down unexpectedly due to an exception          *
    ************************************************************************************************/

  catch(std::logic_error err) {
    *outStream << err.what() << "\n";
    errorFlag = -1000;
  }
  
  // Check if number of thrown exceptions matches the one we expect 
  if (throwCounter != nException) {
    errorFlag++;
    *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
  }
  
  
  if (errorFlag != 0)
    std::cout << "End Result: TEST FAILED\n";
  else
    std::cout << "End Result: TEST PASSED\n";
  
  // reset format state of std::cout
  std::cout.copyfmt(oldFormatState);
  
  return errorFlag;
}
  






