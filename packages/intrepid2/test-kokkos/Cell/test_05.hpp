// @HEADER
// ************************************************************************
//
//                           Intrepid2 Package
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
// Questions? Contact Kyungjoo Kim  (kyukim@sandia.gov), or
//                    Mauro Perego  (mperego@sandia.gov)
//
// ************************************************************************
// @HEADER


/** \file
    \brief  Test of the CellTools class.
    \author Created by P. Bochev, D. Ridzal and K. Peterson
*/
#include "Intrepid2_CellTools.hpp"
#include "Intrepid2_FieldContainer.hpp"
#include "Intrepid2_DefaultCubatureFactory.hpp"
#include "Shards_CellTopology.hpp"

#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_ScalarTraits.hpp"

using namespace std;
using namespace Intrepid2;
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
Kokkos::initialize();
  typedef CellTools<double>       CellTools;
 // typedef shards::CellTopology    CellTopology;
  
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
  supportedTopologies.push_back(shards::getCellTopologyData<Hexahedron<20> >() );
  supportedTopologies.push_back(shards::getCellTopologyData<Hexahedron<27> >() );
  supportedTopologies.push_back(shards::getCellTopologyData<Wedge<6> >() );
  supportedTopologies.push_back(shards::getCellTopologyData<Wedge<15> >() );
  supportedTopologies.push_back(shards::getCellTopologyData<Wedge<18> >() );
  supportedTopologies.push_back(shards::getCellTopologyData<Pyramid<5> >() );
  supportedTopologies.push_back(shards::getCellTopologyData<Pyramid<13> >() );
  
  // Declare iterator to loop over the cell topologies
  std::vector<shards::CellTopology>::iterator topo_iterator;

  *outStream \
    << "\n"
    << "===============================================================================\n"\
    << "| Test 3: Exception testing - only when HAVE_INTREPID2_DEBUG is defined.       |\n"\
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
    
#ifdef HAVE_INTREPID2_DEBUG
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
 Kokkos::finalize(); 
  return errorFlag;
}
  






