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
    \author Created by P. Bochev, D. Ridzal, K. Peterson and Kyungjoo Kim
*/


#include "Intrepid2_config.h"

#ifdef HAVE_INTREPID2_DEBUG
#define INTREPID2_TEST_FOR_DEBUG_ABORT_OVERRIDE_TO_CONTINUE
#endif

#include "Intrepid2_Types.hpp"
#include "Intrepid2_Utils.hpp"

#include "Intrepid2_CellTools.hpp"

#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"


namespace Intrepid2 {

  namespace Test {
#define INTREPID2_TEST_ERROR_EXPECTED( S )                              \
    try {                                                               \
      ++nthrow;                                                         \
      S ;                                                               \
    }                                                                   \
    catch (std::exception &err) {                                      \
      ++ncatch;                                                         \
      *outStream << "Expected Error ----------------------------------------------------------------\n"; \
      *outStream << err.what() << '\n';                                 \
      *outStream << "-------------------------------------------------------------------------------" << "\n\n"; \
    };
#define ConstructWithLabel(obj, ...) obj(#obj, __VA_ARGS__)

    template<typename ValueType, typename DeviceType>
    int CellTools_Test05(const bool verbose) {
      using ExecSpaceType = typename DeviceType::execution_space;

      Teuchos::RCP<std::ostream> outStream;
      Teuchos::oblackholestream bhs; // outputs nothing

      if (verbose)
        outStream = Teuchos::rcp(&std::cout, false);
      else
        outStream = Teuchos::rcp(&bhs, false);

      Teuchos::oblackholestream oldFormatState;
      oldFormatState.copyfmt(std::cout);

      typedef typename
        Kokkos::Impl::is_space<DeviceType>::host_mirror_space::execution_space HostSpaceType ;

      *outStream << "DeviceSpace::  ";   ExecSpaceType::print_configuration(*outStream, false);
      *outStream << "HostSpace::    ";   HostSpaceType::print_configuration(*outStream, false);

      *outStream
        << "===============================================================================\n"
        << "|                                                                             |\n"
        << "|                              Unit Test CellTools                            |\n"
        << "|                                                                             |\n"
        << "|     1) Mapping to and from reference cells with base and extended topologies|\n"
        << "|        using default initial guesses when computing the inverse F^{-1}      |\n"
        << "|     2) Repeat all tests from 1) using user-defined initial guess for F^{-1} |\n"
        << "|     3) Exception testing                                                    |\n"
        << "|                                                                             |\n"
        << "|  Questions? Contact  Pavel Bochev (pbboche@sandia.gov)                      |\n"
        << "|                      Denis Ridzal (dridzal@sandia.gov), or                  |\n"
        << "|                      Kara Peterson (kjpeter@sandia.gov)                     |\n"
        << "|                                                                             |\n"
        << "|  Intrepid's website: http://trilinos.sandia.gov/packages/intrepid           |\n"
        << "|  Trilinos website:   http://trilinos.sandia.gov                             |\n"
        << "|                                                                             |\n"
        << "===============================================================================\n";

      int errorFlag  = 0;

      // Allocate storage and extract all standard cells with base topologies
      std::vector<shards::CellTopology> standardBaseTopologies;
      shards::getTopologies(standardBaseTopologies, 4, shards::STANDARD_CELL, shards::BASE_TOPOLOGY);

      *outStream
        << "\n"
        << "===============================================================================\n"
        << "| Test 3: Exception testing - only when HAVE_INTREPID2_DEBUG is defined.       |\n"
        << "===============================================================================\n\n";
      /*
       *  Test summary:
       *    Calls methods of CellTools class with incorrectly configured arguments. This test is run only
       *    in debug mode because many of the exceptions are checked only in that mode.
       *
       */


      try {
        
#ifdef HAVE_INTREPID2_DEBUG
        typedef Kokkos::DynRankView<ValueType,DeviceType> DynRankView;
        typedef CellTools<DeviceType> ct;

        std::vector<shards::CellTopology> standardBaseTopologies;
        shards::getTopologies(standardBaseTopologies, 4, shards::STANDARD_CELL, shards::BASE_TOPOLOGY);
        const ordinal_type topoSize = standardBaseTopologies.size();
        
        // Loop over admissible topologies
        for (ordinal_type topoOrd=0;topoOrd<topoSize;++topoOrd) {
          const auto cell = standardBaseTopologies[topoOrd];

          // skip cells not supported
          if (!ct::hasReferenceCell(cell))
            continue;

          *outStream << " Testing exceptions for "
                     <<  cell.getName() << "\n";

          // Initialize throw counter for exception testing
          ordinal_type nthrow = 0, ncatch = 0;

          // Some arbitrary dimensions
          ordinal_type C = 10, P = 21, N = cell.getNodeCount(), D = cell.getDimension();// V = cell.getVertexCount();

          // Array arguments
          DynRankView 
            jacobian, jacobianInv, jacobianDet, points, 
            worksetCell, physPoints, refPoints, initGuess,
            refSubcellPoints, paramPoints,
            refEdgeTangent, refFaceTanU, refFaceTanV,
            refSideNormal;

          // for 3D
          ordinal_type subcellDim = 2, subcellOrd = 0;

          /***********************************************************************************************
           *                          Exception tests for setJacobian method                            *
           **********************************************************************************************/
          
          // 1. incorrect jacobian rank
          jacobian = DynRankView("jacobian", C, P, D);
          points = DynRankView("points", P, D);
          worksetCell = DynRankView("worksetCell", C, N, D);
          INTREPID2_TEST_ERROR_EXPECTED( ct::setJacobian(jacobian, points, worksetCell, cell ) );
          
          // 2. Incorrect worksetCell rank
          worksetCell = DynRankView("worksetCell", C, D);
          INTREPID2_TEST_ERROR_EXPECTED( ct::setJacobian(jacobian, points, worksetCell, cell ) );
          
          // 3. Incorrect points rank
          worksetCell = DynRankView("worksetCell", C, N, D);
          points = DynRankView("points", D);
          INTREPID2_TEST_ERROR_EXPECTED( ct::setJacobian(jacobian, points, worksetCell, cell ) );
          
          // 5. Non-matching dim
          jacobian = DynRankView("jacobian", C, P, D, D);
          points = DynRankView("points", C, P, D - 1);
          INTREPID2_TEST_ERROR_EXPECTED( ct::setJacobian(jacobian, points, worksetCell, cell ) );
            
          // 6. Non-matching dim
          jacobian = DynRankView("jacobian", C, P, D, D);
          points = DynRankView("points", C, P - 1, D);
          INTREPID2_TEST_ERROR_EXPECTED( ct::setJacobian(jacobian, points, worksetCell, cell ) );
          
          // 7. Non-matching dim
          jacobian = DynRankView("jacobian", C, P, D, D);
          points = DynRankView("points", C - 1, P, D);
          INTREPID2_TEST_ERROR_EXPECTED( ct::setJacobian(jacobian, points, worksetCell, cell ) );
          
          // 8. Non-matching dim
          jacobian = DynRankView("jacobian", C, P, D, D);
          points = DynRankView("points", C, P, D);
          worksetCell = DynRankView("worksetCell", C, N, D - 1);
          INTREPID2_TEST_ERROR_EXPECTED( ct::setJacobian(jacobian, points, worksetCell, cell ) );
          
          // 9. Non-matching dim
          jacobian = DynRankView("jacobian", C, P, D, D);
          points = DynRankView("points", C, P, D);
          worksetCell = DynRankView("worksetCell", C - 1, N, D);
          INTREPID2_TEST_ERROR_EXPECTED( ct::setJacobian(jacobian, points, worksetCell, cell ) );
          
          // 10. Incompatible ranks
          jacobian = DynRankView("jacobian", C, D, D);
          points = DynRankView("points", C, P, D);
          worksetCell = DynRankView("worksetCell", C, N, D);
          INTREPID2_TEST_ERROR_EXPECTED( ct::setJacobian(jacobian, points, worksetCell, cell ) );
          
          /***********************************************************************************************
           *                          Exception tests for setJacobianInv method                         *
           **********************************************************************************************/
          
          // 11. incompatible ranks
          jacobian = DynRankView("jacobian", C, P, D, D);
          jacobianInv = DynRankView("jacobianInv", P, D, D);
          INTREPID2_TEST_ERROR_EXPECTED( ct::setJacobianInv(jacobianInv, jacobian) );
          
          // 12. incorrect ranks
          jacobian = DynRankView("jacobian", D, D);
          jacobianInv = DynRankView("jacobianInv", D, D);
          INTREPID2_TEST_ERROR_EXPECTED( ct::setJacobianInv(jacobianInv, jacobian) );
          
          // 13. nonmatching dimensions
          jacobian = DynRankView("jacobian", C, P, D, D - 1);
          jacobianInv = DynRankView("jacobianInv", C, P, D, D);
          INTREPID2_TEST_ERROR_EXPECTED( ct::setJacobianInv(jacobianInv, jacobian) );
          
          // 14. nonmatching dimensions
          jacobian = DynRankView("jacobian", C, P, D - 1, D);
          jacobianInv = DynRankView("jacobianInv", C, P, D, D);
          INTREPID2_TEST_ERROR_EXPECTED( ct::setJacobianInv(jacobianInv, jacobian) );
          
          // 15. nonmatching dimensions
          jacobian = DynRankView("jacobian", C, P - 1, D, D);
          jacobianInv = DynRankView("jacobianInv", C, P, D, D);
          INTREPID2_TEST_ERROR_EXPECTED( ct::setJacobianInv(jacobianInv, jacobian) );

          // 16. nonmatching dimensions
          jacobian = DynRankView("jacobian", C - 1, P, D, D);
          jacobianInv = DynRankView("jacobianInv", C, P, D, D);
          INTREPID2_TEST_ERROR_EXPECTED( ct::setJacobianInv(jacobianInv, jacobian) );
          
          /***********************************************************************************************
           *                          Exception tests for setJacobianDet method                         *
           **********************************************************************************************/

          // 17. Incompatible ranks
          jacobian = DynRankView("jacobian", C, P, D, D);
          jacobianDet = DynRankView("jacobianDet", C, P, D);
          INTREPID2_TEST_ERROR_EXPECTED( ct::setJacobianDet(jacobianDet, jacobian) );

          // 18. Incompatible ranks
          jacobian = DynRankView("jacobian", P, D, D);
          jacobianDet = DynRankView("jacobianDet", C, P);
          INTREPID2_TEST_ERROR_EXPECTED( ct::setJacobianDet(jacobianDet, jacobian) );

          // 19. Incorrect rank
          jacobian = DynRankView("jacobian", D, D);
          jacobianDet = DynRankView("jacobianDet", C, P);
          INTREPID2_TEST_ERROR_EXPECTED( ct::setJacobianDet(jacobianDet, jacobian) );

          // 20. Incorrect rank
          jacobian = DynRankView("jacobian", C, P, D, D);
          jacobianDet = DynRankView("jacobianDet", C);
          INTREPID2_TEST_ERROR_EXPECTED( ct::setJacobianDet(jacobianDet, jacobian) );

          // 21. Incorrect dimension
          jacobian = DynRankView("jacobian", C, P, D, D);
          jacobianDet = DynRankView("jacobianDet", C, P-1);
          INTREPID2_TEST_ERROR_EXPECTED( ct::setJacobianDet(jacobianDet, jacobian) );

          // 22. Incorrect dimension
          jacobian = DynRankView("jacobian", C - 1, P, D, D);
          jacobianDet = DynRankView("jacobianDet", C, P);
          INTREPID2_TEST_ERROR_EXPECTED( ct::setJacobianDet(jacobianDet, jacobian) );

          /***********************************************************************************************
           *                        Exception tests for mapToPhysicalFrame method                       *
           **********************************************************************************************/

          // 23. Incorrect refPoint rank
          refPoints = DynRankView("refPoints", P);
          physPoints = DynRankView("physPoints", P, D);
          worksetCell = DynRankView("worksetCell", C, N, D);
          INTREPID2_TEST_ERROR_EXPECTED( ct::mapToPhysicalFrame(physPoints, refPoints, worksetCell, cell ) );

          // 24. Incorrect workset rank
          worksetCell = DynRankView("worksetCell", P, D);
          refPoints = DynRankView("refPoints", P, D);
          physPoints = DynRankView("physPoints", C, P, D);
          INTREPID2_TEST_ERROR_EXPECTED( ct::mapToPhysicalFrame(physPoints, refPoints, worksetCell, cell ) );
          
          // 25. Incompatible ranks
          refPoints = DynRankView("refPoints", C, P, D);
          physPoints = DynRankView("physPoints", P, D);
          worksetCell = DynRankView("worksetCell", C, N, D);
          INTREPID2_TEST_ERROR_EXPECTED( ct::mapToPhysicalFrame(physPoints, refPoints, worksetCell, cell ) );
          
          // 26. Incompatible dimensions
          refPoints = DynRankView("refPoints", C, P, D);
          physPoints = DynRankView("physPoints", C, P, D - 1);
          INTREPID2_TEST_ERROR_EXPECTED( ct::mapToPhysicalFrame(physPoints, refPoints, worksetCell, cell ) );
          
          // 27. Incompatible dimensions
          refPoints = DynRankView("refPoints", C, P, D);
          physPoints = DynRankView("physPoints", C, P - 1, D);
          INTREPID2_TEST_ERROR_EXPECTED( ct::mapToPhysicalFrame(physPoints, refPoints, worksetCell, cell ) );
          
          // 28. Incompatible dimensions
          refPoints = DynRankView("refPoints", C, P, D);
          physPoints = DynRankView("physPoints", C - 1, P, D);
          INTREPID2_TEST_ERROR_EXPECTED( ct::mapToPhysicalFrame(physPoints, refPoints, worksetCell, cell ) );
          
          /***********************************************************************************************
           *          Exception tests for mapToReferenceFrame method (with default initial guesses)     *
           **********************************************************************************************/
          
          // 30. incompatible ranks
          refPoints = DynRankView("refPoints", C, P, D);
          physPoints = DynRankView("physPoints", P, D);
          worksetCell = DynRankView("worksetCell", C, N, D);
          INTREPID2_TEST_ERROR_EXPECTED( ct::mapToReferenceFrame(refPoints, physPoints, worksetCell, cell ) );
          
          // 32. Incompatible ranks with whichCell = -1 (default)
          refPoints = DynRankView("refPoints", P, D);
          physPoints = DynRankView("physPoints", P, D);
          INTREPID2_TEST_ERROR_EXPECTED( ct::mapToReferenceFrame(refPoints, physPoints, worksetCell, cell ) );

          // 33. Nonmatching dimensions
          refPoints = DynRankView("refPoints", C, P, D - 1);
          physPoints = DynRankView("physPoints", C, P, D);
          INTREPID2_TEST_ERROR_EXPECTED( ct::mapToReferenceFrame(refPoints, physPoints, worksetCell, cell ) );

          // 34. Nonmatching dimensions
          refPoints = DynRankView("refPoints", C, P - 1, D);
          physPoints = DynRankView("physPoints", C, P, D);
          INTREPID2_TEST_ERROR_EXPECTED( ct::mapToReferenceFrame(refPoints, physPoints, worksetCell, cell ) );

          // 35. Nonmatching dimensions
          refPoints = DynRankView("refPoints", C - 1, P, D);
          physPoints = DynRankView("physPoints", C, P, D);
          INTREPID2_TEST_ERROR_EXPECTED( ct::mapToReferenceFrame(refPoints, physPoints, worksetCell, cell ) );

          // 36. Incorrect rank for worksetCell
          refPoints = DynRankView("refPoints", C, P, D);
          physPoints = DynRankView("physPoints", C, P, D);
          worksetCell = DynRankView("worksetCell", C, N);
          INTREPID2_TEST_ERROR_EXPECTED( ct::mapToReferenceFrame(refPoints, physPoints, worksetCell, cell ) );
          
          /***********************************************************************************************
           *   Exception tests for mapToReferenceFrameInitGuess method (initial guess is a parameter)   *
           **********************************************************************************************/
          
          // 37. Incompatible ranks
          refPoints = DynRankView("refPoints", C, P, D);
          physPoints = DynRankView("physPoints", C, P, D);
          initGuess = DynRankView("initGuess", P, D);
          worksetCell = DynRankView("worksetCell", C, N, D);
          INTREPID2_TEST_ERROR_EXPECTED( ct::mapToReferenceFrameInitGuess(refPoints, initGuess, physPoints, worksetCell, cell ) );

          // 39. Nonmatching dimensions
          refPoints = DynRankView("refPoints", C, P, D);
          physPoints = DynRankView("physPoints", C, P, D);
          initGuess = DynRankView("initGuess", C, P, D - 1);
          INTREPID2_TEST_ERROR_EXPECTED( ct::mapToReferenceFrameInitGuess(refPoints, initGuess, physPoints, worksetCell, cell ) );

          // 40. Nonmatching dimensions
          initGuess = DynRankView("initGuess", C, P - 1, D);
          INTREPID2_TEST_ERROR_EXPECTED( ct::mapToReferenceFrameInitGuess(refPoints, initGuess, physPoints, worksetCell, cell ) );

          // 41. Nonmatching dimensions
          initGuess = DynRankView("initGuess", C - 1, P, D);
          INTREPID2_TEST_ERROR_EXPECTED( ct::mapToReferenceFrameInitGuess(refPoints, initGuess, physPoints, worksetCell, cell ) );
          
          /***********************************************************************************************
           *                        Exception tests for mapToReferenceSubcell method                    *
           **********************************************************************************************/

          // 42. Incorrect array rank
          refSubcellPoints = DynRankView("refSubcellPoints", P,3);
          paramPoints = DynRankView("paramPoints", P);
          INTREPID2_TEST_ERROR_EXPECTED( ct::mapToReferenceSubcell(refSubcellPoints, paramPoints, subcellDim, subcellOrd, cell ) );

          // 43. Incorrect array rank
          refSubcellPoints = DynRankView("refSubcellPoints", P);
          paramPoints = DynRankView("paramPoints", P, 2);
          INTREPID2_TEST_ERROR_EXPECTED( ct::mapToReferenceSubcell(refSubcellPoints, paramPoints, subcellDim, subcellOrd, cell ) );

          // 44. Incorrect array dimension for face of 3D cell (should be 3)
          refSubcellPoints = DynRankView("refSubcellPoints", P, 2);
          paramPoints = DynRankView("paramPoints", P, 2);
          INTREPID2_TEST_ERROR_EXPECTED( ct::mapToReferenceSubcell(refSubcellPoints, paramPoints, subcellDim, subcellOrd, cell ) );
          
          // 45. Incorrect array dimension for parametrization domain of a face of 3D cell (should be 2)
          refSubcellPoints = DynRankView("refSubcellPoints", P, 3);
          paramPoints = DynRankView("paramPoints", P, 3);
          INTREPID2_TEST_ERROR_EXPECTED( ct::mapToReferenceSubcell(refSubcellPoints, paramPoints, subcellDim, subcellOrd, cell ) );

          /***********************************************************************************************
           *                        Exception tests for getReferenceEdgeTangent method                  *
           **********************************************************************************************/

          // 46. Incorrect rank
          refEdgeTangent = DynRankView("refEdgeTangent", C,P,D);
          INTREPID2_TEST_ERROR_EXPECTED( ct::getReferenceEdgeTangent(refEdgeTangent, 0, cell ) );

          // 47. Incorrect dimension D for Tet<10> cell
          refEdgeTangent = DynRankView("refEdgeTangent", 2);
          INTREPID2_TEST_ERROR_EXPECTED( ct::getReferenceEdgeTangent(refEdgeTangent, 0, cell ) );

          // 48. Invalid edge ordinal for Tet<10>
          refEdgeTangent = DynRankView("refEdgeTangent", C,P,D);
          INTREPID2_TEST_ERROR_EXPECTED( ct::getReferenceEdgeTangent(refEdgeTangent, 10, cell ) );

          /***********************************************************************************************
           *                        Exception tests for getReferenceFaceTangents method                 *
           **********************************************************************************************/

          // 49. Incorrect rank
          refFaceTanU = DynRankView("refFaceTanU", P, D);
          refFaceTanV = DynRankView("refFaceTanV", D);
          INTREPID2_TEST_ERROR_EXPECTED( ct::getReferenceFaceTangents(refFaceTanU, refFaceTanV, 0, cell ) );

          // 50. Incorrect rank
          refFaceTanU = DynRankView("refFaceTanU", D);
          refFaceTanV = DynRankView("refFaceTanV", P, D);
          INTREPID2_TEST_ERROR_EXPECTED( ct::getReferenceFaceTangents(refFaceTanU, refFaceTanV, 0, cell ) );

          // 51. Incorrect dimension for 3D cell
          refFaceTanU = DynRankView("refFaceTanU", D - 1);
          refFaceTanV = DynRankView("refFaceTanV", D);
          INTREPID2_TEST_ERROR_EXPECTED( ct::getReferenceFaceTangents(refFaceTanU, refFaceTanV, 0, cell ) );

          // 52. Incorrect dimension for 3D cell
          refFaceTanU = DynRankView("refFaceTanU", D);
          refFaceTanV = DynRankView("refFaceTanV", D - 1);
          INTREPID2_TEST_ERROR_EXPECTED( ct::getReferenceFaceTangents(refFaceTanU, refFaceTanV, 0, cell ) );

          // 53. Invalid face ordinal
          refFaceTanU = DynRankView("refFaceTanU", D);
          refFaceTanV = DynRankView("refFaceTanV", D);
          INTREPID2_TEST_ERROR_EXPECTED( ct::getReferenceFaceTangents(refFaceTanU, refFaceTanV, 10, cell ) );

          /***********************************************************************************************
           *                        Exception tests for getReferenceSide/FaceNormal methods             *
           **********************************************************************************************/

          // 54-55. Incorrect rank
          refSideNormal = DynRankView("refSideNormal", C,P);
          INTREPID2_TEST_ERROR_EXPECTED( ct::getReferenceSideNormal(refSideNormal, 0, cell ) );

          // 56-57. Incorrect dimension for 3D cell
          refSideNormal = DynRankView("refSideNormal", D - 1);
          INTREPID2_TEST_ERROR_EXPECTED( ct::getReferenceSideNormal(refSideNormal, 0, cell ) );

          // 58-59. Invalid side ordinal for Tet<10>
          refSideNormal = DynRankView("refSideNormal", D);
          INTREPID2_TEST_ERROR_EXPECTED( ct::getReferenceSideNormal(refSideNormal, 10, cell ) );

          // 60. Incorrect dimension for 2D cell
          refSideNormal = DynRankView("refSideNormal", D - 1);
          INTREPID2_TEST_ERROR_EXPECTED( ct::getReferenceSideNormal(refSideNormal, 0, cell ) );

          // 61. Invalid side ordinal for Tri<3>
          refSideNormal = DynRankView("refSideNormal", D);
          INTREPID2_TEST_ERROR_EXPECTED( ct::getReferenceSideNormal(refSideNormal, 10, cell ) );

          // 62. Cannot call the "face" method for 2D cells
          refSideNormal = DynRankView("refSideNormal", D);
          INTREPID2_TEST_ERROR_EXPECTED( ct::getReferenceFaceNormal(refSideNormal, 0, cell ) );

          /***********************************************************************************************
           *          Exception tests for checkPoint/Pointset/PointwiseInclusion methods        *
           **********************************************************************************************/
          // points.resize(2,3,3,4);
          // FieldContainer<int> inCell;

          // 63. Point dimension does not match cell topology
          points = DynRankView("points", D+1);
          INTREPID2_TEST_ERROR_EXPECTED( ct::checkPointInclusion(points, cell ) );

          // // 64. Invalid cell topology
          // CellTopology pentagon_5(shards::getCellTopologyData<shards::Pentagon<> >() );
          // INTREPID2_TEST_ERROR_EXPECTED(ct::checkPointInclusion(point, pentagon_5.getDimension(), pentagon_5 ),
          //                               throwCounter, nException );

          // // 65. Incorrect spatial dimension of points
          // points.resize(10, 10, cell.getDimension() + 1);
          // INTREPID2_TEST_ERROR_EXPECTED(ct::checkPointsetInclusion(points, cell ),
          //                               throwCounter, nException );

          // // 66. Incorrect rank of input array
          // points.resize(10,10,10,3);
          // INTREPID2_TEST_ERROR_EXPECTED(ct::checkPointsetInclusion(points, cell ),
          //                               throwCounter, nException );

          // // 67. Incorrect rank of output array
          // points.resize(10,10,cell.getDimension() );
          // inCell.resize(10);
          // INTREPID2_TEST_ERROR_EXPECTED(ct::checkPointwiseInclusion(inCell, points, cell ),
          //                               throwCounter, nException );

          // // 68. Incorrect rank of output array
          // points.resize(10, cell.getDimension() );
          // inCell.resize(10, 10);
          // INTREPID2_TEST_ERROR_EXPECTED(ct::checkPointwiseInclusion(inCell, points, cell ),
          //                               throwCounter, nException );

          // // 69. Incorrect rank of output array
          // points.resize(cell.getDimension() );
          // inCell.resize(10, 10);
          // INTREPID2_TEST_ERROR_EXPECTED(ct::checkPointwiseInclusion(inCell, points, cell ),
          //                               throwCounter, nException );

          // // 70. Incorrect dimension of output array
          // points.resize(10, 10, cell.getDimension() );
          // inCell.resize(10, 9);
          // INTREPID2_TEST_ERROR_EXPECTED(ct::checkPointwiseInclusion(inCell, points, cell ),
          //                               throwCounter, nException );

          // // 71. Incorrect dimension of output array
          // points.resize(10, 10, cell.getDimension() );
          // inCell.resize(9, 10);
          // INTREPID2_TEST_ERROR_EXPECTED(ct::checkPointwiseInclusion(inCell, points, cell ),
          //                               throwCounter, nException );

          // // 72. Incorrect dimension of output array
          // points.resize(10, cell.getDimension() );
          // inCell.resize(9);
          // INTREPID2_TEST_ERROR_EXPECTED(ct::checkPointwiseInclusion(inCell, points, cell ),
          //                               throwCounter, nException );

          // // 73. Incorrect spatial dimension of input array
          // points.resize(10, 10, cell.getDimension() + 1);
          // inCell.resize(10, 10);
          // INTREPID2_TEST_ERROR_EXPECTED(ct::checkPointwiseInclusion(inCell, points, cell ),
          //                               throwCounter, nException );

          // // 74. Incorrect rank of input array.
          // points.resize(10,10,10,3);
          // INTREPID2_TEST_ERROR_EXPECTED(ct::checkPointwiseInclusion(inCell, points, cell ),
          //                               throwCounter, nException );


          // physPoints = DynRankView("physPoints", C, P, D);
          // inCell.resize(C, P);
          // // 75. Invalid rank of worksetCell
          // worksetCell = DynRankView("worksetCell", C, N, D, D);
          // INTREPID2_TEST_ERROR_EXPECTED(ct::checkPointwiseInclusion(inCell, physPoints, worksetCell, cell ),
          //                               throwCounter, nException );

          // // 76. Invalid dimension 1 (node count) of worksetCell
          // worksetCell = DynRankView("worksetCell", C, N + 1, D);
          // INTREPID2_TEST_ERROR_EXPECTED(ct::checkPointwiseInclusion(inCell, physPoints, worksetCell, cell ),
          //                               throwCounter, nException );

          // // 77. Invalid dimension 2 (spatial dimension) of worksetCell
          // worksetCell = DynRankView("worksetCell", C, N, D + 1);
          // INTREPID2_TEST_ERROR_EXPECTED(ct::checkPointwiseInclusion(inCell, physPoints, worksetCell, cell ),
          //                               throwCounter, nException );

          // // 78. Invalid whichCell value (exceeds cell count in the workset)
          // worksetCell = DynRankView("worksetCell", C, N, D);
          // INTREPID2_TEST_ERROR_EXPECTED(ct::checkPointwiseInclusion(inCell, physPoints, worksetCell, cell, C + 1 ),
          //                               throwCounter, nException );

          // // 79. Invalid whichCell for rank-3 physPoints (must be -1, here it is valid cell ordinal)
          // INTREPID2_TEST_ERROR_EXPECTED(ct::checkPointwiseInclusion(inCell, physPoints, worksetCell, cell, 0 ),
          //                               throwCounter, nException );

          // // 80. Invalid whichCell for rank-2 physPoints (must be a valid cell ordinal, here it is the default -1)
          // physPoints = DynRankView("physPoints", P, D);
          // inCell.resize(P);
          // INTREPID2_TEST_ERROR_EXPECTED(ct::checkPointwiseInclusion(inCell, physPoints, worksetCell, cell ),
          //                               throwCounter, nException );

          // // 81. Incompatible ranks of I/O arrays
          // physPoints = DynRankView("physPoints", C, P, D);
          // inCell.resize(P);
          // INTREPID2_TEST_ERROR_EXPECTED(ct::checkPointwiseInclusion(inCell, physPoints, worksetCell, cell ) );
          //                               throwCounter, nException );

          // // 82. Incompatible ranks of I/O arrays
          // physPoints = DynRankView("physPoints", P, D);
          // inCell.resize(C, P);
          // INTREPID2_TEST_ERROR_EXPECTED(ct::checkPointwiseInclusion(inCell, physPoints, worksetCell, cell, 0),
          //                               throwCounter, nException );

          // // 83. Incompatible dimensions of I/O arrays
          // physPoints = DynRankView("physPoints", C, P, D);
          // inCell.resize(C, P + 1);
          // INTREPID2_TEST_ERROR_EXPECTED(ct::checkPointwiseInclusion(inCell, physPoints, worksetCell, cell ) );
          //                               throwCounter, nException );

          // // 84. Incompatible dimensions of I/O arrays: rank-3 Input
          // physPoints = DynRankView("physPoints", C + 1, P, D);
          // inCell.resize(C, P);
          // INTREPID2_TEST_ERROR_EXPECTED(ct::checkPointwiseInclusion(inCell, physPoints, worksetCell, cell ) );
          //                               throwCounter, nException );

          // // 85. Incompatible dimensions of I/O arrays: rank-2 Input
          // physPoints = DynRankView("physPoints", P, D);
          // inCell.resize(P + 1);
          // INTREPID2_TEST_ERROR_EXPECTED(ct::checkPointwiseInclusion(inCell, physPoints, worksetCell, cell, 0 ),
          //                               throwCounter, nException );


          // /***********************************************************************************************
          //  *               Exception tests for getReferenceVertex/vertices/Node/Nodes methods           *
          //  **********************************************************************************************/

          // FieldContainer<double> subcellNodes;

          // // 86-89. Cell does not have reference cell
          // INTREPID2_TEST_ERROR_EXPECTED(ct::getReferenceVertex(pentagon_5, 0), throwCounter, nException);
          // INTREPID2_TEST_ERROR_EXPECTED(ct::getReferenceNode(pentagon_5, 0), throwCounter, nException);
          // INTREPID2_TEST_ERROR_EXPECTED(ct::getReferenceSubcellVertices(subcellNodes, 0, 0, pentagon_5), throwCounter, nException);
          // INTREPID2_TEST_ERROR_EXPECTED(ct::getReferenceSubcellNodes(subcellNodes, 0, 0, pentagon_5), throwCounter, nException);

          // // Use last cell topology (Wedge<18>) for these tests
          // topo_iterator = supportedTopologies.end() - 1;
          // D = cell.getDimension();
          // int subcDim = D - 1;
          // int S = cell.getSubcellCount(subcDim);
          // V = cell.getVertexCount(subcDim, S - 1);
          // subcellNodes.resize(V, D);
          // // 90. subcell ordinal out of range
          // INTREPID2_TEST_ERROR_EXPECTED(ct::getReferenceSubcellVertices(subcellNodes, subcDim, S + 1, cell ) );
          //                               throwCounter, nException);

          // // 91. subcell dim out of range
          // INTREPID2_TEST_ERROR_EXPECTED(ct::getReferenceSubcellVertices(subcellNodes, D + 1, S, cell ) );
          //                               throwCounter, nException);

          // // 92. Incorrect rank for subcellNodes
          // subcellNodes.resize(V, D, D);
          // INTREPID2_TEST_ERROR_EXPECTED(ct::getReferenceSubcellVertices(subcellNodes, subcDim, S - 1, cell ) );
          //                               throwCounter, nException);

          // // 93. Incorrect dimension for subcellNodes
          // subcellNodes.resize(V - 1, D);
          // INTREPID2_TEST_ERROR_EXPECTED(ct::getReferenceSubcellVertices(subcellNodes, subcDim, S - 1, cell ) );
          //                               throwCounter, nException);

          // // 94. Incorrect dimension for subcellNodes
          // subcellNodes.resize(V, D - 1);
          // INTREPID2_TEST_ERROR_EXPECTED(ct::getReferenceSubcellVertices(subcellNodes, subcDim, S - 1, cell ) );
          //                               throwCounter, nException);


          // N = cell.getNodeCount(subcDim, S - 1);
          // subcellNodes.resize(N, D);
          // // 95. subcell ordinal out of range
          // INTREPID2_TEST_ERROR_EXPECTED(ct::getReferenceSubcellNodes(subcellNodes, subcDim, S + 1, cell ) );
          //                               throwCounter, nException);

          // // 96. subcell dim out of range
          // INTREPID2_TEST_ERROR_EXPECTED(ct::getReferenceSubcellNodes(subcellNodes, D + 1, S, cell ) );
          //                               throwCounter, nException);

          // // 97. Incorrect rank for subcellNodes
          // subcellNodes.resize(N, D, D);
          // INTREPID2_TEST_ERROR_EXPECTED(ct::getReferenceSubcellNodes(subcellNodes, subcDim, S - 1, cell ) );
          //                               throwCounter, nException);

          // // 98. Incorrect dimension for subcellNodes
          // subcellNodes.resize(N - 1, D);
          // INTREPID2_TEST_ERROR_EXPECTED(ct::getReferenceSubcellNodes(subcellNodes, subcDim, S - 1, cell ) );
          //                               throwCounter, nException);

          // // 99. Incorrect dimension for subcellNodes
          // subcellNodes.resize(N, D - 1);
          // INTREPID2_TEST_ERROR_EXPECTED(ct::getReferenceSubcellNodes(subcellNodes, subcDim, S - 1, cell ) );
          //                               throwCounter, nException);

          // something wrong in counting
          // if (nthrow != ncatch) {
          //   errorFlag++;
          //   *outStream << " nthrow = " << nthrow << "  ncatch = " << ncatch << std::endl;
          //   *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
          // }
        }
#endif
      } catch (std::exception &err) {
        *outStream << "Unexpected Error = " << err.what() << "\n";
        errorFlag = -1000;
      }
      
      if (errorFlag != 0)
        std::cout << "End Result: TEST FAILED\n";
      else
        std::cout << "End Result: TEST PASSED\n";
      
      //reset format state of std::cout
      std::cout.copyfmt(oldFormatState);
      
      return errorFlag;
    }
  }
}






