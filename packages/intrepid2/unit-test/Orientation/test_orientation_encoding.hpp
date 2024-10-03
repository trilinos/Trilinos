// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


/** \file
    \brief  Unit test of experimental high order assembly
    \author Created by Kyungjoo Kim
*/

#include "Intrepid2_config.h"

#ifdef HAVE_INTREPID2_DEBUG
#define INTREPID2_TEST_FOR_DEBUG_ABORT_OVERRIDE_TO_CONTINUE
#endif

#include "Intrepid2_Orientation.hpp"

#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"

namespace Intrepid2 {

  namespace Test {

#define INTREPID2_TEST_ERROR_EXPECTED( S )              \
    try {                                                               \
      ++nthrow;                                                         \
      S ;                                                               \
    }                                                                   \
    catch (std::exception &err) {                                        \
      ++ncatch;                                                         \
      *outStream << "Expected Error ----------------------------------------------------------------\n"; \
      *outStream << err.what() << '\n';                                 \
      *outStream << "-------------------------------------------------------------------------------" << "\n\n"; \
    }
    
    template<typename DeviceType>
    int OrientationEncoding(const bool verbose) {

      Teuchos::RCP<std::ostream> outStream;
      Teuchos::oblackholestream bhs; // outputs nothing

      if (verbose)
        outStream = Teuchos::rcp(&std::cout, false);
      else
        outStream = Teuchos::rcp(&bhs,       false);

      Teuchos::oblackholestream oldFormatState;
      oldFormatState.copyfmt(std::cout);

      *outStream
        << "===============================================================================\n"
        << "|                                                                             |\n"
        << "|                 Unit Test (Orientation - encoding/decoding)                 |\n"
        << "|                                                                             |\n"
        << "===============================================================================\n";

      int errorFlag = 0;

      try {
        ordinal_type nthrow = 0, ncatch = 0;
        {
          Orientation ort;
          INTREPID2_TEST_ERROR_EXPECTED( if (ort.isAlignedToReference()) \
                                           throw std::logic_error("Default Orientation is not zero"); );
        }

        if (nthrow != ncatch) {
          errorFlag++;
          *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
          *outStream << "# of catch ("<< ncatch << ") is different from # of throw (" << nthrow << ")\n";
        }

        {
          *outStream << "\n -- Testing Triangle \n\n";

          const auto cellTopo = shards::CellTopology(shards::getCellTopologyData<shards::Triangle<3> >() );
          const ordinal_type elemNodes[6][3] = { { 1, 2, 3 },
                                                 { 1, 3, 2 },
                                                 { 2, 1, 3 },
                                                 { 2, 3, 1 },
                                                 { 3, 1, 2 },
                                                 { 3, 2, 1 } };

          const ordinal_type refEdgeOrts[6][3] = { { 0, 0, 1 },
                                                   { 0, 1, 1 },
                                                   { 1, 0, 1 },
                                                   { 0, 1, 0 },
                                                   { 1, 0, 0 },
                                                   { 1, 1, 0 } };
          
          *outStream << "Triangle element edge reference configuration\n";
          for (auto edgeId=0;edgeId<3;++edgeId) {
            const auto v0 = cellTopo.getNodeMap(1, edgeId, 0);
            const auto v1 = cellTopo.getNodeMap(1, edgeId, 1);
            *outStream << std::setw(3) << edgeId << " edge :: " << v0 << " " << v1 << "\n";
          }
          *outStream << "Triangle element edge orientation with vertex permutations\n";          
          for (auto i=0;i<6;++i) {
            // find orientation
            const auto nodes = Kokkos::View<const ordinal_type[3],Kokkos::HostSpace>(elemNodes[i]);
            const auto ort = Orientation::getOrientation(cellTopo, nodes);
            
            // decode orientation
            ordinal_type edgeOrt[3] = {};
            for (auto edgeId=0;edgeId<3;++edgeId) 
              ort.getEdgeOrientation(edgeOrt, 3);
            
            *outStream << " elemNodes = " 
                       << elemNodes[i][0] << " " 
                       << elemNodes[i][1] << " " 
                       << elemNodes[i][2] << " :: "
                       << " computed edgeOrts = " 
                       << edgeOrt[0] << " "
                       << edgeOrt[1] << " "
                       << edgeOrt[2] << " :: "
                       << " reference edgeOrts = " 
                       << refEdgeOrts[i][0] << " "
                       << refEdgeOrts[i][1] << " "
                       << refEdgeOrts[i][2] << " ::\n";
            
            if (edgeOrt[0] != refEdgeOrts[i][0] || 
                edgeOrt[1] != refEdgeOrts[i][1] || 
                edgeOrt[2] != refEdgeOrts[i][2])  {
              *outStream << "        ^^^^^^^^^^^^^^^^ FAILURE\n";
              ++errorFlag;
            }
          }
        }

        {
          *outStream << "\n -- Testing Quadrilateral \n\n";
          
          const auto cellTopo = shards::CellTopology(shards::getCellTopologyData<shards::Quadrilateral<4> >() );
          const ordinal_type elemNodes[24][4] = { { 1 ,  2 ,  3 ,  4 },
                                                  { 2 ,  1 ,  3 ,  4 },
                                                  { 1 ,  3 ,  2 ,  4 },
                                                  { 2 ,  3 ,  1 ,  4 },
                                                  { 3 ,  1 ,  2 ,  4 },
                                                  { 3 ,  2 ,  1 ,  4 },
                                                  { 1 ,  2 ,  4 ,  3 },
                                                  { 2 ,  1 ,  4 ,  3 },
                                                  { 1 ,  3 ,  4 ,  2 },
                                                  { 2 ,  3 ,  4 ,  1 },
                                                  { 3 ,  1 ,  4 ,  2 },
                                                  { 3 ,  2 ,  4 ,  1 },
                                                  { 1 ,  4 ,  2 ,  3 },
                                                  { 2 ,  4 ,  1 ,  3 },
                                                  { 1 ,  4 ,  3 ,  2 },
                                                  { 2 ,  4 ,  3 ,  1 },
                                                  { 3 ,  4 ,  1 ,  2 },
                                                  { 3 ,  4 ,  2 ,  1 },
                                                  { 4 ,  1 ,  2 ,  3 },
                                                  { 4 ,  2 ,  1 ,  3 },
                                                  { 4 ,  1 ,  3 ,  2 },
                                                  { 4 ,  2 ,  3 ,  1 },
                                                  { 4 ,  3 ,  1 ,  2 },
                                                  { 4 ,  3 ,  2 ,  1 } };
          
          const ordinal_type refEdgeOrts[24][4] = { { 0, 0, 0, 1 },
                                                    { 1, 0, 0, 1 },
                                                    { 0, 1, 0, 1 },
                                                    { 0, 1, 0, 1 },
                                                    { 1, 0, 0, 1 },
                                                    { 1, 1, 0, 1 },
                                                    { 0, 0, 1, 1 },
                                                    { 1, 0, 1, 1 },
                                                    { 0, 0, 1, 1 },
                                                    { 0, 0, 1, 0 },
                                                    { 1, 0, 1, 0 },
                                                    { 1, 0, 1, 0 },
                                                    { 0, 1, 0, 1 },
                                                    { 0, 1, 0, 1 },
                                                    { 0, 1, 1, 1 },
                                                    { 0, 1, 1, 0 },
                                                    { 0, 1, 0, 0 },
                                                    { 0, 1, 1, 0 },
                                                    { 1, 0, 0, 0 },
                                                    { 1, 1, 0, 0 },
                                                    { 1, 0, 1, 0 },
                                                    { 1, 0, 1, 0 },
                                                    { 1, 1, 0, 0 },
                                                    { 1, 1, 1, 0 } };

          *outStream << "Quadrilateral element edge reference configuration\n";
          for (auto edgeId=0;edgeId<4;++edgeId) {
            const auto v0 = cellTopo.getNodeMap(1, edgeId, 0);
            const auto v1 = cellTopo.getNodeMap(1, edgeId, 1);
            *outStream << std::setw(3) << edgeId << " edge :: " << v0 << " " << v1 << "\n";
          }
          *outStream << "Quadrilateral element edge orientation with vertex permutations\n";          
          for (auto i=0;i<24;++i) {
            // find orientation
            const auto nodes = Kokkos::View<const ordinal_type[4],Kokkos::HostSpace>(elemNodes[i]);
            const auto ort = Orientation::getOrientation(cellTopo, nodes);
            
            // decode orientation
            ordinal_type edgeOrt[4] = {};
            for (auto edgeId=0;edgeId<4;++edgeId) 
              ort.getEdgeOrientation(edgeOrt, 4);
            
            *outStream << " elemNodes = " 
                       << elemNodes[i][0] << " " 
                       << elemNodes[i][1] << " " 
                       << elemNodes[i][2] << " " 
                       << elemNodes[i][3] << " :: "
                       << " computed edgeOrts = " 
                       << edgeOrt[0] << " "
                       << edgeOrt[1] << " "
                       << edgeOrt[2] << " "
                       << edgeOrt[3] << " :: "
                       << " reference edgeOrts = " 
                       << refEdgeOrts[i][0] << " "
                       << refEdgeOrts[i][1] << " "
                       << refEdgeOrts[i][2] << " "
                       << refEdgeOrts[i][3] << " ::\n";
            
            if (edgeOrt[0] != refEdgeOrts[i][0] || 
                edgeOrt[1] != refEdgeOrts[i][1] || 
                edgeOrt[2] != refEdgeOrts[i][2] ||
                edgeOrt[3] != refEdgeOrts[i][3])  {
              *outStream << "        ^^^^^^^^^^^^^^^^ FAILURE\n";
              ++errorFlag;
            }
          }
        }

        {
          *outStream << "\n -- Testing hexahedron \n\n";
          // select following permutation order , if one wants to test all possible cases, it will be 40320.
          const auto cellTopo = shards::CellTopology(shards::getCellTopologyData<shards::Hexahedron<8> >() );
          const ordinal_type elemNodes[24][8] = {  {   1,   2,   3,   4,   5,   8,   6,   7  },
                                                   {   2,   1,   3,   4,   6,   8,   5,   7  },
                                                   {   1,   3,   2,   4,   5,   8,   7,   6  },
                                                   {   2,   3,   1,   4,   6,   8,   7,   5  },
                                                                                          
                                                   {   3,   1,   2,   4,   7,   8,   5,   6  },
                                                   {   3,   2,   1,   4,   7,   8,   6,   5  },
                                                   {   1,   2,   4,   3,   8,   5,   6,   7  },
                                                   {   2,   1,   4,   3,   8,   6,   5,   7  },
                                                                                          
                                                   {   1,   3,   4,   2,   8,   5,   7,   6  },
                                                   {   2,   3,   4,   1,   8,   6,   7,   5  },
                                                   {   3,   1,   4,   2,   8,   7,   5,   6  },
                                                   {   3,   2,   4,   1,   8,   7,   6,   5  },
                                                   
                                                   {   1,   4,   2,   3,   5,   6,   7,   8  },
                                                   {   2,   4,   1,   3,   6,   5,   7,   8  },
                                                   {   1,   4,   3,   2,   5,   7,   6,   8  },
                                                   {   2,   4,   3,   1,   6,   7,   5,   8  },
                                                                                           
                                                   {   3,   4,   1,   2,   7,   5,   6,   8  },
                                                   {   3,   4,   2,   1,   7,   6,   5,   8  },
                                                   {   4,   1,   2,   3,   5,   6,   8,   7  },
                                                   {   4,   2,   1,   3,   6,   5,   8,   7  },
                                                                                           
                                                   {   4,   1,   3,   2,   5,   7,   8,   6  },
                                                   {   4,   2,   3,   1,   6,   7,   8,   5  },
                                                   {   4,   3,   1,   2,   7,   5,   8,   6  },
                                                   {   4,   3,   2,   1,   7,   6,   8,   5  } };

          const ordinal_type refEdgeOrts[24][12] = { { 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0 },
                                                     { 1, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0 },
                                                     { 0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0 },
                                                     { 0, 1, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0 },
                                                                                         
                                                     { 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0 },
                                                     { 1, 1, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0 },
                                                     { 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0 },
                                                     { 1, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0 },
                                                                                         
                                                     { 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0 },
                                                     { 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0 },
                                                     { 1, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0 },
                                                     { 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0 },
                                                                                         
                                                     { 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0 },
                                                     { 0, 1, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0 },
                                                     { 0, 1, 1, 1, 0, 1, 0, 1, 0, 0, 0, 0 },
                                                     { 0, 1, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0 },
                                                                                         
                                                     { 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0 },
                                                     { 0, 1, 1, 0, 1, 1, 0, 1, 0, 0, 0, 0 },
                                                     { 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0 },
                                                     { 1, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0 },
                                                                                         
                                                     { 1, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0 },
                                                     { 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0 },
                                                     { 1, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0 },
                                                     { 1, 1, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0 } };

          const ordinal_type refFaceOrts[24][6] = { {  0, 0, 0, 4, 4, 4  },
                                                    {  5, 0, 0, 4, 3, 2  },
                                                    {  0, 5, 0, 4, 4, 4  },
                                                    {  0, 5, 0, 4, 2, 3  },
                                                                     
                                                    {  5, 0, 0, 4, 7, 2  },
                                                    {  5, 5, 0, 4, 2, 7  },
                                                    {  0, 0, 5, 4, 4, 1  },
                                                    {  5, 0, 5, 4, 3, 6  },
                                                                     
                                                    {  0, 0, 5, 4, 0, 1  },
                                                    {  0, 0, 5, 3, 5, 7  },
                                                    {  5, 0, 5, 3, 3, 2  },
                                                    {  5, 0, 5, 3, 5, 7  },
                                                                     
                                                    {  0, 5, 0, 4, 0, 0  },
                                                    {  0, 5, 0, 4, 6, 5  },
                                                    {  0, 5, 5, 4, 0, 0  },
                                                    {  0, 5, 5, 3, 5, 6  },
                                                                     
                                                    {  0, 5, 0, 3, 6, 1  },
                                                    {  0, 5, 5, 3, 1, 6  },
                                                    {  5, 0, 0, 3, 7, 0  },
                                                    {  5, 5, 0, 3, 2, 5  },
                                                                     
                                                    {  5, 0, 5, 3, 7, 4  },
                                                    {  5, 0, 5, 3, 1, 3  },
                                                    {  5, 5, 0, 3, 6, 5  },
                                                    {  5, 5, 5, 3, 1, 3  } };
                                                                             

          *outStream << "Hexahedral element edge reference configuration\n";
          for (auto edgeId=0;edgeId<12;++edgeId) {
            const auto v0 = cellTopo.getNodeMap(1, edgeId, 0);
            const auto v1 = cellTopo.getNodeMap(1, edgeId, 1);
            *outStream << std::setw(3) << edgeId << " edge :: " << v0 << " " << v1 << "\n";
          }
          *outStream << "Hexahedral element edge orientation with vertex permutations\n";          
          for (auto i=0;i<24;++i) {                                          
            // find orientation                                              
            const auto nodes = Kokkos::View<const ordinal_type[8],Kokkos::HostSpace>(elemNodes[i]);
            const auto ort = Orientation::getOrientation(cellTopo, nodes);
            
            // decode orientation
            ordinal_type edgeOrt[12] = {};
            for (auto edgeId=0;edgeId<12;++edgeId) 
              ort.getEdgeOrientation(edgeOrt, 12);
            
            *outStream << " elemNodes = " 
                       << elemNodes[i][0] << " " 
                       << elemNodes[i][1] << " " 
                       << elemNodes[i][2] << " " 
                       << elemNodes[i][3] << " "
                       << elemNodes[i][4] << " "
                       << elemNodes[i][5] << " "
                       << elemNodes[i][6] << " "
                       << elemNodes[i][7] << " :: "
                       << " computed edgeOrts = " 
                       << edgeOrt[0] << " "
                       << edgeOrt[1] << " "
                       << edgeOrt[2] << " "
                       << edgeOrt[3] << " "
                       << edgeOrt[4] << " "
                       << edgeOrt[5] << " "
                       << edgeOrt[6] << " "
                       << edgeOrt[7] << " "
                       << edgeOrt[8] << " "
                       << edgeOrt[9] << " "
                       << edgeOrt[10] << " "
                       << edgeOrt[11] << " :: "
                       << " reference edgeOrts = " 
                       << refEdgeOrts[i][0] << " "
                       << refEdgeOrts[i][1] << " "
                       << refEdgeOrts[i][2] << " "
                       << refEdgeOrts[i][3] << " "
                       << refEdgeOrts[i][4] << " "
                       << refEdgeOrts[i][5] << " "
                       << refEdgeOrts[i][6] << " "
                       << refEdgeOrts[i][7] << " "
                       << refEdgeOrts[i][8] << " "
                       << refEdgeOrts[i][9] << " "
                       << refEdgeOrts[i][10] << " "
                       << refEdgeOrts[i][11] << " ::\n";
            
            if (edgeOrt[0] != refEdgeOrts[i][0] || 
                edgeOrt[1] != refEdgeOrts[i][1] || 
                edgeOrt[2] != refEdgeOrts[i][2] ||
                edgeOrt[3] != refEdgeOrts[i][3] ||
                edgeOrt[4] != refEdgeOrts[i][4] || 
                edgeOrt[5] != refEdgeOrts[i][5] || 
                edgeOrt[6] != refEdgeOrts[i][6] ||
                edgeOrt[7] != refEdgeOrts[i][7] ||
                edgeOrt[8] != refEdgeOrts[i][8] || 
                edgeOrt[9] != refEdgeOrts[i][9] || 
                edgeOrt[10] != refEdgeOrts[i][10] ||
                edgeOrt[11] != refEdgeOrts[i][11])  {
              *outStream << "        ^^^^^^^^^^^^^^^^ FAILURE\n";
              ++errorFlag;
            }
          }
          *outStream << "Hexahedral element face reference configuration\n";
          for (auto faceId=0;faceId<6;++faceId) {
            const auto v0 = cellTopo.getNodeMap(2, faceId, 0);
            const auto v1 = cellTopo.getNodeMap(2, faceId, 1);
            const auto v2 = cellTopo.getNodeMap(2, faceId, 2);
            const auto v3 = cellTopo.getNodeMap(2, faceId, 3);
            *outStream << std::setw(3) << faceId << " face :: " << v0 << " " << v1 << " " << v2 << " " << v3 << "\n";
          }
          *outStream << "Hexahedral element face orientation with vertex permutations\n";          
          for (auto i=0;i<24;++i) {                                          
            // find orientation                                              
            const auto nodes = Kokkos::View<const ordinal_type[8],Kokkos::HostSpace>(elemNodes[i]);
            const auto ort = Orientation::getOrientation(cellTopo, nodes);
            
            // decode orientation
            ordinal_type faceOrt[6] = {};
            for (auto faceId=0;faceId<6;++faceId) 
              ort.getFaceOrientation(faceOrt, 6);
            
            *outStream << " elemNodes = " 
                       << elemNodes[i][0] << " " 
                       << elemNodes[i][1] << " " 
                       << elemNodes[i][2] << " " 
                       << elemNodes[i][3] << " "
                       << elemNodes[i][4] << " "
                       << elemNodes[i][5] << " "
                       << elemNodes[i][6] << " "
                       << elemNodes[i][7] << " :: "
                       << " computed faceOrts = " 
                       << faceOrt[0] << " "
                       << faceOrt[1] << " "
                       << faceOrt[2] << " "
                       << faceOrt[3] << " "
                       << faceOrt[4] << " "
                       << faceOrt[5] << " :: "
                       << " reference faceOrts = " 
                       << refFaceOrts[i][0] << " "
                       << refFaceOrts[i][1] << " "
                       << refFaceOrts[i][2] << " "
                       << refFaceOrts[i][3] << " "
                       << refFaceOrts[i][4] << " "
                       << refFaceOrts[i][5] << " ::\n";
            
            if (faceOrt[0] != refFaceOrts[i][0] || 
                faceOrt[1] != refFaceOrts[i][1] || 
                faceOrt[2] != refFaceOrts[i][2] ||
                faceOrt[3] != refFaceOrts[i][3] ||
                faceOrt[4] != refFaceOrts[i][4] || 
                faceOrt[5] != refFaceOrts[i][5]) {
              *outStream << "        ^^^^^^^^^^^^^^^^ FAILURE\n";
              ++errorFlag;
            }
          }
        }

        {
          *outStream << "\n -- Testing tetrahedron \n\n";
          const auto cellTopo = shards::CellTopology(shards::getCellTopologyData<shards::Tetrahedron<4> >() );
          const ordinal_type elemNodes[24][4] = {  {   1,   2,   3,   4  },
                                                   {   2,   1,   3,   4  },
                                                   {   1,   3,   2,   4  },
                                                   {   2,   3,   1,   4  },
                                                                       
                                                   {   3,   1,   2,   4  },
                                                   {   3,   2,   1,   4  },
                                                   {   1,   2,   4,   3  },
                                                   {   2,   1,   4,   3  },
                                                                       
                                                   {   1,   3,   4,   2  },
                                                   {   2,   3,   4,   1  },
                                                   {   3,   1,   4,   2  },
                                                   {   3,   2,   4,   1  },
                                                   
                                                   {   1,   4,   2,   3  },
                                                   {   2,   4,   1,   3  },
                                                   {   1,   4,   3,   2  },
                                                   {   2,   4,   3,   1  },
                                                                       
                                                   {   3,   4,   1,   2  },
                                                   {   3,   4,   2,   1  },
                                                   {   4,   1,   2,   3  },
                                                   {   4,   2,   1,   3  },
                                                                       
                                                   {   4,   1,   3,   2  },
                                                   {   4,   2,   3,   1  },
                                                   {   4,   3,   1,   2  },
                                                   {   4,   3,   2,   1  } };

          const ordinal_type refEdgeOrts[24][6] = { { 0, 0, 1, 0, 0, 0 },
                                                    { 1, 0, 1, 0, 0, 0 },
                                                    { 0, 1, 1, 0, 0, 0 },
                                                    { 0, 1, 0, 0, 0, 0 },
                                                    
                                                    { 1, 0, 0, 0, 0, 0 },
                                                    { 1, 1, 0, 0, 0, 0 },
                                                    { 0, 0, 1, 0, 0, 1 },
                                                    { 1, 0, 1, 0, 0, 1 },
                                                    
                                                    { 0, 0, 1, 0, 1, 1 },
                                                    { 0, 0, 1, 1, 1, 1 },
                                                    { 1, 0, 1, 1, 0, 1 },
                                                    { 1, 0, 1, 1, 1, 1 },
                                                    
                                                    { 0, 1, 1, 0, 1, 0 },
                                                    { 0, 1, 0, 0, 1, 0 },
                                                    { 0, 1, 1, 0, 1, 1 },
                                                    { 0, 1, 1, 1, 1, 1 },
                                                    
                                                    { 0, 1, 0, 1, 1, 0 },
                                                    { 0, 1, 0, 1, 1, 1 },
                                                    { 1, 0, 0, 1, 0, 0 },
                                                    { 1, 1, 0, 1, 0, 0 },
                                                    
                                                    { 1, 0, 0, 1, 0, 1 },
                                                    { 1, 0, 0, 1, 1, 1 },
                                                    { 1, 1, 0, 1, 1, 0 },
                                                    { 1, 1, 0, 1, 1, 1 } };
          
          const ordinal_type refFaceOrts[24][4] = { { 0, 0, 3, 3 },
                                                    { 4, 0, 3, 2 },
                                                    { 0, 4, 3, 0 },
                                                    { 0, 4, 2, 4 },
                                                                 
                                                    { 4, 0, 2, 5 },
                                                    { 4, 4, 2, 1 },
                                                    { 0, 3, 0, 3 },
                                                    { 4, 3, 0, 2 },
                                                                 
                                                    { 3, 2, 0, 3 },
                                                    { 2, 2, 4, 3 },
                                                    { 1, 3, 4, 2 },
                                                    { 5, 2, 4, 2 },
                                                                 
                                                    { 3, 1, 3, 0 },
                                                    { 3, 1, 2, 4 },
                                                    { 3, 5, 0, 0 },
                                                    { 2, 5, 4, 0 },
                                                                 
                                                    { 2, 1, 5, 4 },
                                                    { 2, 5, 1, 4 },
                                                    { 1, 0, 5, 5 },
                                                    { 1, 4, 5, 1 },
                                                                 
                                                    { 1, 3, 1, 5 },
                                                    { 5, 2, 1, 5 },
                                                    { 5, 1, 5, 1 },
                                                    { 5, 5, 1, 1 } };
          
          *outStream << "Tetrahedral element edge reference configuration\n";
          for (auto edgeId=0;edgeId<6;++edgeId) {
            const auto v0 = cellTopo.getNodeMap(1, edgeId, 0);
            const auto v1 = cellTopo.getNodeMap(1, edgeId, 1);
            *outStream << std::setw(3) << edgeId << " edge :: " << v0 << " " << v1 << "\n";
          }
          *outStream << "Tetrahedral element edge orientation with vertex permutations\n";          
          for (auto i=0;i<24;++i) {                                          
            // find orientation                                              
            const auto nodes = Kokkos::View<const ordinal_type[4],Kokkos::HostSpace>(elemNodes[i]);
            const auto ort = Orientation::getOrientation(cellTopo, nodes);
            
            // decode orientation
            ordinal_type edgeOrt[6] = {};
            for (auto edgeId=0;edgeId<6;++edgeId) 
              ort.getEdgeOrientation(edgeOrt, 6);
            
            *outStream << " elemNodes = " 
                       << elemNodes[i][0] << " " 
                       << elemNodes[i][1] << " " 
                       << elemNodes[i][2] << " " 
                       << elemNodes[i][3] << " "
                       << " computed edgeOrts = " 
                       << edgeOrt[0] << " "
                       << edgeOrt[1] << " "
                       << edgeOrt[2] << " "
                       << edgeOrt[3] << " "
                       << edgeOrt[4] << " "
                       << edgeOrt[5] << " "
                       << " reference edgeOrts = " 
                       << refEdgeOrts[i][0] << " "
                       << refEdgeOrts[i][1] << " "
                       << refEdgeOrts[i][2] << " "
                       << refEdgeOrts[i][3] << " "
                       << refEdgeOrts[i][4] << " "
                       << refEdgeOrts[i][5] << " ::\n";
            
            if (edgeOrt[0] != refEdgeOrts[i][0] || 
                edgeOrt[1] != refEdgeOrts[i][1] || 
                edgeOrt[2] != refEdgeOrts[i][2] ||
                edgeOrt[3] != refEdgeOrts[i][3] ||
                edgeOrt[4] != refEdgeOrts[i][4] || 
                edgeOrt[5] != refEdgeOrts[i][5]) {
              *outStream << "        ^^^^^^^^^^^^^^^^ FAILURE\n";
              ++errorFlag;
            }
          }

          *outStream << "Tetrahedral element face reference configuration\n";
          for (auto faceId=0;faceId<4;++faceId) {
            const auto v0 = cellTopo.getNodeMap(2, faceId, 0);
            const auto v1 = cellTopo.getNodeMap(2, faceId, 1);
            const auto v2 = cellTopo.getNodeMap(2, faceId, 2);
            *outStream << std::setw(3) << faceId << " face :: " << v0 << " " << v1 << " " << v2 << "\n";
          }
          *outStream << "Tetrahedral element face orientation with vertex permutations\n";          
          for (auto i=0;i<24;++i) {                                          
            // find orientation                                              
            const auto nodes = Kokkos::View<const ordinal_type[4],Kokkos::HostSpace>(elemNodes[i]);
            const auto ort = Orientation::getOrientation(cellTopo, nodes);
            
            // decode orientation
            ordinal_type faceOrt[4] = {};
            for (auto faceId=0;faceId<4;++faceId) 
              ort.getFaceOrientation(faceOrt, 4);
            
            *outStream << " elemNodes = " 
                       << elemNodes[i][0] << " " 
                       << elemNodes[i][1] << " " 
                       << elemNodes[i][2] << " " 
                       << elemNodes[i][3] << " "
                       << " computed faceOrts = " 
                       << faceOrt[0] << " "
                       << faceOrt[1] << " "
                       << faceOrt[2] << " "
                       << faceOrt[3] << " "
                       << " reference faceOrts = " 
                       << refFaceOrts[i][0] << " "
                       << refFaceOrts[i][1] << " "
                       << refFaceOrts[i][2] << " "
                       << refFaceOrts[i][3] << " ::\n";
            
            if (faceOrt[0] != refFaceOrts[i][0] || 
                faceOrt[1] != refFaceOrts[i][1] || 
                faceOrt[2] != refFaceOrts[i][2] ||
                faceOrt[3] != refFaceOrts[i][3]) {
              *outStream << "        ^^^^^^^^^^^^^^^^ FAILURE\n";
              ++errorFlag;
            }
          }
        }

      } catch (std::exception &err) {
        std::cout << " Exeption\n";
        *outStream << err.what() << "\n\n";
        errorFlag = -1000;
      }

      if (errorFlag != 0)
        std::cout << "End Result: TEST FAILED = " << errorFlag << "\n";
      else
        std::cout << "End Result: TEST PASSED\n";
      
      // reset format state of std::cout
      std::cout.copyfmt(oldFormatState);
      
      return errorFlag;
    }
  }
}
    
