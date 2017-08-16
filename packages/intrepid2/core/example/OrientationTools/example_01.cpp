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
    \brief  Example of the OrientationTools class (Orientation).
    \author Created by Kyungjoo Kim
*/

#include "Intrepid2_PointTools.hpp"

#include "Intrepid2_OrientationTools.hpp"
#include "Intrepid2_FieldContainer.hpp"
#include "Intrepid2_FieldContainer_Kokkos.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include "Shards_CellTopology.hpp"

#include "Teuchos_RCP.hpp"

using namespace std;
using namespace Intrepid2;
using namespace shards;

int main(int argc, char *argv[]) {

#if defined( INTREPID_USING_EXPERIMENTAL_HIGH_ORDER )
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  
  Kokkos::initialize();

  typedef double value_type;
  typedef OrientationTools<double> OrientationTools;
  typedef shards::CellTopology CellTopology;
  
  std::cout                                                             \
    << "===============================================================================\n" \
    << "|                                                                             |\n" \
    << "|                   Example use of the OrientationTools class                 |\n" \
    << "|                                                                             |\n" \
    << "|  1) Reference edge parametrizations with orientations                       |\n" \
    << "|  2) Reference face parametrizations with orientations                       |\n" \
    << "|                                                                             |\n" \
    << "|  Questions? Contact  Pavel Bochev (pbboche@sandia.gov)                      |\n" \
    << "|                      Denis Ridzal (dridzal@sandia.gov)                      |\n" \
    << "|                      Kara Peterson (kjpeter@sandia.gov), or                 |\n" \
    << "|                      Kyungjoo Kim (kyukim@sandia.gov).                      |\n" \
    << "|                                                                             |\n" \
    << "|  Intrepid's website: http://trilinos.sandia.gov/packages/intrepid           |\n" \
    << "|  Trilinos website:   http://trilinos.sandia.gov                             |\n" \
    << "|                                                                             |\n" \
    << "===============================================================================\n" \
    << "| Summary:                                                                    |\n" \
    << "| Reference parametrizations are modified with their orientations which are   |\n" \
    << "| derived from mesh connectivities. The modified paramterizations are used to |\n" \
    << "| compute high order element dof connectivities.                              |\n" \
    << "===============================================================================\n";
  /*
   * The modified reference points represent the matching conditions accounting for 
   * the orientations. 
   *
   * 1. Populate collocation points according to the order of approximation using PointTools.
   * 
   * 2. Rotate or flip the coordinates of the input point according to the orientation id.
   * 
   *    Vertex enumeration for different orientations
   *    ---------------------------------------------
   *
   *    Reference Line      0-1      o0 = (0,1), o1 = (1,0)
   *
   *    Reference Triangle  2        o0 = (0,1,2), o1 = (1,2,0), o2 = (2,0,1)
   *                        |\       o3 = (0,2,1), o4 = (1,0,2), o5 = (2,1,0)
   *                        0-1 
   *
   *    Reference Quadrilateral 3-2  o0 = (0,1,2,3), o1 = (1,2,3,0), o2 = (2,3,0,1), o3 = (3,0,1,2)
   *                            | |  o4 = (0,3,2,1), o5 = (1,0,3,2), o6 = (2,1,0,3), o7 = (3,2,1,0)
   *                            0-1 
   */

  std::cout << "\n"                                                     \
            << "===============================================================================\n" \
            << "| Line segment parametrizations with orientation                              |\n" \
            << "===============================================================================\n";
  {
    CellTopology cellTopo(shards::getCellTopologyData<shards::Line<2> >() );

    // consider the interior of a line segment
    const int order = 6, offset = 1;
    const int numPts = PointTools::getLatticeSize(cellTopo, order, offset);

    FieldContainer<value_type> refPts(numPts, 1);
    PointTools::getLattice<value_type>(refPts, 
                                       cellTopo, 
                                       order, offset);

    FieldContainer<value_type> ortPts(numPts, 1);
    for (int cellOrt=0;cellOrt<2;++cellOrt) {
      OrientationTools::mapToModifiedReference(ortPts,
                                               refPts,
                                               cellTopo,
                                               cellOrt);
      
      std::stringstream ss;
      for(int pt=0;pt<numPts;++pt) 
        ss << "\t Modified parameter point " 
           << std::setw(12) << std::right << refPts(pt, 0) << std::setw(10) << "  -->  "
           << std::setw(10) << std::right << ortPts(pt, 0) << "\n";
      
      std::cout << "Line segment orientation = " << cellOrt << "\n " 
                << ss.str() << "\n";
      
    }
  }

  std::cout << "\n"                                                  \
            << "===============================================================================\n" \
            << "| Triangle parametrizations with orientation                                  |\n" \
            << "===============================================================================\n";
  {
    CellTopology cellTopo(shards::getCellTopologyData<shards::Triangle<3> >() );

    // consider the interior of a line segment
    const int order = 6, offset = 1;
    const int dimPts = 2, numPts = PointTools::getLatticeSize(cellTopo, order, offset);

    FieldContainer<value_type> refPts(numPts, dimPts);
    PointTools::getLattice<value_type>(refPts, 
                                       cellTopo, 
                                       order, offset);

    FieldContainer<value_type> ortPts(numPts, dimPts);
    for (int cellOrt=0;cellOrt<6;++cellOrt) {
      OrientationTools::mapToModifiedReference(ortPts,
                                               refPts,
                                               cellTopo,
                                               cellOrt);
      
      std::stringstream ss;
      for(int pt=0;pt<numPts;++pt) 
        ss << "\t Modified parameter point " << "  ("
           << std::setw(10) << std::right << refPts(pt, 0) << " , "
           << std::setw(10) << std::right << refPts(pt, 1) << ")"
           << std::setw(10) << "  -->  " << "("           
           << std::setw(10) << std::right << ortPts(pt, 0) << " , "                                         
           << std::setw(10) << std::right << ortPts(pt, 1) << ")\n";                       

      std::cout << "Triangle orientation = " << cellOrt << "\n" 
                << ss.str() << "\n";
    }
  }

  std::cout << "\n"                                                     \
            << "===============================================================================\n" \
            << "| Quadrilateral parametrizations with orientation                             |\n" \
            << "===============================================================================\n";
  {
    CellTopology lineTopo(shards::getCellTopologyData<shards::Line<2> >() );

    // tensor product of line segments
    const int order = 6, offset = 1, numPtsLine = PointTools::getLatticeSize(lineTopo, order, offset);
    const int dimPts = 2, numPts = numPtsLine*numPtsLine;
    
    FieldContainer<value_type> refPtsLine(numPtsLine, 1);
    PointTools::getLattice<value_type>(refPtsLine, 
                                       lineTopo, 
                                       order, offset);
    
    // variable order is not supported in point/cell tools and other facilities
    FieldContainer<value_type> refPts(numPts, dimPts);
    for (int k1=0,idx=0;k1<numPtsLine;++k1) {
      for (int k0=0;k0<numPtsLine;++k0,++idx) {
        refPts(idx, 0) = refPtsLine(k0);
        refPts(idx, 1) = refPtsLine(k1);
      }
    }

    CellTopology cellTopo(shards::getCellTopologyData<shards::Quadrilateral<4> >() );    
    FieldContainer<value_type> ortPts(numPts, dimPts);
    for (int cellOrt=0;cellOrt<8;++cellOrt) {
      OrientationTools::mapToModifiedReference(ortPts,
                                               refPts,
                                               cellTopo,
                                               cellOrt);

      std::stringstream ss;
      for(int pt=0;pt<numPts;++pt) 
        ss << "\t Modified parameter point " << "  ("
           << std::setw(10) << std::right << refPts(pt, 0) << " , "
           << std::setw(10) << std::right << refPts(pt, 1) << ")"
           << std::setw(10) << "  -->  " << "("           
           << std::setw(10) << std::right << ortPts(pt, 0) << " , "                                         
           << std::setw(10) << std::right << ortPts(pt, 1) << ")\n";                       
      
      std::cout << "Quadrilateral orientation = " << cellOrt << "\n" 
                << ss.str() << "\n";

    }
  }

  Kokkos::finalize();
#else
  std::cout << "\t This example is for high order element assembly. \n"
            << "\t Use -D INTREPID_USING_EXPERIMENTAL_HIGH_ORDER in CMAKE_CXX_FLAGS \n";
#endif
  
  return 0;
}

