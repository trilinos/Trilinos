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
    << "|  1) Triangle subcell vertex enumeration                                     |\n" \
    << "|  2) Tetrahedron subcell vertex enumeration                                  |\n" \
    << "|                                                                             |\n" \
    << "|  Questions? Contact  Pavel Bochev (pbboche@sandia.gov)                      |\n" \
    << "|                      Denis Ridzal (dridzal@sandia.gov)                      |\n" \
    << "|                      Kara Peterson (kjpeter@sandia.gov), or                 |\n" \
    << "|                      Kyungjoo Kim (kyukim@sandia.gov).                      |\n" \
    << "|                                                                             |\n" \
    << "|  Intrepid's website: http://trilinos.sandia.gov/packages/intrepid           |\n" \
    << "|  Trilinos website:   http://trilinos.sandia.gov                             |\n" \
    << "|                                                                             |\n" \
    << "===============================================================================\n";

  {
    CellTopology cellTopo(shards::getCellTopologyData<shards::Triangle<3> >() );
    std::cout << cellTopo << std::endl;
  }
  {
    CellTopology cellTopo(shards::getCellTopologyData<shards::Quadrilateral<4> >() );
    std::cout << cellTopo << std::endl;
  }
  {
    CellTopology cellTopo(shards::getCellTopologyData<shards::Tetrahedron<4> >() );
    std::cout << cellTopo << std::endl;
  }
  {
    CellTopology cellTopo(shards::getCellTopologyData<shards::Hexahedron<8> >() );
    std::cout << cellTopo << std::endl;
  }
  {
    CellTopology cellTopo(shards::getCellTopologyData<shards::Pyramid<5> >() );
    std::cout << cellTopo << std::endl;
  }
  {
    CellTopology cellTopo(shards::getCellTopologyData<shards::Wedge<6> >() );
    std::cout << cellTopo << std::endl;
  }

  Kokkos::finalize();
#else
  std::cout << "\t This example is for high order element assembly. \n"
            << "\t Use -D INTREPID_USING_EXPERIMENTAL_HIGH_ORDER in CMAKE_CXX_FLAGS \n";
#endif
  
  return 0;
}

