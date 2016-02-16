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
// Questions? Contact Pavel Bochev  (pbboche@sandia.gov)
//                    Denis Ridzal  (dridzal@sandia.gov), or
//                    Kara Peterson (kjpeter@sandia.gov)
//
// ************************************************************************
// @HEADER


/** \file
    \brief  Unit test of experimental high order assembly
    \author Created by Kyungjoo Kim
*/

#include "Intrepid2_FieldContainer.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Intrepid2_PointTools.hpp"
#if defined( INTREPID_USING_EXPERIMENTAL_HIGH_ORDER )
#include "Intrepid2_OrientationTools.hpp"
#endif
#include "Intrepid2_HGRAD_TET_Cn_FEM.hpp"
#include "Shards_CellTopology.hpp"

#include <iostream>
using namespace Intrepid2;

/** \brief Tests for experimental assembly procedure matching basis values.
    \param argc [in] - number of command-line arguments
    \param argv [in] - command-line arguments
*/
int main(int argc, char *argv[]) {
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  Kokkos::initialize();

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint = argc - 1;

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
    << "|                           Unit Test HGRAD_TET_Cn_FEM                        |\n" \
    << "|                                                                             |\n" \
    << "|     1) High order assembly                                                  |\n" \
    << "|                                                                             |\n" \
    << "|  Questions? Contact  Pavel Bochev (pbboche@sandia.gov) or                   |\n" \
    << "|                      Denis Ridzal (dridzal@sandia.gov) or                   |\n" \
    << "|                      Robert Kirby (robert.c.kirby@ttu.edu) or               |\n" \
    << "|                      Kyungjoo Kim (kyukim@sandia.gov)                       |\n" \
    << "|                                                                             |\n" \
    << "|  Intrepid's website: http://trilinos.sandia.gov/packages/intrepid           |\n" \
    << "|  Trilinos website:   http://trilinos.sandia.gov                             |\n" \
    << "|                                                                             |\n" \
    << "===============================================================================\n";


  int r_val = 0;

#if defined( INTREPID_USING_EXPERIMENTAL_HIGH_ORDER )
  typedef double value_type;

  // Let's instantiate a basis
  try {
    OrientationTools<value_type>::verboseStreamPtr = outStream.get();
    for (int test_order=1;test_order<=10;++test_order) {
      for (int test_face=0;test_face<4;++test_face) {
        for (int test_ort=0;test_ort<6;++test_ort) {
          *outStream << "\n"                                            \
                     << "===============================================================================\n" \
                     << "  Order = " << test_order << " , Face = " << test_face << " , Orientation = " << test_ort << "\n" \
                     << "===============================================================================\n";
          
          // populate points on a triangle
          shards::CellTopology face_topo(shards::getCellTopologyData<shards::Triangle<3> >() );
          
          const int order = test_order, offset = 1;
          const int npts_face = PointTools::getLatticeSize(face_topo, order, offset);

          // reference triangle points
          FieldContainer<value_type> ref_face_pts(npts_face, 2);
          PointTools::getLattice<value_type>(ref_face_pts,
                                             face_topo,
                                             order, offset);
          
          // modified triangle points
          FieldContainer<value_type> ort_face_pts(npts_face, 2);
          OrientationTools<value_type>::mapToModifiedReference(ort_face_pts,
                                                               ref_face_pts,
                                                               face_topo,
                                                               test_ort);

          // create a basis object
          Basis_HGRAD_TET_Cn_FEM<double,FieldContainer<double> >  basis(order , POINTTYPE_WARPBLEND);
          const int nbf = basis.getCardinality();
          
          const shards::CellTopology cell_topo = basis.getBaseCellTopology();
          const int ndim  = cell_topo.getDimension();
          const int nface = cell_topo.getFaceCount();
          
          // create orientation object
          int orts[4] = {};
          orts[test_face] = test_ort;
          
          Orientation ort;
          ort.setFaceOrientation(nface, orts);

          // map triangle points and modified points to reference coordinates
          FieldContainer<value_type> tmp_cell_pts(npts_face, ndim);
          
          const int npts_cell = npts_face*nface;
          FieldContainer<value_type> ref_cell_pts(npts_cell, ndim), ort_cell_pts(npts_cell, ndim);

          for (int face_id=0,offset=0;face_id<nface;++face_id,offset+=npts_face) {
            CellTools<value_type>::mapToReferenceSubcell(tmp_cell_pts,
                                                         ref_face_pts,
                                                         2,
                                                         face_id,
                                                         cell_topo);

            for (int i=0;i<npts_face;++i)
              for (int j=0;j<ndim;++j)
                ref_cell_pts(i + offset, j) = tmp_cell_pts(i, j);
            
            CellTools<value_type>::mapToReferenceSubcell(tmp_cell_pts,
                                                         (orts[face_id] ? ort_face_pts : ref_face_pts),
                                                         2,
                                                         face_id,
                                                         cell_topo);
            
            for (int i=0;i<npts_face;++i)
              for (int j=0;j<ndim;++j)
                ort_cell_pts(i + offset, j) = tmp_cell_pts(i, j);
          }
          
          // temporary cell workspace
          FieldContainer<double> tmp_cell_vals(nbf, npts_cell);
          
          // modified basis values wrt reference coordinates
          FieldContainer<double> ref_cell_vals(nbf, npts_cell);
          basis.getValues(tmp_cell_vals , ref_cell_pts, OPERATOR_VALUE );
          OrientationTools<value_type>::getBasisFunctionsByTopology(ref_cell_vals,
                                                                    tmp_cell_vals,
                                                                    basis);

          // basis values wrt modified coordinates
          FieldContainer<double> ort_cell_vals(nbf, npts_cell);
          basis.getValues(tmp_cell_vals, ort_cell_pts, OPERATOR_VALUE);
          OrientationTools<value_type>::getBasisFunctionsByTopology(ort_cell_vals,
                                                                    tmp_cell_vals,
                                                                    basis);
          
          for (int i=0;i<nbf;++i)
            for (int j=0;j<npts_cell;++j)
              tmp_cell_vals(i, j) = ort_cell_vals(i, j);
          
          // modify the basis wrt orientations
          OrientationTools<value_type>::verbose = true;  
          OrientationTools<value_type>::getModifiedBasisFunctions(ort_cell_vals,
                                                                 tmp_cell_vals,
                                                                 basis,
                                                                 ort);
          OrientationTools<value_type>::verbose = false; 

          // check the basis should be same for edge DOFs
          {
            int ibegin = 0; 
            const int nvert = cell_topo.getVertexCount();
            for (int i=0;i<nvert;++i) {
              const int ord_vert = basis.getDofOrdinal(0, i, 0);
              ibegin += basis.getDofTag(ord_vert)[3];
            }
            if (ibegin < nbf) {
              const int nedge = cell_topo.getEdgeCount();
              for (int i=0;i<nedge;++i) {
                const int ord_edge = basis.getDofOrdinal(1, i, 0);
                ibegin += basis.getDofTag(ord_edge)[3];
              }
            }
            int iend = ibegin; 
            if (iend < nbf) {
              const int nface = cell_topo.getFaceCount();
              for (int i=0;i<nface;++i) {
                const int ord_face = basis.getDofOrdinal(2, i, 0);
                iend += basis.getDofTag(ord_face)[3];
              }
            }
            for (int i=ibegin;i<iend;++i) {
              for (int j=0;j<npts_cell;++j) {
                const value_type diff = std::abs(ref_cell_vals(i,j) - ort_cell_vals(i,j));
                if (diff > INTREPID_TOL) {
                  ++r_val;
                  *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
                  *outStream << " Basis function " << i << " at point " << j << " does not match each other\n";
                }
              }
            }
          }
        }
      }
    }
  }
  catch (std::exception err) {
    std::cout << " Exeption\n";
    *outStream << err.what() << "\n\n";
    r_val = -1000;
  }
#else
  *outStream << "\t This test is for high order element assembly. \n"
             << "\t Use -D INTREPID_USING_EXPERIMENTAL_HIGH_ORDER in CMAKE_CXX_FLAGS \n";
#endif
  if (r_val != 0)
    std::cout << "End Result: TEST FAILED\n";
  else
    std::cout << "End Result: TEST PASSED\n";

  // reset format state of std::cout
  std::cout.copyfmt(oldFormatState);

  Kokkos::finalize();

  return r_val;
}
