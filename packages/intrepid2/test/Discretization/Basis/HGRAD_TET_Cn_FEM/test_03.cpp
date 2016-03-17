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
    \brief  Unit test of experimental high order assembly
    \author Created by Kyungjoo Kim
*/

#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include "Shards_CellTopology.hpp"

#include "Intrepid2_FieldContainer.hpp"
#include "Intrepid2_PointTools.hpp"

#if defined( INTREPID_USING_EXPERIMENTAL_HIGH_ORDER )
#include "Intrepid2_BasisSet.hpp"
#include "Intrepid2_OrientationTools.hpp"
#endif

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

  bool verbose = false;
  int  maxp = INTREPID2_MAX_ORDER;
  for (int i=0;i<argc;++i) {
    if ((strcmp(argv[i],"--verbose")           == 0)) { verbose  = atoi(argv[++i]); continue;}
    if ((strcmp(argv[i],"--maxp")              == 0)) { maxp     = atoi(argv[++i]); continue;}
  }
  
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
    //for (int test_order=1;test_order<=10;++test_order) {
    for (int test_order=1;test_order<=maxp;++test_order) {
      // Step 0 : construct basis function set
      const int order = test_order;
      
      BasisSet_HGRAD_TET_Cn_FEM<double,FieldContainer<double> > basis_set(order , POINTTYPE_EQUISPACED);
      const auto& cell_basis = basis_set.getCellBasis();
      const auto& face_basis = basis_set.getTriangleBasis();
      
      const shards::CellTopology cell_topo = cell_basis.getBaseCellTopology();
      const shards::CellTopology face_topo = face_basis.getBaseCellTopology();
      
      const int nbf_cell = cell_basis.getCardinality();
      const int nbf_face = face_basis.getCardinality();
      
      const int ndim_cell  = cell_topo.getDimension();
      const int ndim_face  = face_topo.getDimension();
      
      const int npts = PointTools::getLatticeSize(face_topo, order, 1);
      
      for (int test_face=0;test_face<4;++test_face) {
        // tricky part
        const bool left_handed = cell_topo.getNodeMap(2, test_face, 1) > cell_topo.getNodeMap(2, test_face, 2);

        for (int test_ort=0;test_ort<6;++test_ort) {
          *outStream << "\n"                                            \
                     << "===============================================================================\n" \
                     << "  Order = " << test_order << " , Face = " << test_face << " , Orientation = " << test_ort << "\n" \
                     << "===============================================================================\n";
          
          // Step 1 : create reference and modified triangle points

          // reference triangle points
          FieldContainer<value_type> ref_face_pts(npts, ndim_face);
          PointTools::getLattice<value_type>(ref_face_pts,
                                             face_topo,
                                             order, 1);

          // modified triangle points
          const int left_ort[] = { 0, 2, 1, 3, 5, 4 };
          FieldContainer<value_type> ort_face_pts(npts, ndim_face);
          OrientationTools<value_type>::mapToModifiedReference(ort_face_pts,
                                                               ref_face_pts,
                                                               face_topo,
                                                               (left_handed ? left_ort[test_ort] : test_ort));


          // Step 2 : map face points to cell points appropriately
          const int nface = cell_topo.getFaceCount();

          // create orientation object
          int orts[4] = {};
          orts[test_face] = test_ort;
          
          Orientation ort;
          ort.setFaceOrientation(nface, orts);

          // map triangle points and modified points to reference coordinates
          FieldContainer<value_type> ref_cell_pts(npts, ndim_cell);
          CellTools<value_type>::mapToReferenceSubcell(ref_cell_pts,
                                                       ref_face_pts,
                                                       ndim_face,
                                                       test_face,
                                                       cell_topo);

          // Step 3 : evaluate modified basis functions with orientation for reference cell points
          FieldContainer<double> ort_cell_vals(nbf_cell, npts);
          {
            // temporary cell workspace
            FieldContainer<double> tmp_cell_vals(nbf_cell, npts);

            cell_basis.getValues(tmp_cell_vals, ref_cell_pts, OPERATOR_VALUE);
            OrientationTools<value_type>::getBasisFunctionsByTopology(ort_cell_vals,
                                                                      tmp_cell_vals,
                                                                      cell_basis);

            for (int i=0;i<nbf_cell;++i)
              for (int j=0;j<npts;++j)
                tmp_cell_vals(i, j) = ort_cell_vals(i, j);

            OrientationTools<value_type>::verbose = verbose;
            OrientationTools<value_type>::reverse = true; // for point matching only
            OrientationTools<value_type>::getModifiedBasisFunctions(ort_cell_vals,
                                                                    tmp_cell_vals,
                                                                    basis_set,
                                                                    ort);
            OrientationTools<value_type>::verbose = false;
          }

          // Step 4 : evaluate reference face basis functions for modified face points
          FieldContainer<double> ref_face_vals(nbf_face, npts);
          {
            // temporary face workspace
            FieldContainer<double> tmp_face_vals(nbf_face, npts);

            face_basis.getValues(tmp_face_vals, ort_face_pts, OPERATOR_VALUE);
            OrientationTools<value_type>::getBasisFunctionsByTopology(ref_face_vals,
                                                                      tmp_face_vals,
                                                                      face_basis);
          }

          // Step 5 : compare the basis functions to face functions
          {
            // strip the range of cell DOFs
            int off_cell = 0;
            {
              const int nvert = cell_topo.getVertexCount();
              for (int i=0;i<nvert;++i) {
                const int ord_vert = cell_basis.getDofOrdinal(0, i, 0);
                off_cell += cell_basis.getDofTag(ord_vert)[3];
              }
              if (off_cell < nbf_cell) {
                const int nedge = cell_topo.getEdgeCount();
                for (int i=0;i<nedge;++i) {
                  const int ord_edge = cell_basis.getDofOrdinal(1, i, 0);
                  off_cell += cell_basis.getDofTag(ord_edge)[3];
                }
              }
              if (off_cell < nbf_cell) {
                for (int i=0;i<test_face;++i) {
                  const int ord_face = cell_basis.getDofOrdinal(2, i, 0);
                  off_cell += cell_basis.getDofTag(ord_face)[3];
                }
              }
            }

            // strip the range of face DOFs
            int off_face = 0;
            {
              const int nvert = face_topo.getVertexCount();
              for (int i=0;i<nvert;++i) {
                const int ord_vert = face_basis.getDofOrdinal(0, i, 0);
                off_face += face_basis.getDofTag(ord_vert)[3];
              }
              if (off_face < nbf_face) {
                const int nedge = face_topo.getEdgeCount();
                for (int i=0;i<nedge;++i) {
                  const int ord_edge = face_basis.getDofOrdinal(1, i, 0);
                  off_face += face_basis.getDofTag(ord_edge)[3];
                }
              }
            }

            const int ndof = nbf_face - off_face;
            for (int i=0;i<ndof;++i) {
              for (int j=0;j<npts;++j) {
                const value_type diff = std::abs(ort_cell_vals(i+off_cell,j) - ref_face_vals(i+off_face,j));
                if (diff > INTREPID_TOL) {
                  ++r_val;
                  *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
                  *outStream << " Basis function " << i << " at point " << j << " does not match each other\n";
                }
              }
            }
          }
        } // test ort
      } // test face
    } // test order
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
    std::cout << "End Result: TEST FAILED, r_val = " << r_val << "\n";
  else
    std::cout << "End Result: TEST PASSED\n";

  // reset format state of std::cout
  std::cout.copyfmt(oldFormatState);

  Kokkos::finalize();

  return r_val;
}
