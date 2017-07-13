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
#include "Intrepid2_HGRAD_TRI_Cn_FEM.hpp"

#if defined( INTREPID_USING_EXPERIMENTAL_HIGH_ORDER )
#include "Intrepid2_BasisSet.hpp"
#include "Intrepid2_OrientationTools.hpp"
#endif

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

#if defined( INTREPID_USING_EXPERIMENTAL_HIGH_ORDER )
  bool verbose = false;
  int  maxp = INTREPID2_MAX_ORDER;
  for (int i=0;i<argc;++i) {
    if ((strcmp(argv[i],"--verbose")           == 0)) { verbose  = atoi(argv[++i]); continue;}
    if ((strcmp(argv[i],"--maxp")              == 0)) { maxp     = atoi(argv[++i]); continue;}
  }
#endif
  
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
    << "|                           Unit Test HGRAD_TRI_Cn_FEM                        |\n" \
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
    for (int test_order=1;test_order<=maxp;++test_order) {
      // Step 0 : construct basis function set
      const int order = test_order;
      
      BasisSet_HGRAD_TRI_Cn_FEM<double,FieldContainer<double> > basis_set(order , POINTTYPE_EQUISPACED);
      const auto& cell_basis = basis_set.getCellBasis();
      const auto& line_basis = basis_set.getLineBasis();
      
      const shards::CellTopology cell_topo = cell_basis.getBaseCellTopology();
      const shards::CellTopology line_topo = line_basis.getBaseCellTopology();
      
      const int nbf_cell = cell_basis.getCardinality();
      const int nbf_line = line_basis.getCardinality();
      
      const int ndim_cell  = cell_topo.getDimension();
      const int ndim_line  = line_topo.getDimension();
      
      const int npts = PointTools::getLatticeSize(line_topo, order, 1);

      for (int test_edge=0;test_edge<3;++test_edge) {
        for (int test_ort=0;test_ort<2;++test_ort) {
          *outStream << "\n"                                              \
                     << "===============================================================================\n" \
                     << "  Order = " << test_order << " , Edge = " << test_edge << " , Orientation = " << test_ort << "\n" \
                     << "===============================================================================\n";

          // Step 1 : create reference and modified line points

          // reference line points
          FieldContainer<value_type> ref_line_pts(npts, ndim_line);
          PointTools::getLattice<value_type>(ref_line_pts,
                                             line_topo,
                                             order, 1);

          // modified line points
          FieldContainer<value_type> ort_line_pts(npts, ndim_line);
          OrientationTools<value_type>::mapToModifiedReference(ort_line_pts,
                                                               ref_line_pts,
                                                               line_topo,
                                                               test_ort);

          // Step 2 : map line points to cell points appropriately
          const int nedge = cell_topo.getEdgeCount();

          // create orientation object
          int orts[3] = {};
          orts[test_edge] = test_ort;

          Orientation ort;
          ort.setEdgeOrientation(nedge, orts);

          // map line points and modified points to reference coordinates
          FieldContainer<value_type> ref_cell_pts(npts, ndim_cell);
          CellTools<value_type>::mapToReferenceSubcell(ref_cell_pts,
                                                       ref_line_pts,
                                                       ndim_line,
                                                       test_edge,
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
            OrientationTools<value_type>::reverse = true;
            OrientationTools<value_type>::getModifiedBasisFunctions(ort_cell_vals,
                                                                    tmp_cell_vals,
                                                                    basis_set,
                                                                    ort);
            OrientationTools<value_type>::verbose = false;
          }

          // Step 4 : evaluate reference line basis functions for modified line points
          FieldContainer<double> ref_line_vals(nbf_line, npts);
          {
            // temporary line workspace
            FieldContainer<double> tmp_line_vals(nbf_line, npts);

            line_basis.getValues(tmp_line_vals, ort_line_pts, OPERATOR_VALUE);
            OrientationTools<value_type>::getBasisFunctionsByTopology(ref_line_vals,
                                                                      tmp_line_vals,
                                                                      line_basis);
          }

          // Step 5 : compare the basis functions to line functions
          {
            // strip the range of edge DOFs
            int off_cell = 0;
            {
              const int nvert = cell_topo.getVertexCount();
              for (int i=0;i<nvert;++i) {
                const int ord_vert = cell_basis.getDofOrdinal(0, i, 0);
                off_cell += cell_basis.getDofTag(ord_vert)[3];
              }
              if (off_cell < nbf_cell) {
                for (int i=0;i<test_edge;++i) {
                  const int ord_edge = cell_basis.getDofOrdinal(1, i, 0);
                  off_cell += cell_basis.getDofTag(ord_edge)[3];
                }
              }
            }

            // strip the range of line DOFs
            int off_line = 0;
            {
              const int nvert = line_topo.getVertexCount();
              for (int i=0;i<nvert;++i) {
                const int ord_vert = line_basis.getDofOrdinal(0, i, 0);
                off_line += line_basis.getDofTag(ord_vert)[3];
              }
            }

            const int ndof = nbf_line - off_line;
            for (int i=0;i<ndof;++i) {
              for (int j=0;j<npts;++j) {
                const value_type diff = std::abs(ort_cell_vals(i+off_cell,j) - ref_line_vals(i+off_line,j));
                if (diff > INTREPID_TOL) {
                  ++r_val;
                  *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
                  *outStream << " Basis function " << i << " at point " << j << " does not match each other\n";
                }
              }
            }
          }
        } // test ort
      } // test edge
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
