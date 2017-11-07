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

/** \file   test_01.hpp
    \brief  Unit tests for the Intrepid2::HDIV_WEDGE_I1_FEM class.
    \author Created by P. Bochev, D. Ridzal, and K. Peterson.
 */

#include "Intrepid2_config.h"

#ifdef HAVE_INTREPID2_DEBUG
#define INTREPID2_TEST_FOR_DEBUG_ABORT_OVERRIDE_TO_CONTINUE
#endif

#include "Intrepid2_HDIV_WEDGE_I1_FEM.hpp"

#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"

namespace Intrepid2 {

namespace Test {
#define INTREPID2_TEST_ERROR_EXPECTED( S )                              \
    try {                                                               \
      ++nthrow;                                                         \
      S ;                                                               \
    }                                                                   \
    catch (std::exception err) {                                        \
      ++ncatch;                                                         \
      *outStream << "Expected Error ----------------------------------------------------------------\n"; \
      *outStream << err.what() << '\n';                                 \
      *outStream << "-------------------------------------------------------------------------------" << "\n\n"; \
    }

template<typename ValueType, typename DeviceSpaceType>
int HDIV_WEDGE_I1_FEM_Test01(const bool verbose) {

  Teuchos::RCP<std::ostream> outStream;
  Teuchos::oblackholestream bhs; // outputs nothing

  if (1)// (verbose)
    outStream = Teuchos::rcp(&std::cout, false);
  else
    outStream = Teuchos::rcp(&bhs,       false);

  Teuchos::oblackholestream oldFormatState;
  oldFormatState.copyfmt(std::cout);

  typedef typename
      Kokkos::Impl::is_space<DeviceSpaceType>::host_mirror_space::execution_space HostSpaceType ;

  *outStream << "DeviceSpace::  "; DeviceSpaceType::print_configuration(*outStream, false);
  *outStream << "HostSpace::    ";   HostSpaceType::print_configuration(*outStream, false);

  *outStream
  << "===============================================================================\n"
  << "|                                                                             |\n"
  << "|                 Unit Test (Basis_HDIV_WEDGE_I1_FEM)                         |\n"
  << "|                                                                             |\n"
  << "|     1) Conversion of Dof tags into Dof ordinals and back                    |\n"
  << "|     2) Basis values for VALUE and DIV operators                             |\n"
  << "|                                                                             |\n"
  << "|  Questions? Contact  Pavel Bochev  (pbboche@sandia.gov),                    |\n"
  << "|                      Denis Ridzal  (dridzal@sandia.gov),                    |\n"
  << "|                      Kara Peterson (kjpeter@sandia.gov).                    |\n"
  << "|                                                                             |\n"
  << "|  Intrepid's website: http://trilinos.sandia.gov/packages/intrepid           |\n"
  << "|  Trilinos website:   http://trilinos.sandia.gov                             |\n"
  << "|                                                                             |\n"
  << "===============================================================================\n";

  typedef Kokkos::DynRankView<ValueType,DeviceSpaceType> DynRankView;
  typedef Kokkos::DynRankView<ValueType,HostSpaceType> DynRankViewHost;
#define ConstructWithLabel(obj, ...) obj(#obj, __VA_ARGS__)

  const ValueType tol = tolerence();
  int errorFlag = 0;

  // for virtual function, value and point types are declared in the class
  typedef ValueType outputValueType;
  typedef ValueType pointValueType;
  Basis_HDIV_WEDGE_I1_FEM<DeviceSpaceType,outputValueType,pointValueType> wedgeBasis;

  *outStream
  << "\n"
  << "===============================================================================\n"
  << "| TEST 1: constructors and exceptions                                         |\n"
  << "===============================================================================\n";

  try {
#ifdef HAVE_INTREPID2_DEBUG
    ordinal_type nthrow = 0, ncatch = 0;
    // Define array containing the 6 vertices of the reference WEDGE and 6 other points.
    DynRankView ConstructWithLabel(wedgeNodes, 12, 3);

    // Generic array for the output values; needs to be properly resized depending on the operator type
    const auto numFields = wedgeBasis.getCardinality();
    const auto numPoints = wedgeNodes.dimension(0);
    const auto spaceDim  = wedgeBasis.getBaseCellTopology().getDimension();

    // exception #1: GRAD cannot be applied to HDIV functions 
    DynRankView ConstructWithLabel(vals, numFields, numPoints, spaceDim );
    INTREPID2_TEST_ERROR_EXPECTED( wedgeBasis.getValues(vals, wedgeNodes, OPERATOR_GRAD));

    // exception #2: CURL cannot be applied to HDIV functions
    INTREPID2_TEST_ERROR_EXPECTED( wedgeBasis.getValues(vals, wedgeNodes, OPERATOR_CURL));

    // Exceptions 3-7: all bf tags/bf Ids below are wrong and should cause getDofOrdinal() and 
    // getDofTag() to access invalid array elements thereby causing bounds check exception
    // exception #3
    INTREPID2_TEST_ERROR_EXPECTED( wedgeBasis.getDofOrdinal(3,0,0));
    // exception #4
    INTREPID2_TEST_ERROR_EXPECTED( wedgeBasis.getDofOrdinal(1,1,1));
    // exception #5
    INTREPID2_TEST_ERROR_EXPECTED( wedgeBasis.getDofOrdinal(0,4,1));
    // exception #6
    INTREPID2_TEST_ERROR_EXPECTED( wedgeBasis.getDofTag(numFields));
    // exception #7
    INTREPID2_TEST_ERROR_EXPECTED( wedgeBasis.getDofTag(-1));

    // Exceptions 8-15 test exception handling with incorrectly dimensioned input/output arrays
    // exception #8: input points array must be of rank-2
    DynRankView ConstructWithLabel(badPoints1, 4, 5, 3);
    INTREPID2_TEST_ERROR_EXPECTED( wedgeBasis.getValues(vals, badPoints1, OPERATOR_VALUE));

    // exception #9 dimension 1 in the input point array must equal space dimension of the cell
    DynRankView ConstructWithLabel(badPoints2, 4, 2);
    INTREPID2_TEST_ERROR_EXPECTED( wedgeBasis.getValues(vals, badPoints2, OPERATOR_VALUE));

    // exception #10 output values must be of rank-3 for OPERATOR_VALUE
    DynRankView ConstructWithLabel(badVals1, 4, 3);
    INTREPID2_TEST_ERROR_EXPECTED( wedgeBasis.getValues(badVals1, wedgeNodes, OPERATOR_VALUE));

    // exception #11 output values must be of rank-2 for OPERATOR_DIV
    DynRankView ConstructWithLabel(badVals2, 4, 3, 1);
    INTREPID2_TEST_ERROR_EXPECTED( wedgeBasis.getValues(badVals2, wedgeNodes, OPERATOR_DIV));

    // exception #12 incorrect 0th dimension of output array (must equal number of basis functions)
    DynRankView ConstructWithLabel(badVals3, wedgeBasis.getCardinality() + 1, wedgeNodes.dimension(0), 3);
    INTREPID2_TEST_ERROR_EXPECTED( wedgeBasis.getValues(badVals3, wedgeNodes, OPERATOR_VALUE));

    // exception #13 incorrect 0th dimension of output array (must equal number of basis functions)
    DynRankView ConstructWithLabel(badVals4, wedgeBasis.getCardinality() + 1, wedgeNodes.dimension(0));
    INTREPID2_TEST_ERROR_EXPECTED( wedgeBasis.getValues(badVals4, wedgeNodes, OPERATOR_DIV));

    // exception #14 incorrect 1st dimension of output array (must equal number of points)
    DynRankView ConstructWithLabel(badVals5, wedgeBasis.getCardinality(), wedgeNodes.dimension(0) + 1, 3);
    INTREPID2_TEST_ERROR_EXPECTED( wedgeBasis.getValues(badVals5, wedgeNodes, OPERATOR_VALUE));

    // exception #15 incorrect 1st dimension of output array (must equal number of points)
    DynRankView ConstructWithLabel(badVals6, wedgeBasis.getCardinality(), wedgeNodes.dimension(0) + 1);
    INTREPID2_TEST_ERROR_EXPECTED( wedgeBasis.getValues(badVals6, wedgeNodes, OPERATOR_DIV));

    // exception #16: incorrect 2nd dimension of output array (must equal the space dimension)
    DynRankView ConstructWithLabel(badVals7, wedgeBasis.getCardinality(), wedgeNodes.dimension(0), 4);
    INTREPID2_TEST_ERROR_EXPECTED( wedgeBasis.getValues(badVals7, wedgeNodes, OPERATOR_VALUE));
#endif

  } catch (std::logic_error err) {
    *outStream << "UNEXPECTED ERROR !!! ----------------------------------------------------------\n";
    *outStream << err.what() << '\n';
    *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    errorFlag = -1000;
  }


  *outStream
  << "\n"
  << "===============================================================================\n"
  << "| TEST 2: correctness of tag to enum and enum to tag lookups                  |\n"
  << "===============================================================================\n";

  try{
    const auto numFields = wedgeBasis.getCardinality();
    const auto allTags = wedgeBasis.getAllDofTags();

    // Loop over all tags, lookup the associated dof enumeration and then lookup the tag again
    const auto dofTagSize = allTags.dimension(0);

    // Loop over all tags, lookup the associated dof enumeration and then lookup the tag again
    for (size_type i = 0; i < dofTagSize; i++) {
      const auto bfOrd  = wedgeBasis.getDofOrdinal(allTags(i, 0), allTags(i, 1), allTags(i, 2));

      const auto myTag = wedgeBasis.getDofTag(bfOrd);
      if( !( (myTag(0) == allTags(i, 0)) &&
          (myTag(1) == allTags(i, 1)) &&
          (myTag(2) == allTags(i, 2)) &&
          (myTag(3) == allTags(i, 3)) ) ) {
        errorFlag++;
        *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
        *outStream << " getDofOrdinal( {" 
            << allTags(i, 0) << ", "
            << allTags(i, 1) << ", "
            << allTags(i, 2) << ", "
            << allTags(i, 3) << "}) = " << bfOrd <<" but \n";
        *outStream << " getDofTag(" << bfOrd << ") = { "
            << myTag(0) << ", "
            << myTag(1) << ", "
            << myTag(2) << ", "
            << myTag(3) << "}\n";
      }
    }

    // Now do the same but loop over basis functions
    for( int bfOrd = 0; bfOrd <numFields; bfOrd++) {
      const auto myTag  = wedgeBasis.getDofTag(bfOrd);
      const auto myBfOrd = wedgeBasis.getDofOrdinal(myTag(0), myTag(1), myTag(2));
      if( bfOrd != myBfOrd) {
        errorFlag++;
        *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
        *outStream << " getDofTag(" << bfOrd << ") = { "
            << myTag(0) << ", "
            << myTag(1) << ", "
            << myTag(2) << ", "
            << myTag(3) << "} but getDofOrdinal({"
            << myTag(0) << ", "
            << myTag(1) << ", "
            << myTag(2) << ", "
            << myTag(3) << "} ) = " << myBfOrd << "\n";
      }
    }
  }
  catch (std::logic_error err){
    *outStream << err.what() << "\n\n";
    errorFlag = -1000;
  };

  *outStream \
  << "\n"
  << "===============================================================================\n"\
  << "| TEST 3: correctness of basis function values                                |\n"\
  << "===============================================================================\n";

  outStream -> precision(20);

  // VALUE: Each row pair gives the 5x3 correct basis set values at an evaluation point
  double basisValues[] = {
      0, -0.500000, 0, 0, 0, 0, -0.500000, 0, 0, 0, 0, -2.00000, 0, 0, 0, \
      0.500000, -0.500000, 0, 0.500000, 0, 0, 0, 0, 0, 0, 0, -2.00000, 0, \
      0, 0, 0, 0, 0, 0, 0.500000, 0, -0.500000, 0.500000, 0, 0, 0, \
      -2.00000, 0, 0, 0, 0, -0.500000, 0, 0, 0, 0, -0.500000, 0, 0, 0, 0, \
      0, 0, 0, 2.00000, 0.500000, -0.500000, 0, 0.500000, 0, 0, 0, 0, 0, 0, \
      0, 0, 0, 0, 2.00000, 0, 0, 0, 0, 0.500000, 0, -0.500000, 0.500000, 0, \
      0, 0, 0, 0, 0, 2.00000, 0.125000, -0.250000, 0, 0.125000, 0.250000, \
      0, -0.375000, 0.250000, 0, 0, 0, -2.00000, 0, 0, 0, 0.250000, \
      -0.375000, 0, 0.250000, 0.125000, 0, -0.250000, 0.125000, 0, 0, 0, \
      -1.00000, 0, 0, 1.00000, 0.125000, -0.375000, 0, 0.125000, 0.125000, \
      0, -0.375000, 0.125000, 0, 0, 0, 0, 0, 0, 2.00000, 0.125000, \
      -0.500000, 0, 0.125000, 0, 0, -0.375000, 0, 0, 0, 0, -0.250000, 0, 0, \
      1.75000, 0, -0.250000, 0, 0, 0.250000, 0, -0.500000, 0.250000, 0, 0, \
      0, -1.25000, 0, 0, 0.750000, 0.250000, -0.250000, 0, 0.250000, \
      0.250000, 0, -0.250000, 0.250000, 0, 0, 0, -1.00000, 0, 0, 1.00000};

  // DIV: each row pair gives the 5 correct values of the divergence of the 5 basis functions
  double basisDivs[] = {   
      // 6 vertices
      1.0, 1.0, 1.0, 1.0, 1.0,
      1.0, 1.0, 1.0, 1.0, 1.0,
      1.0, 1.0, 1.0, 1.0, 1.0,
      1.0, 1.0, 1.0, 1.0, 1.0,
      1.0, 1.0, 1.0, 1.0, 1.0,
      1.0, 1.0, 1.0, 1.0, 1.0,
      // 6 other points
      1.0, 1.0, 1.0, 1.0, 1.0,
      1.0, 1.0, 1.0, 1.0, 1.0,
      1.0, 1.0, 1.0, 1.0, 1.0,
      1.0, 1.0, 1.0, 1.0, 1.0,
      1.0, 1.0, 1.0, 1.0, 1.0,
      1.0, 1.0, 1.0, 1.0, 1.0
  };

  try{

    DynRankViewHost ConstructWithLabel(wedgeNodesHost, 12, 3);
    wedgeNodesHost(0,0) =  0.0;  wedgeNodesHost(0,1) =  0.0;  wedgeNodesHost(0,2) = -1.0;
    wedgeNodesHost(1,0) =  1.0;  wedgeNodesHost(1,1) =  0.0;  wedgeNodesHost(1,2) = -1.0;
    wedgeNodesHost(2,0) =  0.0;  wedgeNodesHost(2,1) =  1.0;  wedgeNodesHost(2,2) = -1.0;
    wedgeNodesHost(3,0) =  0.0;  wedgeNodesHost(3,1) =  0.0;  wedgeNodesHost(3,2) =  1.0;
    wedgeNodesHost(4,0) =  1.0;  wedgeNodesHost(4,1) =  0.0;  wedgeNodesHost(4,2) =  1.0;
    wedgeNodesHost(5,0) =  0.0;  wedgeNodesHost(5,1) =  1.0;  wedgeNodesHost(5,2) =  1.0;

    wedgeNodesHost(6,0) =  0.25; wedgeNodesHost(6,1) =  0.5;  wedgeNodesHost(6,2) = -1.0;
    wedgeNodesHost(7,0) =  0.5;  wedgeNodesHost(7,1) =  0.25; wedgeNodesHost(7,2) =  0.0;
    wedgeNodesHost(8,0) =  0.25; wedgeNodesHost(8,1) =  0.25; wedgeNodesHost(8,2) =  1.0;
    wedgeNodesHost(9,0) =  0.25; wedgeNodesHost(9,1) =  0.0;  wedgeNodesHost(9,2) =  0.75;
    wedgeNodesHost(10,0)=  0.0;  wedgeNodesHost(10,1)=  0.5;  wedgeNodesHost(10,2)= -0.25;
    wedgeNodesHost(11,0)=  0.5;  wedgeNodesHost(11,1)=  0.5;  wedgeNodesHost(11,2)=  0.0;

    const auto wedgeNodes = Kokkos::create_mirror_view(typename DeviceSpaceType::memory_space(), wedgeNodesHost);
    Kokkos::deep_copy(wedgeNodes, wedgeNodesHost);

    // Dimensions for the output arrays:
    const auto numFields = wedgeBasis.getCardinality();
    const auto numPoints = wedgeNodes.dimension(0);
    const auto spaceDim  = wedgeBasis.getBaseCellTopology().getDimension();


    // Check VALUE of basis functions: resize vals to rank-3 container:
    {
      DynRankView ConstructWithLabel(vals, numFields, numPoints, spaceDim);
      wedgeBasis.getValues(vals, wedgeNodes, OPERATOR_VALUE);
      const auto vals_host = Kokkos::create_mirror_view(typename HostSpaceType::memory_space(), vals);
      Kokkos::deep_copy(vals_host, vals);
      for (int i = 0; i < numFields; i++) {
        for (size_type j = 0; j < numPoints; j++) {
          for (size_type k = 0; k < spaceDim; k++) {
            int l = k + i * spaceDim + j * spaceDim * numFields;
            if (std::abs(vals_host(i,j,k) - basisValues[l]) > tol) {
              errorFlag++;
              *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";

              // Output the multi-index of the value where the error is:
              *outStream << " At multi-index { ";
              *outStream << i << " ";*outStream << j << " ";*outStream << k << " ";
              *outStream << "}  computed value: " << vals_host(i,j,k)
                   << " but reference value: " << basisValues[l] << "\n";
            }
          }
        }
      }
    }


    // Check DIV of basis function: resize vals to rank-2 container
    {
      DynRankView ConstructWithLabel(vals, numFields, numPoints);
      wedgeBasis.getValues(vals, wedgeNodes, OPERATOR_DIV);
      const auto vals_host = Kokkos::create_mirror_view(typename HostSpaceType::memory_space(), vals);
      Kokkos::deep_copy(vals_host, vals);
      for (int i = 0; i < numFields; i++) {
        for (size_type j = 0; j < numPoints; j++) {
          int l =  i + j * numFields;
          if (std::abs(vals_host(i,j) - basisDivs[l]) > tol) {
            errorFlag++;
            *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";

            // Output the multi-index of the value where the error is:
            *outStream << " At multi-index { ";
            *outStream << i << " ";*outStream << j << " ";
            *outStream << "}  computed divergence component: " << vals_host(i,j)
                   << " but reference divergence component: " << basisDivs[l] << "\n";
          }
        }
      }
    }
  } catch (std::logic_error err) {
    *outStream << err.what() << "\n\n";
    errorFlag = -1000;
  }

  *outStream
  << "\n"
  << "===============================================================================\n"
  << "| TEST 4: correctness of DoF locations                                        |\n"
  << "===============================================================================\n";

  try {
    const auto numFields = wedgeBasis.getCardinality();
    const auto spaceDim  = wedgeBasis.getBaseCellTopology().getDimension();

    // Check exceptions.
    ordinal_type nthrow = 0, ncatch = 0;
#ifdef HAVE_INTREPID2_DEBUG
    {
      DynRankView ConstructWithLabel(badVals, 1,2,3);
      INTREPID2_TEST_ERROR_EXPECTED( wedgeBasis.getDofCoords(badVals) );
    }
    {
      DynRankView ConstructWithLabel(badVals, 4,3);
      INTREPID2_TEST_ERROR_EXPECTED( wedgeBasis.getDofCoords(badVals) );
    }
    {
      DynRankView ConstructWithLabel(badVals, 5,2);
      INTREPID2_TEST_ERROR_EXPECTED( wedgeBasis.getDofCoords(badVals) );
    }
#endif
    if (nthrow != ncatch) {
      errorFlag++;
      *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
      *outStream << "# of catch ("<< ncatch << ") is different from # of throw (" << ncatch << ")\n";
    }

    DynRankView ConstructWithLabel(bvals, numFields, numFields, spaceDim);
    DynRankView ConstructWithLabel(cvals, numFields, spaceDim);

    // Check mathematical correctness.
    wedgeBasis.getDofCoords(cvals);
    wedgeBasis.getValues(bvals, cvals, OPERATOR_VALUE);

    // Check mathematical correctness
    DynRankViewHost ConstructWithLabel(normals, numFields,spaceDim); // normals at each point basis point
    normals(0,0)  =  0.0; normals(0,1)  = -2.0; normals(0,2)  =  0.0;
    normals(1,0)  =  2.0; normals(1,1)  =  2.0; normals(1,2)  =  0.0;
    normals(2,0)  = -2.0; normals(2,1)  =  0.0; normals(2,2)  =  0.0;
    normals(3,0)  =  0.0; normals(3,1)  =  0.0; normals(3,2)  = -0.5;
    normals(4,0)  =  0.0; normals(4,1)  =  0.0; normals(4,2)  =  0.5;

    auto cvals_host = Kokkos::create_mirror_view(typename HostSpaceType::memory_space(), cvals);
    Kokkos::deep_copy(cvals_host, cvals);

    auto bvals_host = Kokkos::create_mirror_view(typename HostSpaceType::memory_space(), bvals);
    Kokkos::deep_copy(bvals_host, bvals);

    for (ordinal_type i=0;i<numFields;++i) {
      for (ordinal_type j=0;j<numFields;++j) {

        ValueType normal = 0.0;
        for(size_type d=0;d<spaceDim;++d) {
          normal += bvals_host(i,j,d)*normals(j,d);
        }

        const ValueType expected_normal = (i == j);
        if (std::abs(normal - expected_normal) > tol || std::isnan(normal)) {
          errorFlag++;
          std::stringstream ss;
          ss << "\nNormal component of basis function " << i << " at (" << cvals_host(j,0) << ", " << cvals_host(j,1)<< ", " << cvals_host(j,2) << ") is " << normal << " but should be " << expected_normal << "\n";
          *outStream << ss.str();
        }
      }
    }
  } catch (std::logic_error err) {
    *outStream << err.what() << "\n\n";
    errorFlag = -1000;
  }

  if (errorFlag != 0)
    std::cout << "End Result: TEST FAILED\n";
  else
    std::cout << "End Result: TEST PASSED\n";

  // reset format state of std::cout
  std::cout.copyfmt(oldFormatState);

  return errorFlag;
}
}
}
