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

/** \file test_01.cpp
    \brief  Unit tests for the Intrepid2::Basis_HGRAD_LINE_Cn_FEM class.
    \author Created by P. Bochev, D. Ridzal, K. Peterson, Kyungjoo Kim
*/

#include "Intrepid2_config.h"

#ifdef HAVE_INTREPID2_DEBUG
#define INTREPID2_TEST_FOR_DEBUG_ABORT_OVERRIDE_TO_CONTINUE
#endif

#include "Intrepid2_Types.hpp"
#include "Intrepid2_Utils.hpp"

//#include "Intrepid2_PointTools.hpp"
#include "Intrepid2_HGRAD_LINE_Cn_FEM.hpp"

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

    template<typename OutValueType, typename PointValueType, typename DeviceSpaceType>
    int HGRAD_LINE_Cn_FEM_Test01(const bool verbose) {

      Teuchos::RCP<std::ostream> outStream;
      Teuchos::oblackholestream bhs; // outputs nothing

      if (verbose)
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
        << "|               Unit Test (Basis_HGRAD_LINE_Cn_FEM)                           |\n"
        << "|                                                                             |\n"
        << "|     1) Conversion of Dof tags into Dof ordinals and back                    |\n"
        << "|     2) Basis values for VALUE, GRAD, CURL, and Dk operators                 |\n"
        << "|                                                                             |\n"
        << "|  Questions? Contact  Pavel Bochev  (pbboche@sandia.gov),                    |\n"
        << "|                      Denis Ridzal  (dridzal@sandia.gov),                    |\n"
        << "|                      Kara Peterson (kjpeter@sandia.gov),                    |\n"
        << "|                      Kyungjoo Kim  (kyukim@sandia.gov).                     |\n"
        << "|                                                                             |\n"
        << "|  Intrepid's website: http://trilinos.sandia.gov/packages/intrepid           |\n"
        << "|  Trilinos website:   http://trilinos.sandia.gov                             |\n"
        << "|                                                                             |\n"
        << "===============================================================================\n";

      typedef Kokkos::DynRankView<PointValueType,DeviceSpaceType> DynRankViewPointValueType;
      typedef Kokkos::DynRankView<OutValueType,DeviceSpaceType> DynRankViewOutValueType;
      typedef typename ScalarTraits<OutValueType>::scalar_type scalar_type;
      typedef Kokkos::DynRankView<scalar_type, DeviceSpaceType> DynRankViewScalarValueType;
      //typedef Kokkos::DynRankView<PointValueType,HostSpaceType> DynRankViewHostPointValueType;

#define ConstructWithLabelScalar(obj, ...) obj(#obj, __VA_ARGS__)

      constexpr int spaceDim = 1;
      const scalar_type tol = tolerence();
      int errorFlag = 0;

      // for virtual function, value and point types are declared in the class
      typedef OutValueType outputValueType;
      typedef PointValueType pointValueType;
      //typedef ValueType weightValueType;

      typedef Basis_HGRAD_LINE_Cn_FEM<DeviceSpaceType,outputValueType,pointValueType> LineBasisType;
      constexpr ordinal_type maxOrder = Parameters::MaxOrder;

      *outStream
        << "\n"
        << "===============================================================================\n"
        << "| TEST 1: Basis creation, exceptions tests                                    |\n"
        << "===============================================================================\n";

      try{
#ifdef HAVE_INTREPID2_DEBUG
        ordinal_type nthrow = 0, ncatch = 0;
        constexpr ordinal_type order = 5;
          if(order <= maxOrder) {

          LineBasisType lineBasis(order);

          // Define array containing array of nodes to evaluate
          DynRankViewPointValueType ConstructWithLabelPointView(lineNodes, 10, 1);

          // Generic array for the output values; needs to be properly resized depending on the operator type
          const auto numFields = lineBasis.getCardinality();
          const auto numPoints = lineNodes.extent(0);


          // Exceptions 1-5: all bf tags/bf Ids below are wrong and should cause getDofOrdinal() and
          // getDofTag() to access invalid array elements thereby causing bounds check exception
          {
            INTREPID2_TEST_ERROR_EXPECTED( lineBasis.getDofOrdinal(2,0,0) );
            INTREPID2_TEST_ERROR_EXPECTED( lineBasis.getDofOrdinal(1,1,1) );
            INTREPID2_TEST_ERROR_EXPECTED( lineBasis.getDofOrdinal(1,0,7) );
            INTREPID2_TEST_ERROR_EXPECTED( lineBasis.getDofTag(numFields) );
            INTREPID2_TEST_ERROR_EXPECTED( lineBasis.getDofTag(-1) );

            // no exception; if happens, it is unexpected;
            lineBasis.getDofOrdinal(1,0,3);
            lineBasis.getDofTag(5);
          }

          // Exceptions 6-16 test exception handling with incorrectly dimensioned input/output arrays
          {
            DynRankViewOutValueType ConstructWithLabelOutView(vals, numFields, numPoints);
            {
              // exception #6: input points array must be of rank-2
              DynRankViewPointValueType ConstructWithLabelPointView(badPoints, 4, 5, 3);
              INTREPID2_TEST_ERROR_EXPECTED( lineBasis.getValues(vals, badPoints, OPERATOR_VALUE) );
            }
            {
              // exception #7: dimension 1 in the input point array must equal space dimension of the cell
              DynRankViewPointValueType ConstructWithLabelPointView(badPoints, 4, 3);
              INTREPID2_TEST_ERROR_EXPECTED( lineBasis.getValues(vals, badPoints, OPERATOR_VALUE) );
            }
            {
              // exception #8: output values must be of rank-2 for OPERATOR_VALUE
              DynRankViewOutValueType ConstructWithLabelOutView(badVals, 4, 3, 1);
              INTREPID2_TEST_ERROR_EXPECTED( lineBasis.getValues(badVals, lineNodes, OPERATOR_VALUE) );
            }
            lineBasis.getValues(vals, lineNodes, OPERATOR_VALUE);
          }
          {
            DynRankViewOutValueType ConstructWithLabelOutView(badVals, 4, 3);
            
            // exception #9: output values must be of rank-3 for OPERATOR_GRAD
            INTREPID2_TEST_ERROR_EXPECTED( lineBasis.getValues(badVals, lineNodes, OPERATOR_GRAD) );

            // exception #10: output values must be of rank-3 for OPERATOR_CURL
            INTREPID2_TEST_ERROR_EXPECTED( lineBasis.getValues(badVals, lineNodes, OPERATOR_CURL) );

            // exception #11: output values must be of rank-2 for OPERATOR_DIV
            INTREPID2_TEST_ERROR_EXPECTED( lineBasis.getValues(badVals, lineNodes, OPERATOR_DIV) );

            // exception #12: output values must be of rank-2 for OPERATOR_D1
            INTREPID2_TEST_ERROR_EXPECTED( lineBasis.getValues(badVals, lineNodes, OPERATOR_D1) );
          }
          {
            // exception #13: incorrect 0th dimension of output array (must equal number of basis functions)
            DynRankViewOutValueType ConstructWithLabelOutView(badVals, lineBasis.getCardinality() + 1, lineNodes.extent(0));
            INTREPID2_TEST_ERROR_EXPECTED( lineBasis.getValues(badVals, lineNodes, OPERATOR_VALUE) );
          }
          {
            // exception #14: incorrect 1st dimension of output array (must equal number of points)
            DynRankViewOutValueType ConstructWithLabelOutView(badVals, lineBasis.getCardinality(), lineNodes.extent(0) + 1);
            INTREPID2_TEST_ERROR_EXPECTED( lineBasis.getValues(badVals, lineNodes, OPERATOR_VALUE) );
          }
          {
            // exception #15: incorrect 2nd dimension of output array (must equal spatial dimension)
            DynRankViewOutValueType ConstructWithLabelOutView(badVals, lineBasis.getCardinality(), lineNodes.extent(0), 2);
            INTREPID2_TEST_ERROR_EXPECTED( lineBasis.getValues(badVals, lineNodes, OPERATOR_GRAD) );
          }
        }
        if (nthrow != ncatch) {
          errorFlag++;
          *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
          *outStream << "# of catch ("<< ncatch << ") is different from # of throw (" << nthrow << ")\n";
        }
#endif
      } catch (std::exception err) {
        *outStream << "UNEXPECTED ERROR !!! ----------------------------------------------------------\n";
        *outStream << err.what() << '\n';
        *outStream << "-------------------------------------------------------------------------------" << "\n\n";
        errorFlag = -1000;
      };
      
      *outStream
        << "\n"
        << "===============================================================================\n"
        << "| TEST 2: correctness of basis function values                                |\n"
        << "===============================================================================\n";
      outStream->precision(20);

      // check Kronecker property for Lagrange polynomials.      
      try {
        const EPointType pts[3] = { POINTTYPE_EQUISPACED, POINTTYPE_WARPBLEND, POINTTYPE_GAUSS };
        for (auto idx=0;idx<3;++idx) {
          *outStream << " -- Testing " << EPointTypeToString(pts[idx]) << " -- \n";
          for (auto ip=1;ip<maxOrder;++ip) {

            LineBasisType lineBasis(ip, pts[idx]);

            const auto numDofs   = lineBasis.getCardinality();
            const auto numPoints = lineBasis.getCardinality();
            
            DynRankViewScalarValueType ConstructWithLabelScalar(dofCoords_scalar, numPoints, spaceDim);
            lineBasis.getDofCoords(dofCoords_scalar);
            
            DynRankViewPointValueType ConstructWithLabelPointView(dofCoords, numDofs , spaceDim);
            RealSpaceTools<DeviceSpaceType>::clone(dofCoords,dofCoords_scalar);

            DynRankViewOutValueType ConstructWithLabelOutView(vals, numDofs, numPoints);
            lineBasis.getValues(vals, dofCoords, OPERATOR_VALUE);

            // host mirror for comparison
            auto valsHost = Kokkos::create_mirror_view(vals);
            Kokkos::deep_copy(valsHost, vals);
            
            for (auto i=0;i<numDofs;++i) 
              for (int j=0;j<numPoints;++j) {
                const scalar_type  exactVal = (i == j);
                const auto val = get_scalar_value(valsHost(i,j));
                if (std::isnan(val) || std::abs(val-exactVal) > tol) {
                  errorFlag++;
                  *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
                  *outStream << " Basis function at i= " << i << ", j=" << j << ": "
                             << valsHost(i,j) << " != " << exactVal << "\n\n";
                }
              }
          }
        }
      } catch (std::exception err) {
        *outStream << err.what() << "\n\n";
        errorFlag = -1000;
      };

      if (errorFlag != 0)
        std::cout << "End Result: TEST FAILED\n";
      else
        std::cout << "End Result: TEST PASSED\n";
      
      std::cout.copyfmt(oldFormatState);

      return errorFlag;
    }
  }
}
