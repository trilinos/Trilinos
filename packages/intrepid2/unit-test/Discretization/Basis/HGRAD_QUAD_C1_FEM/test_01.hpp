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
    \brief  Unit tests for the Intrepid2::HGRAD_QUAD_C1_FEM class.
    \author Created by P. Bochev, D. Ridzal, and K. Peterson.
*/
#include "Intrepid2_config.h"

#ifdef HAVE_INTREPID2_DEBUG
#define INTREPID2_TEST_FOR_DEBUG_ABORT_OVERRIDE_TO_CONTINUE
#endif

#include "Intrepid2_Basis.hpp"
#include "Intrepid2_HGRAD_QUAD_C1_FEM.hpp"

#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"

namespace Intrepid2 {

  namespace Test {

#define INTREPID2_TEST_ERROR_EXPECTED( S )              \
    try {                                                               \
      ++nthrow;                                                         \
      S ;                                                               \
    }                                                                   \
    catch (std::exception err) {                                      \
      ++ncatch;                                                         \
      *outStream << "Expected Error ----------------------------------------------------------------\n"; \
      *outStream << err.what() << '\n';                                 \
      *outStream << "-------------------------------------------------------------------------------" << "\n\n"; \
    }

    template<typename ValueType, typename DeviceSpaceType>
    int HGRAD_QUAD_C1_FEM_Test01(const bool verbose) {

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
        << "|                 Unit Test (Basis_HGRAD_QUAD_C1_FEM)                         |\n"
        << "|                                                                             |\n"
        << "|     1) Conversion of Dof tags into Dof ordinals and back                    |\n"
        << "|     2) Basis values for VALUE, GRAD, CURL, and Dk operators                 |\n"
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
      typedef Kokkos::DynRankView<ValueType,HostSpaceType>   DynRankViewHost;
#define ConstructWithLabel(obj, ...) obj(#obj, __VA_ARGS__)

      const ValueType tol = tolerence();
      int errorFlag = 0;

      // for virtual function, value and point types are declared in the class
      typedef ValueType outputValueType;
      typedef ValueType pointValueType;
      Basis_HGRAD_QUAD_C1_FEM<DeviceSpaceType,outputValueType,pointValueType> quadBasis;
      //typedef typename decltype(quadBasis)::outputViewType outputViewType;
      //typedef typename decltype(quadBasis)::pointViewType  pointViewType;

      *outStream
        << "\n"
        << "===============================================================================\n"
        << "| TEST 1: Basis creation, exceptions tests                                    |\n"
        << "===============================================================================\n";


      try{
        ordinal_type nthrow = 0, ncatch = 0;
#ifdef HAVE_INTREPID2_DEBUG

        // Define array containing the 4 vertices of the reference QUAD and its center.
        DynRankView ConstructWithLabel(quadNodes, 5, 2);

        // Generic array for the output values; needs to be properly resized depending on the operator type
        const ordinal_type numFields = quadBasis.getCardinality();
        const ordinal_type numPoints = quadNodes.extent(0);
        const ordinal_type spaceDim  = quadBasis.getBaseCellTopology().getDimension();
        const ordinal_type D2Cardin  = getDkCardinality(OPERATOR_D2, spaceDim);

        const ordinal_type workSize  = numFields*numPoints*D2Cardin;
        DynRankView ConstructWithLabel(work, workSize);

        // resize vals to rank-2 container with dimensions
        DynRankView vals = DynRankView(work.data(), numFields, numPoints);

        {
          // exception #1: DIV cannot be applied to scalar functions
          INTREPID2_TEST_ERROR_EXPECTED( quadBasis.getValues(vals, quadNodes, OPERATOR_DIV) );
        }

        // Exceptions 2-6: all bf tags/bf Ids below are wrong and should cause getDofOrdinal() and
        // getDofTag() to access invalid array elements thereby causing bounds check exception
        {
          INTREPID2_TEST_ERROR_EXPECTED( quadBasis.getDofOrdinal(3,0,0) ); // #2
          INTREPID2_TEST_ERROR_EXPECTED( quadBasis.getDofOrdinal(1,1,1) ); // #3
          INTREPID2_TEST_ERROR_EXPECTED( quadBasis.getDofOrdinal(0,4,0) ); // #4
          INTREPID2_TEST_ERROR_EXPECTED( quadBasis.getDofTag(numFields) ); // #5
          INTREPID2_TEST_ERROR_EXPECTED( quadBasis.getDofTag(-1)        ); // #6
        }

        // Exceptions 7-17 test exception handling with incorrectly dimensioned input/output arrays
        {
          // exception #7: input points array must be of rank-2
          DynRankView ConstructWithLabel(badPoints, 4, 5, 3);
          INTREPID2_TEST_ERROR_EXPECTED( quadBasis.getValues(vals, badPoints, OPERATOR_VALUE) );
        }
        {
          // exception #8 dimension 1 in the input point array must equal space dimension of the cell
          DynRankView ConstructWithLabel(badPoints, 4, 3);
          INTREPID2_TEST_ERROR_EXPECTED( quadBasis.getValues(vals, badPoints, OPERATOR_VALUE) );
        }
        {
          // exception #9 output values must be of rank-2 for OPERATOR_VALUE
          DynRankView ConstructWithLabel(badVals, 4, 3, 1);
          INTREPID2_TEST_ERROR_EXPECTED( quadBasis.getValues(badVals, quadNodes, OPERATOR_VALUE) );
        }
        {
          // exception #10 output values must be of rank-3 for OPERATOR_GRAD
          DynRankView ConstructWithLabel(badVals, 4, 3);
          INTREPID2_TEST_ERROR_EXPECTED( quadBasis.getValues(badVals, quadNodes, OPERATOR_GRAD) );

          // exception #11 output values must be of rank-3 for OPERATOR_CURL
          INTREPID2_TEST_ERROR_EXPECTED( quadBasis.getValues(badVals, quadNodes, OPERATOR_CURL) );

          // exception #12 output values must be of rank-3 for OPERATOR_D2
          INTREPID2_TEST_ERROR_EXPECTED( quadBasis.getValues(badVals, quadNodes, OPERATOR_D2) );
        }
        {
          // exception #13 incorrect 0th dimension of output array (must equal number of basis functions)
          DynRankView ConstructWithLabel(badVals, quadBasis.getCardinality() + 1, quadNodes.extent(0));
          INTREPID2_TEST_ERROR_EXPECTED( quadBasis.getValues(badVals, quadNodes, OPERATOR_VALUE) );
        }
        {
          // exception #14 incorrect 1st dimension of output array (must equal number of points)
          DynRankView ConstructWithLabel(badVals, quadBasis.getCardinality(), quadNodes.extent(0) + 1);
          INTREPID2_TEST_ERROR_EXPECTED( quadBasis.getValues(badVals, quadNodes, OPERATOR_VALUE) );
        }
        {
          // exception #15: incorrect 2nd dimension of output array (must equal the space dimension)
          DynRankView ConstructWithLabel(badVals, quadBasis.getCardinality(), quadNodes.extent(0), 4);
          INTREPID2_TEST_ERROR_EXPECTED( quadBasis.getValues(badVals, quadNodes, OPERATOR_GRAD) );
        }
        {
          // exception #16: incorrect 2nd dimension of output array (must equal D2 cardinality in 2D)
          DynRankView ConstructWithLabel(badVals, quadBasis.getCardinality(), quadNodes.extent(0), 40);
          INTREPID2_TEST_ERROR_EXPECTED( quadBasis.getValues(badVals, quadNodes, OPERATOR_D2) );

          // exception #17: incorrect 2nd dimension of output array (must equal D3 cardinality in 2D)
          INTREPID2_TEST_ERROR_EXPECTED( quadBasis.getValues(badVals, quadNodes, OPERATOR_D3) );
        }
#endif
        if (nthrow != ncatch) {
          errorFlag++;
          *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
          *outStream << "# of catch ("<< ncatch << ") is different from # of throw (" << nthrow << ")\n";
        }

      } catch (std::logic_error err) {
        *outStream << "UNEXPECTED ERROR !!! ----------------------------------------------------------\n";
        *outStream << err.what() << '\n';
        *outStream << "-------------------------------------------------------------------------------" << "\n\n";
        errorFlag = -1000;
      };

      *outStream
        << "\n"
        << "===============================================================================\n"
        << "| TEST 2: correctness of tag to enum and enum to tag lookups                  |\n"
        << "===============================================================================\n";

      try{
        const ordinal_type numFields = quadBasis.getCardinality();
        const auto allTags = quadBasis.getAllDofTags();

        // Loop over all tags, lookup the associated dof enumeration and then lookup the tag again
        const ordinal_type dofTagSize = allTags.extent(0);
        for (ordinal_type i=0;i<dofTagSize;++i) {
          const auto bfOrd = quadBasis.getDofOrdinal(allTags(i,0), allTags(i,1), allTags(i,2));

          const auto myTag = quadBasis.getDofTag(bfOrd);
          if( !( (myTag(0) == allTags(i,0)) &&
                 (myTag(1) == allTags(i,1)) &&
                 (myTag(2) == allTags(i,2)) &&
                 (myTag(3) == allTags(i,3)) ) ) {
            errorFlag++;
            *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
            *outStream << " getDofOrdinal( {"
                       << allTags(i,0) << ", "
                       << allTags(i,1) << ", "
                       << allTags(i,2) << ", "
                       << allTags(i,3) << "}) = " << bfOrd <<" but \n";
            *outStream << " getDofTag(" << bfOrd << ") = { "
                       << myTag(0) << ", "
                       << myTag(1) << ", "
                       << myTag(2) << ", "
                       << myTag(3) << "}\n";
          }
        }

        // Now do the same but loop over basis functions
        for(ordinal_type bfOrd=0;bfOrd<numFields;++bfOrd) {
          const auto myTag  = quadBasis.getDofTag(bfOrd);
          const auto myBfOrd = quadBasis.getDofOrdinal(myTag(0), myTag(1), myTag(2));
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
      } catch (std::logic_error err){
        *outStream << "UNEXPECTED ERROR !!! ----------------------------------------------------------\n";
        *outStream << err.what() << '\n';
        *outStream << "-------------------------------------------------------------------------------" << "\n\n";
        errorFlag = -1000;
      };

      *outStream
        << "\n"
        << "===============================================================================\n"
        << "| TEST 3: correctness of basis function values                                |\n"
        << "===============================================================================\n";

      outStream -> precision(20);

      try {
        // VALUE: Each row gives the 4 correct basis set values at an evaluation point
        const ValueType basisValues[][4] = {
          { 1.0, 0.0, 0.0, 0.0  },
          { 0.0, 1.0, 0.0, 0.0  },
          { 0.0, 0.0, 1.0, 0.0  },
          { 0.0, 0.0, 0.0, 1.0  },
          { 0.25,0.25,0.25,0.25 }
        };

        // GRAD and D1: each row gives the 8 correct values of the gradients of the 4 basis functions
        const ValueType basisGrads[][4][2] = {
          { { -0.5, -0.5  },    { 0.5,  0.0  },     { 0.0,  0.0 },    {  0.0,  0.5 } },
          { { -0.5,  0.0  },    { 0.5, -0.5  },     { 0.0,  0.5 },    {  0.0,  0.0 } },
          { {  0.0,  0.0  },    { 0.0, -0.5  },     { 0.5,  0.5 },    { -0.5,  0.0 } },
          { {  0.0, -0.5  },    { 0.0,  0.0  },     { 0.5,  0.0 },    { -0.5,  0.5 } },
          { { -0.25,-0.25 },    { 0.25,-0.25 },     { 0.25, 0.25},    { -0.25, 0.25} }
        };

        // CURL: each row gives the 8 correct values of the curls of the 4 basis functions
        const ValueType basisCurls[][4][2] = {
          { { -0.5,  0.5 },    { 0.0, -0.5 },    { 0.0,  0.0 },    { 0.5,  0.0 } },
          { {  0.0,  0.5 },    {-0.5, -0.5 },    { 0.5,  0.0 },    { 0.0,  0.0 } },
          { {  0.0,  0.0 },    {-0.5,  0.0 },    { 0.5, -0.5 },    { 0.0,  0.5 } },
          { { -0.5,  0.0 },    { 0.0,  0.0 },    { 0.0, -0.5 },    { 0.5,  0.5 } },
          { {-0.25,  0.25},    {-0.25,-0.25},    { 0.25,-0.25},    { 0.25, 0.25} }
        };

        //D2: each row gives the 12 correct values of all 2nd derivatives of the 4 basis functions
        const ValueType basisD2[][4][3] = {
          { { 0.0, 0.25, 0.0 },   { 0.0,-0.25, 0.0 },   { 0.0, 0.25, 0.0 },   { 0.0,-0.25, 0.0 } },
          { { 0.0, 0.25, 0.0 },   { 0.0,-0.25, 0.0 },   { 0.0, 0.25, 0.0 },   { 0.0,-0.25, 0.0 } },
          { { 0.0, 0.25, 0.0 },   { 0.0,-0.25, 0.0 },   { 0.0, 0.25, 0.0 },   { 0.0,-0.25, 0.0 } },
          { { 0.0, 0.25, 0.0 },   { 0.0,-0.25, 0.0 },   { 0.0, 0.25, 0.0 },   { 0.0,-0.25, 0.0 } },
          { { 0.0, 0.25, 0.0 },   { 0.0,-0.25, 0.0 },   { 0.0, 0.25, 0.0 },   { 0.0,-0.25, 0.0 } }
        };


        // Define array containing the 4 vertices of the reference QUAD and its center.
        DynRankViewHost ConstructWithLabel(quadNodesHost, 5, 2);

        quadNodesHost(0,0) = -1.0;  quadNodesHost(0,1) = -1.0;
        quadNodesHost(1,0) =  1.0;  quadNodesHost(1,1) = -1.0;
        quadNodesHost(2,0) =  1.0;  quadNodesHost(2,1) =  1.0;
        quadNodesHost(3,0) = -1.0;  quadNodesHost(3,1) =  1.0;
        quadNodesHost(4,0) =  0.0;  quadNodesHost(4,1) =  0.0;

        auto quadNodes = Kokkos::create_mirror_view(typename DeviceSpaceType::memory_space(), quadNodesHost);
        Kokkos::deep_copy(quadNodes, quadNodesHost);

        // Generic array for the output values; needs to be properly resized depending on the operator type
        const ordinal_type numFields = quadBasis.getCardinality();
        const ordinal_type numPoints = quadNodes.extent(0);
        const ordinal_type spaceDim  = quadBasis.getBaseCellTopology().getDimension();
        const ordinal_type D2Cardin  = getDkCardinality(OPERATOR_D2, spaceDim);

        // Check VALUE of basis functions: resize vals to rank-2 container:
        {
          DynRankView ConstructWithLabel(vals, numFields, numPoints);
          quadBasis.getValues(vals, quadNodes, OPERATOR_VALUE);
          auto vals_host = Kokkos::create_mirror_view(typename HostSpaceType::memory_space(), vals);
          Kokkos::deep_copy(vals_host, vals);
          for (ordinal_type i=0;i<numFields;++i)
            for (ordinal_type j=0;j<numPoints;++j)
              if (std::abs(vals_host(i,j) - basisValues[j][i]) > tol) {
                errorFlag++;
                *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";

                // Output the multi-index of the value where the error is:
                *outStream << " At multi-index { ";
                *outStream << i << " ";*outStream << j << " ";
                *outStream << "}  computed value: " << vals_host(i,j)
                           << " but reference value: " << basisValues[j][i] << "\n";
              }
        }

        // Check GRAD of basis function: resize vals to rank-3 container
        {
          const EOperator ops[] = { OPERATOR_GRAD,
                                    OPERATOR_D1,
                                    OPERATOR_MAX };
          for (auto h=0;ops[h]!=OPERATOR_MAX;++h) {
            const auto op = ops[h];
            DynRankView ConstructWithLabel(vals, numFields, numPoints, spaceDim);
            quadBasis.getValues(vals, quadNodes, op);
            auto vals_host = Kokkos::create_mirror_view(typename HostSpaceType::memory_space(), vals);
            Kokkos::deep_copy(vals_host, vals);
            for (ordinal_type i=0;i<numFields;++i)
              for (ordinal_type j=0;j<numPoints;++j)
                for (ordinal_type k=0;k<spaceDim;++k)
                  if (std::abs(vals_host(i,j,k) - basisGrads[j][i][k]) > tol) {
                    errorFlag++;
                    *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";

                    // Output the multi-index of the value where the error is:
                    *outStream << " At multi-index { ";
                    *outStream << i << " ";*outStream << j << " ";*outStream << k << " ";
                    *outStream << "}  computed grad component: " << vals_host(i,j,k)
                               << " but reference grad component: " << basisGrads[j][i][k] << "\n";
                  }
          }
        }

        // Check CURL of basis function: resize vals just for illustration!
        {
          DynRankView ConstructWithLabel(vals, numFields, numPoints, spaceDim);
          quadBasis.getValues(vals, quadNodes, OPERATOR_CURL);
          auto vals_host = Kokkos::create_mirror_view(typename HostSpaceType::memory_space(), vals);
          Kokkos::deep_copy(vals_host, vals);
          for (ordinal_type i=0;i<numFields;++i)
            for (ordinal_type j=0;j<numPoints;++j)
              for (ordinal_type k=0;k<spaceDim;++k)
                if (std::abs(vals_host(i,j,k) - basisCurls[j][i][k]) > tol) {
                  errorFlag++;
                  *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";

                  // Output the multi-index of the value where the error is:
                  *outStream << " At multi-index { ";
                  *outStream << i << " ";*outStream << j << " ";*outStream << k << " ";
                  *outStream << "}  computed curl component: " << vals_host(i,j,k)
                             << " but reference curl component: " << basisCurls[j][i][k] << "\n";
                }
        }


        // Check D2 of basis function
        {
          DynRankView ConstructWithLabel(vals, numFields, numPoints, D2Cardin);
          quadBasis.getValues(vals, quadNodes, OPERATOR_D2);
          auto vals_host = Kokkos::create_mirror_view(typename HostSpaceType::memory_space(), vals);
          Kokkos::deep_copy(vals_host, vals);
          for (ordinal_type i=0;i<numFields;++i)
            for (ordinal_type j=0;j<numPoints;++j)
              for (ordinal_type k=0;k<D2Cardin;++k)
                if (std::abs(vals_host(i,j,k) - basisD2[j][i][k]) > tol) {
                  errorFlag++;
                  *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";

                  // Output the multi-index of the value where the error is:
                  *outStream << " At multi-index { ";
                  *outStream << i << " ";*outStream << j << " ";*outStream << k << " ";
                  *outStream << "}  computed D2 component: " << vals_host(i,j,k)
                             << " but reference D2 component: " << basisD2[j][i][k] << "\n";
                }
        }


        // Check all higher derivatives - must be zero.
        {
          const EOperator ops[] = { OPERATOR_D3,
                                    OPERATOR_D4,
                                    OPERATOR_D5,
                                    OPERATOR_D6,
                                    OPERATOR_D7,
                                    OPERATOR_D8,
                                    OPERATOR_D9,
                                    OPERATOR_D10,
                                    OPERATOR_MAX };
          for (auto h=0;ops[h]!=OPERATOR_MAX;++h) {
            const auto op = ops[h];
            const ordinal_type DkCardin  = getDkCardinality(op, spaceDim);
            DynRankView ConstructWithLabel(vals, numFields, numPoints, DkCardin);
            quadBasis.getValues(vals, quadNodes, op);
            auto vals_host = Kokkos::create_mirror_view(typename HostSpaceType::memory_space(), vals);
            Kokkos::deep_copy(vals_host, vals);
            for (ordinal_type i1=0;i1<numFields;++i1)
              for (ordinal_type i2=0;i2<numPoints;++i2)
                for (ordinal_type i3=0;i3<DkCardin;++i3)
                  if (std::abs(vals_host(i1,i2,i3)) > tol) {
                    errorFlag++;
                    *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";

                    // Get the multi-index of the value where the error is and the operator order
                    const auto ord = Intrepid2::getOperatorOrder(op);
                    *outStream << " At multi-index { "<<i1<<" "<<i2 <<" "<<i3;
                    *outStream << "}  computed D"<< ord <<" component: " << vals_host(i1,i2,i3)
                               << " but reference D" << ord << " component:  0 \n";
                  }
          }
        }
      } catch (std::logic_error err) {
        *outStream << "UNEXPECTED ERROR !!! ----------------------------------------------------------\n";
        *outStream << err.what() << '\n';
        *outStream << "-------------------------------------------------------------------------------" << "\n\n";
        errorFlag = -1000;
      };


      *outStream
        << "\n"
        << "===============================================================================\n"
        << "| TEST 4: correctness of DoF locations                                        |\n"
        << "===============================================================================\n";

      try{
        Basis_HGRAD_QUAD_C1_FEM<DeviceSpaceType> quadBasis;
        const ordinal_type numFields = quadBasis.getCardinality();
        const ordinal_type spaceDim  = quadBasis.getBaseCellTopology().getDimension();

        // Check exceptions.
        ordinal_type nthrow = 0, ncatch = 0;
#ifdef HAVE_INTREPID2_DEBUG
        {
          DynRankView ConstructWithLabel(badVals, 1,2,3);
          INTREPID2_TEST_ERROR_EXPECTED( quadBasis.getDofCoords(badVals) );
        }
        {
          DynRankView ConstructWithLabel(badVals, 3,2);
          INTREPID2_TEST_ERROR_EXPECTED( quadBasis.getDofCoords(badVals) );
        }
        {
          DynRankView ConstructWithLabel(badVals, 4,3);
          INTREPID2_TEST_ERROR_EXPECTED( quadBasis.getDofCoords(badVals) );
        }
#endif
        if (nthrow != ncatch) {
          errorFlag++;
          *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
          *outStream << "# of catch ("<< ncatch << ") is different from # of throw (" << ncatch << ")\n";
        }

        DynRankView ConstructWithLabel(bvals, numFields, numFields);
        DynRankView ConstructWithLabel(cvals, numFields, spaceDim);

        // Check mathematical correctness.
        quadBasis.getDofCoords(cvals);
        quadBasis.getValues(bvals, cvals, OPERATOR_VALUE);

        auto bvals_host = Kokkos::create_mirror_view(typename HostSpaceType::memory_space(), bvals);
        Kokkos::deep_copy(bvals_host, bvals);

        auto cvals_host = Kokkos::create_mirror_view(typename HostSpaceType::memory_space(), cvals);
        Kokkos::deep_copy(cvals_host, cvals);

        for (ordinal_type i=0;i<numFields;++i) {
          for (ordinal_type j=0;j<numFields;++j) {
            if (i != j && (std::abs(bvals_host(i,j) - 0.0) > tol)) {
              errorFlag++;
              std::stringstream ss;
              ss << "\n Value of basis function " << i << " at (" << cvals_host(i,0) << ", " << cvals_host(i,1) << ") is " << bvals_host(i,j) << " but should be 0.0\n";
              *outStream << ss.str();
            }
            else if ((i == j) && (std::abs(bvals_host(i,j) - 1.0) > tol)) {
              errorFlag++;
              std::stringstream ss;
              ss << "\n Value of basis function " << i << " at (" << cvals_host(i,0) << ", " << cvals_host(i,1) << ") is " << bvals_host(i,j) << " but should be 1.0\n";
              *outStream << ss.str();
            }
          }
        }

      } catch (std::logic_error err){
        *outStream << "UNEXPECTED ERROR !!! ----------------------------------------------------------\n";
        *outStream << err.what() << '\n';
        *outStream << "-------------------------------------------------------------------------------" << "\n\n";
        errorFlag = -1000;
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
  }
}
