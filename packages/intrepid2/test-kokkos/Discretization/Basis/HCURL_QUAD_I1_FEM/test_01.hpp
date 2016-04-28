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

/** \file   test_01.cpp
    \brief  Unit tests for the Intrepid2::HCURL_QUAD_I1_FEM class.
    \author Created by P. Bochev, D. Ridzal and K. Peterson.
*/
//#include "Teuchos_GlobalMPISession.hpp"

#include "Intrepid2_config.h"

#ifdef HAVE_INTREPID2_DEBUG
#define INTREPID2_TEST_FOR_DEBUG_ABORT_OVERRIDE_TO_CONTINUE
#endif

#include "Intrepid2_HCURL_QUAD_I1_FEM.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"

using namespace std;

namespace Intrepid2 {
namespace Test {

#define INTREPID2_TEST_ERROR_EXPECTED ( S , nthrows, ncatch )                                                              \
{                                                                                                                          \
  try {                                                                                                                    \
    ++nthrows ;                                                                                                                    \
    S ;                                                                                                                    \
  }                                                                                                                        \
  catch (std::logic_error err) {                                                                                           \
      ++ncatch;                                                                                                      \
      *outStream << "Expected Error " << nException << " -------------------------------------------------------------\n"; \
      *outStream << err.what() << '\n';                                                                                    \
      *outStream << "-------------------------------------------------------------------------------" << "\n\n";           \
  };                                                                                                                       \
}

  template<typename ValueType, typename DeviceSpaceType>
  int HCURL_QUAD_I1_FEM_Test01(const bool verbose) {

  typedef ValueType value_type;

  Teuchos::RCP<std::ostream> outStream;
  Teuchos::oblackholestream bhs; // outputs nothing

  if (verbose)
    outStream = Teuchos::rcp(&std::cout, false);
  else
    outStream = Teuchos::rcp(&bhs, false);
  
  // Save the format state of the original std::cout.
  Teuchos::oblackholestream oldFormatState;
  oldFormatState.copyfmt(std::cout);
  
  *outStream \
    << "===============================================================================\n" \
    << "|                                                                             |\n" \
    << "|                 Unit Test (Basis_HCURL_QUAD_I1_FEM)                          |\n" \
    << "|                                                                             |\n" \
    << "|     1) Conversion of Dof tags into Dof ordinals and back                    |\n" \
    << "|     2) Basis values for VALUE and CURL operators                            |\n" \
    << "|                                                                             |\n" \
    << "|  Questions? Contact  Pavel Bochev  (pbboche@sandia.gov),                    |\n" \
    << "|                      Denis Ridzal  (dridzal@sandia.gov),                    |\n" \
    << "|                      Kara Peterson (kjpeter@sandia.gov).                    |\n" \
    << "|                                                                             |\n" \
    << "|  Intrepid's website: http://trilinos.sandia.gov/packages/intrepid           |\n" \
    << "|  Trilinos website:   http://trilinos.sandia.gov                             |\n" \
    << "|                                                                             |\n" \
    << "===============================================================================\n"\
    << "| TEST 1: Basis creation, exception testing                                   |\n"\
    << "===============================================================================\n";

    typedef Kokkos::DynRankView<value_type,DeviceSpaceType> DynRankView;
#define ConstructWithLabel(obj, ...) obj(#obj, __VA_ARGS__)

    const value_type tol = Parameters::Tolerence;
    int errorFlag = 0;
  
// ------------------------------------------------------------
  ordinal_type nthrow = 0, ncatch = 0;
  try{
#ifdef HAVE_INTREPID2_DEBUG
    Basis_HCURL_QUAD_I1_FEM<DeviceSpaceType> quadBasis;

    // Array with the 4 vertices of the reference Quadrilateral, its center and 4 more points
    DynRankView ConstructWithLabel(quadNodes, 9, 2);

    // Generic array for the output values; needs to be properly resized depending on the operator type
    const auto numFields = quadBasis.getCardinality();
    const auto numPoints = quadNodes.dimension(0);
    const auto spaceDim  = quadBasis.getBaseCellTopology().getDimension();
    const auto D2Cardin  = getDkCardinality(OPERATOR_D2, spaceDim); // is this needed?

    const auto workSize  = numFields*numPoints*D2Cardin;
    DynRankView ConstructWithLabel(work, workSize);

    // resize vals to rank-3 container with dimensions (num. basis functions, num. points, arbitrary)
//    vals.resize(quadBasis.getCardinality(), quadNodes.dimension(0), 4 );
    DynRankView vals = DynRankView(work.data(), numFields, numPoints, 4);

{
    // exception #1: GRAD cannot be applied to HCURL functions 
    INTREPID2_TEST_COMMAND( quadBasis.getValues(vals, quadNodes, OPERATOR_GRAD) );
//    INTREPID2_TEST_ERROR_EXPECTED( quadBasis.getValues(vals, quadNodes, OPERATOR_GRAD), throwCounter, nException );
}

{
    // exception #2: DIV cannot be applied to HCURL functions
    // resize vals to rank-2 container with dimensions (num. points, num. basis functions)
//    vals.resize(quadBasis.getCardinality(), quadNodes.dimension(0) );
    vals = DynRankView(work.data(), numFields, numPoints);
    INTREPID2_TEST_ERROR_EXPECTED( quadBasis.getValues(vals, quadNodes, OPERATOR_DIV) );
}
{
    // Exceptions 3-7: all bf tags/bf Ids below are wrong and should cause getDofOrdinal() and 
    // getDofTag() to access invalid array elements thereby causing bounds check exception
    // exception #3
    INTREPID2_TEST_ERROR_EXPECTED( quadBasis.getDofOrdinal(3,0,0) );
    // exception #4
    INTREPID2_TEST_ERROR_EXPECTED( quadBasis.getDofOrdinal(1,1,1) );
    // exception #5
    INTREPID2_TEST_ERROR_EXPECTED( quadBasis.getDofOrdinal(0,4,1) );
    // exception #6
    INTREPID2_TEST_ERROR_EXPECTED( quadBasis.getDofTag(12) );
    // exception #7
    INTREPID2_TEST_ERROR_EXPECTED( quadBasis.getDofTag(-1) );
}

    // Exceptions 8-15 test exception handling with incorrectly dimensioned input/output arrays
    // exception #8: input points array must be of rank-2
{
    DynRankView ConstructWithLabel(badPoints1, 4, 5, 3);
    INTREPID2_TEST_ERROR_EXPECTED( quadBasis.getValues(vals, badPoints1, OPERATOR_VALUE) );
}
{
    // exception #9 dimension 1 in the input point array must equal space dimension of the cell
    DynRankView ConstructWithLabel(badPoints2, 4, spaceDim+1);
    INTREPID2_TEST_ERROR_EXPECTED( quadBasis.getValues(vals, badPoints2, OPERATOR_VALUE) );
}
{
    // exception #10 output values must be of rank-3 for OPERATOR_VALUE in 2D
    DynRankView ConstructWithLabel(badVals1, 4, 3);
    INTREPID2_TEST_ERROR_EXPECTED( quadBasis.getValues(badVals1, quadNodes, OPERATOR_VALUE) );
}
{
    DynRankView ConstructWithLabel(badCurls1, 4, 3, 2);
    // exception #11 output values must be of rank-2 for OPERATOR_CURL
    INTREPID2_TEST_ERROR_EXPECTED( quadBasis.getValues(badCurls1, quadNodes, OPERATOR_CURL) );
}
{
    // exception #12 incorrect 0th dimension of output array (must equal number of basis functions)
//    FieldContainer<double> badVals2(quadBasis.getCardinality() + 1, quadNodes.dimension(0), quadBasis.getBaseCellTopology().getDimension());
    DynRankView ConstructWithLabel(badVals2, numFields + 1, numPoints, spaceDim);
    INTREPID2_TEST_ERROR_EXPECTED( quadBasis.getValues(badVals2, quadNodes, OPERATOR_VALUE) ) ;
}
{
    // exception #13 incorrect 1st  dimension of output array (must equal number of points)
    DynRankView ConstructWithLabels(badVals3, numFields, numPoints + 1, spaceDim);
    INTREPID2_TEST_ERROR_EXPECTED( quadBasis.getValues(badVals3, quadNodes, OPERATOR_VALUE) ) ;
}
{
    // exception #14: incorrect 2nd dimension of output array for VALUE (must equal the space dimension)
    DynRankView ConstructWithLabel(badVals4, numFields, numPoints, spaceDim - 1);
    INTREPID2_TEST_ERROR_EXPECTED( quadBasis.getValues(badVals4, quadNodes, OPERATOR_VALUE) ) ;
}
{
    // exception #15: D2 cannot be applied to HCURL functions 
    // resize vals to rank-3 container with dimensions (num. basis functions, num. points, arbitrary)
    vals = DynRankView(work.data(), numFields, numPoints, D2Cardin);
    INTREPID2_TEST_ERROR_EXPECTED( quadBasis.getValues(vals, quadNodes, OPERATOR_D2) );
}
#endif
    
  }
  catch (std::logic_error err) {
    *outStream << "UNEXPECTED ERROR !!! ----------------------------------------------------------\n";
    *outStream << err.what() << '\n';
    *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    errorFlag = -1000;
  };
  
  // Check if number of thrown exceptions matches the one we expect 
  // Note Teuchos throw number will not pick up exceptions 3-7 and therefore will not match.
  if (nthrow != ncatch) {
    errorFlag++;
    *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
  }
//#endif
  
/*
  *outStream \
    << "\n"
    << "===============================================================================\n"\
    << "| TEST 2: correctness of tag to enum and enum to tag lookups                  |\n"\
    << "===============================================================================\n";
  
  try{
    std::vector<std::vector<int> > allTags = quadBasis.getAllDofTags();
    
    // Loop over all tags, lookup the associated dof enumeration and then lookup the tag again
    for (unsigned i = 0; i < allTags.size(); i++) {
      int bfOrd  = quadBasis.getDofOrdinal(allTags[i][0], allTags[i][1], allTags[i][2]);
      
      std::vector<int> myTag = quadBasis.getDofTag(bfOrd);
       if( !( (myTag[0] == allTags[i][0]) &&
              (myTag[1] == allTags[i][1]) &&
              (myTag[2] == allTags[i][2]) &&
              (myTag[3] == allTags[i][3]) ) ) {
        errorFlag++;
        *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
        *outStream << " getDofOrdinal( {" 
          << allTags[i][0] << ", " 
          << allTags[i][1] << ", " 
          << allTags[i][2] << ", " 
          << allTags[i][3] << "}) = " << bfOrd <<" but \n";   
        *outStream << " getDofTag(" << bfOrd << ") = { "
          << myTag[0] << ", " 
          << myTag[1] << ", " 
          << myTag[2] << ", " 
          << myTag[3] << "}\n";        
      }
    }
    
    // Now do the same but loop over basis functions
    for( int bfOrd = 0; bfOrd < quadBasis.getCardinality(); bfOrd++) {
      std::vector<int> myTag  = quadBasis.getDofTag(bfOrd);
      int myBfOrd = quadBasis.getDofOrdinal(myTag[0], myTag[1], myTag[2]);
      if( bfOrd != myBfOrd) {
        errorFlag++;
        *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
        *outStream << " getDofTag(" << bfOrd << ") = { "
          << myTag[0] << ", " 
          << myTag[1] << ", " 
          << myTag[2] << ", " 
          << myTag[3] << "} but getDofOrdinal({" 
          << myTag[0] << ", " 
          << myTag[1] << ", " 
          << myTag[2] << ", " 
          << myTag[3] << "} ) = " << myBfOrd << "\n";
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
  
  // VALUE: Each row pair gives the 4x2 correct basis set values at an evaluation point: (P,F,D) layout
  double basisValues[] = {
    0.500000, 0, 0, 0, 0, 0, 0, -0.500000, 0.500000, 0, 0, 0.500000, 0, \
    0, 0, 0, 0, 0, 0, 0.500000, -0.500000, 0, 0, 0, 0, 0, 0, 0, \
    -0.500000, 0, 0, -0.500000, 0.250000, 0, 0, 0.250000, -0.250000, 0, \
    0, -0.250000, 0.375000, 0, 0, 0.250000, -0.125000, 0, 0, -0.250000, \
    0.125000, 0, 0, 0.250000, -0.375000, 0, 0, -0.250000, 0.250000, 0, 0, \
    0.125000, -0.250000, 0, 0, -0.375000, 0.250000, 0, 0, 0.375000, \
    -0.250000, 0, 0, -0.125000
  };
  
  // CURL: correct values in (F,P) format
  double basisCurls[] = {
    0.25, 0.25, 0.25, 0.25,
    0.25, 0.25, 0.25, 0.25,
    0.25, 0.25, 0.25, 0.25,
    0.25, 0.25, 0.25, 0.25,
    0.25, 0.25, 0.25, 0.25,
    0.25, 0.25, 0.25, 0.25,
    0.25, 0.25, 0.25, 0.25,
    0.25, 0.25, 0.25, 0.25,
    0.25, 0.25, 0.25, 0.25
  };

  
  try{
        
    // Dimensions for the output arrays:
    int numFields = quadBasis.getCardinality();
    int numPoints = quadNodes.dimension(0);
    int spaceDim  = quadBasis.getBaseCellTopology().getDimension();
    
    // Generic array for values and curls that will be properly sized before each call
    FieldContainer<double> vals;
    
    // Check VALUE of basis functions: resize vals to rank-3 container:
    vals.resize(numFields, numPoints, spaceDim);
    quadBasis.getValues(vals, quadNodes, OPERATOR_VALUE);
    for (int i = 0; i < numFields; i++) {
      for (int j = 0; j < numPoints; j++) {
        for (int k = 0; k < spaceDim; k++) {
           
          // compute offset for (P,F,D) data layout: indices are P->j, F->i, D->k
          int l = k + i * spaceDim + j * spaceDim * numFields;
           if (std::abs(vals(i,j,k) - basisValues[l]) > INTREPID_TOL) {
             errorFlag++;
             *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
        
             // Output the multi-index of the value where the error is: 
             *outStream << " At multi-index { ";
             *outStream << i << " ";*outStream << j << " ";*outStream << k << " ";
             *outStream << "}  computed value: " << vals(i,j,k)
               << " but reference value: " << basisValues[l] << "\n";
            }
         }
      }
    }

    // Check CURL of basis function: resize vals to rank-2 container
    vals.resize(numFields, numPoints);
    quadBasis.getValues(vals, quadNodes, OPERATOR_CURL);
    for (int i = 0; i < numFields; i++) {
      for (int j = 0; j < numPoints; j++) {
        int l =  i + j * numFields;
        if (std::abs(vals(i,j) - basisCurls[l]) > INTREPID_TOL) {
          errorFlag++;
          *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
          
          // Output the multi-index of the value where the error is:
          *outStream << " At multi-index { ";
          *outStream << i << " ";*outStream << j << " ";
          *outStream << "}  computed curl component: " << vals(i,j)
            << " but reference curl component: " << basisCurls[l] << "\n";
        }
      }
    }  
   }    
  
  // Catch unexpected errors
  catch (std::logic_error err) {
    *outStream << err.what() << "\n\n";
    errorFlag = -1000;
  };

  *outStream \
    << "\n"
    << "===============================================================================\n"\
    << "| TEST 4: correctness of DoF locations                                        |\n"\
    << "===============================================================================\n";

  try{
    Teuchos::RCP<Basis<double, FieldContainer<double> > > basis =
      Teuchos::rcp(new Basis_HCURL_QUAD_I1_FEM<double, FieldContainer<double> >);
    Teuchos::RCP<DofCoordsInterface<FieldContainer<double> > > coord_iface =
      Teuchos::rcp_dynamic_cast<DofCoordsInterface<FieldContainer<double> > >(basis);

    int spaceDim = 2;
    FieldContainer<double> cvals;
    FieldContainer<double> bvals(basis->getCardinality(), basis->getCardinality(),2); // last dimension is spatial dim

    // Check exceptions.
#ifdef HAVE_INTREPID2_DEBUG
    cvals.resize(1,2,3);
    INTREPID2_TEST_ERROR_EXPECTED( coord_iface->getDofCoords(cvals), throwCounter, nException );
    cvals.resize(3,2);
    INTREPID2_TEST_ERROR_EXPECTED( coord_iface->getDofCoords(cvals), throwCounter, nException );
    cvals.resize(4,3);
    INTREPID2_TEST_ERROR_EXPECTED( coord_iface->getDofCoords(cvals), throwCounter, nException );
#endif
    cvals.resize(4,spaceDim);
    INTREPID2_TEST_ERROR_EXPECTED( coord_iface->getDofCoords(cvals), throwCounter, nException ); nException--;
    // Check if number of thrown exceptions matches the one we expect
    if (throwCounter != nException) {
      errorFlag++;
      *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
    }

    // Check mathematical correctness
    FieldContainer<double> tangents(basis->getCardinality(),spaceDim); // tangents at each point basis point
    tangents(0,0) =  2.0; tangents(0,1) =  0.0;
    tangents(1,0) =  0.0; tangents(1,1) =  2.0;
    tangents(2,0) = -2.0; tangents(2,1) =  0.0;
    tangents(3,0) =  0.0; tangents(3,1) = -2.0;

    basis->getValues(bvals, cvals, OPERATOR_VALUE);
    char buffer[120];
    for (int i=0; i<bvals.dimension(0); i++) {
      for (int j=0; j<bvals.dimension(1); j++) {

        double tangent = 0.0;
        for(int d=0;d<spaceDim;d++) 
           tangent += bvals(i,j,d)*tangents(j,d);

        if ((i != j) && (std::abs(tangent - 0.0) > INTREPID_TOL)) {
          errorFlag++;
          sprintf(buffer, "\nValue of basis function %d at (%6.4e, %6.4e) is %6.4e but should be %6.4e!\n", i, cvals(i,0), cvals(i,1), tangent, 0.0);
          *outStream << buffer;
        }
        else if ((i == j) && (std::abs(tangent - 1.0) > INTREPID_TOL)) {
          errorFlag++;
          sprintf(buffer, "\nValue of basis function %d at (%6.4e, %6.4e) is %6.4e but should be %6.4e!\n", i, cvals(i,0), cvals(i,1), tangent, 1.0);
          *outStream << buffer;
        }
      }
    }

  }
  catch (std::logic_error err){
    *outStream << err.what() << "\n\n";
    errorFlag = -1000;
  };
*/
  
  if (errorFlag != 0)
    std::cout << "End Result: TEST FAILED\n";
  else
    std::cout << "End Result: TEST PASSED\n";
  
  // reset format state of std::cout
  std::cout.copyfmt(oldFormatState);
  return errorFlag;
}

} //end Test
} //end Intrepid2
