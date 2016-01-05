// @HEADER
// ************************************************************************
//
//                           Intrepid Package
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


/** \file   test_03.cpp
    \brief  Unit test of FieldContainer class.
    \author Created by P. Bochev and D. Ridzal.
*/

#include "Intrepid_FieldContainer.hpp"
#include "Intrepid_RealSpaceTools.hpp"
#include "Shards_Array.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_GlobalMPISession.hpp"


using namespace Intrepid;

SHARDS_ARRAY_DIM_TAG_SIMPLE_DECLARATION( Cell )
SHARDS_ARRAY_DIM_TAG_SIMPLE_IMPLEMENTATION( Cell )

SHARDS_ARRAY_DIM_TAG_SIMPLE_DECLARATION( Field )
SHARDS_ARRAY_DIM_TAG_SIMPLE_IMPLEMENTATION( Field )

SHARDS_ARRAY_DIM_TAG_SIMPLE_DECLARATION( Point )
SHARDS_ARRAY_DIM_TAG_SIMPLE_IMPLEMENTATION( Point )

SHARDS_ARRAY_DIM_TAG_SIMPLE_DECLARATION( Dim )
SHARDS_ARRAY_DIM_TAG_SIMPLE_IMPLEMENTATION( Dim )

#define INTREPID_TEST_COMMAND( S )                                                                                  \
{                                                                                                                   \
  try {                                                                                                             \
    S ;                                                                                                             \
  }                                                                                                                 \
  catch (std::logic_error err) {                                                                                    \
      *outStream << "Expected Error ----------------------------------------------------------------\n";            \
      *outStream << err.what() << '\n';                                                                             \
      *outStream << "-------------------------------------------------------------------------------" << "\n\n";    \
  };                                                                                                                \
}


int main(int argc, char *argv[]) {

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  
  // This little trick lets us print to cout only if a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  
  Teuchos::RCP<std::ostream> outStream;
  Teuchos::oblackholestream bhs; // outputs nothing
  
  if (iprint > 0)
    outStream = Teuchos::rcp(&std::cout, false);
  else
    outStream = Teuchos::rcp(&bhs, false);
  
  // Save the format state of the original cout .
  Teuchos::oblackholestream oldFormatState;
  oldFormatState.copyfmt(std::cout);
  
  *outStream  \
    << "===============================================================================\n" \
    << "|                                                                             |\n" \
    << "|                           Unit Test FieldContainer                          |\n" \
    << "|                                                                             |\n" \
    << "|     1) Testing usage of various constructors / wrappers                     |\n" \
    << "|     2) Testing usage of resize                                              |\n" \
    << "|                                                                             |\n" \
    << "|  Questions? Contact  Pavel Bochev (pbboche@sandia.gov) or                   |\n" \
    << "|                      Denis Ridzal (dridzal@sandia.gov).                     |\n" \
    << "|                                                                             |\n" \
    << "|  Intrepid's website: http://trilinos.sandia.gov/packages/intrepid           |\n" \
    << "|  Trilinos website:   http://trilinos.sandia.gov                             |\n" \
    << "|                                                                             |\n" \
    << "===============================================================================\n";
  

  // Set error flag.
  int errorFlag  = 0;

  double zero = INTREPID_TOL;

  try {
  
    // Dimensions array.
    Teuchos::Array<int> dimensions;
  
    *outStream << "\n" \
      << "===============================================================================\n"\
      << "| TEST 1: Constructors / Wrappers for a particular rank-4 container           |\n"\
      << "===============================================================================\n\n";

    { // start variable scope

      // Initialize dimensions for rank-4 multi-index value
      dimensions.resize(4);
      dimensions[0] = 5;
      dimensions[1] = 3;
      dimensions[2] = 2;
      dimensions[3] = 7;

      // Allocate data through a Teuchos::Array, initialized to 0.
      Teuchos::Array<double> data(dimensions[0]*dimensions[1]*dimensions[2]*dimensions[3]);
    
      // Create a FieldContainer using a deep copy via Teuchos::Array.
      FieldContainer<double> fc_array(dimensions, data);

      // modify the (1,1,1,1) entry
      fc_array(1,1,1,1) = 1.0;
      // verify that the data array has NOT changed
      if (std::abs(data[dimensions[1]*dimensions[2]*dimensions[3] + dimensions[2]*dimensions[3] + dimensions[3] + 1]) > zero) {
        *outStream << "\n\nError in constructor using Array (ArrayView) and deep copy.\n\n";
        errorFlag = -1000;
      }

      // test getData access function
      if (std::abs((fc_array.getData())[dimensions[1]*dimensions[2]*dimensions[3] +
                                        dimensions[2]*dimensions[3] + dimensions[3] + 1] - fc_array(1,1,1,1)) > zero) {
        *outStream << "\n\nError in getData() member of FieldContainer.\n\n";
        errorFlag = -1000;
      }

      // Create a FieldContainer using a deep copy via Teuchos::ArrayRCP.
      Teuchos::RCP< Teuchos::Array<double>  > rcp_to_data = rcpFromRef(data); // first create RCP
      Teuchos::ArrayRCP<double> arrayrcp = Teuchos::arcp(rcp_to_data);        // now create ArrayRCP
      FieldContainer<double> fc_arrayrcp_deep(dimensions, arrayrcp());        // IMPORTANT!!!: use '()' after arrayrcp
                                                                              //    for direct conversion to ArrayView
      // modify the (1,1,1,1) entry
      fc_arrayrcp_deep(1,1,1,1) = 1.0;
      // verify that the data array has NOT changed
      if (std::abs(data[dimensions[1]*dimensions[2]*dimensions[3] + dimensions[2]*dimensions[3] + dimensions[3] + 1]) > zero) {
        *outStream << "\n\nError in constructor using ArrayRCP (ArrayView) and deep copy.\n\n";
        errorFlag = -1000;
      }

      // Create a FieldContainer using a shallow copy via Teuchos::ArrayRCP.
      FieldContainer<double> fc_arrayrcp_shallow(dimensions, arrayrcp);
      // modify the (1,1,1,1) entry
      fc_arrayrcp_shallow(1,1,1,1) = 1.0;
      // verify that the data array HAS changed
      if (std::abs(data[dimensions[1]*dimensions[2]*dimensions[3] + dimensions[2]*dimensions[3] + dimensions[3] + 1] - fc_arrayrcp_shallow(1,1,1,1)) > zero) {
        *outStream << "\n\nError in constructor using ArrayRCP and shallow copy.\n\n";
        errorFlag = -1000;
      }

      // Create a FieldContainer using a deep copy via Scalar*.
      FieldContainer<double> fc_scalarptr_deep(dimensions, data.getRawPtr(), true);
      // modify the (1,1,1,1) entry
      fc_scalarptr_deep(1,1,1,1) = 2.0;
      // verify that the data array has NOT changed
      if (std::abs(data[dimensions[1]*dimensions[2]*dimensions[3] + dimensions[2]*dimensions[3] + dimensions[3] + 1] - 1.0) > zero) {
        *outStream << "\n\nError in constructor using Scalar* and deep copy.\n\n";
        errorFlag = -1000;
      }

      // Create a FieldContainer using a shallow copy via Scalar*.
      FieldContainer<double> fc_scalarptr_shallow(dimensions, data.getRawPtr());
      // modify the (1,1,1,1) entry
      fc_scalarptr_shallow(1,1,1,1) = 2.0;
      // verify that the data array HAS changed
      if (std::abs(data[dimensions[1]*dimensions[2]*dimensions[3] + dimensions[2]*dimensions[3] + dimensions[3] + 1] - fc_scalarptr_shallow(1,1,1,1)) > zero) {
        *outStream << "\n\nError in constructor using Scalar* and shallow copy.\n\n";
        errorFlag = -1000;
      }

      // Create a FieldContainer using a deep copy via shards::Array.
      shards::Array<double,shards::NaturalOrder,Cell,Field,Point,Dim> shards_array(data.getRawPtr(),dimensions[0],dimensions[1],dimensions[2],dimensions[3]);
      FieldContainer<double> fc_shards_deep(shards_array, true);
      // modify the (1,1,1,1) entry
      fc_shards_deep(1,1,1,1) = 3.0;
      // verify that the data array has NOT changed
      if (std::abs(data[dimensions[1]*dimensions[2]*dimensions[3] + dimensions[2]*dimensions[3] + dimensions[3] + 1] - 2.0) > zero) {
        *outStream << "\n\nError in constructor using shards::Array and deep copy.\n\n";
        errorFlag = -1000;
      }

      // Create a FieldContainer using a shallow copy via shards::Array.
      FieldContainer<double> fc_shards_shallow(shards_array);
      // modify the (1,1,1,1) entry
      fc_shards_shallow(1,1,1,1) = 3.0;
      // verify that the data array HAS changed
      if (std::abs(data[dimensions[1]*dimensions[2]*dimensions[3] + dimensions[2]*dimensions[3] + dimensions[3] + 1] - fc_shards_shallow(1,1,1,1)) > zero) {
        *outStream << "\n\nError in constructor using shards::Array and shallow copy.\n\n";
        errorFlag = -1000;
      }

    } // end variable scope


    *outStream << "\n" \
      << "===============================================================================\n"\
      << "| TEST 1 cont'd: Run through constructors / wrappers of various ranks         |\n"\
      << "===============================================================================\n\n";

    for (int rank=0; rank<10; rank++) {
      dimensions.resize(rank);
      int total_size = 1;
      if (rank == 0) {
        total_size = 0;
      }
      for (int dim=0; dim<rank; dim++) {
        dimensions[dim] = 2;
        total_size *= dimensions[dim];
      }

      // Allocate data through a Teuchos::Array, initialized to 0.
      Teuchos::Array<double> data(total_size);

      // Create a FieldContainer using a deep copy via Teuchos::Array.
      FieldContainer<double> fc_array(dimensions, data);
      // Create a FieldContainer using a deep copy via Teuchos::ArrayRCP.
      Teuchos::RCP< Teuchos::Array<double>  > rcp_to_data = rcpFromRef(data); // first create RCP
      Teuchos::ArrayRCP<double> arrayrcp = Teuchos::arcp(rcp_to_data);        // now create ArrayRCP
      FieldContainer<double> fc_arrayrcp_deep(dimensions, arrayrcp());        // IMPORTANT!!!: use '()' after arrayrcp
                                                                              //    for direct conversion to ArrayView
      // Create a FieldContainer using a shallow copy via Teuchos::ArrayRCP.
      FieldContainer<double> fc_arrayrcp_shallow(dimensions, arrayrcp);
      // Create a FieldContainer using a deep copy via Scalar*.
      FieldContainer<double> fc_scalarptr_deep(dimensions, data.getRawPtr(), true);
      // Create a FieldContainer using a shallow copy via Scalar*.
      FieldContainer<double> fc_scalarptr_shallow(dimensions, data.getRawPtr());
    }

    { // start variable scope
      // Create FieldContainers using a deep copy via shards::Array.
      Teuchos::Array<double> data(2*2*2*2*2*2);
      shards::Array<double,shards::NaturalOrder,Cell> shards_array_c(&data[0],2);
      shards::Array<double,shards::NaturalOrder,Cell,Field> shards_array_cf(&data[0],2,2);
      shards::Array<double,shards::NaturalOrder,Cell,Field,Point> shards_array_cfp(&data[0],2,2,2);
      shards::Array<double,shards::NaturalOrder,Cell,Field,Point,Dim> shards_array_cfpd(&data[0],2,2,2,2);
      shards::Array<double,shards::NaturalOrder,Cell,Field,Point,Dim,Dim> shards_array_cfpdd(&data[0],2,2,2,2,2);
      shards::Array<double,shards::NaturalOrder,Cell,Field,Point,Dim,Dim,Dim> shards_array_cfpddd(&data[0],2,2,2,2,2,2);
      FieldContainer<double> fc_shards_c_deep(shards_array_c, true);
      FieldContainer<double> fc_shards_cf_deep(shards_array_cf, true);
      FieldContainer<double> fc_shards_cfp_deep(shards_array_cfp, true);
      FieldContainer<double> fc_shards_cfpd_deep(shards_array_cfpd, true);
      FieldContainer<double> fc_shards_cfpdd_deep(shards_array_cfpdd, true);
      FieldContainer<double> fc_shards_cfpddd_deep(shards_array_cfpddd, true);
      // Create FieldContainers using a shallow copy via shards::Array.
      FieldContainer<double> fc_shards_c_shallow(shards_array_c);
      FieldContainer<double> fc_shards_cf_shallow(shards_array_cf);
      FieldContainer<double> fc_shards_cfp_shallow(shards_array_cfp);
      FieldContainer<double> fc_shards_cfpd_shallow(shards_array_cfpd);
      FieldContainer<double> fc_shards_cfpdd_shallow(shards_array_cfpdd);
      FieldContainer<double> fc_shards_cfpddd_shallow(shards_array_cfpddd);
    } // end variable scope


    *outStream << "\n" \
      << "===============================================================================\n"\
      << "| TEST 2: Usage of resize                                                     |\n"\
      << "===============================================================================\n\n";

    { // start variable scope
      // Initialize dimensions for rank-4 multi-index value
      dimensions.resize(5);
      dimensions[0] = 5;
      dimensions[1] = 3;
      dimensions[2] = 2;
      dimensions[3] = 7;
      dimensions[4] = 2;

      // temporary ints
      int d0, d1, d2, d3;

      // Allocate data through a Teuchos::Array, initialized to 0.
      Teuchos::Array<double> data(dimensions[0]*dimensions[1]*dimensions[2]*dimensions[3]*dimensions[4]);

      // Create a FieldContainer using a deep copy via Teuchos::Array.
      FieldContainer<double> fc_array(dimensions, data);
      // modify the (1,1,1,1) entry
      double mod_entry    = 1.0;
      fc_array(1,1,1,1,1) = mod_entry;
      int enumeration = fc_array.getEnumeration(1,1,1,1,1);

      // Resize into a (dimensions[0]*dimensions[1]) x (dimensions[2]) x (dimensions[3]) x (dimensions[4]) rank-4 array.
      // This is effectively a reshape, the data should not be touched.
      fc_array.resize(dimensions[0]*dimensions[1], dimensions[2], dimensions[3], dimensions[4]);
      // verify that the data array has NOT changed
      fc_array.getMultiIndex(d0,d1,d2,d3, enumeration);
      if (std::abs(fc_array(d0,d1,d2,d3) - mod_entry) > zero) {
        *outStream << "\n\nError in resize.\n\n";
        errorFlag = -1000;
      }

      // Resize into a (dimensions[0]*dimensions[1]*dimensions[2]) x (dimensions[3]) x (dimensions[4]) rank-3 array.
      // This is effectively a reshape, the data should not be touched.
      fc_array.resize(dimensions[0]*dimensions[1]*dimensions[2], dimensions[3], dimensions[4]);
      // verify that the data array has NOT changed
      fc_array.getMultiIndex(d0,d1,d2, enumeration);
      if (std::abs(fc_array(d0,d1,d2) - mod_entry) > zero) {
        *outStream << "\n\nError in resize.\n\n";
        errorFlag = -1000;
      }

      // Resize into a (dimensions[0]*dimensions[1]*dimensions[2]*dimensions[3]) x (dimensions[4]) rank-2 array.
      // This is effectively a reshape, the data should not be touched.
      fc_array.resize(dimensions[0]*dimensions[1]*dimensions[2]*dimensions[3], dimensions[4]);
      // verify that the data array has NOT changed
      fc_array.getMultiIndex(d0,d1, enumeration);
      if (std::abs(fc_array(d0,d1) - mod_entry) > zero) {
        *outStream << "\n\nError in resize.\n\n";
        errorFlag = -1000;
      }

      // Resize into a (dimensions[0]*dimensions[1]*dimensions[2]*dimensions[3]*dimensions[4]) rank-1 array.
      // This is effectively a reshape, the data should not be touched.
      fc_array.resize(dimensions[0]*dimensions[1]*dimensions[2]*dimensions[3]*dimensions[4]);
      // verify that the data array has NOT changed
      fc_array.getMultiIndex(d0, enumeration);
      if (std::abs(fc_array(d0) - mod_entry) > zero) {
        *outStream << "\n\nError in resize.\n\n";
        errorFlag = -1000;
      }

      // More advanced example allocating new memory, with shards::Array.
      data.assign(dimensions[0]*dimensions[1]*dimensions[2]*dimensions[3]*dimensions[4], 3.0);
      shards::Array<double,shards::NaturalOrder,Cell,Field,Point,Dim,Dim>
        shards_array(data.getRawPtr(),dimensions[0],dimensions[1],dimensions[2],dimensions[3],dimensions[4]);
      // Create a FieldContainer using a shallow copy via shards::Array.
      FieldContainer<double> fc_shards_shallow(shards_array);
      fc_shards_shallow.resize(4,4,4,4,4);  // keep old entries + allocate new memory and fill with zeros
      fc_shards_shallow.resize(4*4*4*4*4);  // reshape into rank-1 vector
      if (std::abs(RealSpaceTools<double>::vectorNorm(data.getRawPtr(), fc_array.size(), NORM_ONE) -
                   RealSpaceTools<double>::vectorNorm(fc_shards_shallow, NORM_ONE)) > zero) {
        *outStream << "\n\nError in resize.\n\n";
        errorFlag = -1000;
      }


    } // end variable scope


  } // outer try block
  catch (std::logic_error err) {
    *outStream << "UNEXPECTED ERROR !!! ----------------------------------------------------------\n";
    *outStream  << err.what() << "\n";
    *outStream << "-------------------------------------------------------------------------------" << "\n\n";
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
