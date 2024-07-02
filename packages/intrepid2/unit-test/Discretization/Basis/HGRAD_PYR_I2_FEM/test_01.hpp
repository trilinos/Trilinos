// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file test_01.cpp
    \brief  Unit tests for the Intrepid2::HGRAD_PYR_I2_FEM class.
    \author Created by P. Bochev, D. Ridzal, K. Peterson, M. Perego and Kyungjoo Kim.
*/

#include "Intrepid2_config.h"

#ifdef HAVE_INTREPID2_DEBUG
#define INTREPID2_TEST_FOR_DEBUG_ABORT_OVERRIDE_TO_CONTINUE
#endif

#include "Intrepid2_HGRAD_PYR_I2_FEM.hpp"

#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"

namespace Intrepid2 {

  namespace Test {
#define INTREPID2_TEST_ERROR_EXPECTED( S )                              \
    try {                                                               \
      ++nthrow;                                                         \
      S ;                                                               \
    }                                                                   \
    catch (std::exception &err) {                                        \
      ++ncatch;                                                         \
      *outStream << "Expected Error ----------------------------------------------------------------\n"; \
      *outStream << err.what() << '\n';                                 \
      *outStream << "-------------------------------------------------------------------------------" << "\n\n"; \
    }

    template<typename ValueType, typename DeviceType>
    int HGRAD_PYR_I2_FEM_Test01(const bool verbose) {

      Teuchos::RCP<std::ostream> outStream;
      Teuchos::oblackholestream bhs; // outputs nothing
      
      if (verbose)
        outStream = Teuchos::rcp(&std::cout, false);
      else
        outStream = Teuchos::rcp(&bhs,       false);

      Teuchos::oblackholestream oldFormatState;
      oldFormatState.copyfmt(std::cout);

      using DeviceSpaceType = typename DeviceType::execution_space;
      typedef typename
        Kokkos::DefaultHostExecutionSpace HostSpaceType ;

      *outStream << "DeviceSpace::  "; DeviceSpaceType().print_configuration(*outStream, false);
      *outStream << "HostSpace::    ";   HostSpaceType().print_configuration(*outStream, false);

      *outStream
        << "===============================================================================\n"
        << "|                                                                             |\n"
        << "|              Unit Test (Basis_HGRAD_PYR_I2_FEM)                 |\n"
        << "|                                                                             |\n"
        << "|     1) Conversion of Dof tags into Dof ordinals and back                    |\n"
        << "|     2) Basis values for VALUE, GRAD, and Dk operators                       |\n"
        << "|                                                                             |\n"
        << "|  Questions? Contact  Pavel Bochev  (pbboche@sandia.gov),                    |\n"
        << "|                      Denis Ridzal  (dridzal@sandia.gov),                    |\n"
        << "|                      Kara Peterson (kjpeter@sandia.gov).                    |\n"
        << "|                      Kyungjoo Kim  (kyukim@sandia.gov).                     |\n"
        << "|                                                                             |\n"
        << "|  Intrepid's website: http://trilinos.sandia.gov/packages/intrepid           |\n"
        << "|  Trilinos website:   http://trilinos.sandia.gov                             |\n"
        << "|                                                                             |\n"
        << "===============================================================================\n";

      typedef Kokkos::DynRankView<ValueType,DeviceType> DynRankView;
      typedef Kokkos::DynRankView<ValueType,HostSpaceType>   DynRankViewHost;
#define ConstructWithLabel(obj, ...) obj(#obj, __VA_ARGS__)

      const ValueType tol = tolerence();
      int errorFlag = 0;

      // for virtual function, value and point types are declared in the class
      typedef ValueType outputValueType;
      typedef ValueType pointValueType;
      Basis_HGRAD_PYR_I2_FEM<DeviceType,outputValueType,pointValueType> pyrBasis;

     *outStream
       << "\n"
       << "===============================================================================\n"
       << "| TEST 1: constructors and exceptions                                         |\n"
       << "===============================================================================\n";

     try {
       ordinal_type nthrow = 0, ncatch = 0;
#ifdef HAVE_INTREPID2_DEBUG
       // Define array containing the 4 vertices of the reference PYR and its center.
       DynRankView ConstructWithLabel(pyrNodes, 10, 3);

       // Generic array for the output values; needs to be properly resized depending on the operator type
       const auto numFields = pyrBasis.getCardinality();
       const auto numPoints = pyrNodes.extent(0);
       const auto spaceDim  = pyrBasis.getBaseCellTopology().getDimension();

       DynRankView vals("vals", numFields, numPoints);
       DynRankView vals_vec("vals", numFields, numPoints, spaceDim);

       {
         // exception #1: CURL cannot be applied to scalar functions
         INTREPID2_TEST_ERROR_EXPECTED( pyrBasis.getValues(vals_vec, pyrNodes, OPERATOR_CURL) );
         // exception #2: DIV cannot be applied to scalar functions
         INTREPID2_TEST_ERROR_EXPECTED( pyrBasis.getValues(vals_vec, pyrNodes, OPERATOR_DIV) );
       }

       // Exceptions 3-7: all bf tags/bf Ids below are wrong and should cause getDofOrdinal() and
       // getDofTag() to access invalid array elements thereby causing bounds check exception
       {
         // exception #3
         INTREPID2_TEST_ERROR_EXPECTED( pyrBasis.getDofOrdinal(3,0,0) );
         // exception #4
         INTREPID2_TEST_ERROR_EXPECTED( pyrBasis.getDofOrdinal(1,1,1) );
         // exception #5
         INTREPID2_TEST_ERROR_EXPECTED( pyrBasis.getDofOrdinal(0,6,0) );
         // exception #6
         INTREPID2_TEST_ERROR_EXPECTED( pyrBasis.getDofTag(numFields) );
         // exception #7
         INTREPID2_TEST_ERROR_EXPECTED( pyrBasis.getDofTag(-1) );
       }

       // Exceptions 8-18 test exception handling with incorrectly dimensioned input/output arrays
       {
         // exception #8: input points array must be of rank-2
         DynRankView ConstructWithLabel(badPoints1, 4, 5, 3);
         INTREPID2_TEST_ERROR_EXPECTED( pyrBasis.getValues(vals, badPoints1, OPERATOR_VALUE) );
       }
       {
         // exception #9 dimension 1 in the input point array must equal space dimension of the cell
         DynRankView ConstructWithLabel(badPoints2, 4, 2);
         INTREPID2_TEST_ERROR_EXPECTED( pyrBasis.getValues(vals, badPoints2, OPERATOR_VALUE) );
       }
       {
         // exception #10 output values must be of rank-2 for OPERATOR_VALUE
         DynRankView ConstructWithLabel(badVals1, 4, 3, 1);
         INTREPID2_TEST_ERROR_EXPECTED( pyrBasis.getValues(badVals1, pyrNodes, OPERATOR_VALUE) );
       }
       {
         // exception #11 output values must be of rank-3 for OPERATOR_GRAD
         DynRankView ConstructWithLabel(badVals2, 4, 3);
         INTREPID2_TEST_ERROR_EXPECTED( pyrBasis.getValues(badVals2, pyrNodes, OPERATOR_GRAD) );

         // exception #12 output values must be of rank-3 for OPERATOR_D1
         INTREPID2_TEST_ERROR_EXPECTED( pyrBasis.getValues(badVals2, pyrNodes, OPERATOR_D1) );

         // exception #13 output values must be of rank-3 for OPERATOR_D2
         INTREPID2_TEST_ERROR_EXPECTED( pyrBasis.getValues(badVals2, pyrNodes, OPERATOR_D2) );
       }
       {
         // exception #14 incorrect 0th dimension of output array (must equal number of basis functions)
         DynRankView ConstructWithLabel(badVals3, numFields + 1, numPoints);
         INTREPID2_TEST_ERROR_EXPECTED( pyrBasis.getValues(badVals3, pyrNodes, OPERATOR_VALUE) );
       }
       {
         // exception #15 incorrect 1st dimension of output array (must equal number of points)
         DynRankView ConstructWithLabel(badVals4, numFields, numPoints + 1);
         INTREPID2_TEST_ERROR_EXPECTED( pyrBasis.getValues(badVals4, pyrNodes, OPERATOR_VALUE) );
       }
       {
         // exception #16: incorrect 2nd dimension of output array (must equal the space dimension)
         DynRankView ConstructWithLabel(badVals5, numFields, numPoints, 4);
         INTREPID2_TEST_ERROR_EXPECTED( pyrBasis.getValues(badVals5, pyrNodes, OPERATOR_GRAD) );
       }
       {
         // exception #17: incorrect 2nd dimension of output array (must equal D2 cardinality in 3D)
         DynRankView ConstructWithLabel(badVals6, numFields, numPoints, 40);
         INTREPID2_TEST_ERROR_EXPECTED( pyrBasis.getValues(badVals6, pyrNodes, OPERATOR_D2) );
         //exception #18: incorrect 2nd dimension of output array (must equal D3 cardinality in 3D)
         INTREPID2_TEST_ERROR_EXPECTED( pyrBasis.getValues(badVals6, pyrNodes, OPERATOR_D3) );
       }
#endif
       if (nthrow != ncatch) {
         errorFlag++;
         *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
         *outStream << "# of catch ("<< ncatch << ") is different from # of throw (" << nthrow << ")\n";
       }
     } catch (std::exception &err) {
       *outStream << "UNEXPECTED ERROR !!! ----------------------------------------------------------\n";
       *outStream << err.what() << '\n';
       *outStream << "-------------------------------------------------------------------------------" << "\n\n";
       errorFlag = -1000;
     }

     *outStream                                 \
       << "\n"
       << "===============================================================================\n" \
       << "| TEST 2: correctness of tag to enum and enum to tag lookups                  |\n" \
       << "===============================================================================\n";

     try {
       const ordinal_type numFields = pyrBasis.getCardinality();
       const auto allTags = pyrBasis.getAllDofTags();

       // Loop over all tags, lookup the associated dof enumeration and then lookup the tag again
       const ordinal_type dofTagSize = allTags.extent(0);
       for (ordinal_type i=0;i<dofTagSize;++i) {
         const auto bfOrd = pyrBasis.getDofOrdinal(allTags(i,0), allTags(i,1), allTags(i,2));

         const auto myTag = pyrBasis.getDofTag(bfOrd);
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
       for (ordinal_type bfOrd=0;bfOrd<numFields;++bfOrd) {
         const auto myTag = pyrBasis.getDofTag(bfOrd);

         const auto myBfOrd = pyrBasis.getDofOrdinal(myTag(0), myTag(1), myTag(2));
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
     } catch (std::logic_error &err) {
       *outStream << err.what() << "\n\n";
       errorFlag = -1000;
     }

     *outStream                                 \
       << "\n"
       << "===============================================================================\n" \
       << "| TEST 3: correctness of basis function values                                |\n" \
       << "===============================================================================\n";

     outStream -> precision(20);
      // Basis values are stored in (F,P) format in a data file. Read file and do the test     
      std::vector<ValueType> basisValues;           // Flat array for the gradient values.
      { 
        std::ifstream dataFile;
        dataFile.open("./testdata/PYR_I2_Vals.dat");
        INTREPID2_TEST_FOR_EXCEPTION( !dataFile.good(), std::logic_error,
                                      ">>> ERROR (HGRAD_PYR_I2_Serendipity/test01): could not open values data file, test aborted.");
        while (!dataFile.eof() ){
          ValueType temp;
          std::string line;                            // string for one line of input file
          std::getline(dataFile, line);           // get next line from file
          std::stringstream data_line(line);           // convert to stringstream
          while (data_line >> temp)               // extract value from line
            basisValues.push_back(temp);           // push into vector
        }
      }

      // GRAD and D1 values are stored in (F,P,D) format in a data file. Read file and do the test     
      std::vector<ValueType> basisGrads;           // Flat array for the gradient values.
      { 
        std::ifstream dataFile;
        dataFile.open("./testdata/PYR_I2_GradVals.dat");
        INTREPID2_TEST_FOR_EXCEPTION( !dataFile.good(), std::logic_error,
                                      ">>> ERROR (HGRAD_PYR_I2_Serendipity/test01): could not open GRAD values data file, test aborted.");
        while (!dataFile.eof() ){
          ValueType temp;
          std::string line;                            // string for one line of input file
          std::getline(dataFile, line);           // get next line from file
          std::stringstream data_line(line);           // convert to stringstream
          while (data_line >> temp)               // extract value from line
            basisGrads.push_back(temp);           // push into vector
        }
      }
  
      //D2: flat array with the values of D2 applied to basis functions. Multi-index is (F,P,D2cardinality)
      std::vector<ValueType> basisD2;
      { 
        std::ifstream dataFile;

        dataFile.open("./testdata/PYR_I2_D2Vals.dat");        
        INTREPID2_TEST_FOR_EXCEPTION( !dataFile.good(), std::logic_error,
                                      ">>> ERROR (HGRAD_PYR_I2_Serendipity/test01): could not open D2 values data file, test aborted.");
        while (!dataFile.eof() ){
          ValueType temp;
          std::string line;                            // string for one line of input file
          std::getline(dataFile, line);           // get next line from file
          std::stringstream data_line(line);           // convert to stringstream
          while (data_line >> temp)               // extract value from line
            basisD2.push_back(temp);           // push into vector
        }
      }

     try {
      DynRankViewHost ConstructWithLabel(pyrNodesHost, 18, 3);
      pyrNodesHost(0,0)  = -1.0;  pyrNodesHost(0,1)  = -1.0;  pyrNodesHost(0,2)  =  0;
      pyrNodesHost(1,0)  =  1.0;  pyrNodesHost(1,1)  = -1.0;  pyrNodesHost(1,2)  =  0;
      pyrNodesHost(2,0)  =  1.0;  pyrNodesHost(2,1)  =  1.0;  pyrNodesHost(2,2)  =  0;
      pyrNodesHost(3,0)  = -1.0;  pyrNodesHost(3,1)  =  1.0;  pyrNodesHost(3,2)  =  0;
      pyrNodesHost(4,0)  =  0.0;  pyrNodesHost(4,1)  =  0.0;  pyrNodesHost(4,2)  =  1.0;
      pyrNodesHost(5,0)  =  0.0;  pyrNodesHost(5,1)  = -1.0;  pyrNodesHost(5,2)  =  0.0;
      pyrNodesHost(6,0)  =  1.0;  pyrNodesHost(6,1)  =  0.0;  pyrNodesHost(6,2)  =  0.0;
      pyrNodesHost(7,0)  =  0.0;  pyrNodesHost(7,1)  =  1.0;  pyrNodesHost(7,2)  =  0.0;
      pyrNodesHost(8,0)  = -1.0;  pyrNodesHost(8,1)  =  0.0;  pyrNodesHost(8,2)  =  0.0;
      pyrNodesHost(9,0)  = -0.5;  pyrNodesHost(9,1)  = -0.5;  pyrNodesHost(9,2)  =  0.5;
      pyrNodesHost(10,0) =  0.5;  pyrNodesHost(10,1) = -0.5;  pyrNodesHost(10,2) =  0.5;
      pyrNodesHost(11,0) =  0.5;  pyrNodesHost(11,1) =  0.5;  pyrNodesHost(11,2) =  0.5;
      pyrNodesHost(12,0) = -0.5;  pyrNodesHost(12,1) =  0.5;  pyrNodesHost(12,2) =  0.5;

      pyrNodesHost(13,0) =  0.25; pyrNodesHost(13,1) =  0.5;  pyrNodesHost(13,2) = 0.2;
      pyrNodesHost(14,0) = -0.7 ; pyrNodesHost(14,1) =  0.3;  pyrNodesHost(14,2) = 0.3;
      pyrNodesHost(15,0) =  0.;   pyrNodesHost(15,1) = -0.05; pyrNodesHost(15,2) = 0.95;
      pyrNodesHost(16,0) = -0.15; pyrNodesHost(16,1) = -0.2;  pyrNodesHost(16,2) = 0.75;
      pyrNodesHost(17,0) = -0.4;  pyrNodesHost(17,1) =  0.9;  pyrNodesHost(17,2) = 0.0;


       auto pyrNodes = Kokkos::create_mirror_view(typename DeviceType::memory_space(), pyrNodesHost);
       Kokkos::deep_copy(pyrNodes, pyrNodesHost);

       // Dimensions for the output arrays:
       const ordinal_type numFields = pyrBasis.getCardinality();
       const ordinal_type numPoints = pyrNodes.extent(0);
       const ordinal_type spaceDim  = pyrBasis.getBaseCellTopology().getDimension();
       const ordinal_type D2Cardin  = getDkCardinality(OPERATOR_D2, spaceDim);

       // Check VALUE of basis functions: resize vals to rank-2 container:
       {
         DynRankView vals = DynRankView("vals", numFields, numPoints);
         pyrBasis.getValues(vals, pyrNodes, OPERATOR_VALUE);
         auto vals_host = Kokkos::create_mirror_view(typename HostSpaceType::memory_space(), vals);
         Kokkos::deep_copy(vals_host, vals);
         for (ordinal_type i=0;i<numFields;++i) {
           for (ordinal_type j=0;j<numPoints;++j) {
             const ordinal_type l =  i*numPoints + j;
             if (std::abs(vals_host(i,j) - basisValues[l]) > tol) {
               errorFlag++;
               *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";

               // Output the multi-index of the value where the error is:
               *outStream << " At multi-index { ";
               *outStream << i << " ";*outStream << j << " ";
               *outStream << "}  computed value: " << vals_host(i,j)
                          << " but reference value: " << basisValues[l] << "\n";
             }
           }
         }
       }

       // Check GRAD of basis function: resize vals to rank-3 container
       {
         DynRankView vals = DynRankView("vals", numFields, numPoints, spaceDim);
         pyrBasis.getValues(vals, pyrNodes, OPERATOR_GRAD);
         auto vals_host = Kokkos::create_mirror_view(typename HostSpaceType::memory_space(), vals);
         Kokkos::deep_copy(vals_host, vals);
         for (ordinal_type i=0;i<numFields;++i) {
           for (ordinal_type j=0;j<numPoints;++j) {
             for (ordinal_type k=0;k<spaceDim;++k) {
               const ordinal_type l = i + j * numFields + k * numFields * numPoints;
               if (std::abs(vals_host(i,j,k) - basisGrads[l]) > tol) {
                 errorFlag++;
                 *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";

                 // Output the multi-index of the value where the error is:
                 *outStream << " At multi-index { ";
                 *outStream << i << " ";*outStream << j << " ";*outStream << k << " ";
                 *outStream << "}  computed grad component: " << vals_host(i,j,k)
                            << " but reference grad component: " << basisGrads[l] << "\n";
               }
             }
           }
         }
       }

       // Check D1 of basis function
       {
         DynRankView vals = DynRankView("vals", numFields, numPoints, spaceDim);
         pyrBasis.getValues(vals, pyrNodes, OPERATOR_D1);
         auto vals_host = Kokkos::create_mirror_view(typename HostSpaceType::memory_space(), vals);
         Kokkos::deep_copy(vals_host, vals);
         for (ordinal_type i=0;i<numFields;++i) {
           for (ordinal_type j=0;j<numPoints;++j) {
             for (ordinal_type k=0;k<spaceDim;++k) {
               const ordinal_type l = i + j * numFields + k * numFields * numPoints;
               if (std::abs(vals_host(i,j,k) - basisGrads[l]) > tol) {
                 errorFlag++;
                 *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";

                 // Output the multi-index of the value where the error is:
                 *outStream << " At multi-index { ";
                 *outStream << i << " ";*outStream << j << " ";*outStream << k << " ";
                 *outStream << "}  computed D1 component: " << vals_host(i,j,k)
                            << " but reference D1 component: " << basisGrads[l] << "\n";
               }
             }
           }
         }
       }

       // Check D2 of basis function
       {
         DynRankView vals = DynRankView("vals", numFields, numPoints, D2Cardin);
         pyrBasis.getValues(vals, pyrNodes, OPERATOR_D2);
         auto vals_host = Kokkos::create_mirror_view(typename HostSpaceType::memory_space(), vals);
         Kokkos::deep_copy(vals_host, vals);
         for (ordinal_type i=0;i<numFields;++i) {
           for (ordinal_type j=0;j<numPoints;++j) {
             // derivatives are singular when z = 1; using the same eps, it can be comparable
             //if (j == 4) continue; 
             for (ordinal_type k=0;k<D2Cardin;++k) {
               const ordinal_type l = i + j * numFields + k * numFields * numPoints;
               if (std::abs(vals_host(i,j,k) - basisD2[l]) > tol) {
                 errorFlag++;
                 *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";

                 // Output the multi-index of the value where the error is:
                 *outStream << " At multi-index { ";
                 *outStream << i << " ";*outStream << j << " ";*outStream << k << " ";
                 *outStream << "}  computed D2 component: " << vals_host(i,j,k)
                            << " but reference D2 component: " << basisD2[l] << "\n";
               }
             }
           }
         }
       }
     } catch (std::logic_error &err) {
       *outStream << err.what() << "\n\n";
       errorFlag = -1000;
     }

     *outStream
     << "\n"
     << "===============================================================================\n"
     << "| TEST 4: Function Space is Correct                                           |\n"
     << "===============================================================================\n";
     
     try {
       const EFunctionSpace fs = pyrBasis.getFunctionSpace();
       
       if (fs != FUNCTION_SPACE_HGRAD)
       {
         *outStream << std::setw(70) << "------------- TEST FAILURE! -------------" << "\n";
         
         // Output the multi-index of the value where the error is:
         *outStream << " Expected a function space of FUNCTION_SPACE_HGRAD (enum value " << FUNCTION_SPACE_HGRAD << "),";
         *outStream << " but got " << fs << "\n";
         if (fs == FUNCTION_SPACE_MAX)
         {
           *outStream << "Note that this matches the default value defined by superclass, FUNCTION_SPACE_MAX.  Likely the subclass has failed to set the superclass functionSpace_ field.\n";
         }
         errorFlag++;
       }
     } catch (std::logic_error &err){
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

     return errorFlag;
    }
  }
}






