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

/** \file test_01.cpp
\brief  Unit test for the PointTools class.
\author Created by R. Kirby.
*/

#include "Intrepid_PointTools.hpp"
#include "Intrepid_FieldContainer.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Shards_CellTopology.hpp"

using namespace std;
using namespace Intrepid;


#define INTREPID_TEST_COMMAND( S , throwCounter, nException )                                                              \
{                                                                                                                          \
  ++nException;                                                                                                            \
  try {                                                                                                                    \
    S ;                                                                                                                    \
  }                                                                                                                        \
  catch (std::logic_error err) {                                                                                           \
      ++throwCounter;                                                                                                      \
      *outStream << "Expected Error " << nException << " -------------------------------------------------------------\n"; \
      *outStream << err.what() << '\n';                                                                                    \
      *outStream << "-------------------------------------------------------------------------------" << "\n\n";           \
  };                                                                                                                       \
}


int main(int argc, char *argv[]) {

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  // This little trick lets us print to std::cout only if
  // a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
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
  << "|                       Unit Test (PointTools)                                |\n" \
  << "|                                                                             |\n" \
  << "|     1) Construction of equispaced and warped lattices on simplices          |\n" \
  << "|                                                                             |\n" \
  << "|  Questions? Contact  Pavel Bochev (pbboche@sandia.gov) or                   |\n" \
  << "|                      Denis Ridzal (dridzal@sandia.gov) or                   |\n" \
  << "|                      Robert Kirby (robert.c.kirby@ttu.edu)                  |\n" \
  << "|                                                                             |\n" \
  << "|  Intrepid's website: http://trilinos.sandia.gov/packages/intrepid           |\n" \
  << "|  Trilinos website:   http://trilinos.sandia.gov                             |\n" \
  << "|                                                                             |\n" \
  << "===============================================================================\n";



  int errorFlag = 0;

  shards::CellTopology myLine_2( shards::getCellTopologyData< shards::Line<2> >() );
  shards::CellTopology myTri_3( shards::getCellTopologyData< shards::Triangle<3> >() );
  shards::CellTopology myTet_4( shards::getCellTopologyData< shards::Tetrahedron<4> >() );


  *outStream \
  << "\n"
  << "===============================================================================\n"\
  << "| TEST 1: size of lattices                                                   |\n"\
  << "===============================================================================\n";

  try{
    // first try the lattices with offset = 0.  This is a spot-check

    if (PointTools::getLatticeSize( myLine_2 , 4 , 0 ) != 5) {
      errorFlag++;
      *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
      *outStream << " size of 4th order lattice on a line with no offset: " << PointTools::getLatticeSize( myLine_2 , 4 , 0 )  << "\n";
      *outStream << " should be 5\n";
    }


    if (PointTools::getLatticeSize( myTri_3 , 3 , 0 ) != 10) {
      errorFlag++;
      *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
      *outStream << " size of 3rd order lattice on a line with no offset: " << PointTools::getLatticeSize( myTri_3 , 3 , 0 )  << "\n";
      *outStream << " should be 10\n";    
    }


    if (PointTools::getLatticeSize( myTet_4 , 3 , 0 ) != 20) {
      errorFlag++;
      *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
      *outStream << " size of 3rd order lattice on a tet with no offset: " << PointTools::getLatticeSize( myTet_4 , 3 , 0 )  << "\n";
      *outStream << " should be 20\n"; 
    }

                        
    // check with the offset = 1
    if (PointTools::getLatticeSize( myLine_2 , 3 , 1 ) != 2) {
      errorFlag++;
      *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
      *outStream << " size of 3rd order lattice on a line with 1 offset: " << PointTools::getLatticeSize( myLine_2 , 3 , 1 )  << "\n";
      *outStream << " should be 2\n";       
    }

    if (PointTools::getLatticeSize( myTri_3 , 4 , 1 ) != 3) {
      errorFlag++;
      *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
      *outStream << " size of 4th order lattice on a triangle with 1 offset: " << PointTools::getLatticeSize( myTri_3 , 4 , 1 )  << "\n";
      *outStream << " should be 3\n";           
    }

    if (PointTools::getLatticeSize( myTet_4 , 5 , 1 ) != 4) {
      errorFlag++;
      *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
      *outStream << " size of 5th order lattice on a tetrahedron with 1 offset: " << PointTools::getLatticeSize( myTet_4 , 5 , 1 )  << "\n";
      *outStream << " should be 4\n";   
    }

  }
  catch (std::exception &err) {
    *outStream << "UNEXPECTED ERROR !!! ----------------------------------------------------------\n";
    *outStream << err.what() << '\n';
    *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    errorFlag = -1000;
  };

  // Now verify that we throw an exception on some of the non-supported cell types.

  *outStream \
  << "\n"
  << "===============================================================================\n"\
  << "| TEST 2: check for unsupported cell types                                     \n"\
  << "===============================================================================\n";
  try{
    try {
      PointTools::getLatticeSize( shards::getCellTopologyData< shards::Quadrilateral<4> >() , 3 , 0 );
    }
    catch (std::invalid_argument err) {
      *outStream << err.what() << "\n";
    }

    try {
      PointTools::getLatticeSize( shards::getCellTopologyData< shards::Hexahedron<8> >() , 3 , 0 );
    }
    catch (std::invalid_argument err) {
      *outStream << err.what() << "\n";
    }
  }
  catch (std::exception &err) {
    *outStream << "UNEXPECTED ERROR !!! ----------------------------------------------------------\n";
    *outStream << err.what() << '\n';
    *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    errorFlag = -1000;
  };

  // Check for malformed point arrays


  *outStream \
  << "\n"
  << "===============================================================================\n"\
  << "| TEST 2: malformed point arrays                                               \n"\
  << "===============================================================================\n";
  try{
    // line: not enough points allocated
    try {
      FieldContainer<double> pts(3,1);
      PointTools::getLatticeSize( shards::getCellTopologyData< shards::Line<2> >() , 5 , 0 );
    }
    catch (std::invalid_argument err) {
      *outStream << err.what() << "\n";
    }
    // line: wrong dimension for points
    try {
      FieldContainer<double> pts(6,2);
      PointTools::getLatticeSize( shards::getCellTopologyData< shards::Line<2> >() , 5 , 0 );
    }
    catch (std::invalid_argument err) {
      *outStream << err.what() << "\n";
    }
    // triangle: too many points allocated
    try {
      FieldContainer<double> pts(4,2);
      PointTools::getLatticeSize( shards::getCellTopologyData< shards::Triangle<3> >() , 3 , 1 );
    }
    catch (std::invalid_argument err) {
      *outStream << err.what() << "\n";
    }
    // triangle: wrong dimension for points
    try {
      FieldContainer<double> pts(6,1);
      PointTools::getLatticeSize( shards::getCellTopologyData< shards::Triangle<3> >() , 3 , 0 );
    }
    catch (std::invalid_argument err) {
      *outStream << err.what() << "\n";
    }
    // tetrahedron: not enough points allocated
    try {
      FieldContainer<double> pts(4,2);
      PointTools::getLatticeSize( shards::getCellTopologyData< shards::Tetrahedron<4> >() , 2 , 0 );
    }
    catch (std::invalid_argument err) {
      *outStream << err.what() << "\n";
    }
    // tetrahedron: wrong dimension for points
    try {
      FieldContainer<double> pts(4,2);
      PointTools::getLatticeSize( shards::getCellTopologyData< shards::Tetrahedron<4> >() , 1 , 0 );
    }
    catch (std::invalid_argument err) {
      *outStream << err.what() << "\n";
    }

  }
  catch (std::exception &err) {
    *outStream << "UNEXPECTED ERROR !!! ----------------------------------------------------------\n";
    *outStream << err.what() << '\n';
    *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    errorFlag = -1000;
  };


  *outStream \
  << "\n"
  << "===============================================================================\n"\
  << "| TEST 3: check values of triangular lattice compared to Warburton's code      \n"\
  << "===============================================================================\n";
  try {
    // triangle case
    const int order = 4;
    const int offset = 0;
    int numPts = PointTools::getLatticeSize( myTri_3 , order , offset );
    int ptDim = 2;
    FieldContainer<double> warpBlendPts( numPts , ptDim );
    PointTools::getLattice<double,FieldContainer<double> >( warpBlendPts ,
                                                            myTri_3 ,
                                                            order ,
                                                            offset ,
                                                            POINTTYPE_WARPBLEND );
    FieldContainer<double> verts( 1, 3 , 2 );
    verts(0,0,0) = -1.0;
    verts(0,0,1) = -1.0/sqrt(3.0);
    verts(0,1,0) = 1.0;
    verts(0,1,1) = -1.0/sqrt(3.0);
    verts(0,2,0) = 0.0;
    verts(0,2,1) = 2.0/sqrt(3.0);

    // holds points on the equilateral triangle
    FieldContainer<double> warpBlendMappedPts( numPts , ptDim );

    CellTools<double>::mapToPhysicalFrame(warpBlendMappedPts ,
                                          warpBlendPts ,
                                          verts ,
                                          myTri_3 ,
                                          0 );

    // Values from MATLAB code
    double points[] = { -1.000000000000000 , -0.577350269189626 ,
                        -0.654653670707977 , -0.577350269189626 ,
                        -0.000000000000000 , -0.577350269189626 ,
                        0.654653670707977  , -0.577350269189626 ,
                        1.000000000000000  , -0.577350269189626 ,
                        -0.827326835353989 , -0.278271574919028 ,
                        -0.327375261332958 , -0.189010195256608 ,
                        0.327375261332958 , -0.189010195256608 ,
                        0.827326835353989,  -0.278271574919028,
                        -0.500000000000000,   0.288675134594813,
                        0.000000000000000,   0.378020390513215,
                        0.500000000000000,   0.288675134594813,
                        -0.172673164646011,   0.855621844108654,
                        0.172673164646011,   0.855621844108654,
                        0,   1.154700538379252 };

    // compare
    for (int i=0;i<numPts;i++) {
      for (int j=0;j<2;j++) {
        int l = 2*i+j;
        if (std::abs(warpBlendMappedPts(i,j) - points[l]) > INTREPID_TOL) {
          errorFlag++;
          *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
          
          // Output the multi-index of the value where the error is:
          *outStream << " At multi-index { ";
          *outStream << i << " ";*outStream << j << " ";
          *outStream << "}  computed value: " << warpBlendMappedPts(i,j)
                     << " but correct value: " << points[l] << "\n";
        }
      }
    }


  }
  catch (std::exception &err) {
    *outStream << "UNEXPECTED ERROR !!! ----------------------------------------------------------\n";
    *outStream << err.what() << '\n';
    *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    errorFlag = -1000;
  };


  *outStream \
  << "\n"
  << "===============================================================================\n"\
  << "| TEST 4: check values of tetrahedral lattice compared to Warburton's code     \n"\
  << "===============================================================================\n";
  try {
    // triangle case
    const int order = 6;
    const int offset = 0;
    int numPts = PointTools::getLatticeSize( myTet_4 , order , offset );
    int ptDim = 3;
    FieldContainer<double> warpBlendPts( numPts , ptDim );
    PointTools::getLattice<double,FieldContainer<double> >( warpBlendPts ,
                                                            myTet_4 ,
                                                            order ,
                                                            offset ,
                                                            POINTTYPE_WARPBLEND );

    FieldContainer<double> verts(1,4,3);
    verts(0,0,0) = -1.0;
    verts(0,0,1) = -1.0/sqrt(3.0);
    verts(0,0,2) = -1.0/sqrt(6.0);
    verts(0,1,0) = 1.0;
    verts(0,1,1) = -1.0/sqrt(3.0);
    verts(0,1,2) = -1.0/sqrt(6.0);
    verts(0,2,0) = 0.0;
    verts(0,2,1) = 2.0 / sqrt(3.0);
    verts(0,2,2) = -1.0/sqrt(6.0);
    verts(0,3,0) = 0.0;
    verts(0,3,1) = 0.0;
    verts(0,3,2) = 3.0 / sqrt(6.0);


    // points on the equilateral tet
    FieldContainer<double> warpBlendMappedPts( numPts , ptDim );

    CellTools<double>::mapToPhysicalFrame(warpBlendMappedPts ,
                                          warpBlendPts ,
                                          verts ,
                                          myTet_4 ,
                                          0 );

    // Values from MATLAB code
    double points[] = {  -1.000000000000000,  -0.577350269189626,  -0.408248290463863,
                        -0.830223896278567,  -0.577350269189626,  -0.408248290463863,
                        -0.468848793470714,  -0.577350269189626,  -0.408248290463863,
                        -0.000000000000000,  -0.577350269189626,  -0.408248290463863,
                        0.468848793470714,  -0.577350269189626, -0.408248290463863,
                        0.830223896278567,  -0.577350269189626,  -0.408248290463863,
                        1.000000000000000,  -0.577350269189626,  -0.408248290463863,
                        -0.915111948139283,  -0.430319850411323,  -0.408248290463863,
                        -0.660434383303427,  -0.381301968982318,  -0.408248290463863,
                        -0.239932664820086,  -0.368405260495326,  -0.408248290463863,
                        0.239932664820086,  -0.368405260495326,  -0.408248290463863,
                        0.660434383303426, -0.381301968982318,  -0.408248290463863,
                        0.915111948139283,  -0.430319850411323,  -0.408248290463863,
                        -0.734424396735357,  -0.117359831084509,  -0.408248290463863,
                        -0.439014646886819,  -0.023585152684228,  -0.408248290463863,
                        -0.000000000000000,  -0.000000000000000,  -0.408248290463863,
                        0.439014646886819,  -0.023585152684228,  -0.408248290463863,
                        0.734424396735357,  -0.117359831084509,  -0.408248290463863,
                        -0.500000000000000,   0.288675134594813,  -0.408248290463863,
                        -0.199081982066733,   0.391990413179555,  -0.408248290463863,
                        0.199081982066733,   0.391990413179555,  -0.408248290463863,
                        0.500000000000000,   0.288675134594813,  -0.408248290463863,
                        -0.265575603264643,   0.694710100274135,  -0.408248290463863,
                        -0.000000000000000,  0.762603937964635,  -0.408248290463863,
                        0.265575603264643,   0.694710100274135,  -0.408248290463863,
                        -0.084888051860716,   1.007670119600949,  -0.408248290463863,
                        0.084888051860716,   1.007670119600949,  -0.408248290463863,
                        0,   1.154700538379252,  -0.408248290463863,
                        -0.915111948139284,  -0.528340129596858,  -0.269626682252082,
                        -0.660434383303427,  -0.512000835787190,  -0.223412180441618,
                        -0.239932664820086,  -0.507701932958193,  -0.211253047073435,
                        0.239932664820086,  -0.507701932958193,  -0.211253047073435,
                        0.660434383303426,  -0.512000835787190,  -0.223412180441618,
                        0.915111948139284,  -0.528340129596858,  -0.269626682252082,
                        -0.773622922202284,  -0.315952535579882,  -0.223412180441618,
                        -0.421605613935553,  -0.243414114697549,  -0.172119771139157,
                        -0.000000000000000,  -0.224211101329670,  -0.158541190167514,
                        0.421605613935553,  -0.243414114697549,  -0.172119771139157,
                        0.773622922202284,  -0.315952535579882,  -0.223412180441618,
                        -0.559649103902302,   0.046063183547205,  -0.211253047073435,
                        -0.194172509561981,   0.112105550664835,  -0.158541190167514,
                        0.194172509561981,   0.112105550664835,  -0.158541190167514,
                        0.559649103902302,   0.046063183547205,  -0.211253047073435,
                        -0.319716439082216,   0.461638749410988,  -0.211253047073435,
                        -0.000000000000000,   0.486828229395098,  -0.172119771139157,
                        0.319716439082216,   0.461638749410988,  -0.211253047073435,
                        -0.113188538898858,   0.827953371367071,  -0.223412180441618,
                        0.113188538898858,   0.827953371367071,  -0.223412180441618,
                        -0.000000000000000,   1.056680259193716,  -0.269626682252082,
                        -0.734424396735357,  -0.424020123154587,   0.025434853622935,
                        -0.439014646886819,  -0.392761897021160,   0.113846468290170,
                        -0.000000000000000,  -0.384900179459751,   0.136082763487954,
                        0.439014646886819,  -0.392761897021160,   0.113846468290170,
                        0.734424396735357,  -0.424020123154587,   0.025434853622935,
                        -0.559649103902302,  -0.183816888326860,   0.113846468290170,
                        -0.194172509561981,  -0.112105550664835,   0.158541190167514,
                        0.194172509561981,  -0.112105550664835,   0.158541190167514,
                        0.559649103902302,  -0.183816888326860,   0.113846468290170,
                        -0.333333333333333,   0.192450089729875,   0.136082763487954,
                        -0.000000000000000,   0.224211101329670,   0.158541190167514,
                        0.333333333333333,   0.192450089729875,   0.136082763487954,
                        -0.120634457015483,   0.576578785348020,   0.113846468290170,
                        0.120634457015482,   0.576578785348020,   0.113846468290170,
                        -0.000000000000000,   0.848040246309174,   0.025434853622935,
                        -0.500000000000000,  -0.288675134594813,   0.408248290463863,
                        -0.199081982066733,  -0.254236708399899,   0.505654869247127,
                        0.199081982066733,  -0.254236708399899,   0.505654869247127,
                        0.500000000000000,  -0.288675134594813,   0.408248290463863,
                        -0.319716439082216,  -0.045291699705599,   0.505654869247127,
                        -0.000000000000000,  -0.000000000000000,   0.516359313417471,
                        0.319716439082216,  -0.045291699705599,   0.505654869247127,
                        -0.120634457015483,   0.299528408105498,   0.505654869247127,
                        0.120634457015483,  0.299528408105498,   0.505654869247127,
                        -0.000000000000000,   0.577350269189626,   0.408248290463863,
                        -0.265575603264643,  -0.153330146035039,   0.791061727304791,
                        -0.000000000000000,  -0.130698866804872,   0.855072651347100,
                        0.265575603264643,  -0.153330146035039,   0.791061727304791,
                        -0.113188538898858,   0.065349433402436,   0.855072651347100,
                        0.113188538898858,   0.065349433402436,   0.855072651347099,
                        0,   0.306660292070078,   0.791061727304791,
                        -0.084888051860717,  -0.049010139592768,   1.086123263179808,
                        0.084888051860717,  -0.049010139592768,   1.086123263179808,
                        0.000000000000000,   0.098020279185535,   1.086123263179808,
                        0,                   0,   1.224744871391589
    };

    // compare
    for (int i=0;i<numPts;i++) {
      for (int j=0;j<ptDim;j++) {
        int l = ptDim*i+j;
        if (std::abs(warpBlendMappedPts(i,j) - points[l]) > INTREPID_TOL) {
          errorFlag++;
          *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
          
          // Output the multi-index of the value where the error is:
          *outStream << " At multi-index { ";
          *outStream << i << " ";*outStream << j << " ";
          *outStream << "}  computed value: " << warpBlendMappedPts(i,j)
                     << " but correct value: " << points[l] << "\n";
        }
      }
    }


  }
  catch (std::exception &err) {
    *outStream << "UNEXPECTED ERROR !!! ----------------------------------------------------------\n";
    *outStream << err.what() << '\n';
    *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    errorFlag = -1000;
  };



  if (errorFlag != 0)
    std::cout << "End Result: TEST FAILED\n";
  else
    std::cout << "End Result: TEST PASSED\n";

  // reset format state of std::cout
  std::cout.copyfmt(oldFormatState);

  return errorFlag;
}
