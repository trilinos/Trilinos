// @HEADER
// ************************************************************************
//
//                           Intrepid Package
//                 Copyright (2007) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Pavel Bochev (pbboche@sandia.gov) or
//                    Denis Ridzal (dridzal@sandia.gov) or
//                    Robert C. Kirby (robert.c.kirby@sandia.gov)
//
// ************************************************************************
// @HEADER

/** \file test_01.cpp
\brief  Unit tests for the Intrepid::Basis_F0_TRI_C2_FEM_FIAT class.
\author Created by R. Kirby via FIAT.
*/

#include "Intrepid_DefaultBasisFactory.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"

using namespace std;
using namespace Intrepid;

int main(int argc, char *argv[]) {

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
    << "|                 Unit Test (Basis_F0_TRI_C2_FEM_FIAT)                                   |\n" \
    << "|                                                                             |\n" \
    << "|     1) Basis creation, conversion of Dof tags into enumeration and back     |\n" \
    << "|     2) Basis values for VALUE, GRAD, CURL, and Dk operators                 |\n" \
    << "|                                                                             |\n" \
    << "|  Questions? Contact  Pavel Bochev (pbboche@sandia.gov) or                   |\n" \
    << "|                      Denis Ridzal (dridzal@sandia.gov) or                   |\n" \
    << "|                      Robert C. Kirby (robert.c.kirby@ttu.edu)               |\n" \
    << "|                                                                             |\n" \
    << "|  Intrepid's website: http://trilinos.sandia.gov/packages/intrepid           |\n" \
    << "|  Trilinos website:   http://trilinos.sandia.gov                             |\n" \
    << "|  FIAT website:       http://www.fenics.org/fiat                             |\n" \
    << "|                                                                             |\n" \
    << "===============================================================================\n"\
    << "| TEST 1: Basis creation, exception testing                                   |\n"\
    << "===============================================================================\n";

  int errorFlag = 0;
  int beginThrowNumber = TestForException_getThrowNumber();
  int endThrowNumber = beginThrowNumber + 1;
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  endThrowNumber += 5;
#endif

  // Reference element points are the standard equispaced lattice of degree 2
  Teuchos::Array< Point<double> > elNodes;
  Point<double> tempPt(2, FRAME_REFERENCE);
  elNodes.assign(6, tempPt);  
#ifndef DOXYGEN_SHOULD_SKIP_THIS
  elNodes[0] = Point<double>( 0.0000000000000000e+00 , 0.0000000000000000e+00 , FRAME_REFERENCE);
  elNodes[1] = Point<double>( 5.0000000000000000e-01 , 0.0000000000000000e+00 , FRAME_REFERENCE);
  elNodes[2] = Point<double>( 1.0000000000000000e+00 , 0.0000000000000000e+00 , FRAME_REFERENCE);
  elNodes[3] = Point<double>( 0.0000000000000000e+00 , 5.0000000000000000e-01 , FRAME_REFERENCE);
  elNodes[4] = Point<double>( 5.0000000000000000e-01 , 5.0000000000000000e-01 , FRAME_REFERENCE);
  elNodes[5] = Point<double>( 0.0000000000000000e+00 , 1.0000000000000000e+00 , FRAME_REFERENCE);
#endif


  try{
    FieldContainer<double> vals;
    DefaultBasisFactory<double> BFactory;
    Teuchos::RCP<Basis<double> > triBasis = BFactory.create(FIELD_FORM_0, 
                                                            CELL_TRI, 
                                                            RECONSTRUCTION_SPACE_COMPLETE, 
                                                            2, 
                                                            BASIS_FEM_FIAT, 
                                                            COORDINATES_CARTESIAN);
    // Exception 1: DIV cannot be applied to scalar functions
    try {
      triBasis->getValues(vals, elNodes, OPERATOR_DIV);
    }
    catch (std::logic_error err) {
      *outStream << "Expected Error 1----------------------------------------------------------------\n";
      *outStream << err.what() << '\n';
      *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    };

#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
    // Exceptions 2-6: all tags/bd Ids below are wrong and should cause getLocalDofEnumeration() and 
    // getLocalDofTag() to access invalid array elements thereby causing Teuchos bounds check exception
    try {      
      LocalDofTag myTag = {{3,0,0,1}};
      triBasis -> getLocalDofEnumeration(myTag);
    }
    catch (std::logic_error err) {
      *outStream << "Expected Error 2----------------------------------------------------------------\n";
      *outStream << err.what() << '\n';
      *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    };

    try {
      LocalDofTag myTag = {{1,1,2,0}};
      triBasis -> getLocalDofEnumeration(myTag);
    }
    catch (std::logic_error err) {
      *outStream << "Expected Error 3----------------------------------------------------------------\n";
      *outStream << err.what() << '\n';
      *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    };
    try {
      LocalDofTag myTag = {{0,3,0,0}};
      triBasis -> getLocalDofEnumeration(myTag);
    }
    catch (std::logic_error err) {
      *outStream << "Expected Error 4----------------------------------------------------------------\n";
      *outStream << err.what() << '\n';
      *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    }; 

    try {
      int bfId = 6;
      triBasis -> getLocalDofTag(bfId);
    }
    catch (std::logic_error err) {
      *outStream << "Expected Error 5----------------------------------------------------------------\n";
      *outStream << err.what() << '\n';
      *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    };

    try {
      int bfId = -1;
      triBasis -> getLocalDofTag(bfId);
    }
    catch (std::logic_error err) {
      *outStream << "Expected Error 6----------------------------------------------------------------\n";
      *outStream << err.what() << '\n';
      *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    };
#endif
  }
  catch (std::logic_error err) {
    *outStream << "UNEXPECTED ERROR !!! ----------------------------------------------------------\n";
    *outStream << err.what() << '\n';
    *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    errorFlag = -1000;
  };

  // Check if number of thrown exceptions matches the one we expect (1 + 5)
  if (TestForException_getThrowNumber() != endThrowNumber) {
    errorFlag++;
    *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
  }

  *outStream \
    << "\n"
    << "===============================================================================\n"\
    << "| TEST 2: correctness of tag to enum and enum to tag lookups                  |\n"\
    << "===============================================================================\n";

  try{
    DefaultBasisFactory<double> BFactory;
    Teuchos::RCP<Basis<double> > triBasis = BFactory.create(FIELD_FORM_0, 
                                                            CELL_TRI, 
                                                            RECONSTRUCTION_SPACE_COMPLETE, 
                                                            2, 
                                                            BASIS_FEM_FIAT, 
                                                            COORDINATES_CARTESIAN);
    Teuchos::Array<LocalDofTag> allTags;
    allTags = triBasis -> getAllLocalDofTags();

    // Loop over all tags, lookup the associated dof enumeration and then lookup the tag again
    for (unsigned i = 0; i < allTags.size(); i++) {
      int        bfId  = triBasis -> getLocalDofEnumeration(allTags[i]);
      LocalDofTag myTag = triBasis -> getLocalDofTag(bfId);
       if( !( (myTag.tag_[0] == allTags[i].tag_[0]) &&
             (myTag.tag_[1] == allTags[i].tag_[1]) &&
             (myTag.tag_[2] == allTags[i].tag_[2]) &&
             (myTag.tag_[3] == allTags[i].tag_[3]) ) ) {
        errorFlag++;
        *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
        *outStream << " getLocalDofEnumeration( {" 
          << allTags[i].tag_[0] << ", " 
          << allTags[i].tag_[1] << ", " 
          << allTags[i].tag_[2] << ", " 
          << allTags[i].tag_[3] << "}) = " << bfId <<" but \n";   
        *outStream << " getLocalDofTag(" << bfId << ") = { "
          << myTag.tag_[0] << ", " 
          << myTag.tag_[1] << ", " 
          << myTag.tag_[2] << ", " 
          << myTag.tag_[3] << "}\n";        
      }
    }

    // Now do the same but loop over basis functions
    for( int bfId = 0; bfId < triBasis -> getNumLocalDof(); bfId++) {
      LocalDofTag myTag  = triBasis -> getLocalDofTag(bfId);
      int myBfId = triBasis -> getLocalDofEnumeration(myTag);
      if( bfId != myBfId) {
        errorFlag++;
        *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
        *outStream << " getLocalDofTag(" << bfId << ") = { "
          << myTag.tag_[0] << ", " 
          << myTag.tag_[1] << ", " 
          << myTag.tag_[2] << ", " 
          << myTag.tag_[3] << "} but getLocalDofEnumeration({" 
          << myTag.tag_[0] << ", " 
          << myTag.tag_[1] << ", " 
          << myTag.tag_[2] << ", " 
          << myTag.tag_[3] << "} ) = " << myBfId << "\n";
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

  // VALUE: Each correct basis function at each correct point (point increasing fastest)
#ifndef DOXYGEN_SHOULD_SKIP_THIS
  double basisValues[] = { 1.0000000000000000e+00 , -2.7755575615628914e-17 , 0.0000000000000000e+00 , 1.3877787807814457e-17 , 2.7755575615628914e-17 , 5.5511151231257827e-17 , 1.2490009027033011e-16 , 1.3877787807814457e-17 , 0.0000000000000000e+00 , 1.0000000000000000e+00 , 2.7755575615628914e-17 , 2.7755575615628914e-17 , 8.3266726846886741e-17 , 1.0000000000000000e+00 , 0.0000000000000000e+00 , 1.3877787807814457e-17 , 5.5511151231257827e-17 , 2.7755575615628914e-17 , 9.7144514654701197e-17 , 2.7755575615628914e-17 , -4.1633363423443370e-17 , -6.9388939039072284e-18 , -2.7755575615628914e-17 , 1.0000000000000002e+00 , 1.3877787807814457e-16 , -1.3877787807814457e-17 , -4.1633363423443370e-17 , -6.9388939039072284e-18 , 1.0000000000000002e+00 , -2.7755575615628914e-17 , 1.6653345369377348e-16 , 5.5511151231257827e-17 , 1.0000000000000000e+00 , 0.0000000000000000e+00 , 0.0000000000000000e+00 , 0.0000000000000000e+00 };
#endif

  // DERIVS: Components of derivatives increasing fastest, followed by points, then basis functions
#ifndef DOXYGEN_SHOULD_SKIP_THIS
double basisD1s[] = { -3.0000000000000000e+00 , -2.9999999999999991e+00 , -1.0000000000000009e+00 , 3.3306690738754696e-16 , 0.0000000000000000e+00 , -1.0000000000000009e+00 , 4.0000000000000009e+00 , 1.4432899320127035e-15 , -1.7763568394002505e-15 , -1.7763568394002505e-15 , 1.7763568394002505e-15 , 4.0000000000000000e+00 , -1.0000000000000000e+00 , -1.0000000000000002e+00 , 1.0000000000000000e+00 , -1.1102230246251565e-16 , 0.0000000000000000e+00 , -1.0000000000000004e+00 , 5.1810407815841009e-17 , -2.0000000000000004e+00 , 5.5511151231257827e-16 , 2.0000000000000009e+00 , -5.5511151231257827e-16 , 2.0000000000000004e+00 , 1.0000000000000009e+00 , 9.9999999999999978e-01 , 3.0000000000000009e+00 , -3.3306690738754696e-16 , 0.0000000000000000e+00 , -9.9999999999999978e-01 , -4.0000000000000018e+00 , -4.0000000000000000e+00 , 1.7763568394002505e-15 , 4.0000000000000018e+00 , -1.7763568394002505e-15 , -8.8817841970012523e-16 , -9.9999999999999967e-01 , -9.9999999999999978e-01 , -9.9999999999999967e-01 , -5.5511151231257827e-17 , 0.0000000000000000e+00 , 1.0000000000000007e+00 , 1.9999999999999993e+00 , -4.1633363423443370e-16 , 2.0000000000000009e+00 , 4.4408920985006262e-16 , -2.0000000000000009e+00 , -6.1062266354383610e-16 , 1.0000000000000009e+00 , 1.0000000000000004e+00 , 9.9999999999999956e-01 , 0.0000000000000000e+00 , 0.0000000000000000e+00 , 9.9999999999999978e-01 , -2.0000000000000004e+00 , -1.9999999999999984e+00 , 1.9999999999999982e+00 , 1.9999999999999980e+00 , -1.9999999999999982e+00 , -2.0000000000000000e+00 , 1.0000000000000013e+00 , 1.0000000000000009e+00 , -1.0000000000000004e+00 , 3.3306690738754696e-16 , 0.0000000000000000e+00 , 2.9999999999999991e+00 , -9.8439774850097100e-16 , 2.2204460492503131e-16 , 4.0000000000000018e+00 , 4.4408920985006262e-16 , -4.0000000000000018e+00 , -4.0000000000000009e+00 };

double basisD2s[] = { 3.9999999999999911e+00 , 3.9999999999999925e+00 , 3.9999999999999898e+00 , 3.9999999999999956e+00 , 6.2172489379008766e-15 , -2.2204460492503131e-16 , 0.0000000000000000e+00 , 0.0000000000000000e+00 , 4.0000000000000000e+00 , -7.9999999999999867e+00 , -3.9999999999999987e+00 , -6.6613381477509392e-15 , 5.6843418860808042e-15 , 3.9999999999999991e+00 , 5.3290705182007514e-15 , -5.6843418860808042e-15 , -3.9999999999999991e+00 , -7.9999999999999876e+00 , 4.0000000000000027e+00 , 4.0000000000000018e+00 , 4.0000000000000018e+00 , 4.0000000000000044e+00 , 1.5543122344752192e-15 , 2.2204460492503131e-16 , 0.0000000000000000e+00 , 0.0000000000000000e+00 , 4.0000000000000062e+00 , -8.0000000000000071e+00 , -4.0000000000000036e+00 , 1.9984014443252818e-15 , 3.5527136788004997e-15 , 3.9999999999999973e+00 , -4.4408920985006262e-15 , -3.5527136788004997e-15 , -3.9999999999999973e+00 , -8.0000000000000053e+00 , 4.0000000000000080e+00 , 4.0000000000000053e+00 , 4.0000000000000089e+00 , 4.0000000000000062e+00 , -2.2204460492503131e-15 , 2.2204460492503131e-16 , 0.0000000000000000e+00 , 0.0000000000000000e+00 , 4.0000000000000044e+00 , -8.0000000000000142e+00 , -4.0000000000000036e+00 , 9.5479180117763462e-15 , 1.4210854715201945e-15 , 3.9999999999999889e+00 , -1.2878587085651816e-14 , -1.4210854715201945e-15 , -3.9999999999999889e+00 , -8.0000000000000107e+00 , 4.0000000000000009e+00 , 4.0000000000000018e+00 , 4.0000000000000018e+00 , 3.9999999999999991e+00 , 1.3322676295501878e-15 , 2.2204460492503131e-16 , 0.0000000000000000e+00 , 0.0000000000000000e+00 , 4.0000000000000018e+00 , -8.0000000000000000e+00 , -4.0000000000000036e+00 , -8.8817841970012523e-16 , -4.2632564145606025e-15 , 4.0000000000000036e+00 , 1.7763568394002505e-15 , 4.2632564145606025e-15 , -4.0000000000000036e+00 , -8.0000000000000036e+00 , 3.9999999999999960e+00 , 3.9999999999999964e+00 , 3.9999999999999987e+00 , 3.9999999999999916e+00 , -2.8865798640254070e-15 , 0.0000000000000000e+00 , 0.0000000000000000e+00 , 0.0000000000000000e+00 , 3.9999999999999916e+00 , -7.9999999999999876e+00 , -3.9999999999999933e+00 , 7.5495165674510645e-15 , -6.3948846218408979e-15 , 3.9999999999999862e+00 , -6.6613381477509392e-15 , 6.3948846218408979e-15 , -3.9999999999999862e+00 , -7.9999999999999911e+00 , 4.0000000000000009e+00 , 4.0000000000000000e+00 , 4.0000000000000009e+00 , 3.9999999999999933e+00 , 0.0000000000000000e+00 , 6.6613381477509392e-16 , 0.0000000000000000e+00 , 0.0000000000000000e+00 , 3.9999999999999880e+00 , -7.9999999999999938e+00 , -4.0000000000000000e+00 , 1.5543122344752192e-15 , -1.4210854715202010e-14 , 3.9999999999999933e+00 , 1.7763568394002505e-15 , 1.4210854715202010e-14 , -3.9999999999999933e+00 , -7.9999999999999929e+00 };
#endif


  // CURL: each correct values of the curls of the basis functions at the points (point increasing fastest
#ifndef DOXYGEN_SHOULD_SKIP_THIS
  double basisCurls[] = { -2.9999999999999991e+00 , 3.0000000000000000e+00 , 3.3306690738754696e-16 , 1.0000000000000009e+00 , -1.0000000000000009e+00 , -0.0000000000000000e+00 , 1.4432899320127035e-15 , -4.0000000000000009e+00 , -1.7763568394002505e-15 , 1.7763568394002505e-15 , 4.0000000000000000e+00 , -1.7763568394002505e-15 , -1.0000000000000002e+00 , 1.0000000000000000e+00 , -1.1102230246251565e-16 , -1.0000000000000000e+00 , -1.0000000000000004e+00 , -0.0000000000000000e+00 , -2.0000000000000004e+00 , -5.1810407815841009e-17 , 2.0000000000000009e+00 , -5.5511151231257827e-16 , 2.0000000000000004e+00 , 5.5511151231257827e-16 , 9.9999999999999978e-01 , -1.0000000000000009e+00 , -3.3306690738754696e-16 , -3.0000000000000009e+00 , -9.9999999999999978e-01 , -0.0000000000000000e+00 , -4.0000000000000000e+00 , 4.0000000000000018e+00 , 4.0000000000000018e+00 , -1.7763568394002505e-15 , -8.8817841970012523e-16 , 1.7763568394002505e-15 , -9.9999999999999978e-01 , 9.9999999999999967e-01 , -5.5511151231257827e-17 , 9.9999999999999967e-01 , 1.0000000000000007e+00 , -0.0000000000000000e+00 , -4.1633363423443370e-16 , -1.9999999999999993e+00 , 4.4408920985006262e-16 , -2.0000000000000009e+00 , -6.1062266354383610e-16 , 2.0000000000000009e+00 , 1.0000000000000004e+00 , -1.0000000000000009e+00 , 0.0000000000000000e+00 , -9.9999999999999956e-01 , 9.9999999999999978e-01 , -0.0000000000000000e+00 , -1.9999999999999984e+00 , 2.0000000000000004e+00 , 1.9999999999999980e+00 , -1.9999999999999982e+00 , -2.0000000000000000e+00 , 1.9999999999999982e+00 , 1.0000000000000009e+00 , -1.0000000000000013e+00 , 3.3306690738754696e-16 , 1.0000000000000004e+00 , 2.9999999999999991e+00 , -0.0000000000000000e+00 , 2.2204460492503131e-16 , 9.8439774850097100e-16 , 4.4408920985006262e-16 , -4.0000000000000018e+00 , -4.0000000000000009e+00 , 4.0000000000000018e+00 };
#endif


  try{
    FieldContainer<double> vals;
    DefaultBasisFactory<double> BFactory;
    Teuchos::RCP<Basis<double> > triBasis = BFactory.create(FIELD_FORM_0, 
                                                            CELL_TRI, 
                                                            RECONSTRUCTION_SPACE_COMPLETE, 
                                                            2, 
                                                            BASIS_FEM_FIAT, 
                                                            COORDINATES_CARTESIAN);

    // Check VALUE of basis functions
    triBasis -> getValues(vals, elNodes, OPERATOR_VALUE);
    for (int i=0; i < vals.getSize(); i++) {
      bool fail = false;
      if (std::abs(basisValues[i]) < INTREPID_FIAT_TOL ) {
        if (std::abs(vals[i]-basisValues[i]) > INTREPID_FIAT_TOL ) {
          fail = true;
        }
      }
      else if (std::abs(vals[i] - basisValues[i]) / std::abs(basisValues[i] ) > INTREPID_FIAT_TOL ) {
        fail = true;
      }
      if (fail) {
        errorFlag++;
        *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
        // Get the multi-index of the value where the error is:
        Teuchos::Array<int> myIndex;
        vals.getMultiIndex(myIndex,i);
        *outStream << " At multi-index { ";
        for(int j = 0; j < vals.getRank(); j++) {
          *outStream << myIndex[j] << " ";
        }
        *outStream << "}  computed value: " << vals[i] 
          << " but reference value: " << basisValues[i] << "\n";
      }
    }

    // Check GRAD of basis function
    triBasis -> getValues(vals, elNodes, OPERATOR_GRAD);
    for (int i=0; i < vals.getSize(); i++) {
      bool fail = false;
      if (std::abs(basisD1s[i]) < 10.0 * INTREPID_FIAT_TOL ) {
        if (std::abs(vals[i]-basisD1s[i]) > 10.0 * INTREPID_FIAT_TOL ) {
          fail = true;
        }
      }
      else if (std::abs(vals[i] - basisD1s[i]) / std::abs(basisD1s[i] ) > 10.0 * INTREPID_FIAT_TOL ) {
        fail = true;
      }
      if (fail) {
        errorFlag++;
        *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
        // Get the multi-index of the value where the error is:
        Teuchos::Array<int> myIndex;
        vals.getMultiIndex(myIndex,i);
        *outStream << " At multi-index { ";
        for(int j = 0; j < vals.getRank(); j++) {
          *outStream << myIndex[j] << " ";
        }
        *outStream << "}  computed grad component: " << vals[i] 
          << " but reference grad component: " << basisD1s[i] << "\n";
      }
    }

    // Check CURL of basis function
    triBasis -> getValues(vals, elNodes, OPERATOR_CURL);
    for (int i=0; i < vals.getSize(); i++) {
      bool fail = false;
      if (std::abs(basisCurls[i]) < 10.0 * INTREPID_FIAT_TOL ) {
        if (std::abs(vals[i]-basisCurls[i]) > 10.0 * INTREPID_FIAT_TOL ) {
          fail = true;
        }
      }
      else if (std::abs(vals[i] - basisCurls[i]) / std::abs(basisCurls[i] ) > 10.0 * INTREPID_FIAT_TOL ) {
        fail = true;
      }
      if (fail) {
        errorFlag++;
        *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
        // Get the multi-index of the value where the error is:
        Teuchos::Array<int> myIndex;
        vals.getMultiIndex(myIndex,i);
        *outStream << " At multi-index { ";
        for(int j = 0; j < vals.getRank(); j++) {
          *outStream << myIndex[j] << " ";
        }
        *outStream << "}  computed curl component: " << vals[i] 
          << " but reference curl component: " << basisCurls[i] << "\n";
      }
    }


    // Check D1 of basis function
    triBasis -> getValues(vals, elNodes, OPERATOR_D1);
    for (int i=0; i < vals.getSize(); i++) {
      bool fail = false;
      if (std::abs(basisD1s[i]) < pow(15.0,1) * INTREPID_FIAT_TOL ) {
        if (std::abs(vals[i]-basisD1s[i]) > pow(15.0,1) * INTREPID_FIAT_TOL ) {
          fail = true;
        }
      }
      else if (std::abs(vals[i] - basisD1s[i]) / std::abs(basisD1s[i] ) > pow(15.0,1) * INTREPID_FIAT_TOL ) {
        fail = true;
      }
      if (fail) {
        errorFlag++;
        *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
        // Get the multi-index of the value where the error is:
        Teuchos::Array<int> myIndex;
        vals.getMultiIndex(myIndex,i);
        *outStream << " At multi-index { ";
        for(int j = 0; j < vals.getRank(); j++) {
          *outStream << myIndex[j] << " ";
        }
        *outStream << "}  computed D1 component: " << vals[i] 
          << " but reference D1 component: " << basisD1s[i] << "\n";
      }
    }

    // Check D2 of basis function
    triBasis -> getValues(vals, elNodes, OPERATOR_D2);
    for (int i=0; i < vals.getSize(); i++) {
      bool fail = false;
      if (std::abs(basisD2s[i]) < pow(15.0,2) * INTREPID_FIAT_TOL ) {
        if (std::abs(vals[i]-basisD2s[i]) > pow(15.0,2) * INTREPID_FIAT_TOL ) {
          fail = true;
        }
      }
      else if (std::abs(vals[i] - basisD2s[i]) / std::abs(basisD2s[i] ) > pow(15.0,2) * INTREPID_FIAT_TOL ) {
        fail = true;
      }
      if (fail) {
        errorFlag++;
        *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
        // Get the multi-index of the value where the error is:
        Teuchos::Array<int> myIndex;
        vals.getMultiIndex(myIndex,i);
        *outStream << " At multi-index { ";
        for(int j = 0; j < vals.getRank(); j++) {
          *outStream << myIndex[j] << " ";
        }
        *outStream << "}  computed D2 component: " << vals[i] 
          << " but reference D2 component: " << basisD2s[i] << "\n";
      }
    }

      // Check D3 of basis function: should be zero 
      triBasis -> getValues(vals, elNodes, OPERATOR_D3);
      for (int i=0; i < vals.getSize(); i++) {
        if ( std::abs(vals[i])  > pow( 15.0 , 3 ) * INTREPID_FIAT_TOL) {
          errorFlag++;
          *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
          // Get the multi-index of the value where the error is and the operator order
          Teuchos::Array<int> myIndex;
          vals.getMultiIndex(myIndex,i);
          *outStream << " At multi-index { ";
          for(int j = 0; j < vals.getRank(); j++) {
            *outStream << myIndex[j] << " ";
          }
          *outStream << "}  computed D3 component: " << vals[i] 
            << " but reference D3 component:  0 \n";
        }
      }


  }
  
  // Catch unexpected errors
  catch (std::logic_error err) {
    *outStream << err.what() << "\n\n";
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

