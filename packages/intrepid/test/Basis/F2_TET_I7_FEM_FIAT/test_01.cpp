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
\brief  Unit tests for the Intrepid::Basis_F2_TET_I7_FEM_FIAT class.
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
    << "|                 Unit Test (Basis_F2_TET_I7_FEM_FIAT)                                   |\n" \
    << "|                                                                             |\n" \
    << "|     1) Basis creation, conversion of Dof tags into enumeration and back     |\n" \
    << "|     2) Basis values for VALUE and DIV operators                             |\n" \
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

#ifdef HAVE_INTREPID_FIAT_VERY_HIGH_ORDER

  int beginThrowNumber = TestForException_getThrowNumber();
  int endThrowNumber = beginThrowNumber + 1;
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  endThrowNumber += 5;
#endif

  // Reference element points are the standard equispaced lattice of degree 7
  Teuchos::Array< Point<double> > elNodes;
  Point<double> tempPt(3, FRAME_REFERENCE);
  elNodes.assign(120, tempPt);  
#ifndef DOXYGEN_SHOULD_SKIP_THIS
  elNodes[0] = Point<double>( 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , FRAME_REFERENCE);
  elNodes[1] = Point<double>( 1.428571428571428e-01 , 0.000000000000000e+00 , 0.000000000000000e+00 , FRAME_REFERENCE);
  elNodes[2] = Point<double>( 2.857142857142857e-01 , 0.000000000000000e+00 , 0.000000000000000e+00 , FRAME_REFERENCE);
  elNodes[3] = Point<double>( 4.285714285714285e-01 , 0.000000000000000e+00 , 0.000000000000000e+00 , FRAME_REFERENCE);
  elNodes[4] = Point<double>( 5.714285714285714e-01 , 0.000000000000000e+00 , 0.000000000000000e+00 , FRAME_REFERENCE);
  elNodes[5] = Point<double>( 7.142857142857142e-01 , 0.000000000000000e+00 , 0.000000000000000e+00 , FRAME_REFERENCE);
  elNodes[6] = Point<double>( 8.571428571428571e-01 , 0.000000000000000e+00 , 0.000000000000000e+00 , FRAME_REFERENCE);
  elNodes[7] = Point<double>( 1.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , FRAME_REFERENCE);
  elNodes[8] = Point<double>( 0.000000000000000e+00 , 1.428571428571428e-01 , 0.000000000000000e+00 , FRAME_REFERENCE);
  elNodes[9] = Point<double>( 1.428571428571428e-01 , 1.428571428571428e-01 , 0.000000000000000e+00 , FRAME_REFERENCE);
  elNodes[10] = Point<double>( 2.857142857142857e-01 , 1.428571428571428e-01 , 0.000000000000000e+00 , FRAME_REFERENCE);
  elNodes[11] = Point<double>( 4.285714285714285e-01 , 1.428571428571428e-01 , 0.000000000000000e+00 , FRAME_REFERENCE);
  elNodes[12] = Point<double>( 5.714285714285714e-01 , 1.428571428571428e-01 , 0.000000000000000e+00 , FRAME_REFERENCE);
  elNodes[13] = Point<double>( 7.142857142857142e-01 , 1.428571428571428e-01 , 0.000000000000000e+00 , FRAME_REFERENCE);
  elNodes[14] = Point<double>( 8.571428571428571e-01 , 1.428571428571428e-01 , 0.000000000000000e+00 , FRAME_REFERENCE);
  elNodes[15] = Point<double>( 0.000000000000000e+00 , 2.857142857142857e-01 , 0.000000000000000e+00 , FRAME_REFERENCE);
  elNodes[16] = Point<double>( 1.428571428571428e-01 , 2.857142857142857e-01 , 0.000000000000000e+00 , FRAME_REFERENCE);
  elNodes[17] = Point<double>( 2.857142857142857e-01 , 2.857142857142857e-01 , 0.000000000000000e+00 , FRAME_REFERENCE);
  elNodes[18] = Point<double>( 4.285714285714285e-01 , 2.857142857142857e-01 , 0.000000000000000e+00 , FRAME_REFERENCE);
  elNodes[19] = Point<double>( 5.714285714285714e-01 , 2.857142857142857e-01 , 0.000000000000000e+00 , FRAME_REFERENCE);
  elNodes[20] = Point<double>( 7.142857142857142e-01 , 2.857142857142857e-01 , 0.000000000000000e+00 , FRAME_REFERENCE);
  elNodes[21] = Point<double>( 0.000000000000000e+00 , 4.285714285714285e-01 , 0.000000000000000e+00 , FRAME_REFERENCE);
  elNodes[22] = Point<double>( 1.428571428571428e-01 , 4.285714285714285e-01 , 0.000000000000000e+00 , FRAME_REFERENCE);
  elNodes[23] = Point<double>( 2.857142857142857e-01 , 4.285714285714285e-01 , 0.000000000000000e+00 , FRAME_REFERENCE);
  elNodes[24] = Point<double>( 4.285714285714285e-01 , 4.285714285714285e-01 , 0.000000000000000e+00 , FRAME_REFERENCE);
  elNodes[25] = Point<double>( 5.714285714285714e-01 , 4.285714285714285e-01 , 0.000000000000000e+00 , FRAME_REFERENCE);
  elNodes[26] = Point<double>( 0.000000000000000e+00 , 5.714285714285714e-01 , 0.000000000000000e+00 , FRAME_REFERENCE);
  elNodes[27] = Point<double>( 1.428571428571428e-01 , 5.714285714285714e-01 , 0.000000000000000e+00 , FRAME_REFERENCE);
  elNodes[28] = Point<double>( 2.857142857142857e-01 , 5.714285714285714e-01 , 0.000000000000000e+00 , FRAME_REFERENCE);
  elNodes[29] = Point<double>( 4.285714285714285e-01 , 5.714285714285714e-01 , 0.000000000000000e+00 , FRAME_REFERENCE);
  elNodes[30] = Point<double>( 0.000000000000000e+00 , 7.142857142857142e-01 , 0.000000000000000e+00 , FRAME_REFERENCE);
  elNodes[31] = Point<double>( 1.428571428571428e-01 , 7.142857142857142e-01 , 0.000000000000000e+00 , FRAME_REFERENCE);
  elNodes[32] = Point<double>( 2.857142857142857e-01 , 7.142857142857142e-01 , 0.000000000000000e+00 , FRAME_REFERENCE);
  elNodes[33] = Point<double>( 0.000000000000000e+00 , 8.571428571428571e-01 , 0.000000000000000e+00 , FRAME_REFERENCE);
  elNodes[34] = Point<double>( 1.428571428571428e-01 , 8.571428571428571e-01 , 0.000000000000000e+00 , FRAME_REFERENCE);
  elNodes[35] = Point<double>( 0.000000000000000e+00 , 1.000000000000000e+00 , 0.000000000000000e+00 , FRAME_REFERENCE);
  elNodes[36] = Point<double>( 0.000000000000000e+00 , 0.000000000000000e+00 , 1.428571428571428e-01 , FRAME_REFERENCE);
  elNodes[37] = Point<double>( 1.428571428571428e-01 , 0.000000000000000e+00 , 1.428571428571428e-01 , FRAME_REFERENCE);
  elNodes[38] = Point<double>( 2.857142857142857e-01 , 0.000000000000000e+00 , 1.428571428571428e-01 , FRAME_REFERENCE);
  elNodes[39] = Point<double>( 4.285714285714285e-01 , 0.000000000000000e+00 , 1.428571428571428e-01 , FRAME_REFERENCE);
  elNodes[40] = Point<double>( 5.714285714285714e-01 , 0.000000000000000e+00 , 1.428571428571428e-01 , FRAME_REFERENCE);
  elNodes[41] = Point<double>( 7.142857142857142e-01 , 0.000000000000000e+00 , 1.428571428571428e-01 , FRAME_REFERENCE);
  elNodes[42] = Point<double>( 8.571428571428571e-01 , 0.000000000000000e+00 , 1.428571428571428e-01 , FRAME_REFERENCE);
  elNodes[43] = Point<double>( 0.000000000000000e+00 , 1.428571428571428e-01 , 1.428571428571428e-01 , FRAME_REFERENCE);
  elNodes[44] = Point<double>( 1.428571428571428e-01 , 1.428571428571428e-01 , 1.428571428571428e-01 , FRAME_REFERENCE);
  elNodes[45] = Point<double>( 2.857142857142857e-01 , 1.428571428571428e-01 , 1.428571428571428e-01 , FRAME_REFERENCE);
  elNodes[46] = Point<double>( 4.285714285714285e-01 , 1.428571428571428e-01 , 1.428571428571428e-01 , FRAME_REFERENCE);
  elNodes[47] = Point<double>( 5.714285714285714e-01 , 1.428571428571428e-01 , 1.428571428571428e-01 , FRAME_REFERENCE);
  elNodes[48] = Point<double>( 7.142857142857142e-01 , 1.428571428571428e-01 , 1.428571428571428e-01 , FRAME_REFERENCE);
  elNodes[49] = Point<double>( 0.000000000000000e+00 , 2.857142857142857e-01 , 1.428571428571428e-01 , FRAME_REFERENCE);
  elNodes[50] = Point<double>( 1.428571428571428e-01 , 2.857142857142857e-01 , 1.428571428571428e-01 , FRAME_REFERENCE);
  elNodes[51] = Point<double>( 2.857142857142857e-01 , 2.857142857142857e-01 , 1.428571428571428e-01 , FRAME_REFERENCE);
  elNodes[52] = Point<double>( 4.285714285714285e-01 , 2.857142857142857e-01 , 1.428571428571428e-01 , FRAME_REFERENCE);
  elNodes[53] = Point<double>( 5.714285714285714e-01 , 2.857142857142857e-01 , 1.428571428571428e-01 , FRAME_REFERENCE);
  elNodes[54] = Point<double>( 0.000000000000000e+00 , 4.285714285714285e-01 , 1.428571428571428e-01 , FRAME_REFERENCE);
  elNodes[55] = Point<double>( 1.428571428571428e-01 , 4.285714285714285e-01 , 1.428571428571428e-01 , FRAME_REFERENCE);
  elNodes[56] = Point<double>( 2.857142857142857e-01 , 4.285714285714285e-01 , 1.428571428571428e-01 , FRAME_REFERENCE);
  elNodes[57] = Point<double>( 4.285714285714285e-01 , 4.285714285714285e-01 , 1.428571428571428e-01 , FRAME_REFERENCE);
  elNodes[58] = Point<double>( 0.000000000000000e+00 , 5.714285714285714e-01 , 1.428571428571428e-01 , FRAME_REFERENCE);
  elNodes[59] = Point<double>( 1.428571428571428e-01 , 5.714285714285714e-01 , 1.428571428571428e-01 , FRAME_REFERENCE);
  elNodes[60] = Point<double>( 2.857142857142857e-01 , 5.714285714285714e-01 , 1.428571428571428e-01 , FRAME_REFERENCE);
  elNodes[61] = Point<double>( 0.000000000000000e+00 , 7.142857142857142e-01 , 1.428571428571428e-01 , FRAME_REFERENCE);
  elNodes[62] = Point<double>( 1.428571428571428e-01 , 7.142857142857142e-01 , 1.428571428571428e-01 , FRAME_REFERENCE);
  elNodes[63] = Point<double>( 0.000000000000000e+00 , 8.571428571428571e-01 , 1.428571428571428e-01 , FRAME_REFERENCE);
  elNodes[64] = Point<double>( 0.000000000000000e+00 , 0.000000000000000e+00 , 2.857142857142857e-01 , FRAME_REFERENCE);
  elNodes[65] = Point<double>( 1.428571428571428e-01 , 0.000000000000000e+00 , 2.857142857142857e-01 , FRAME_REFERENCE);
  elNodes[66] = Point<double>( 2.857142857142857e-01 , 0.000000000000000e+00 , 2.857142857142857e-01 , FRAME_REFERENCE);
  elNodes[67] = Point<double>( 4.285714285714285e-01 , 0.000000000000000e+00 , 2.857142857142857e-01 , FRAME_REFERENCE);
  elNodes[68] = Point<double>( 5.714285714285714e-01 , 0.000000000000000e+00 , 2.857142857142857e-01 , FRAME_REFERENCE);
  elNodes[69] = Point<double>( 7.142857142857142e-01 , 0.000000000000000e+00 , 2.857142857142857e-01 , FRAME_REFERENCE);
  elNodes[70] = Point<double>( 0.000000000000000e+00 , 1.428571428571428e-01 , 2.857142857142857e-01 , FRAME_REFERENCE);
  elNodes[71] = Point<double>( 1.428571428571428e-01 , 1.428571428571428e-01 , 2.857142857142857e-01 , FRAME_REFERENCE);
  elNodes[72] = Point<double>( 2.857142857142857e-01 , 1.428571428571428e-01 , 2.857142857142857e-01 , FRAME_REFERENCE);
  elNodes[73] = Point<double>( 4.285714285714285e-01 , 1.428571428571428e-01 , 2.857142857142857e-01 , FRAME_REFERENCE);
  elNodes[74] = Point<double>( 5.714285714285714e-01 , 1.428571428571428e-01 , 2.857142857142857e-01 , FRAME_REFERENCE);
  elNodes[75] = Point<double>( 0.000000000000000e+00 , 2.857142857142857e-01 , 2.857142857142857e-01 , FRAME_REFERENCE);
  elNodes[76] = Point<double>( 1.428571428571428e-01 , 2.857142857142857e-01 , 2.857142857142857e-01 , FRAME_REFERENCE);
  elNodes[77] = Point<double>( 2.857142857142857e-01 , 2.857142857142857e-01 , 2.857142857142857e-01 , FRAME_REFERENCE);
  elNodes[78] = Point<double>( 4.285714285714285e-01 , 2.857142857142857e-01 , 2.857142857142857e-01 , FRAME_REFERENCE);
  elNodes[79] = Point<double>( 0.000000000000000e+00 , 4.285714285714285e-01 , 2.857142857142857e-01 , FRAME_REFERENCE);
  elNodes[80] = Point<double>( 1.428571428571428e-01 , 4.285714285714285e-01 , 2.857142857142857e-01 , FRAME_REFERENCE);
  elNodes[81] = Point<double>( 2.857142857142857e-01 , 4.285714285714285e-01 , 2.857142857142857e-01 , FRAME_REFERENCE);
  elNodes[82] = Point<double>( 0.000000000000000e+00 , 5.714285714285714e-01 , 2.857142857142857e-01 , FRAME_REFERENCE);
  elNodes[83] = Point<double>( 1.428571428571428e-01 , 5.714285714285714e-01 , 2.857142857142857e-01 , FRAME_REFERENCE);
  elNodes[84] = Point<double>( 0.000000000000000e+00 , 7.142857142857142e-01 , 2.857142857142857e-01 , FRAME_REFERENCE);
  elNodes[85] = Point<double>( 0.000000000000000e+00 , 0.000000000000000e+00 , 4.285714285714285e-01 , FRAME_REFERENCE);
  elNodes[86] = Point<double>( 1.428571428571428e-01 , 0.000000000000000e+00 , 4.285714285714285e-01 , FRAME_REFERENCE);
  elNodes[87] = Point<double>( 2.857142857142857e-01 , 0.000000000000000e+00 , 4.285714285714285e-01 , FRAME_REFERENCE);
  elNodes[88] = Point<double>( 4.285714285714285e-01 , 0.000000000000000e+00 , 4.285714285714285e-01 , FRAME_REFERENCE);
  elNodes[89] = Point<double>( 5.714285714285714e-01 , 0.000000000000000e+00 , 4.285714285714285e-01 , FRAME_REFERENCE);
  elNodes[90] = Point<double>( 0.000000000000000e+00 , 1.428571428571428e-01 , 4.285714285714285e-01 , FRAME_REFERENCE);
  elNodes[91] = Point<double>( 1.428571428571428e-01 , 1.428571428571428e-01 , 4.285714285714285e-01 , FRAME_REFERENCE);
  elNodes[92] = Point<double>( 2.857142857142857e-01 , 1.428571428571428e-01 , 4.285714285714285e-01 , FRAME_REFERENCE);
  elNodes[93] = Point<double>( 4.285714285714285e-01 , 1.428571428571428e-01 , 4.285714285714285e-01 , FRAME_REFERENCE);
  elNodes[94] = Point<double>( 0.000000000000000e+00 , 2.857142857142857e-01 , 4.285714285714285e-01 , FRAME_REFERENCE);
  elNodes[95] = Point<double>( 1.428571428571428e-01 , 2.857142857142857e-01 , 4.285714285714285e-01 , FRAME_REFERENCE);
  elNodes[96] = Point<double>( 2.857142857142857e-01 , 2.857142857142857e-01 , 4.285714285714285e-01 , FRAME_REFERENCE);
  elNodes[97] = Point<double>( 0.000000000000000e+00 , 4.285714285714285e-01 , 4.285714285714285e-01 , FRAME_REFERENCE);
  elNodes[98] = Point<double>( 1.428571428571428e-01 , 4.285714285714285e-01 , 4.285714285714285e-01 , FRAME_REFERENCE);
  elNodes[99] = Point<double>( 0.000000000000000e+00 , 5.714285714285714e-01 , 4.285714285714285e-01 , FRAME_REFERENCE);
  elNodes[100] = Point<double>( 0.000000000000000e+00 , 0.000000000000000e+00 , 5.714285714285714e-01 , FRAME_REFERENCE);
  elNodes[101] = Point<double>( 1.428571428571428e-01 , 0.000000000000000e+00 , 5.714285714285714e-01 , FRAME_REFERENCE);
  elNodes[102] = Point<double>( 2.857142857142857e-01 , 0.000000000000000e+00 , 5.714285714285714e-01 , FRAME_REFERENCE);
  elNodes[103] = Point<double>( 4.285714285714285e-01 , 0.000000000000000e+00 , 5.714285714285714e-01 , FRAME_REFERENCE);
  elNodes[104] = Point<double>( 0.000000000000000e+00 , 1.428571428571428e-01 , 5.714285714285714e-01 , FRAME_REFERENCE);
  elNodes[105] = Point<double>( 1.428571428571428e-01 , 1.428571428571428e-01 , 5.714285714285714e-01 , FRAME_REFERENCE);
  elNodes[106] = Point<double>( 2.857142857142857e-01 , 1.428571428571428e-01 , 5.714285714285714e-01 , FRAME_REFERENCE);
  elNodes[107] = Point<double>( 0.000000000000000e+00 , 2.857142857142857e-01 , 5.714285714285714e-01 , FRAME_REFERENCE);
  elNodes[108] = Point<double>( 1.428571428571428e-01 , 2.857142857142857e-01 , 5.714285714285714e-01 , FRAME_REFERENCE);
  elNodes[109] = Point<double>( 0.000000000000000e+00 , 4.285714285714285e-01 , 5.714285714285714e-01 , FRAME_REFERENCE);
  elNodes[110] = Point<double>( 0.000000000000000e+00 , 0.000000000000000e+00 , 7.142857142857142e-01 , FRAME_REFERENCE);
  elNodes[111] = Point<double>( 1.428571428571428e-01 , 0.000000000000000e+00 , 7.142857142857142e-01 , FRAME_REFERENCE);
  elNodes[112] = Point<double>( 2.857142857142857e-01 , 0.000000000000000e+00 , 7.142857142857142e-01 , FRAME_REFERENCE);
  elNodes[113] = Point<double>( 0.000000000000000e+00 , 1.428571428571428e-01 , 7.142857142857142e-01 , FRAME_REFERENCE);
  elNodes[114] = Point<double>( 1.428571428571428e-01 , 1.428571428571428e-01 , 7.142857142857142e-01 , FRAME_REFERENCE);
  elNodes[115] = Point<double>( 0.000000000000000e+00 , 2.857142857142857e-01 , 7.142857142857142e-01 , FRAME_REFERENCE);
  elNodes[116] = Point<double>( 0.000000000000000e+00 , 0.000000000000000e+00 , 8.571428571428571e-01 , FRAME_REFERENCE);
  elNodes[117] = Point<double>( 1.428571428571428e-01 , 0.000000000000000e+00 , 8.571428571428571e-01 , FRAME_REFERENCE);
  elNodes[118] = Point<double>( 0.000000000000000e+00 , 1.428571428571428e-01 , 8.571428571428571e-01 , FRAME_REFERENCE);
  elNodes[119] = Point<double>( 0.000000000000000e+00 , 0.000000000000000e+00 , 1.000000000000000e+00 , FRAME_REFERENCE);
#endif


  try{
    FieldContainer<double> vals;
    DefaultBasisFactory<double> BFactory;
    Teuchos::RCP<Basis<double> > triBasis = BFactory.create(FIELD_FORM_2, 
                                                            CELL_TET, 
                                                            RECONSTRUCTION_SPACE_INCOMPLETE, 
                                                            7, 
                                                            BASIS_FEM_FIAT, 
                                                            COORDINATES_CARTESIAN);
    // Exception 1: CURL cannot be applied to H(div) functions
    try {
      triBasis->getValues(vals, elNodes, OPERATOR_CURL);
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
      LocalDofTag myTag = {{4,0,0,1}};
      triBasis -> getLocalDofEnumeration(myTag);
    }
    catch (std::logic_error err) {
      *outStream << "Expected Error 2----------------------------------------------------------------\n";
      *outStream << err.what() << '\n';
      *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    };

    try {
      LocalDofTag myTag = {{1,1,21953,0}};
      triBasis -> getLocalDofEnumeration(myTag);
    }
    catch (std::logic_error err) {
      *outStream << "Expected Error 3----------------------------------------------------------------\n";
      *outStream << err.what() << '\n';
      *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    };
    try {
      LocalDofTag myTag = {{0,64,0,0}};
      triBasis -> getLocalDofEnumeration(myTag);
    }
    catch (std::logic_error err) {
      *outStream << "Expected Error 4----------------------------------------------------------------\n";
      *outStream << err.what() << '\n';
      *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    }; 

    try {
      int bfId = 280;
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
    Teuchos::RCP<Basis<double> > triBasis = BFactory.create(FIELD_FORM_2, 
                                                            CELL_TET, 
                                                            RECONSTRUCTION_SPACE_INCOMPLETE, 
                                                            7, 
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

  try{
    FieldContainer<double> vals;
    DefaultBasisFactory<double> BFactory;
    Teuchos::RCP<Basis<double> > triBasis = BFactory.create(FIELD_FORM_2, 
                                                            CELL_TET, 
                                                            RECONSTRUCTION_SPACE_INCOMPLETE, 
                                                            7, 
                                                            BASIS_FEM_FIAT, 
                                                            COORDINATES_CARTESIAN);

    // Check VALUE of basis functions
    triBasis -> getValues(vals, elNodes, OPERATOR_VALUE);
    ifstream VALfile("data/0.dat");
    if (VALfile.is_open()) {
      *outStream << "  Basis value: \n";
      for (int i=0;i<vals.getSize();i++) {
	double cur, fabscur;
	VALfile >> cur;
	fabscur = std::fabs(cur);
	if ( ( (fabscur < 1.885582e+05 * INTREPID_FIAT_TOL)
	       && (std::fabs(vals[i]) > 1.885582e+05 * INTREPID_FIAT_TOL ) ) 
	     || ( ( fabscur >= 1.885582e+05 * INTREPID_FIAT_TOL )
		  && ( std::fabs( vals[i] - cur ) / fabscur ) >= 1.885582e+05 * INTREPID_FIAT_TOL  ) ) {
	  errorFlag++;
	  *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
	  *outStream << std::setw(70) << " Wrong basis values" << "\n";
	  *outStream << i << " " << vals[i] << " " << cur << "\n";
	}
      }
      VALfile.close();
    }
    else {
      errorFlag = -999;
    }

    // Check DIV of basis functions
    triBasis -> getValues(vals, elNodes, OPERATOR_DIV);
    ifstream DIVfile("data/DIV.dat");
    if (DIVfile.is_open()) {
      *outStream << "  Basis divergences: \n";
      for (int i=0;i<vals.getSize();i++) {
	double cur, fabscur;
	DIVfile >> cur;
	fabscur = std::fabs(cur);
	if ( ( (fabscur < 3.640341e+06 * INTREPID_FIAT_TOL)
	       && (std::fabs(vals[i]) > 3.640341e+06 * INTREPID_FIAT_TOL ) ) 
	     || ( ( fabscur >= 3.640341e+06 * INTREPID_FIAT_TOL )
		  && ( std::fabs( vals[i] - cur ) / fabscur ) >= 3.640341e+06 * INTREPID_FIAT_TOL  ) ) {
	  errorFlag++;
	  *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
	  *outStream << std::setw(70) << " Wrong basis values" << "\n";
	  *outStream << i << " " << vals[i] << " " << cur << "\n";
	}
      }
      DIVfile.close();
    }
    else {
      errorFlag = -999;
    }
  }
  
  // Catch unexpected errors
  catch (std::logic_error err) {
    *outStream << err.what() << "\n\n";
    errorFlag = -1000;
  };

#endif

#ifndef HAVE_INTREPID_FIAT_VERY_HIGH_ORDER
  *outStream << "\nWarning: FIAT very high order disabled!\n";
#endif

  if (errorFlag != 0)
    std::cout << "End Result: TEST FAILED\n";
  else
    std::cout << "End Result: TEST PASSED\n";

  // reset format state of std::cout
  std::cout.copyfmt(oldFormatState);

  return errorFlag;
}

