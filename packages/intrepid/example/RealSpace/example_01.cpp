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
//                    Denis Ridzal (dridzal@sandia.gov).
//
// ************************************************************************
// @HEADER

/** \file
\brief  Illustrates use of the Point class.
\author Created by P. Bochev and D. Ridzal
*/

#include "Intrepid_RealSpace.hpp"

using namespace std;
using namespace Intrepid;

int main(int argc, char *argv[]) {
  cout \
  << "===============================================================================\n" \
  << "|                                                                             |\n" \
  << "|                       Example use of the Point class                        |\n" \
  << "|                                                                             |\n" \
  << "|  Questions? Contact  Pavel Bochev (pbboche@sandia.gov) or                   |\n" \
  << "|                      Denis Ridzal (dridzal@sandia.gov).                     |\n" \
  << "|                                                                             |\n" \
  << "|  Intrepid's website: http://trilinos.sandia.gov/packages/intrepid           |\n" \
  << "|  Trilinos website:   http://trilinos.sandia.gov                             |\n" \
  << "|                                                                             |\n" \
  << "===============================================================================\n\n";
  
  cout.precision(16);
  cout \
  << "\tThis machine's epsilon from Teuchos =" << INTREPID_EPSILON <<endl \
  << "\tINTREPID_THRESHOLD = " << INTREPID_THRESHOLD << endl \
  << "\tINTREPID_MINUS_ONE = " << INTREPID_MINUS_ONE << endl \
  << "\tINTREPID_PLUS_ONE  = " << INTREPID_PLUS_ONE  << endl << endl \
  << "===============================================================================\n"\
  << "| EXAMPLE 1: class Point in 3D                                                |\n"\
  << "===============================================================================\n\n";
  // Create arrays of coefficients
  double vec[]  = {1.0, 2.0, 3.0};
  double vec2[] = {-4.0, -1.0, 10.0};
  
  // Using constructor that takes pointer, dimension and uses default for the frame kind.
  Point<double> p1(vec, 3);
  cout << "\tCreated point:  p1 = " << p1 << endl;
  Point<double> p2(vec2, 3);
  cout << "\tCreated point:  p2 = " << p2 << endl;
  
  cout << "Examples of operations between points in 3D:\n";
  cout << "\t p1 + p2           = " << p1+p2 << endl;
  cout << "\t p1 - p2           = " << p1-p2 << endl;
  cout << "\t p1 . p2           = " << p1*p2 << endl;
  cout << "\t p1 x p2           = " << (p1^p2) << endl;
  cout << "\t 5.0 * p1          = " << 5.0*p1 << endl;
  cout << "\t 1/(p1.p2)*(p1xp2) = " << (1.0 / (p1*p2)) * (p1 ^ p2) << endl;
  cout << "Computing distance between points:"<< endl;
  cout << "\t |p1-p2| = " << p1.distance(p2) << endl;
  cout << "\t |p2-p1| = " << p2.distance(p1) << endl;
  cout << "\t |p1-p1| = " << p1.distance(p1) << "\n";
  cout << "Computing norms of points:"<< endl;
  cout << "\t || p1 ||_2  = " << p1.norm(NORM_TWO) \
    << "\t || p2 ||_2  = " << p2.norm(NORM_TWO)<< endl;
  cout << "\t || p1 ||_1  = " << p1.norm(NORM_ONE) \
    << "\t || p2 ||_1  = " << p2.norm(NORM_ONE)<< endl;
  cout << "\t || p1 ||_I  = " << p1.norm(NORM_INF) \
    << "\t || p2 ||_I  = " << p2.norm(NORM_INF)<< endl;
  //
  cout << "\n" \
  << "===============================================================================\n"\
  << "| EXAMPLE 2: class Point in 2D                                                |\n"\
  << "===============================================================================\n\n";
  
  // Using constructor that takes pointer, dimension and uses default for the frame kind.
  Point<double> p4(vec, 2);
  cout << "\tCreated point:  p4 = " << p4 << endl;
  Point<double> p5(vec2, 2);
  cout << "\tCreated point:  p5 = " << p5 << endl;
  
  cout << "Examples of operations between points in 2D:\n";
  cout << "\t p4 + p5           = " << p4+p5 << endl;
  cout << "\t p4 - p5           = " << p4-p5 << endl;
  cout << "\t p4 . p5           = " << p4*p5 << endl;
  cout << "\t 5.0 * p4          = " << 5.0*p4 << endl;
  cout << "\t 1/(p4.p5)*(p4+p5) = " << (1.0 / (p4*p5)) * (p4 + p5) << endl;
  cout << "Computing distance between points:"<< endl;
  cout << "\t |p4-p5| = " << p4.distance(p5) << endl;
  cout << "\t |p5-p4| = " << p5.distance(p4) << endl;
  cout << "\t |p4-p4| = " << p4.distance(p4) << endl;
  cout << "Computing norms of 2D points:"<< endl;
  cout << "\t || p4 ||_2  = " << p4.norm(NORM_TWO) \
    << "\t || p5 ||_2  = " << p5.norm(NORM_TWO)<< endl;
  cout << "\t || p4 ||_1  = " << p4.norm(NORM_ONE) \
    << "\t || p5 ||_1  = " << p5.norm(NORM_ONE)<< endl;
  cout << "\t || p4 ||_I  = " << p4.norm(NORM_INF) \
    << "\t || p5 ||_I  = " << p5.norm(NORM_INF)<< endl;
  //
  cout << "\n" \
  << "===============================================================================\n"\
  << "| EXAMPLE 3: class Point in 2D                                                |\n"\
  << "===============================================================================\n\n";
  
  // Using constructor that takes pointer, dimension and uses default for the frame kind.
  Point<double> p7(vec, 1);
  cout << "\tCreated point:  p7 = " << p7 << endl;
  Point<double> p8(vec2, 1);
  cout << "\tCreated point:  p8 = " << p8 << endl;
  cout << "Examples of operations between points in 2D:\n";
  cout << "\t p7 + p8           = " << p7+p8 << endl;
  cout << "\t p7 - p8           = " << p7-p8 << endl;
  cout << "\t p7 . p8           = " << p7*p8 << endl;
  cout << "\t 5.0 * p7          = " << 5.0*p7 << endl;
  cout << "\t 1/(p7.p8)*(p7+p8) = " << (1.0 / (p7*p8)) * (p7 + p8) << endl;
  cout << "Computing distance between points in 1D:"<< endl;
  cout << "\t |p7-p8| = " << p7.distance(p8) << endl;
  cout << "\t |p8-p7| = " << p8.distance(p7) << endl;
  cout << "\t |p7-p7| = " << p7.distance(p7) << endl;
  cout << "Computing norms of 1D points:"<< endl;
  cout << "\t || p7 ||_2  = " << p7.norm(NORM_TWO) \
    << "\t || p8 ||_2  = " << p8.norm(NORM_TWO)<< endl;
  cout << "\t || p7 ||_1  = " << p7.norm(NORM_ONE) \
    << "\t || p8 ||_1  = " << p8.norm(NORM_ONE)<< endl;
  cout << "\t || p7 ||_I  = " << p7.norm(NORM_INF) \
    << "\t || p8 ||_I  = " << p8.norm(NORM_INF)<< endl;  
  //
  cout << "\n" \
  << "===============================================================================\n"\
  << "| EXAMPLE 4: other operations on Point objects                                |\n"\
  << "===============================================================================\n\n";
  //
  cout << "Changing the frame kind of a point: " << endl;
  cout << "\t Original Point p1 is : " << p1 << endl;
  p1.setFrameKind(FRAME_REFERENCE);
  cout << "\t      Now Point p1 is : " << p1 << endl;
  
  cout << "Checking if a Point belongs to a reference cell: " << endl;
  
  // Points below illustrate using constructors that take 1,2 or 3 point coordinates 
  // and the user sets the frame kind explicitely to override the default (FRAME_PHYSICAL).

  // 1D point that is close to the right endpoint of the reference edge cell
  Point<double> p_in_edge(1.0-INTREPID_EPSILON,FRAME_REFERENCE);
  
  // 2D point that is close to the top right corner of the reference quad cell
  Point<double> p_in_quad(1.0,1.0-INTREPID_EPSILON,FRAME_REFERENCE);
  
  // 2D point that is close to the midpoint of the slanted edge of the reference tri cell
  Point<double> p_in_tri(0.5-INTREPID_EPSILON,0.5-INTREPID_EPSILON,FRAME_REFERENCE);
  
  // 3D point that is close to 0th vertex of the reference hex cell
  Point<double> p_in_hex(1.0-INTREPID_EPSILON,1.0-INTREPID_EPSILON,1.0-INTREPID_EPSILON,FRAME_REFERENCE);
  
  // 3D point that is close to the slanted face of the reference tet cell
  Point<double> p_in_tet(0.5-INTREPID_EPSILON,0.5-INTREPID_EPSILON,0.5-INTREPID_EPSILON,FRAME_REFERENCE);
  
  // 3D point close to the top face of the reference prism 
  Point<double> p_in_prism(0.5,0.25,1.0-INTREPID_EPSILON,FRAME_REFERENCE);
  
  // 3D point close to the top of the reference pyramid
  Point<double> p_in_pyramid(-INTREPID_EPSILON,INTREPID_EPSILON,(1.0-INTREPID_EPSILON),FRAME_REFERENCE);
  
  // Check if the points are in their respective reference cells
  EFailCode in_edge    = p_in_edge.inclusion(CELL_EDGE);
  EFailCode in_tri     = p_in_tri.inclusion(CELL_TRI);
  EFailCode in_quad    = p_in_quad.inclusion(CELL_QUAD);
  EFailCode in_tet     = p_in_tet.inclusion(CELL_TET);
  EFailCode in_hex     = p_in_hex.inclusion(CELL_HEX);
  EFailCode in_prism   = p_in_prism.inclusion(CELL_TRIPRISM);
  EFailCode in_pyramid = p_in_pyramid.inclusion(CELL_PYRAMID);
  //
  if(in_edge == FAIL_CODE_SUCCESS) {
    cout <<  p_in_edge << " is inside reference edge " << endl;
  }
  if(in_tri == FAIL_CODE_SUCCESS) {
    cout << p_in_tri << " is inside reference triangle " << endl;
  }
  if(in_quad == FAIL_CODE_SUCCESS) {
    cout << p_in_quad << " is inside reference quad " << endl;
  }
  if(in_tet == FAIL_CODE_SUCCESS) {
    cout << p_in_tet << " is inside reference tet " << endl;
  }
  if(in_hex == FAIL_CODE_SUCCESS) {
    cout << p_in_hex << " is inside reference hex " << endl;
  }
  if(in_prism == FAIL_CODE_SUCCESS) {
    cout << p_in_prism << " is inside reference prism " << endl;
  }
  if(in_pyramid == FAIL_CODE_SUCCESS) {
    cout << p_in_pyramid << " is inside reference pyramid " << endl;
  }
  
  // Now make 1,2 and 3D points with very small coefficients, but larger than threshold
  double small = 2.0*INTREPID_THRESHOLD;
  Point<double> p_eps_1D(small,FRAME_REFERENCE);
  Point<double> p_eps_2D(small,small,FRAME_REFERENCE);
  Point<double> p_eps_3D(small,small,small,FRAME_REFERENCE);
  
  // Add these points to the good reference points above:
  cout << "Adding small perturbations to these points..." << endl;
  p_in_edge    += p_eps_1D;
  p_in_tri     += p_eps_2D;
  p_in_quad    += p_eps_2D;
  p_in_tet     += p_eps_3D;
  p_in_hex     += p_eps_3D;
  p_in_prism   += p_eps_3D;
  p_in_pyramid += p_eps_3D;
  
  // Now check again if the points are in their respective reference cells.
  cout << "Checking if the perturbed Points belongs to reference cell: " << endl;
  in_edge    = p_in_edge.inclusion(CELL_EDGE);
  in_tri     = p_in_tri.inclusion(CELL_TRI);
  in_quad    = p_in_quad.inclusion(CELL_QUAD);
  in_tet     = p_in_tet.inclusion(CELL_TET);
  in_hex     = p_in_hex.inclusion(CELL_HEX);
  in_prism   = p_in_prism.inclusion(CELL_TRIPRISM);
  in_pyramid = p_in_pyramid.inclusion(CELL_PYRAMID);
  
  //
  if(in_edge == FAIL_CODE_NOT_IN_REF_CELL) {
    cout <<  p_in_edge << " is NOT inside reference edge " << endl;
  }
  if(in_tri == FAIL_CODE_NOT_IN_REF_CELL) {
    cout << p_in_tri << " is NOT inside reference triangle " << endl;
  }
  if(in_quad == FAIL_CODE_NOT_IN_REF_CELL) {
    cout << p_in_quad << " is NOT inside reference quad " << endl;
  }
  if(in_tet == FAIL_CODE_NOT_IN_REF_CELL) {
    cout << p_in_tet << " is NOT inside reference tet " << endl;
  }
  if(in_hex == FAIL_CODE_NOT_IN_REF_CELL) {
    cout << p_in_hex << " is NOT inside reference hex " << endl;
  }
  if(in_prism == FAIL_CODE_NOT_IN_REF_CELL) {
    cout << p_in_prism << " is NOT inside reference prism " << endl;
  }
  if(in_pyramid == FAIL_CODE_NOT_IN_REF_CELL) {
    cout << p_in_pyramid << " is NOT inside reference pyramid " << endl;
  }
  cout << "You need at least 16 digits to see the added pertirbation!" << endl;
  return 0;
}
