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
\brief  Example of the CubatureTensor class.
\author Created by P. Bochev and D. Ridzal.
*/

#include "Intrepid_CubatureTensor.hpp"

using namespace std;
using namespace Intrepid;


/*
  Prints points and weights of a tensor-product cubature rule to std::cout.
*/
void printInfo(ECell cellType, int cubDegree) {

  CubatureTensor<double> tensorCub(cellType, cubDegree);

  int ambientDim =  MultiCell<double>::getCellDim(cellType);

  int numCubPoints = tensorCub.getNumPoints();

  Teuchos::Array< Point<double> > cubPoints;
  Teuchos::Array<double> cubWeights;

  Point<double> tempPoint(ambientDim);
  cubPoints.assign(numCubPoints,tempPoint);
  cubWeights.assign(numCubPoints,0.0);

  tensorCub.getCubature(numCubPoints, cubPoints, cubWeights);

  cout << "Cell type:      " << MultiCell<double>::getCellName(cellType) << "\n";
  cout << "Degree:         " << cubDegree << "\n";
  cout << "Num. of points: " << numCubPoints << "\n";
  cout << "Cubature points:\n";
  for (int i=0; i<numCubPoints; i++) {
    cout << cubPoints[i] << endl;
  }
  cout << "Cubature weights:\n";
  for (int i=0; i<numCubPoints; i++) {
    cout << "  " << cubWeights[i] << endl;
  }

  cout << "\n";
}



/*
  Computes volume of a given reference cell (compatible w/ tensor-product cubature).
*/
double computeRefCellVolume(ECell cellType, int cubDegree) {

  CubatureTensor<double> tensorCub(cellType, cubDegree);

  int ambientDim =  MultiCell<double>::getCellDim(cellType);

  int numCubPoints = tensorCub.getNumPoints();

  Teuchos::Array< Point<double> > cubPoints;
  Teuchos::Array<double> cubWeights;

  Point<double> tempPoint(ambientDim);
  cubPoints.assign(numCubPoints,tempPoint);
  cubWeights.assign(numCubPoints,0.0);

  tensorCub.getCubature(numCubPoints, cubPoints, cubWeights);

  double vol = 0.0;
  for (int i=0; i<numCubPoints; i++)
    vol += cubWeights[i];

  return vol;
}


/*
  'Inverted cup' function evaluation.
    in 1D, for point p(x)    : 1-x^2
    in 2D, for point p(x,y)  : 1-x^2-y^2
    in 3D, for point p(x,y,z): 1-x^2-y^2-z^2
*/
double computeCupFunction(Point<double> p) {
  double val = 1.0;
  for (int i=0; i<p.getDim(); i++) {
    val -= p.getCoordinates()[i]*p.getCoordinates()[i];
  } 
  return val;
}


/*
  Compute integral of the inverted cup function over reference cells.
*/
double computeIntegral(ECell cellType, int cubDegree) {

  CubatureTensor<double> tensorCub(cellType, cubDegree);

  int ambientDim =  MultiCell<double>::getCellDim(cellType);

  int numCubPoints = tensorCub.getNumPoints();

  Teuchos::Array< Point<double> > cubPoints;
  Teuchos::Array<double> cubWeights;

  Point<double> tempPoint(ambientDim);
  cubPoints.assign(numCubPoints,tempPoint);
  cubWeights.assign(numCubPoints,0.0);

  tensorCub.getCubature(numCubPoints, cubPoints, cubWeights);

  double vol = 0.0;
  for (int i=0; i<numCubPoints; i++)
    vol += computeCupFunction(cubPoints[i])*cubWeights[i];

  return vol;
}



int main(int argc, char *argv[]) {
  cout \
  << "===============================================================================\n" \
  << "|                                                                             |\n" \
  << "|               Example use of the CubatureTensor class                       |\n" \
  << "|                                                                             |\n" \
  << "|     1) Creating tensor cubature rules and extracting point and weight data  |\n" \
  << "|     2) Integrating functions over quads, hexes, and triprisms               |\n" \
  << "|                                                                             |\n" \
  << "|  Questions? Contact  Pavel Bochev (pbboche@sandia.gov) or                   |\n" \
  << "|                      Denis Ridzal (dridzal@sandia.gov).                     |\n" \
  << "|                                                                             |\n" \
  << "|  Intrepid's website: http://trilinos.sandia.gov/packages/intrepid           |\n" \
  << "|  Trilinos website:   http://trilinos.sandia.gov                             |\n" \
  << "|                                                                             |\n" \
  << "===============================================================================\n"\
  << "| EXAMPLE 1: object creation, point & weight data extraction                  |\n"\
  << "===============================================================================\n";

  try {
    ECell cellType = CELL_QUAD;
    int cubDegree  = 0;
    printInfo(cellType, cubDegree);

    cellType = CELL_QUAD;
    cubDegree  = 5;
    printInfo(cellType, cubDegree);

    cellType = CELL_HEX;
    cubDegree  = 0;
    printInfo(cellType, cubDegree);

    cellType = CELL_HEX;
    cubDegree  = 5;
    printInfo(cellType, cubDegree);

    cellType = CELL_TRIPRISM;
    cubDegree  = 0;
    printInfo(cellType, cubDegree);

    cellType = CELL_TRIPRISM;
    cubDegree  = 5;
    printInfo(cellType, cubDegree);
  }
  catch (std::logic_error err) {
    cout << err.what() << "\n";
  };

  cout \
  << "===============================================================================\n"\
  << "| EXAMPLE 2: integration of functions over quads, hexes, and triprisms        |\n"\
  << "===============================================================================\n";

  cout << "\nReference cell volumes:\n";

  try {
    cout << "QUAD     volume w/ lowest  degree rule --> " << computeRefCellVolume(CELL_QUAD, 0) << '\n';
    cout << "QUAD     volume w/ highest degree rule --> " << computeRefCellVolume(CELL_QUAD, INTREPID_MAX_CUBATURE_DEGREE_EDGE) << '\n';
    cout << "HEX      volume w/ lowest  degree rule --> " << computeRefCellVolume(CELL_HEX, 0) << '\n';
    cout << "HEX      volume w/ highest degree rule --> " << computeRefCellVolume(CELL_HEX, INTREPID_MAX_CUBATURE_DEGREE_EDGE) << '\n';
    cout << "TRIPRISM volume w/ lowest  degree rule --> " << computeRefCellVolume(CELL_TRIPRISM, 0) << '\n';
    cout << "TRIPRISM volume w/ highest degree rule --> " << computeRefCellVolume(CELL_TRIPRISM, INTREPID_MAX_CUBATURE_DEGREE_EDGE) << '\n';
  }
  catch (std::logic_error err) {
    cout << err.what() << "\n";
  };

  cout << "\n3D and 4D integrals of functions 1-x^2-y^2 and  1-x^2-y^2-y^2, respectively:\n";

  try {
    cout << "QUAD     integral w/ 1st     degree rule --> " << computeIntegral(CELL_QUAD, 1) << '\n';
    cout << "QUAD     integral w/ 2nd     degree rule --> " << computeIntegral(CELL_QUAD, 2) << '\n';
    cout << "QUAD     integral w/ highest degree rule --> " << computeIntegral(CELL_QUAD, INTREPID_MAX_CUBATURE_DEGREE_EDGE) << '\n';
    cout << "HEX      integral w/ 1st     degree rule --> " << computeIntegral(CELL_HEX, 1) << '\n';
    cout << "HEX      integral w/ 2nd     degree rule --> " << computeIntegral(CELL_HEX, 2) << '\n';
    cout << "HEX      integral w/ highest degree rule --> " << computeIntegral(CELL_HEX, INTREPID_MAX_CUBATURE_DEGREE_EDGE) << '\n';
    cout << "TRIPRISM integral w/ 1st     degree rule --> " << computeIntegral(CELL_TRIPRISM, 1) << '\n';
    cout << "TRIPRISM integral w/ 2nd     degree rule --> " << computeIntegral(CELL_TRIPRISM, 2) << '\n';
    cout << "TRIPRISM integral w/ highest degree rule --> " <<
             computeIntegral(CELL_TRIPRISM, std::min(INTREPID_MAX_CUBATURE_DEGREE_TRI,INTREPID_MAX_CUBATURE_DEGREE_EDGE)) << '\n';
  }
  catch (std::logic_error err) {
    cout << err.what() << "\n";
  };

  return 0;
}
