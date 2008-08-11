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
\brief  Example of the CubatureDirect class.
\author Created by P. Bochev and D. Ridzal.
*/

#include "Intrepid_CubatureDirect.hpp"

using namespace std;
using namespace Intrepid;


/*
  Prints points and weights of a direct cubature rule to std::cout.
*/
void printInfo(ECell cellType, int cubDegree) {

  CubatureDirect<double> dCub(cellType, cubDegree);

  int ambientDim =  MultiCell<double>::getCellDim(cellType);

  int numCubPoints = dCub.getNumPoints();

  Teuchos::Array< Point<double> > cubPoints;
  Teuchos::Array<double> cubWeights;

  Point<double> tempPoint(ambientDim);
  cubPoints.assign(numCubPoints,tempPoint);
  cubWeights.assign(numCubPoints,0.0);

  dCub.getCubature(numCubPoints, cubPoints, cubWeights);

  cout << "Cell type:      " << MultiCell<double>::getCellName(cellType) << "\n";
  cout << "Degree:         " << cubDegree << "\n";
  cout << "Cubature name:  " << dCub.getName() << "\n";
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
  Computes volume of a given reference cell (compatible w/ direct cubature).
*/
double computeSimplexVolume(ECell cellType, int cubDegree) {

  CubatureDirect<double> dCub(cellType, cubDegree);

  int ambientDim =  MultiCell<double>::getCellDim(cellType);

  int numCubPoints = dCub.getNumPoints();

  Teuchos::Array< Point<double> > cubPoints;
  Teuchos::Array<double> cubWeights;

  Point<double> tempPoint(ambientDim);
  cubPoints.assign(numCubPoints,tempPoint);
  cubWeights.assign(numCubPoints,0.0);

  dCub.getCubature(numCubPoints, cubPoints, cubWeights);

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

  CubatureDirect<double> dCub(cellType, cubDegree);

  int ambientDim =  MultiCell<double>::getCellDim(cellType);

  int numCubPoints = dCub.getNumPoints();

  Teuchos::Array< Point<double> > cubPoints;
  Teuchos::Array<double> cubWeights;

  Point<double> tempPoint(ambientDim);
  cubPoints.assign(numCubPoints,tempPoint);
  cubWeights.assign(numCubPoints,0.0);

  dCub.getCubature(numCubPoints, cubPoints, cubWeights);

  double vol = 0.0;
  for (int i=0; i<numCubPoints; i++)
    vol += computeCupFunction(cubPoints[i])*cubWeights[i];

  return vol;
}



int main(int argc, char *argv[]) {
  cout \
  << "===============================================================================\n" \
  << "|                                                                             |\n" \
  << "|               Example use of the CubatureDirect class                       |\n" \
  << "|                                                                             |\n" \
  << "|     1) Creating direct cubature rules and extracting point and weight data  |\n" \
  << "|     2) Integrating functions over simplices                                 |\n" \
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
    ECell cellType = CELL_EDGE;
    int cubDegree  = 0;
    printInfo(cellType, cubDegree);

    cellType = CELL_EDGE;
    cubDegree  = INTREPID_MAX_CUBATURE_DEGREE_EDGE;
    printInfo(cellType, cubDegree);

    cellType = CELL_TRI;
    cubDegree  = 13;
    printInfo(cellType, cubDegree);

    cellType = CELL_TET;
    cubDegree  = 5;
    printInfo(cellType, cubDegree);
  }
  catch (std::logic_error err) {
    cout << err.what() << "\n";
  };

  cout \
  << "===============================================================================\n"\
  << "| EXAMPLE 2: integration of functions over simplices                          |\n"\
  << "===============================================================================\n";

  cout << "\nSimplex volumes:\n";

  try {
    cout << "EDGE volume w/ lowest  degree rule --> " << computeSimplexVolume(CELL_EDGE, 0) << '\n';
    cout << "EDGE volume w/ highest degree rule --> " << computeSimplexVolume(CELL_EDGE, INTREPID_MAX_CUBATURE_DEGREE_EDGE) << '\n';
    cout << "TRI  volume w/ lowest  degree rule --> " << computeSimplexVolume(CELL_TRI, 0) << '\n';
    cout << "TRI  volume w/ highest degree rule --> " << computeSimplexVolume(CELL_TRI, INTREPID_MAX_CUBATURE_DEGREE_TRI) << '\n';
    cout << "TET  volume w/ lowest  degree rule --> " << computeSimplexVolume(CELL_TET, 0) << '\n';
    cout << "TET  volume w/ highest degree rule --> " << computeSimplexVolume(CELL_TET, INTREPID_MAX_CUBATURE_DEGREE_TET) << '\n';
  }
  catch (std::logic_error err) {
    cout << err.what() << "\n";
  };

  cout << "\n1D, 2D, and 3D integrals of functions 1-x^2, 1-x^2-y^2, and 1-x^2-y^2-y^2, respectively:\n";

  try {
    cout << "EDGE integral w/ 1st     degree rule --> " << computeIntegral(CELL_EDGE, 1) << '\n';
    cout << "EDGE integral w/ 2nd     degree rule --> " << computeIntegral(CELL_EDGE, 2) << '\n';
    cout << "EDGE integral w/ highest degree rule --> " << computeIntegral(CELL_EDGE, INTREPID_MAX_CUBATURE_DEGREE_EDGE) << '\n';
    cout << "TRI  integral w/ 1st     degree rule --> " << computeIntegral(CELL_TRI, 1) << '\n';
    cout << "TRI  integral w/ 2nd     degree rule --> " << computeIntegral(CELL_TRI, 2) << '\n';
    cout << "TRI  integral w/ highest degree rule --> " << computeIntegral(CELL_TRI, INTREPID_MAX_CUBATURE_DEGREE_TRI) << '\n';
    cout << "TET  integral w/ 1st     degree rule --> " << computeIntegral(CELL_TET, 1) << '\n';
    cout << "TET  integral w/ 2nd     degree rule --> " << computeIntegral(CELL_TET, 2) << '\n';
    cout << "TET  integral w/ highest degree rule --> " << computeIntegral(CELL_TET, INTREPID_MAX_CUBATURE_DEGREE_TET) << '\n';
  }
  catch (std::logic_error err) {
    cout << err.what() << "\n";
  };

  return 0;
}
