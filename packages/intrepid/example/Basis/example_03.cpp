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
//                    Robert Kirby (robert.c.kirby@ttu.edu)
//
// ************************************************************************
// @HEADER

/** \file
\brief  Illustrates use of the new Lagrange class on quads.
\author Created by R. Kirby
*/
#include <iostream>

#include "Intrepid_F0_HEX_DD.hpp"
#include "Intrepid_RealSpace.hpp"
#include "Intrepid_FieldContainer.hpp"
#include "Intrepid_Types.hpp"
#include "Lagrange1d.hpp"

using namespace std;
using namespace Intrepid;

int main(int argc, char *argv[]) {
  cout \
  << "===============================================================================\n" \
  << "|                                                                             |\n" \
  << "|                       Example use of the Lagrange class                     |\n" \
  << "|                                                                             |\n" \
  << "|    1) Creating HEX Lagrange elements                                        |\n" \
  << "|                                                                             |\n" \
  << "|  Questions? Contact  Pavel Bochev (pbboche@sandia.gov) or                   |\n" \
  << "|                      Denis Ridzal (dridzal@sandia.gov) or                   |\n" \
  << "|                      Robert Kirby (robert.c.kirby@ttu.edu)                  |\n" \
  << "|                                                                             |\n" \
  << "|  Intrepid's website: http://trilinos.sandia.gov/packages/intrepid           |\n" \
  << "|  Trilinos website:   http://trilinos.sandia.gov                             |\n" \
  << "|                                                                             |\n" \
  << "===============================================================================\n\n";
  
  cout.precision(16);


  Basis_F0_HEX_DD<double> myBasis( 1 );

  // need to initialize before evaluation
  myBasis.initialize();

  vector<double> pts_vec(3);
  double vec[3];
  Lagrange::equispacedPoints<double>(2,-1.0,1.0,pts_vec);

  Teuchos::Array<Point<double> > pts;
  for (unsigned i=0;i<pts_vec.size();i++) {
    for (unsigned j=0;j<pts_vec.size();j++) {
      for (unsigned k=0;k<pts_vec.size();k++) {
	vec[0] = pts_vec[i];
	vec[1] = pts_vec[j];
	vec[2] = pts_vec[k];
	Point<double> pt(vec,3,FRAME_REFERENCE);
	pts.push_back( pt );
      }
    }
  }
  
  FieldContainer<double> myContainer;


  myBasis.getValues(myContainer,pts,OPERATOR_VALUE);
  cout << "Values:" << endl;
  cout << pts << endl;
  cout << myContainer << endl;
  
  myBasis.getValues(myContainer,pts,OPERATOR_GRAD); 
  cout << "Gradients" << endl;
  cout << pts << endl;
  cout << myContainer << endl;


  return 0;
}
