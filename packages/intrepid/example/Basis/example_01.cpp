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
#include <iostream>

#include "Lagrange1d.hpp"

using namespace std;

int main(int argc, char *argv[]) {
  cout \
  << "===============================================================================\n" \
  << "|                                                                             |\n" \
  << "|                       Example use of the Lagrange class                     |\n" \
  << "|                                                                             |\n" \
  << "|    1) Creating Lagrange object in 1D                                        |\n" \
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

  /* cubic Lagrange polynomials on [0,1] */
  const int n = 3;
  vector<double> pts(n+1);
  vector<double> ys(n+1);
  vector<double> coeffs(n+1);

  /* get the equispaced points */
  Lagrange::equispacedPoints( n , 0.0 , 1.0 , pts );

  for (unsigned i=0;i<n+1;i++) {
    ys[i] = 0.0;
  }
  ys[2] = 1.0;

  Lagrange::dividedDifferences( pts , ys , coeffs );

  for (unsigned i=0;i<n+1;i++) {
    cout << coeffs[i] << " ";
  }
  cout << endl;

  cout << Lagrange::evaluateDividedDifferencePoly( pts , coeffs , pts[2] ) << endl;


  /* construct the divided difference table */
  Lagrange::Lagrange<double> U(pts);

  /* confirm the Kronecker delta */
  for (int bf = 0;bf<U.getDegree()+1;bf++) {
    cout << "Bf " << bf << endl;
    for (unsigned pt =0;pt<pts.size();pt++) {
      cout << "pt " << pt << " val = " << U.eval(bf,pts[pt]) << endl;
    }
  }


  return 0;
}
