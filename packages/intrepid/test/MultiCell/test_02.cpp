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
\brief  Example of the Multicell class.
\author Created by P. Bochev and D. Ridzal
*/
#include "Intrepid_MultiCell.hpp"

#define MY_HEXAGON   CELL_POLY0
#define HONEYCOMB    CELL_POLY1

using namespace std;
using namespace Intrepid;


int main(int argc, char *argv[]) {

  cout \
  << "===============================================================================\n" \
  << "|                                                                             |\n" \
  << "|                   Example use of the MultiCell class                        |\n" \
  << "|                                                                             |\n" \
  << "|  1) Creating and saving custom cell templates                                |\n" \
  << "|                                                                             |\n" \
  << "|  Questions? Contact  Pavel Bochev (pbboche@sandia.gov) or                   |\n" \
  << "|                      Denis Ridzal (dridzal@sandia.gov).                     |\n" \
  << "|                                                                             |\n" \
  << "|  Intrepid's website: http://trilinos.sandia.gov/packages/intrepid           |\n" \
  << "|  Trilinos website:   http://trilinos.sandia.gov                             |\n" \
  << "|                                                                             |\n" \
  << "===============================================================================\n"\
  << "| EXAMPLE 1: class MultiCell in 3D / custom cell shapes                       |\n"\
  << "|       We are creating a 3D beehive multicell consisting                     |\n"\
  << "|       of two cells (prisms with hexagonal base).                            |\n"\
  << "===============================================================================\n";
 
  // first define a hexagonal cell, which is a face of the honeycomb cell
  ConnMapTemplate hexagon[3] =
  {   // hexagon
    { // hexagon->1cell
      2,
      6,
      {2,2,2,2,2,2},
      {CELL_EDGE,CELL_EDGE,CELL_EDGE,CELL_EDGE,CELL_EDGE,CELL_EDGE},
      {{0,1}, {1,2}, {2,3}, {3,4}, {4,5}, {5,0}}
    },
    { // hexagon->2cell
      2,
      1,
      {6},
      {MY_HEXAGON},
      {{0,1,2,3,4,5}}
    },
    { // hexagon->3cell
      2,
      0,
      {0},
      {},
      {{0}}
    }
  };  // end hexagon

  // call static setConnMapCustom function to store the hexagon template
  MultiCell<double>::setConnMapCustom(MY_HEXAGON, hexagon);

  // define a honeycomb cell
  ConnMapTemplate honeycomb[3] =
  {   // honeycomb cell
    { // honeycomb->1cell
      3,
      18,
      {2,2,2,2,2,2},
      {CELL_EDGE,CELL_EDGE,CELL_EDGE,CELL_EDGE,CELL_EDGE,CELL_EDGE,
       CELL_EDGE,CELL_EDGE,CELL_EDGE,CELL_EDGE,CELL_EDGE,CELL_EDGE,
       CELL_EDGE,CELL_EDGE,CELL_EDGE,CELL_EDGE,CELL_EDGE,CELL_EDGE},
      {{0,1}, {1,2}, {2,3}, {3,4}, {4,5}, {5,0},
       {0,6}, {1,7}, {2,8}, {3,9}, {4,10}, {5,11},
       {6,7}, {7,8}, {8,9}, {9,10}, {10,11}, {11,6}}
    },
    { // honeycomb->2cell
      3,
      8,
      {4,4,4,4,4,4,6,6},
      // note that the previously defined hexagon template is used below
      {CELL_QUAD,CELL_QUAD,CELL_QUAD,CELL_QUAD,CELL_QUAD,CELL_QUAD,CELL_HEXAGON,CELL_HEXAGON},
      {{0,1,7,6},{1,2,8,7},{2,3,9,8},{3,4,10,9},{4,5,11,10},{5,0,6,11},{0,1,2,3,4,5},{6,7,8,9,10,11}}
    },
    { // honeycomb->3cell
      3,
      1,
      {12},
      {HONEYCOMB},
      {{0,1,2,3,4,5,6,7,8,9,10,11}}
    }
  };  // end honeycomb

  // call static setConnMapCustom function to store the honeycomb template
  MultiCell<double>::setConnMapCustom(HONEYCOMB, honeycomb);

  double hivenodes[] = {1.1, 1.2, 1.3,
                      2.1, 2.2, 2.3,
                      3.1, 3.2, 3.3,
                      4.1, 4.2, 4.3,
                      5.1, 5.2, 5.3,
                      6.1, 6.2, 6.3,
                      7.1, 7.2, 7.3,
                      8.1, 8.2, 8.3,
                      9.1, 9.2, 9.3,
                      10.1, 10.2, 10.3,
                      11.1, 11.2, 11.3,
                      12.1, 12.2, 12.3, //
                      3.1, 3.2, 3.3,
                      13.1, 13.2, 13.3,
                      14.1, 14.2, 14.3,
                      15.1, 15.2, 15.3,
                      16.1, 16.2, 16.3,
                      4.1, 4.2, 4.3,
                      9.1, 9.2, 9.3,
                      17.1, 17.2, 17.3,
                      18.1, 18.2, 18.3,
                      19.1, 19.2, 19.3,
                      20.1, 20.2, 20.3,
                      10.1, 10.2, 10.3};

  MultiCell<double> beehive(2, HONEYCOMB, hivenodes);

  cout << beehive << endl;

  cout << "Accessing node coordinates ...\n";
  cout << "beehive[0][10] =\n" << beehive[0][10] << "\n";
  cout << "beehive[1][5] =\n" << beehive[1][5] << "\n\n";

  cout << "END TEST 1: class MultiCell in 3D / custom cell shapes\n";

  return 0;
}
