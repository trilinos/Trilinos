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

/** \file   Intrepid_CubatureSparseDef.hpp
    \brief  Definition file for the Intrepid::CubatureSparse class.
    \author Created by P. Bochev, D. Ridzal, and M. Keegan.
*/


namespace Intrepid {

/**************************************************************************
**  Function Definitions for Class CubatureSparse
***************************************************************************/

template <class Scalar, int dimension_, class ArrayType>
CubatureSparse<Scalar,dimension_,ArrayType>::CubatureSparse(const int degree) :
  degree_(degree) {

  if(dimension_ == 2)
  {
    if(degree == 1)
      level_ = 1;
    else if(degree <= 3)
      level_ = 2;
    else if(degree <= 7)
      level_ = 3;
    else if(degree <= 11)
      level_ = 4;
    else if(degree <= 15)
      level_ = 5;
    else if(degree <= 19)
      level_ = 6;
    else if(degree <= 23)
      level_ = 7;
    else if(degree <= 27)
      level_ = 8;
    else if(degree <= 31)
      level_ = 9;
    else if(degree <= 35)
      level_ = 10;
    else if(degree <= 39)
      level_ = 11;
    else if(degree <= 43)
      level_ = 12;
    else if(degree <= 47)
      level_ = 13;
    else if(degree <= 51)
      level_ = 14;
    else if(degree <= 55)
      level_ = 15;
    else if(degree <= 59)
      level_ = 16;
  }
  else if(dimension_ == 3)
  {
    if(degree == 1)
      level_ = 1;
    else if(degree <= 3)
      level_ = 2;
    else if(degree <= 5)
      level_ = 3;
    else if(degree <= 9)
      level_ = 4;
    else if(degree <= 13)
      level_ = 5;
    else if(degree <= 17)
      level_ = 6;
    else if(degree <= 21)
      level_ = 7;
    else if(degree <= 25)
      level_ = 8;
    else if(degree <= 29)
      level_ = 9;
    else if(degree <= 33)
      level_ = 10;
    else if(degree <= 37)
      level_ = 11;
    else if(degree <= 41)
      level_ = 12;
    else if(degree <= 45)
      level_ = 13;
    else if(degree <= 49)
      level_ = 14;
    else if(degree <= 53)
      level_ = 15;
    else if(degree <= 57)
      level_ = 16;
  }

  numPoints_ = calculateNumPoints(dimension_,level_);
}



template <class Scalar, int dimension_, class ArrayType>
void CubatureSparse<Scalar,dimension_,ArrayType>::getCubature(ArrayType & cubPoints,
                                                              ArrayType & cubWeights) const{
  Teuchos::Array<Scalar> dummy_point(1);
  dummy_point[0] = 0.0;
  Scalar dummy_weight = 1.0;
  SGNodes<Scalar, dimension_> grid;

  iterateThroughDimensions(level_, dimension_, grid, dummy_point, dummy_weight);
  
  grid.copyToArrays(cubPoints, cubWeights);
} // end getCubature



template<class Scalar, int dimension_, class ArrayType>
int CubatureSparse<Scalar,dimension_,ArrayType>::getNumPoints() const {
  return numPoints_;
} // end getNumPoints



template<class Scalar, int dimension_, class ArrayType>
int CubatureSparse<Scalar,dimension_,ArrayType>::getDimension() const {
  return dimension_;
} // end dimension



template <class Scalar, int dimension_, class ArrayType>
void CubatureSparse<Scalar,dimension_,ArrayType>::getAccuracy(std::vector<int> & accuracy) const {
  accuracy.assign(1, degree_);
} //end getAccuracy



/************************************************************************
**     Function Definition for iterateThroughDimensions() 
**                 and its helper functions
*************************************************************************/

template< class Scalar, int DIM>
void iterateThroughDimensions(int level,
                              int dims_left,
                              SGNodes<Scalar,DIM> & cubPointsND,
                              Teuchos::Array<Scalar> & partial_node,
                              Scalar partial_weight)
{
  int l = level;
  int d = DIM;
  int add_on = d - dims_left;
  int start = dims_left > 1 ? 1 : (int)std::max(1, l);
  int end = l + add_on;

  for(int k_i = start; k_i <= end; k_i++)
  {
    /*******************
    **  Slow-Gauss
    ********************/
    int order1D = 2*k_i-1;

    /*******************
    **  Fast-Gauss
    ********************/
    //int order1D = (int)pow(2,k_i) - 1;

    int cubDegree1D = 2*order1D - 1;
    CubatureDirectLineGauss<Scalar> Cub1D(cubDegree1D);
    FieldContainer<Scalar> cubPoints1D(order1D, 1);
    FieldContainer<Scalar> cubWeights1D(order1D);

    Cub1D.getCubature(cubPoints1D, cubWeights1D);

    for(int node1D = 0; node1D < order1D; node1D++)
    {
      Teuchos::Array<Scalar> node(d-dims_left+1);
      node[d-dims_left] = cubPoints1D(node1D,0);
      for(int m = 0; m < d-dims_left; m++)
        node[m] = partial_node[m];

      Scalar weight = cubWeights1D(node1D)*partial_weight;

      if(dims_left != 1)
      {
        iterateThroughDimensions(l - k_i, dims_left-1, cubPointsND, node, weight);
      }
      else
      {
        weight = pow(-1.0, end - k_i)*combination(d-1, k_i - l)*weight;
        cubPointsND.addNode(&node[0], weight);
      }
    }
  }
}

} // end namespace Intrepid
