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

/** \file   Intrepid_CubatureSparse.hpp
    \brief  Header file for the Intrepid::CubatureSparse class.
    \author Created by P. Bochev, D. Ridzal, and M. Keegan.
*/

#ifndef INTREPID_CUBATURE_SPARSE_HPP
#define INTREPID_CUBATURE_SPARSE_HPP

#include "Intrepid_ConfigDefs.hpp"
#include "Intrepid_Cubature.hpp"
#include "Intrepid_CubatureDirectLineGauss.hpp"
#include "Intrepid_CubatureSparseHelper.hpp"
#include "Teuchos_TestForException.hpp"


/** \def INTREPID_CUBATURE_SPARSE2D_GAUSS_MAX
  \brief The maximum degree of the polynomial that can be integrated exactly by
         a sparse rule of the Gauss(-Legendre) type in 2D.
*/
#define INTREPID_CUBATURE_SPARSE2D_GAUSS_MAX 59

/** \def INTREPID_CUBATURE_SPARSE3D_GAUSS_MAX
  \brief The maximum degree of the polynomial that can be integrated exactly by
         a sparse of the Gauss(-Legendre) type in 3D.
*/
#define INTREPID_CUBATURE_SPARSE3D_GAUSS_MAX 57


namespace Intrepid{

template<class Scalar, int dimension_, class ArrayPoint = FieldContainer<Scalar>, class ArrayWeight = ArrayPoint>
class CubatureSparse : public Intrepid::Cubature<Scalar,ArrayPoint,ArrayWeight> {
  private:

  int level_;

  int numPoints_;

  const int degree_;

  
  public:

  ~CubatureSparse() {}


  CubatureSparse(const int degree);

  /** \brief Returns cubature points and weights
             (return arrays must be pre-sized/pre-allocated).

      \param cubPoints       [out]     - Array containing the cubature points.
      \param cubWeights      [out]     - Array of corresponding cubature weights.
  */
  virtual void getCubature(ArrayPoint  & cubPoints,
                           ArrayWeight & cubWeights) const;

  /** \brief Returns the number of cubature points.
  */
  virtual int getNumPoints() const;

  /** \brief Returns dimension of the integration domain.
  */
  virtual int getDimension() const;

  /** \brief Returns algebraic accuracy (e.g. max. degree of polynomial
             that is integrated exactly).
  */
  virtual void getAccuracy(std::vector<int> & accuracy) const;

}; // end class CubatureSparse 

// helper functions

template<class Scalar, int DIM>
void iterateThroughDimensions(int level,
                              int dims_left,
                              SGNodes<Scalar,DIM> & cubPointsND,
                              Teuchos::Array<Scalar> & partial_node,
                              Scalar partial_weight);

inline int factorial(int num)
{
  int answer = 1;
  if(num >= 1)
  {
    while(num > 0)
    {
      answer = answer*num;
      num--;
    }
  }
  else if(num == 0)
    answer = 1;
  else
    answer = -1;

  return answer;
}

inline double combination(int top, int bot)
{
  double answer = factorial(top)/(factorial(bot) * factorial(top-bot));
  return answer;
}

inline int iterateThroughDimensionsForNumCalc(int dims_left,
                                              int level,
                                              int levels_left,
                                              int level_so_far,
                                              Teuchos::Array<int> & nodes,
                                              int product,
                                              bool no_uni_quad)
{
  int numNodes = 0;
  for(int j = 1; j <= levels_left; j++)
  {
    bool temp_bool = no_uni_quad;
    int temp_knots = nodes[j-1]*product;
    int temp_lsf = level_so_far + j;

    if(j==1)
      temp_bool = false;

    if(dims_left == 1)
    {
      if(temp_lsf < level && temp_bool == true)
        numNodes += 0;
      else
      {
        numNodes += temp_knots;
      }
    }
    else
    {
      numNodes += iterateThroughDimensionsForNumCalc(dims_left-1,level, levels_left-j+1, temp_lsf, nodes, temp_knots, temp_bool);
    }
  }
  return numNodes;
}

inline int calculateNumPoints(int dim, int level)
{
  //int* uninum = new int[level];
  Teuchos::Array<int> uninum(level);
  uninum[0] = 1;
  for(int i = 1; i <= level-1; i++)
  {
    uninum[i] = 2*i;
  }

  int numOfNodes = iterateThroughDimensionsForNumCalc(dim, level, level, 0, uninum, 1, true);
  return numOfNodes;
}


} // end namespace Intrepid


// include templated definitions
#include <Intrepid_CubatureSparseDef.hpp>

#endif
