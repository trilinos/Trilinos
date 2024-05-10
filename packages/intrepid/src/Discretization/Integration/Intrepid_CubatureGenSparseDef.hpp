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

/** \file   Intrepid_CubatureGenSparseDef.hpp
    \brief  Definition file for the Intrepid::CubatureGenSparse class.
    \author Created by P. Bochev, D. Ridzal, and M. Keegan.
*/


namespace Intrepid {

/**************************************************************************
**  Function Definitions for Class CubatureGenSparse
***************************************************************************/

template <class Scalar, int dimension_, class ArrayPoint, class ArrayWeight>
CubatureGenSparse<Scalar,dimension_,ArrayPoint,ArrayWeight>::CubatureGenSparse(const int degree) :
    degree_(degree) {

  SGNodes<int, dimension_> list;
  SGNodes<int,dimension_> bigger_rules;

  bool continue_making_first_list = true;
  bool more_bigger_rules = true;

  int poly_exp[dimension_];
  int level[dimension_];
  int temp_big_rule[dimension_];
  
  for(int i = 0; i<dimension_; i++){
    poly_exp[i] = 0;
    temp_big_rule[i] = 0;
  }

  while(continue_making_first_list){
    for(int i = 0; i < dimension_; i++)
    {
      int max_exp = 0;
      if(i == 0)
        max_exp = std::max(degree_,1) - Sum(poly_exp,1,dimension_-1);
      else if(i == dimension_ -1)
        max_exp = std::max(degree_,1) - Sum(poly_exp,0,dimension_-2);
      else
        max_exp = std::max(degree_,1) - Sum(poly_exp,0,dimension_-1) + poly_exp[i];

      if(poly_exp[i] < max_exp)
      {
        poly_exp[i]++;
        break;
      }
      else
      {
        if(i == dimension_-1)
          continue_making_first_list = false;
        else
          poly_exp[i] = 0;
          
      }
    }

    if(continue_making_first_list)
    {
      for(int j = 0; j < dimension_;j++)
      {
        /*******************
        **  Slow-Gauss
        ********************/
        level[j] = (int)std::ceil((((Scalar)poly_exp[j])+3.0)/4.0);
        /*******************
        **  Fast-Gauss
        ********************/
        //level[j] = intstd::ceil(std::log(poly_exp[j]+3)/std::log(2) - 1);
      }
      list.addNode(level,1);
      
    }
  }


  while(more_bigger_rules)
  {
    bigger_rules.addNode(temp_big_rule,1);

    for(int i = 0; i < dimension_; i++)
    {
      if(temp_big_rule[i] == 0){
        temp_big_rule[i] = 1;
        break;
      }
      else{
        if(i == dimension_-1)
          more_bigger_rules = false;
        else
          temp_big_rule[i] = 0;
      }
    }
  } 

  for(int x = 0; x < list.size(); x++){
    for(int y = 0; y < bigger_rules.size(); y++)
    { 
      SGPoint<int, dimension_> next_rule;
      for(int t = 0; t < dimension_; t++)
        next_rule.coords[t] = list.nodes[x].coords[t] + bigger_rules.nodes[y].coords[t];

      bool is_in_set = false;
      for(int z = 0; z < list.size(); z++)
      {
        if(next_rule == list.nodes[z]){
          is_in_set = true;
          break;
        }
      }

      if(is_in_set)
      {
        int big_sum[dimension_];
        for(int i = 0; i < dimension_; i++)
          big_sum[i] = bigger_rules.nodes[y].coords[i];
        Scalar coeff = std::pow(-1.0, Sum(big_sum, 0, dimension_-1));
        
        Scalar point[dimension_];
        int point_record[dimension_];

        for(int j = 0; j<dimension_; j++)
          point_record[j] = 1;

        bool more_points = true;

        while(more_points)
        {
          Scalar weight = 1.0;
        
          for(int w = 0; w < dimension_; w++){
            /*******************
            **  Slow-Gauss
            ********************/
            int order1D = 2*list.nodes[x].coords[w]-1;
            /*******************
            **  Fast-Gauss
            ********************/
            //int order1D = (int)std::pow(2.0,next_rule.coords[w]) - 1;

            int cubDegree1D = 2*order1D - 1;
            CubatureDirectLineGauss<Scalar> Cub1D(cubDegree1D);
            FieldContainer<Scalar> cubPoints1D(order1D, 1);
            FieldContainer<Scalar> cubWeights1D(order1D);

            Cub1D.getCubature(cubPoints1D, cubWeights1D);

            point[w] =  cubPoints1D(point_record[w]-1, 0);
            weight = weight * cubWeights1D(point_record[w]-1);
          }     
          weight = weight*coeff;
          grid.addNode(point, weight);

          for(int v = 0; v < dimension_; v++)
          {
            if(point_record[v] < 2*list.nodes[x].coords[v]-1){
              (point_record[v])++;
              break;
            }
            else{
              if(v == dimension_-1)
                more_points = false;
              else
                point_record[v] = 1;
            }
          }
        }
      }
    }
  }

  numPoints_ = grid.size();
}



template <class Scalar, int dimension_, class ArrayPoint, class ArrayWeight>
void CubatureGenSparse<Scalar,dimension_,ArrayPoint,ArrayWeight>::getCubature(ArrayPoint  & cubPoints,
                                                                              ArrayWeight & cubWeights) const{
  grid.copyToArrays(cubPoints, cubWeights);
} // end getCubature

template<class Scalar, int dimension_, class ArrayPoint, class ArrayWeight>
void CubatureGenSparse<Scalar,dimension_, ArrayPoint,ArrayWeight>::getCubature(ArrayPoint& cubPoints,
                                                                               ArrayWeight& cubWeights,
                                                                               ArrayPoint& cellCoords) const
{
    TEUCHOS_TEST_FOR_EXCEPTION( (true), std::logic_error,
                      ">>> ERROR (CubatureGenSparse): Cubature defined in reference space calling method for physical space cubature.");
}


template <class Scalar, int dimension_, class ArrayPoint, class ArrayWeight>
int CubatureGenSparse<Scalar,dimension_,ArrayPoint,ArrayWeight>::getNumPoints() const {
  return numPoints_;
} // end getNumPoints



template <class Scalar, int dimension_, class ArrayPoint, class ArrayWeight>
int CubatureGenSparse<Scalar,dimension_,ArrayPoint,ArrayWeight>::getDimension() const {
  return dimension_;
} // end dimension



template <class Scalar, int dimension_, class ArrayPoint, class ArrayWeight>
void CubatureGenSparse<Scalar,dimension_,ArrayPoint,ArrayWeight>::getAccuracy(std::vector<int> & accuracy) const {
  accuracy.assign(1, degree_);
} //end getAccuracy


} // end namespace Intrepid

#if defined(Intrepid_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Intrepid package is deprecated"
#endif
#endif

