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
