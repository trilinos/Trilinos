#ifndef INTREPID_CUBATURE_SPARSE_HELPER_HPP
#define INTREPID_CUBATURE_SPARSE_HELPER_HPP

#include "Intrepid_ConfigDefs.hpp"
#include "Intrepid_Types.hpp"
#include "Intrepid_FieldContainer.hpp"
#include "Teuchos_TestForException.hpp"
#include "Teuchos_Array.hpp"

namespace Intrepid{

/************************************************************************
**  Class Definition for class SGPoint
**  Function: Helper Class with cosntruction of Sparse Grid
*************************************************************************/
template<class Scalar, int D>
class SGPoint{
public:
  Scalar coords[D];

  SGPoint();
  SGPoint(Scalar p[D]);
  bool const operator==(const SGPoint<Scalar, D> & right);
  bool const operator<(const SGPoint<Scalar, D> & right);
  bool const operator>(const SGPoint<Scalar, D> & right);
  //friend ostream & operator<<(ostream & o, const SGPoint<D> & p);
};

/************************************************************************
**  Class Definition for class SGNodes
**  function: Helper Class with constrution of Sparse Grid
************************************************************************/
template<class Scalar, int D, class ArrayType=FieldContainer<Scalar> >
class SGNodes{
public:
  Teuchos::Array< SGPoint<Scalar, D> > nodes;
  Teuchos::Array< Scalar > weights;
  bool addNode(Scalar new_node[D], Scalar weight);
  void copyToArrays(ArrayType & cubPoints, ArrayType & cubWeights) const;
  //void copyToTeuchos(Teuchos::Array< Scalar* > & cubPoints, Teuchos::Array<Scalar> & cubWeights) const;
  int size() const {return nodes.size();}
};

/**************************************************************************
**  Function Definitions for Class SGPoint
***************************************************************************/
template<class Scalar, int D>
SGPoint<Scalar,D>::SGPoint()
{
  for(int i = 0; i < D; i++)
  {
    coords[i] = 0;
  }
}

template<class Scalar, int D>
SGPoint<Scalar, D>::SGPoint(Scalar p[D])
{
  for(int i = 0; i < D; i++)
  {
    coords[i] = p[i];
  }
}

template<class Scalar, int D>
bool const SGPoint<Scalar, D>::operator==(const SGPoint<Scalar, D> & right)
{
  bool equal = true;

  for(int i = 0; i < D; i++)
  {
    if(coords[i] != right.coords[i])
      return false;
  }

  return equal;
}

template<class Scalar, int D>
bool const SGPoint<Scalar, D>::operator<(const SGPoint<Scalar, D> & right)
{
  for(int i = 0; i < D; i++)
  {
    if(coords[i] < right.coords[i])
      return true;
    else if(coords[i] > right.coords[i])
      return false;
  }
  
  return false;
}

template<class Scalar, int D>
bool const SGPoint<Scalar, D>::operator>(const SGPoint<Scalar, D> & right)
{
  if(this < right || this == right)
    return false;

  return true;
}

template<class Scalar, int D>
std::ostream & operator<<(std::ostream & o, SGPoint<Scalar, D> & p)
{
  o << "(";
  for(int i = 0; i<D;i++)
    o<< p.coords[i] << " ";
  o << ")";
  return o;
}


/**************************************************************************
**  Function Definitions for Class SGNodes
***************************************************************************/

template<class Scalar, int D, class ArrayType>
bool SGNodes<Scalar,D,ArrayType>::addNode(Scalar new_node[D], Scalar weight)
{
  SGPoint<Scalar, D> new_point(new_node);
  bool new_and_added = true;

  if(nodes.size() == 0)
  {
    nodes.push_back(new_point);
    weights.push_back(weight);
  }
  else
  {   
    int left = -1;
    int right = (int)nodes.size();
    int mid_node = (int)ceil(nodes.size()/2.0)-1;

    bool iterate_continue = true;

    while(iterate_continue)
    {
      if(new_point == nodes[mid_node]){
        weights[mid_node] += weight;
        iterate_continue = false;
        new_and_added = false;  
      }
      else if(new_point < nodes[mid_node]){
        if(right - left <= 2)
        {
          //insert the node into the vector to the left of mid_node
          nodes.insert(nodes.begin()+mid_node, new_point);
          weights.insert(weights.begin()+mid_node,weight);
          iterate_continue = false;
        }
        else 
        {
          right = mid_node;
          mid_node += (int)ceil((left-mid_node)/2.0);
        }
      }
      else{ //new_point > nodes[mid_node];

        if(mid_node == (int)nodes.size()-1)
        {
          nodes.push_back(new_point);
          weights.push_back(weight);
          iterate_continue = false;
        }
        else if(right - left <= 2)
        {
          //insert the node into the vector to the right of mid_node
          nodes.insert(nodes.begin()+mid_node+1, new_point);
          weights.insert(weights.begin()+mid_node+1,weight);
          iterate_continue = false;
        }
        else 
        {
          left = mid_node;
          mid_node += (int)ceil((right-mid_node)/2.0);
        }
      }
    }
  }

  return new_and_added;
}

template<class Scalar, int D, class ArrayType>
void SGNodes<Scalar,D,ArrayType>::copyToArrays(ArrayType & cubPoints, ArrayType & cubWeights) const
{
  int numPoints = size();

  for(int i = 0; i < numPoints; i++)
  {
    for (int j=0; j<D; j++) {
      cubPoints(i,j) = nodes[i].coords[j];
    }
    cubWeights(i) = weights[i];
  }
}

/*
template< class Scalar, int D>
void SGNodes<Scalar,D>::copyToTeuchos(Teuchos::Array< Scalar* > & cubPoints, Teuchos::Array<Scalar> & cubWeights) const
{
  int numPoints = size();

  Scalar tempPoint[D];
  cubPoints.assign(numPoints,tempPoint);
  cubWeights.assign(numPoints, 0.0);

  for(int i = 0; i < numPoints; i++)
  {
    cubPoints[i] = nodes[i].coords;
    cubWeights[i] = weights[i];
  }
}
*/

}
#endif
