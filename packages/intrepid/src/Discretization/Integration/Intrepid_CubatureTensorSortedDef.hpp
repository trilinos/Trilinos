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

/** \file   Intrepid_CubatureTensorSortedDef.hpp
    \brief  Definition file for the Intrepid::CubatureTensorSorted class.
    \author Created by D. Kouri and D. Ridzal.
*/

namespace Intrepid {

template <class Scalar, class ArrayPoint, class ArrayWeight> 
CubatureTensorSorted<Scalar,ArrayPoint,ArrayWeight>::CubatureTensorSorted(
                                               int numPoints, int dimension) {
  /*
    This constructor initializes the nodes and weights for an Ndim quadrature 
    rule and sets the nodes and weights lists to zero.
  */
  points_.clear(); weights_.clear();
  numPoints_ = numPoints;
  dimension_ = dimension;
}

template <class Scalar, class ArrayPoint, class ArrayWeight> 
CubatureTensorSorted<Scalar,ArrayPoint,ArrayWeight>::CubatureTensorSorted(
		     CubatureLineSorted<Scalar> & cubLine) {
  /*
    This constructor takes a 1D rule and casts it as a tensor product rule.
  */
  dimension_ = 1;
  numPoints_ = cubLine.getNumPoints();
  degree_.resize(1);
  cubLine.getAccuracy(degree_);

  int loc = 0;
  std::vector<Scalar> node(1,0.0);
  typename std::map<Scalar,int>::iterator it;
  points_.clear(); weights_.clear(); 
  for (it = cubLine.begin(); it != cubLine.end(); it++) {
    node[0] = cubLine.getNode(it);
    points_.insert(std::pair<std::vector<Scalar>,int>(node,loc));
    weights_.push_back(cubLine.getWeight(it->second));
    loc++;
  }
}

template <class Scalar, class ArrayPoint, class ArrayWeight> 
CubatureTensorSorted<Scalar,ArrayPoint,ArrayWeight>::CubatureTensorSorted(
                 int dimension, std::vector<int> numPoints1D, 
		 std::vector<EIntrepidBurkardt> rule1D, bool isNormalized) {
  /*
    This constructor builds a tensor product rule according to quadInfo myRule.
  */  
  TEST_FOR_EXCEPTION((dimension!=(int)numPoints1D.size()||
		      dimension!=(int)rule1D.size()),std::out_of_range,
           ">>> ERROR (CubatureTensorSorted): Dimension mismatch for inputs.");

  dimension_ = dimension;  
  degree_.resize(dimension);
  std::vector<int> degree(1,0);
  CubatureTensorSorted<Scalar> newRule(0,1);
  for (int i=0; i<dimension; i++) {
    // Compute 1D rules
    CubatureLineSorted<Scalar> rule1(rule1D[i],numPoints1D[i],isNormalized);
    rule1.getAccuracy(degree);
    degree_[i] = degree[0];
    // Build Tensor Rule
    newRule = kron_prod<Scalar>(newRule,rule1);
  }
  numPoints_ = newRule.getNumPoints();
  typename std::map<std::vector<Scalar>,int>::iterator it;
  points_.clear(); weights_.clear();
  int loc = 0;
  std::vector<Scalar> node(dimension_,0.0);
  for (it=newRule.begin(); it!=newRule.end(); it++) {
    node = it->first;
    points_.insert(std::pair<std::vector<Scalar>,int>(node,loc));
    weights_.push_back(newRule.getWeight(node));
    loc++;
  }
} 

template <class Scalar, class ArrayPoint, class ArrayWeight> 
CubatureTensorSorted<Scalar,ArrayPoint,ArrayWeight>::CubatureTensorSorted(
                        int dimension, std::vector<int> numPoints1D, 
			std::vector<EIntrepidBurkardt> rule1D, 
			std::vector<EIntrepidGrowth> growth1D, 
			bool isNormalized) {
  /*
    This constructor builds a tensor product rule according to quadInfo myRule.
  */  
  TEST_FOR_EXCEPTION((dimension!=(int)numPoints1D.size()||
		      dimension!=(int)rule1D.size()||
		      dimension!=(int)growth1D.size()),std::out_of_range,
           ">>> ERROR (CubatureTensorSorted): Dimension mismatch for inputs.");
  dimension_ = dimension;  
  degree_.resize(dimension);
  std::vector<int> degree(1);
  CubatureTensorSorted<Scalar> newRule(0,1);
  for (int i=0; i<dimension; i++) {
    // Compute 1D rules
    int numPoints = growthRule1D(numPoints1D[i],growth1D[i],rule1D[i]);
    CubatureLineSorted<Scalar> rule1(rule1D[i],numPoints,isNormalized);
    rule1.getAccuracy(degree);
    degree_[i] = degree[0];
    // Build Tensor Rule
    newRule = kron_prod<Scalar>(newRule,rule1);
  }
  numPoints_ = newRule.getNumPoints();

  typename std::map<std::vector<Scalar>,int>::iterator it;
  points_.clear(); weights_.clear();
  int loc = 0;
  std::vector<Scalar> node;
  for (it=newRule.begin(); it!=newRule.end(); it++) {
    node = it->first;
    points_.insert(std::pair<std::vector<Scalar>,int>(node,loc));
    weights_.push_back(newRule.getWeight(node));
    loc++;
  }
} 

template <class Scalar, class ArrayPoint, class ArrayWeight> 
CubatureTensorSorted<Scalar,ArrayPoint,ArrayWeight>::CubatureTensorSorted(
			 int dimension, int maxNumPoints, 
			 std::vector<EIntrepidBurkardt> rule1D, 
			 std::vector<EIntrepidGrowth> growth1D, 
			 bool isNormalized) {
  /*
    This constructor builds a tensor product rule according to quadInfo myRule.
  */  
  TEST_FOR_EXCEPTION((dimension!=(int)rule1D.size()||
		      dimension!=(int)growth1D.size()),std::out_of_range,
            ">>> ERROR (CubatureTensorSorted): Dimension mismatch for inputs.");
  dimension_ = dimension;
  degree_.resize(dimension);
  std::vector<int> degree(1);
  CubatureTensorSorted<Scalar> newRule(0,1);
  for (int i=0; i<dimension; i++) {
    // Compute 1D rules   
    int numPoints = growthRule1D(maxNumPoints,growth1D[i],rule1D[i]);
    CubatureLineSorted<Scalar> rule1(rule1D[i],numPoints,isNormalized);
    rule1.getAccuracy(degree);
    degree_[i] = degree[0];
    // Build Tensor Rule
    newRule = kron_prod<Scalar>(newRule,rule1);
  }
  numPoints_ = newRule.getNumPoints();
 
  typename std::map<std::vector<Scalar>,int>::iterator it;
  points_.clear(); weights_.clear();
  int loc = 0;
  std::vector<Scalar> node;
  for (it=newRule.begin(); it!=newRule.end(); it++) {
    node = it->first;
    points_.insert(std::pair<std::vector<Scalar>,int>(node,loc));
    weights_.push_back(newRule.getWeight(node));
    loc++;
  }
} 

/* =========================================================================
                     Access Operator - ruleTP
   ========================================================================= */
template <class Scalar, class ArrayPoint, class ArrayWeight>
int CubatureTensorSorted<Scalar,ArrayPoint,ArrayWeight>::getNumPoints() const {
  return numPoints_;
} // end getNumPoints

template <class Scalar, class ArrayPoint, class ArrayWeight>
void CubatureTensorSorted<Scalar,ArrayPoint,ArrayWeight>::getAccuracy(
				  	   std::vector<int> & accuracy) const {
  accuracy = degree_;
} // end getAccuracy

template <class Scalar, class ArrayPoint, class ArrayWeight>
void CubatureTensorSorted<Scalar,ArrayPoint,ArrayWeight>::getCubature(
		      ArrayPoint & cubPoints, ArrayWeight & cubWeights) const {

  typename std::map<std::vector<Scalar>,int>::const_iterator it;
  for (it=points_.begin(); it!=points_.end();it++) {
    for (int j=0; j<dimension_; j++) {
      cubPoints(it->second,j)  = it->first[j];
    }
    cubWeights(it->second) = weights_[it->second];
  }

  /*
  typename std::map<std::vector<Scalar>,int>::const_iterator it = 
    points_.begin();
  for (int i=0; i<numPoints_; i++) {
    for (int j=0; j<dimension_; j++) {
      cubPoints(i,j)  = it->first[j];
    }
    cubWeights(i) = weights_[it->second];
    it++;
  }
  */
} // end getCubature

template <class Scalar, class ArrayPoint, class ArrayWeight>
int CubatureTensorSorted<Scalar,ArrayPoint,ArrayWeight>::getDimension() const {
  return dimension_;
} // end getDimension

template <class Scalar, class ArrayPoint, class ArrayWeight> 
typename std::map<std::vector<Scalar>,int>::iterator CubatureTensorSorted<Scalar,ArrayPoint,ArrayWeight>::begin() {
  return points_.begin();
}

template <class Scalar, class ArrayPoint, class ArrayWeight> 
typename std::map<std::vector<Scalar>,int>::iterator CubatureTensorSorted<Scalar,ArrayPoint,ArrayWeight>::end() {
  return points_.end();
}

template <class Scalar, class ArrayPoint, class ArrayWeight> 
void CubatureTensorSorted<Scalar,ArrayPoint,ArrayWeight>::insert(
                    typename std::map<std::vector<Scalar>,int>::iterator it, 
		    std::vector<Scalar> point, 
		    Scalar weight) {
  points_.insert(it,std::pair<std::vector<Scalar>,int>(point,
						       (int)points_.size()));
  weights_.push_back(weight);
  numPoints_++;
  return;
}

template <class Scalar, class ArrayPoint, class ArrayWeight> 
std::vector<Scalar> CubatureTensorSorted<Scalar,ArrayPoint,ArrayWeight>::getNode(typename std::map<std::vector<Scalar>,int>::iterator it) {
  /*
    Access node for ruleTP
  */ 
  return it->first;
}

template <class Scalar, class ArrayPoint, class ArrayWeight> 
Scalar CubatureTensorSorted<Scalar,ArrayPoint,ArrayWeight>::getWeight(
								   int node) { 
  /*
    Access weight for ruleTP
  */   
  return weights_[node];
}

template <class Scalar, class ArrayPoint, class ArrayWeight> 
Scalar CubatureTensorSorted<Scalar,ArrayPoint,ArrayWeight>::getWeight(
						   std::vector<Scalar> point) {
  //
  //  Access weight for ruleTP
  //   
  return weights_[points_[point]];
}

template <class Scalar, class ArrayPoint, class ArrayWeight>
void CubatureTensorSorted<Scalar,ArrayPoint,ArrayWeight>::update(
			     Scalar alpha2, 
			     CubatureTensorSorted<Scalar> & cubRule2, 
			     Scalar alpha1) {
  
  // Initialize an iterator on std::map<std::vector<Scalar>,Scalar>
  typename std::map<std::vector<Scalar>,int>::iterator it;
  
  // Temporary Container for updated rule
  typename std::map<std::vector<Scalar>,int> newPoints;
  std::vector<Scalar> newWeights(0,0.0);
  std::vector<Scalar> node(dimension_,0.0);
   int loc = 0;

  // Intersection of rule1 and rule2
  typename std::map<std::vector<Scalar>,int> inter; 
  std::set_intersection(points_.begin(),points_.end(),
			cubRule2.begin(),cubRule2.end(),
			inserter(inter,inter.begin()),inter.value_comp());
  for (it=inter.begin(); it!=inter.end(); it++) { 
    node = it->first;
    newWeights.push_back( alpha1*weights_[it->second]
			 +alpha2*cubRule2.getWeight(node));
    newPoints.insert(std::pair<std::vector<Scalar>,int>(node,loc));
    //points_.erase(node); cubRule2.erase(node);
    loc++;    
  }
  int isize = inter.size(); 

  // Set Difference rule1 \ rule2
  int size = points_.size();
  if (isize!=size) {
    typename std::map<std::vector<Scalar>,int> diff1; 
    std::set_difference(points_.begin(),points_.end(),
			cubRule2.begin(),cubRule2.end(),
			inserter(diff1,diff1.begin()),diff1.value_comp());
    for (it=diff1.begin(); it!=diff1.end(); it++) {      
      node = it->first;
      newWeights.push_back(alpha1*weights_[it->second]);
      newPoints.insert(std::pair<std::vector<Scalar>,int>(node,loc));    
      loc++;
    }  
  }

  // Set Difference rule2 \ rule1
  size = cubRule2.getNumPoints();
  if (isize!=size) {
    typename std::map<std::vector<Scalar>,int> diff2; 
    std::set_difference(cubRule2.begin(),cubRule2.end(),
			points_.begin(),points_.end(),
			inserter(diff2,diff2.begin()),diff2.value_comp());
    for (it=diff2.begin(); it!=diff2.end(); it++) {      
      node = it->first;
      newWeights.push_back(alpha2*cubRule2.getWeight(it->second));
      newPoints.insert(std::pair<std::vector<Scalar>,int>(node,loc));  
      loc++;
    }        
  }  
 
  points_.clear();  points_.insert(newPoints.begin(),newPoints.end());
  weights_.clear(); weights_.assign(newWeights.begin(),newWeights.end());
  numPoints_ = (int)points_.size(); 
}

template <class Scalar, class ArrayPoint, class ArrayWeight>
void CubatureTensorSorted<Scalar,ArrayPoint,ArrayWeight>::normalize(){
  Scalar sum = 0.0;
  
  typename std::vector<Scalar>::iterator it;
  for (it=weights_.begin(); it!=weights_.end(); it++)
    sum += *it;

  for (it=weights_.begin(); it!=weights_.end(); it++)
    *it /= sum;
}


/* ======================================================================
                             Kronecker Products
   ====================================================================== */ 
template <class Scalar>
CubatureTensorSorted<Scalar> kron_prod(CubatureTensorSorted<Scalar> & rule1, 
				       CubatureLineSorted<Scalar> & rule2) {
  /* 
    Compute the Kronecker Product of a Tensor Product rule and a 1D rule.
  */
  int s1   = rule1.getNumPoints();
  // int s2   = rule2.getNumPoints();
  int Ndim = rule1.getDimension();

  if (s1==0) {
    CubatureTensorSorted<Scalar> TPrule(rule2);
    return TPrule;
  }
  else {
    // Initialize Arrays Containing Updated Nodes and Weights
    CubatureTensorSorted<Scalar> TPrule(0,Ndim+1); 

    Scalar weight = 0.0;
    Scalar node2  = 0.0;

    // Perform Kronecker Products
    // Compute Kronecker Product of Nodes
    typename std::map<std::vector<Scalar>,int>::iterator it = TPrule.begin();
    typename std::map<std::vector<Scalar>,int>::iterator it_i;
    typename std::map<Scalar,int>::iterator it_j;
    for (it_i=rule1.begin(); it_i!=rule1.end(); it_i++) {
      for (it_j=rule2.begin(); it_j!=rule2.end(); it_j++) {
	std::vector<Scalar> node = rule1.getNode(it_i);
	node2   = rule2.getNode(it_j);
	weight  = rule1.getWeight(node)*rule2.getWeight(node2);
	node.push_back(node2);
	TPrule.insert(it,node,weight);
	it = TPrule.end();
      }
    }
    return TPrule;
  }
}
} // end Intrepid namespace
