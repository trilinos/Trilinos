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

/** \file   Intrepid_CubatureLineSortedDef.hpp
    \brief  Definition file for the Intrepid::CubatureLineSorted class.
    \author Created by D. Kouri and D. Ridzal.
*/

namespace Intrepid {

template <class Scalar, class ArrayPoint, class ArrayWeight> 
CubatureLineSorted<Scalar,ArrayPoint,ArrayWeight>::CubatureLineSorted(
		      int degree, EIntrepidBurkardt rule, bool isNormalized ) {
  TEUCHOS_TEST_FOR_EXCEPTION((degree < 0),std::out_of_range,
    ">>> ERROR (CubatureLineSorted): No rule implemented for desired polynomial degree.");
  degree_    = degree;
  rule_type_ = rule;
 
  if (rule==BURK_CHEBYSHEV1||rule==BURK_CHEBYSHEV2||
      rule==BURK_LAGUERRE  ||rule==BURK_LEGENDRE  ||
      rule==BURK_HERMITE) {
    numPoints_ = (degree+1)/2+1;
  }
  else if (rule==BURK_CLENSHAWCURTIS||rule==BURK_FEJER2) {
    numPoints_ = degree+1;
  }
  else if (rule==BURK_TRAPEZOIDAL) {
    numPoints_ = 2;
  }
  else if (rule==BURK_PATTERSON) {
    int l = 0, o = (degree-0.5)/1.5;
    for (int i=0; i<8; i++) {
      l = (int)pow(2.0,(double)i+1.0)-1;
      if (l>=o) {
	numPoints_ = l;
	break;
      }
    }
  }
  else if (rule==BURK_GENZKEISTER) {
    int o_ghk[8] = {1,3,9,19,35,37,41,43}; 
    int o = (degree-0.5)/1.5;
    for (int i=0; i<8; i++) {
      if (o_ghk[i]>=o) {
	numPoints_ = o_ghk[i];
	break;
      }
    }
  }

  Teuchos::Array<Scalar> nodes(numPoints_), weights(numPoints_);

  if (rule==BURK_CHEBYSHEV1) { // Gauss-Chebyshev Type 1
    IntrepidBurkardtRules::chebyshev1_compute<Scalar>(numPoints_,
				    nodes.getRawPtr(),weights.getRawPtr());
  }
  else if (rule==BURK_CHEBYSHEV2) { // Gauss-Chebyshev Type 2
    IntrepidBurkardtRules::chebyshev2_compute<Scalar>(numPoints_,
                                    nodes.getRawPtr(),weights.getRawPtr());
  }
  else if (rule==BURK_CLENSHAWCURTIS) { // Clenshaw-Curtis    
    IntrepidBurkardtRules::clenshaw_curtis_compute<Scalar>(numPoints_,
                                    nodes.getRawPtr(),weights.getRawPtr());
  }
  else if (rule==BURK_FEJER2) { // Fejer Type 2
    IntrepidBurkardtRules::fejer2_compute<Scalar>(numPoints_,
                                    nodes.getRawPtr(),weights.getRawPtr());
  }
  else if (rule==BURK_LEGENDRE) { // Gauss-Legendre
    IntrepidBurkardtRules::legendre_compute<Scalar>(numPoints_,
                                    nodes.getRawPtr(),weights.getRawPtr());
  }
  else if (rule==BURK_PATTERSON) { // Gauss-Patterson
    IntrepidBurkardtRules::patterson_lookup<Scalar>(numPoints_,
                                    nodes.getRawPtr(),weights.getRawPtr());
  }
  else if (rule==BURK_TRAPEZOIDAL) { // Trapezoidal Rule
    IntrepidBurkardtRules::trapezoidal_compute<Scalar>(numPoints_,
                                    nodes.getRawPtr(),weights.getRawPtr());
  }
  else if (rule==BURK_HERMITE) { // Gauss-Hermite
    IntrepidBurkardtRules::hermite_compute<Scalar>(numPoints_, 
                                    nodes.getRawPtr(),weights.getRawPtr());
  }
  else if (rule==BURK_GENZKEISTER) { // Hermite-Genz-Keister
    IntrepidBurkardtRules::hermite_genz_keister_lookup<Scalar>(numPoints_,
                                    nodes.getRawPtr(),weights.getRawPtr());
    if (numPoints_>=37) {
      for (int i=0; i<numPoints_; i++) {
	weights[i] *= sqrt(M_PI);
      }
    }
  }
  else if (rule==BURK_LAGUERRE) { // Gauss-Laguerre
    IntrepidBurkardtRules::laguerre_compute<Scalar>(numPoints_,
                                    nodes.getRawPtr(),weights.getRawPtr());
  }

  if (isNormalized) {
    Scalar sum = 0.0;
    for (int i=0; i<numPoints_; i++) {
      sum += weights[i];
    }
    for (int i=0; i<numPoints_; i++) {
      weights[i] /= sum;
    }
  }

  points_.clear(); weights_.clear();
  typename std::map<Scalar,int>::iterator it(points_.begin());
  for (int i=0; i<numPoints_; i++) {
    points_.insert(it,std::pair<Scalar,int>(nodes[i],i));
    weights_.push_back(weights[i]);
    it = points_.end();
  }
} // end constructor
  
template <class Scalar, class ArrayPoint, class ArrayWeight> 
CubatureLineSorted<Scalar,ArrayPoint,ArrayWeight>::CubatureLineSorted(
                   EIntrepidBurkardt rule, int numPoints, bool isNormalized ) {
  TEUCHOS_TEST_FOR_EXCEPTION((numPoints < 0),std::out_of_range, 
     ">>> ERROR (CubatureLineSorted): No rule implemented for desired number of points.");
  numPoints_ = numPoints;
  rule_type_ = rule;

  Teuchos::Array<Scalar> nodes(numPoints_), weights(numPoints_);
  if (rule==BURK_CHEBYSHEV1) { // Gauss-Chebyshev Type 1
    degree_ = 2*numPoints-1;
    IntrepidBurkardtRules::chebyshev1_compute<Scalar>(numPoints_,
                                        nodes.getRawPtr(),weights.getRawPtr());
  }
  else if (rule==BURK_CHEBYSHEV2) { // Gauss-Chebyshev Type 2
    degree_ = 2*numPoints-1;
    IntrepidBurkardtRules::chebyshev2_compute<Scalar>(numPoints_,
                                        nodes.getRawPtr(),weights.getRawPtr());
  }
  else if (rule==BURK_CLENSHAWCURTIS) { // Clenshaw-Curtis    
    degree_ = numPoints-1;
    IntrepidBurkardtRules::clenshaw_curtis_compute<Scalar>(numPoints_,
                                        nodes.getRawPtr(),weights.getRawPtr());
  }
  else if (rule==BURK_FEJER2) { // Fejer Type 2
    degree_ = numPoints-1;
    IntrepidBurkardtRules::fejer2_compute<Scalar>(numPoints_,
                                        nodes.getRawPtr(),weights.getRawPtr());
  }
  else if (rule==BURK_LEGENDRE) { // Gauss-Legendre
    degree_ = 2*numPoints-1;
    IntrepidBurkardtRules::legendre_compute<Scalar>(numPoints_,
                                        nodes.getRawPtr(),weights.getRawPtr());
  }
  else if (rule==BURK_PATTERSON) { // Gauss-Patterson
    bool correctNumPoints = false;
    for (int i=0; i<8; i++) {
      int l = (int)pow(2.0,(double)i+1.0)-1;
      if (numPoints==l) {
	correctNumPoints = true;
	break;
      }
    }
    TEUCHOS_TEST_FOR_EXCEPTION((correctNumPoints==false),std::out_of_range,
	">>> ERROR (CubatureLineSorted): Number of points must be numPoints = 1, 3, 7, 15, 31, 63, 127, 255.");
    Scalar degree = 1.5*(double)numPoints+0.5;
    degree_ = (int)degree;
    IntrepidBurkardtRules::patterson_lookup<Scalar>(numPoints_,
                                        nodes.getRawPtr(),weights.getRawPtr());
  }
  else if (rule==BURK_TRAPEZOIDAL) { // Trapezoidal Rule
    degree_ = 2;
    IntrepidBurkardtRules::trapezoidal_compute<Scalar>(numPoints_,
					nodes.getRawPtr(),weights.getRawPtr());
  }
  else if (rule==BURK_HERMITE) { // Gauss-Hermite
    degree_ = 2*numPoints-1;
    IntrepidBurkardtRules::hermite_compute<Scalar>(numPoints_,
                                        nodes.getRawPtr(),weights.getRawPtr());
  }
  else if (rule==BURK_GENZKEISTER) { // Hermite-Genz-Keister
    bool correctNumPoints = false;
    int o_ghk[8] = {1,3,9,19,35,37,41,43};
    for (int i=0; i<8; i++) {
      if (o_ghk[i]==numPoints) {
	correctNumPoints = true;
	break;
      }
    }
    TEUCHOS_TEST_FOR_EXCEPTION((correctNumPoints==false),std::out_of_range,
       ">>> ERROR (CubatureLineSorted): Number of points must be numPoints = 1, 3, 9, 35, 37, 41, 43.");
    Scalar degree = 1.5*(double)numPoints+0.5;
    degree_ = (int)degree;
    IntrepidBurkardtRules::hermite_genz_keister_lookup<Scalar>(numPoints_,
					nodes.getRawPtr(),weights.getRawPtr());
    if (numPoints_>=37) {
      for (int i=0; i<numPoints_; i++) {
	weights[i] *= sqrt(M_PI);
      }
    }
  }
  else if (rule==BURK_LAGUERRE) { // Gauss-Laguerre
    degree_ = 2*numPoints-1;
    IntrepidBurkardtRules::laguerre_compute<Scalar>(numPoints_,
					nodes.getRawPtr(),weights.getRawPtr());
  }
  
  if (isNormalized) {
    Scalar sum = 0.0;
    for (int i=0; i<numPoints_; i++) {
      sum += weights[i];
    }
    for (int i=0; i<numPoints_; i++) {
      weights[i] /= sum;
    }
  }
  points_.clear(); weights_.clear();
  typename std::map<Scalar,int>::iterator it(points_.begin());
  for (int i=0; i<numPoints; i++) {
    points_.insert(it,std::pair<Scalar,int>(nodes[i],i));
    weights_.push_back(weights[i]);
    it = points_.end(); 
  }
} // end constructor
  
template <class Scalar, class ArrayPoint, class ArrayWeight>
CubatureLineSorted<Scalar,ArrayPoint,ArrayWeight>::CubatureLineSorted(
                 std::vector<Scalar> & points, std::vector<Scalar> & weights) {

  int size = (int)weights.size();
  TEUCHOS_TEST_FOR_EXCEPTION(((int)points.size()!=size),std::out_of_range,
	     ">>> ERROR (CubatureLineSorted): Input dimension mismatch.");
  points_.clear(); weights.clear();
  for (int loc=0; loc<size; loc++) {
    points_.insert(std::pair<Scalar,int>(points[loc],loc));
    weights_.push_back(weights[loc]);
  }
  numPoints_ = size;
}

template <class Scalar, class ArrayPoint, class ArrayWeight>
const char* CubatureLineSorted<Scalar,ArrayPoint,ArrayWeight>::getName() const {
  return cubature_name_;
} // end getName

template <class Scalar, class ArrayPoint, class ArrayWeight>
int CubatureLineSorted<Scalar,ArrayPoint,ArrayWeight>::getNumPoints() const {
  return numPoints_;
} // end getNumPoints

template <class Scalar, class ArrayPoint, class ArrayWeight>
void CubatureLineSorted<Scalar,ArrayPoint,ArrayWeight>::getAccuracy(
                                           std::vector<int> & accuracy) const {
  accuracy.assign(1, degree_);
} // end getAccuracy

template <class Scalar, class ArrayPoint, class ArrayWeight>
const char* CubatureLineSorted<Scalar,ArrayPoint,ArrayWeight>::cubature_name_ = "INTREPID_CUBATURE_LINESORTED";

template <class Scalar, class ArrayPoint, class ArrayWeight>
void CubatureLineSorted<Scalar,ArrayPoint,ArrayWeight>::getCubature(
                      ArrayPoint & cubPoints, ArrayWeight & cubWeights) const {
  typename std::map<Scalar,int>::const_iterator it;
  int i = 0;
  for (it = points_.begin(); it!=points_.end(); it++) {
    cubPoints(i)  = it->first;
    cubWeights(i) = weights_[it->second];
    i++;
  }
} // end getCubature

template <class Scalar, class ArrayPoint, class ArrayWeight>
void CubatureLineSorted<Scalar,ArrayPoint,ArrayWeight>::getCubature(ArrayPoint& cubPoints,
                                                                    ArrayWeight& cubWeights,
                                                                    ArrayPoint& cellCoords) const
{
    TEUCHOS_TEST_FOR_EXCEPTION( (true), std::logic_error,
                      ">>> ERROR (CubatureLineSorted): Cubature defined in reference space calling method for physical space cubature.");
}

template <class Scalar, class ArrayPoint, class ArrayWeight>
int CubatureLineSorted<Scalar,ArrayPoint,ArrayWeight>::getDimension() const {
  return 1;
} // end getDimension

template <class Scalar, class ArrayPoint, class ArrayWeight>
Scalar CubatureLineSorted<Scalar,ArrayPoint,ArrayWeight>::getNode(
                                  typename std::map<Scalar,int>::iterator it) { 
  return it->first;
} // end getNode

template <class Scalar, class ArrayPoint, class ArrayWeight>
Scalar CubatureLineSorted<Scalar,ArrayPoint,ArrayWeight>::getWeight(int node) {
  return weights_[node];
} // end getWeight

template <class Scalar, class ArrayPoint, class ArrayWeight>
Scalar CubatureLineSorted<Scalar,ArrayPoint,ArrayWeight>::getWeight(
                                                                Scalar point) {
  return weights_[points_[point]];
} // end getWeight


template <class Scalar, class ArrayPoint, class ArrayWeight>
typename std::map<Scalar,int>::iterator CubatureLineSorted<Scalar,ArrayPoint,ArrayWeight>::begin(void) {
  return points_.begin();
} // end begin

template <class Scalar, class ArrayPoint, class ArrayWeight> 
typename std::map<Scalar,int>::iterator CubatureLineSorted<Scalar,ArrayPoint,ArrayWeight>::end(void) {
  return points_.end();
} // end end

template <class Scalar, class ArrayPoint, class ArrayWeight>
void CubatureLineSorted<Scalar,ArrayPoint,ArrayWeight>::update(
         Scalar alpha2, CubatureLineSorted<Scalar> & cubRule2, Scalar alpha1) {

  // Initialize an iterator on std::map<Scalar,Scalar>
  typename std::map<Scalar,int>::iterator it;

  // Temporary Container for updated rule
  typename std::map<Scalar,int> newPoints;
  std::vector<Scalar> newWeights;
 
  int loc = 0;
  Scalar node = 0.0;

  // Set Intersection rule1 and rule2
  typename std::map<Scalar,int> inter; 
  std::set_intersection(points_.begin(),points_.end(),
			cubRule2.begin(),cubRule2.end(),
			inserter(inter,inter.begin()),inter.value_comp());
  for (it=inter.begin(); it!=inter.end(); it++) {
    node = it->first;
    newWeights.push_back( alpha1*weights_[it->second]
			 +alpha2*cubRule2.getWeight(node));
    newPoints.insert(std::pair<Scalar,int>(node,loc));
    loc++;
  }
  int isize = inter.size();

  // Set Difference rule1 \ rule2
  int size = weights_.size();
  if (isize!=size) {
    typename std::map<Scalar,int> diff1; 
    std::set_difference(points_.begin(),points_.end(),
			cubRule2.begin(),cubRule2.end(),
			inserter(diff1,diff1.begin()),diff1.value_comp());
    for (it=diff1.begin(); it!=diff1.end(); it++) {
      node = it->first;
      newWeights.push_back(alpha1*weights_[it->second]);
      newPoints.insert(std::pair<Scalar,int>(node,loc));
      loc++;
    }
  }

  // Set Difference rule2 \ rule1
  size = cubRule2.getNumPoints();
  if(isize!=size) {    
    typename std::map<Scalar,int> diff2; 
    std::set_difference(cubRule2.begin(),cubRule2.end(),
			points_.begin(),points_.end(),
			inserter(diff2,diff2.begin()),diff2.value_comp());
    for (it=diff2.begin(); it!=diff2.end(); it++) {
      node = it->first;
      newWeights.push_back(alpha2*cubRule2.getWeight(it->second));
      newPoints.insert(std::pair<Scalar,int>(node,loc));
      loc++;
    }
  }

  points_.clear();  points_.insert(newPoints.begin(),newPoints.end());
  weights_.clear(); weights_.assign(newWeights.begin(),newWeights.end());
  numPoints_ = (int)points_.size(); 
}

int growthRule1D(int index, EIntrepidGrowth growth, EIntrepidBurkardt rule) {
  //
  //  Compute the growth sequence for 1D quadrature rules according to growth.
  //  For more information on growth rules, see 
  //  
  //  J. Burkardt. 1D Quadrature Rules For Sparse Grids.
  //  http://people.sc.fsu.edu/~jburkardt/presentations/sgmga_1d_rules.pdf.
  //  
  //  J. Burkardt. SGMGA: Sparse Grid Mixed Growth Anisotropic Rules.
  //  http://people.sc.fsu.edu/~jburkardt/cpp_src/sgmga/sgmga.html.
  //  
  //  Drew P. Kouri
  //  Sandia National Laboratories - CSRI
  //  May 27, 2011
  //

  int level = index-1;
  //int level = index;
  if (rule==BURK_CLENSHAWCURTIS) { // Clenshaw-Curtis
    if (growth==GROWTH_SLOWLIN) {
      return level+1;
    }
    else if (growth==GROWTH_SLOWLINODD) {
      return 2*((level+1)/2)+1;
    }
    else if (growth==GROWTH_MODLIN) {
      return 2*level+1;
    }
    else if (growth==GROWTH_SLOWEXP) {
      if (level==0) {
	return 1;
      }
      else { 
	int o = 2;
	while(o<2*level+1) {
	  o = 2*(o-1)+1;
	}
	return o;
      }
    }
    else if (growth==GROWTH_MODEXP||growth==GROWTH_DEFAULT) {
      if (level==0) {
	return 1;
      }
      else {
	int o = 2;
	while (o<4*level+1) {
	  o = 2*(o-1)+1;
	}
	return o;
      }
    }
    else if (growth==GROWTH_FULLEXP) {
      if (level==0) {
	return 1;
      }
      else {
	return (int)pow(2.0,(double)level)+1;
      }
    }
  }
  else if (rule==BURK_FEJER2) { // Fejer Type 2
    if (growth==GROWTH_SLOWLIN) {
      return level+1;
    }
    else if (growth==GROWTH_SLOWLINODD) {
      return 2*((level+1)/2)+1;
    }
    else if (growth==GROWTH_MODLIN) {
      return 2*level+1;
    }
    else if (growth==GROWTH_SLOWEXP) {
      int o = 1;
      while (o<2*level+1) {
	o = 2*o+1;
      }
      return o;
    }
    else if (growth==GROWTH_MODEXP||growth==GROWTH_DEFAULT) {
      int o = 1;
      while (o<4*level+1) {
	o = 2*o+1;
      }
      return o;
    }
    else if (growth==GROWTH_FULLEXP) {
      return (int)pow(2.0,(double)level+1.0)-1;
    }
  }

  else if (rule==BURK_PATTERSON) { // Gauss-Patterson
    if (growth==GROWTH_SLOWLIN||
	growth==GROWTH_SLOWLINODD||
	growth==GROWTH_MODLIN) {
      std::cout << "Specified Growth Rule Not Allowed!\n";
      return 0;
    }
    else if (growth==GROWTH_SLOWEXP) {
      if (level==0) {
	return 1;
      }
      else {
	int p = 5;
	int o = 3;
	while (p<2*level+1) {
	  p = 2*p+1;
	  o = 2*o+1;
	}
	return o;
      }
    }
    else if (growth==GROWTH_MODEXP||growth==GROWTH_DEFAULT) {
      if (level==0) {
	return 1;
      }
      else {
	int p = 5;
	int o = 3;
	while (p<4*level+1) {
	  p = 2*p+1;
	  o = 2*o+1;
	}
	return o;
      }
    }
    else if (growth==GROWTH_FULLEXP) {
      return (int)pow(2.0,(double)level+1.0)-1;
    }
  }

  else if (rule==BURK_LEGENDRE) { // Gauss-Legendre
    if (growth==GROWTH_SLOWLIN) {
      return level+1;
    }
    else if (growth==GROWTH_SLOWLINODD) {
      return 2*((level+1)/2)+1;
    }
    else if (growth==GROWTH_MODLIN||growth==GROWTH_DEFAULT) {
      return 2*level+1;
    }
    else if (growth==GROWTH_SLOWEXP) {
      int o = 1;
      while (2*o-1<2*level+1) {
	o = 2*o+1;
      }
      return o;
    }
    else if (growth==GROWTH_MODEXP) {
      int o = 1;
      while (2*o-1<4*level+1) {
	o = 2*o+1;
      }
      return o;
    }
    else if (growth==GROWTH_FULLEXP) {
      return (int)pow(2.0,(double)level+1.0)-1;
    }
  }

  else if (rule==BURK_HERMITE) { // Gauss-Hermite
    if (growth==GROWTH_SLOWLIN) {
      return level+1;
    }
    else if (growth==GROWTH_SLOWLINODD) {
      return 2*((level+1)/2)+1;
    }
    else if (growth==GROWTH_MODLIN||growth==GROWTH_DEFAULT) {
      return 2*level+1;
    }
    else if (growth==GROWTH_SLOWEXP) {
      int o = 1;
      while (2*o-1<2*level+1) {
	o = 2*o+1;
      }
      return o;
    }
    else if (growth==GROWTH_MODEXP) {
      int o = 1;
      while (2*o-1<4*level+1) {
	o = 2*o+1;
      }
      return o;
    }
    else if (growth==GROWTH_FULLEXP) {
      return (int)pow(2.0,(double)level+1.0)-1;
    }
  }
  
  else if (rule==BURK_LAGUERRE) { // Gauss-Laguerre
    if (growth==GROWTH_SLOWLIN) {
      return level+1;
    }
    else if (growth==GROWTH_SLOWLINODD) {
      return 2*((level+1)/2)+1;
    }
    else if (growth==GROWTH_MODLIN||growth==GROWTH_DEFAULT) {
      return 2*level+1;
    }
    else if (growth==GROWTH_SLOWEXP) {
      int o = 1;
      while (2*o-1<2*level+1) {
	o = 2*o+1;
      }
      return o;
    }
    else if (growth==GROWTH_MODEXP) {
      int o = 1;
      while (2*o-1<4*level+1) {
	o = 2*o+1;
      }
      return o;
    }
    else if (growth==GROWTH_FULLEXP) {
      return (int)pow(2.0,(double)level+1.0)-1;
    }
  }

  else if (rule==BURK_CHEBYSHEV1) { // Gauss-Chebyshev Type 1
    if (growth==GROWTH_SLOWLIN) {
      return level+1;
    }
    else if (growth==GROWTH_SLOWLINODD) {
      return 2*((level+1)/2)+1;
    }
    else if (growth==GROWTH_MODLIN||growth==GROWTH_DEFAULT) {
      return 2*level+1;
    }
    else if (growth==GROWTH_SLOWEXP) {
      int o = 1;
      while (2*o-1<2*level+1) {
	o = 2*o+1;
      }
      return o;
    }
    else if (growth==GROWTH_MODEXP) {
      int o = 1;
      while (2*o-1<4*level+1) {
	o = 2*o+1;
      }
      return o;
    }
    else if (growth==GROWTH_FULLEXP) {
      return (int)pow(2.0,(double)level+1.0)-1;
    }
  }


  else if (rule==BURK_CHEBYSHEV2) { // Gauss-Chebyshev Type 2
    if (growth==GROWTH_SLOWLIN) {
      return level+1;
    }
    else if (growth==GROWTH_SLOWLINODD) {
      return 2*((level+1)/2)+1;
    }
    else if (growth==GROWTH_MODLIN||growth==GROWTH_DEFAULT) {
      return 2*level+1;
    }
    else if (growth==GROWTH_SLOWEXP) {
      int o = 1;
      while (2*o-1<2*level+1) {
	o = 2*o+1;
      }
      return o;
    }
    else if (growth==GROWTH_MODEXP) {
      int o = 1;
      while (2*o-1<4*level+1) {
	o = 2*o+1;
      }
      return o;
    }
    else if (growth==GROWTH_FULLEXP) {
      return (int)pow(2.0,(double)level+1.0)-1;
    }
  }
  
  else if (rule==BURK_GENZKEISTER) { // Hermite-Genz-Keister  
    static int o_hgk[5] = { 1, 3, 9, 19, 35 };
    static int p_hgk[5] = { 1, 5, 15, 29, 51 };
    if (growth==GROWTH_SLOWLIN||
	growth==GROWTH_SLOWLINODD||
	growth==GROWTH_MODLIN) {
      std::cout << "Specified Growth Rule Not Allowed!\n";
      return 0;
    }
    else if (growth==GROWTH_SLOWEXP) { 
      int l = 0, p = p_hgk[l], o = o_hgk[l];
      while (p<2*level+1 && l<4) {
	l++;
	p = p_hgk[l];
	o = o_hgk[l];
      }
      return o;
    }
    else if (growth==GROWTH_MODEXP||growth==GROWTH_DEFAULT) {
      int l = 0, p = p_hgk[l], o = o_hgk[l];
      while (p<4*level+1 && l<4) {
	l++;
	p = p_hgk[l];
	o = o_hgk[l];
      }
      return o;
    }
    else if (growth==GROWTH_FULLEXP) {
      int l = level; l = std::max(l,0); l = std::min(l,4);
      return o_hgk[l];
    }
  }  

  else if (rule==BURK_TRAPEZOIDAL) { // Trapezoidal
    if (growth==GROWTH_SLOWLIN) {
      return level+1;
    }
    else if (growth==GROWTH_SLOWLINODD) {
      return 2*((level+1)/2)+1;
    }
    else if (growth==GROWTH_MODLIN) {
      return 2*level+1;
    }
    else if (growth==GROWTH_SLOWEXP) {
      if (level==0) {
	return 1;
      }
      else { 
	int o = 2;
	while(o<2*level+1) {
	  o = 2*(o-1)+1;
	}
	return o;
      }
    }
    else if (growth==GROWTH_MODEXP||growth==GROWTH_DEFAULT) {
      if (level==0) {
	return 1;
      }
      else {
	int o = 2;
	while (o<4*level+1) {
	  o = 2*(o-1)+1;
	}
	return o;
      }
    }
    else if (growth==GROWTH_FULLEXP) {
      if (level==0) {
	return 1;
      }
      else {
	return (int)pow(2.0,(double)level)+1;
      }
    }
  }
  return 0;
} // end growthRule1D

} // Intrepid namespace

#if defined(Intrepid_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Intrepid package is deprecated"
#endif
#endif

