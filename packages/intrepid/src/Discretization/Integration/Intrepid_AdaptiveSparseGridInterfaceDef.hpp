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

/** \file   Intrepid_AdaptiveSparseGridInterfaceDef.hpp
    \brief  Definition file for the Intrepid::AdaptiveSparseGridInterface class.
    \author Created by D. Kouri and D. Ridzal.
*/

namespace Intrepid {

template<class Scalar, class UserVector>
AdaptiveSparseGridInterface<Scalar,UserVector>::AdaptiveSparseGridInterface(
	      int dimension,
	      std::vector<EIntrepidBurkardt> rule1D,
	      std::vector<EIntrepidGrowth> growth1D,
	      int maxLevel,
	      bool isNormalized) {

  TEST_FOR_EXCEPTION((dimension!=(int)rule1D.size()||
		      dimension!=(int)growth1D.size()),std::out_of_range,
    ">>> ERROR (AdaptiveSparseGridInterface): Dimension mismatch for inputs.");

  dimension_    = dimension;
  rule1D_       = rule1D;
  growth1D_     = growth1D;
  maxLevel_     = maxLevel;
  isNormalized_ = isNormalized;
}

template<class Scalar, class UserVector>
void AdaptiveSparseGridInterface<Scalar,UserVector>::init(UserVector & output) {
  std::vector<int> index(dimension_,1);
  CubatureTensorSorted<Scalar> cubRule(
            dimension_,index,rule1D_,growth1D_,isNormalized_);
  
  // Evaluate the initial contribution to the integral
  initialDiff_ = 1.0;
  output.Update(-1.0,output);
  eval_cubature(output,cubRule);

  // Compute the initial error indicator
  initialDiff_ = error_indicator(output);
  if (fabs(initialDiff_)<INTREPID_TOL) 
    initialDiff_ = 1.0;
}

template<class Scalar, class UserVector>
bool AdaptiveSparseGridInterface<Scalar,UserVector>::max_level(
                                             std::vector<int> index) {
  int dimension = (int)index.size();
  int sum = 0;
  for (int i=0; i<dimension; i++) {
    sum += index[i];
  }
  if (sum <= maxLevel_ + dimension - 1) 
    return true;
  return false;
}

template<class Scalar, class UserVector>
void AdaptiveSparseGridInterface<Scalar,UserVector>::eval_cubature(
	       UserVector & output, 
	       CubatureTensorSorted<Scalar> & cubRule) {

  //int dimf      = 0;                      // Dimension of the integrand
  Scalar weight = 0.0;
  std::vector<Scalar> point(dimension_,(Scalar)0.0);
  //std::vector<Scalar> f(1,0.0);
  Teuchos::RCP<UserVector> f = output.Create(); output.Update(-1.0,output);

  typename std::map<std::vector<Scalar>,int>::iterator it;
  for (it=cubRule.begin(); it!=cubRule.end(); it++) {
    // Evaluate Function
    point.assign((it->first).begin(),(it->first).end()); // Extract point
    f->Update(-1.0,*f);
    eval_integrand(*f,point);     // Evaluate Integrand at point
 
    // Update integral
    weight = cubRule.getWeight(it->second);
    output.Update(weight,*f);
  }
}

template<class Scalar, class UserVector>
void AdaptiveSparseGridInterface<Scalar,UserVector>::getRule(
        std::vector<EIntrepidBurkardt> & rule1D) {
  rule1D.clear();
  rule1D.resize(rule1D_.size());
  rule1D = rule1D_;
}

template<class Scalar, class UserVector>
void AdaptiveSparseGridInterface<Scalar,UserVector>::getGrowth(
	std::vector<EIntrepidGrowth> & growth1D) {
  growth1D.clear();
  growth1D.resize(growth1D_.size());
  growth1D = growth1D_;
}

template<class Scalar, class UserVector>
int AdaptiveSparseGridInterface<Scalar,UserVector>::getDimension(void) {
  return dimension_;
}

template<class Scalar, class UserVector>
bool AdaptiveSparseGridInterface<Scalar,UserVector>::isNormalized() {
  return isNormalized_;
}

template<class Scalar, class UserVector>
Scalar AdaptiveSparseGridInterface<Scalar,UserVector>::getInitialDiff() {
  return initialDiff_;
}

} // end Intrepid namespace



