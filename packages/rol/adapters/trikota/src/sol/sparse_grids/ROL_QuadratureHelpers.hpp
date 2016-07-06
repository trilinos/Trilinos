// @HEADER
// ************************************************************************
//
//               Rapid Optimization Library (ROL) Package
//                 Copyright (2014) Sandia Corporation
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
// Questions? Contact lead developers:
//              Drew Kouri   (dpkouri@sandia.gov) and
//              Denis Ridzal (dridzal@sandia.gov)
//
// ************************************************************************
// @HEADER

#ifndef ROL_QUADRATUREHELPERS_HPP
#define ROL_QUADRATUREHELPERS_HPP

#include "sandia_rules.hpp"

namespace ROL {

template<class Real>
class Quadrature;

template<class Real>
Quadrature<Real> kron_prod(Quadrature<Real> & rule1, Quadrature<Real> & rule2) {
  // Compute the Kronecker Product of a Tensor Product rule and a 1D rule.
  if (rule1.getNumPoints()==0) {
    Quadrature<Real> rule(rule2);
    return rule;
  }
  else {
    // Initialize Arrays Containing Updated Nodes and Weights
    int dim1 = rule1.getDimension();
    int dim2 = rule2.getDimension();
    Quadrature<Real> rule(dim1+dim2); 

    Real weight(0);
    std::vector<Real> node2(dim2);

    // Perform Kronecker Products
    // Compute Kronecker Product of Nodes
    typename std::map<std::vector<Real>,int>::iterator it = rule.begin();
    typename std::map<std::vector<Real>,int>::iterator it_i;
    typename std::map<std::vector<Real>,int>::iterator it_j;
    for (it_i=rule1.begin(); it_i!=rule1.end(); it_i++) {
      for (it_j=rule2.begin(); it_j!=rule2.end(); it_j++) {
	std::vector<Real> node = rule1.getNode(it_i);
	node2   = rule2.getNode(it_j);
	weight  = rule1.getWeight(node)*rule2.getWeight(node2);
	//node.push_back(node2[0]);
        node.insert(node.end(),node2.begin(),node2.end());
	rule.insert(it,node,weight);
	it = rule.end();
      }
    }
    return rule;
  }
}

} // end ROL namespace

#endif
