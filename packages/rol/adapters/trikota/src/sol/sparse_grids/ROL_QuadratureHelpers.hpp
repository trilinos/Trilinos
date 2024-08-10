// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_QUADRATUREHELPERS_HPP
#define ROL_QUADRATUREHELPERS_HPP

#include <set>
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

template<class Real>
class SparseGridIndexSet {
private:
  std::multimap<Real, std::vector<int> > activeIndices_;
  std::set<std::vector<int> > oldIndices_;
  const int dim_;
  const int maxLevel_;
  const bool useMax_;

public:
  SparseGridIndexSet(const int dim, const int maxLevel, const bool useMax = false)
    : dim_(dim), maxLevel_(maxLevel), useMax_(useMax) {}

  void reset(void) {
    activeIndices_.clear();
    oldIndices_.clear();
  }

  bool isOldEmpty(void) const {
    return (oldIndices_.begin() == oldIndices_.end());
  }

  bool isActiveEmpty(void) const {
    return (activeIndices_.begin() == activeIndices_.end());
  }

  bool isEmpty(void) const {
    return isActiveEmpty() && isOldEmpty();
  }

  bool isMember(const std::vector<int> &index) const {
    bool flagOld = (oldIndices_.count(index)>0);
    bool flagAct = false;
    typename std::multimap<Real,std::vector<int> >::const_iterator it;
    for (it = activeIndices_.begin(); it != activeIndices_.end(); ++it) {
      if (it->second==index) {
        flagAct = true;
        break;
      }
    }
    return flagOld || flagAct;
  }

  bool isAdmissible(const std::vector<int> &index) const {
    std::vector<int> ind = index;
    for ( int i = 0; i < dim_; ++i ) {
      if ( ind[i] > 1 ) {
        ind[i]--;
        if ( !(oldIndices_.count(ind)) ) {
          return false;
        }
        ind[i]++;
      }
    }
    return true;
  }

  bool isMaxLevelExceeded(const std::vector<int> &index) const {
    int  maxLevel = (useMax_ ? maxLevel_ : maxLevel_+dim_-1);
    int  level    = 0;
    for ( int l = 0; l < dim_; ++l ) {
      if ( useMax_ ) {
        level = std::max(level,index[l]);
      }
      else {
        level += index[l];
      }
    }
    return (level > maxLevel ? true : false);
  }

  void add(const Real &error, const std::vector<int> &index) {
    // Added index to active set
    if ( isActiveEmpty() ) {
      activeIndices_.insert(std::pair<Real, std::vector<int> >(error, index));
    }
    else {
      activeIndices_.insert(activeIndices_.end()--,
        std::pair<Real, std::vector<int> >(error, index));
    }
  }

  bool get(Real &error, std::vector<int> &index) {
    bool isEmpty = isActiveEmpty();
    if (!isEmpty) {
      // Select index corresponding to largest error
      typename std::multimap<Real,std::vector<int> >::iterator it = activeIndices_.end();
      it--;
      error = it->first;
      index = it->second;
      activeIndices_.erase(it);
      oldIndices_.insert(oldIndices_.end(),index);
    }
    else {
      error = 0;
      index.assign(dim_,1);
    }
    return isEmpty;
  }

  void print(const std::string &name="", const int batchID=0) const {
    if (batchID == 0) {
      // Print Old Indices
      std::stringstream nameOld;
      nameOld << name << "_oldIndexSet.txt";
      std::ofstream fileOld(nameOld.str());
      std::set<std::vector<int> >::const_iterator itOld;
      for (itOld = oldIndices_.begin(); itOld != oldIndices_.end(); ++itOld) {
        std::vector<int> index = *itOld;
        for (int i = 0; i < dim_; ++i) {
          fileOld << std::setw(10) << std::left << index[i];
        }
        fileOld << std::endl;
      }
      fileOld.close();
      // Print Active Indices
      std::stringstream nameAct;
      nameAct << name << "_activeIndexSet.txt";
      std::ofstream fileAct(nameAct.str());
      typename std::multimap<Real, std::vector<int> >::const_iterator itAct;
      for (itAct = activeIndices_.begin(); itAct != activeIndices_.end(); ++itAct) {
        std::vector<int> index = itAct->second;
        for (int i = 0; i < dim_; ++i) {
          fileAct << std::setw(10) << std::left << index[i];
        }
        fileAct << std::endl;
      }
      fileAct.close();
    }
  }
};

} // end ROL namespace

#endif
