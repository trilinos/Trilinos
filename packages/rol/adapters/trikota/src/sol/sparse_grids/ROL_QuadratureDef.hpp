// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_QUADRATUREDEF_HPP
#define ROL_QUADRATUREDEF_HPP

/** \file   ROL_QuadratureDef.hpp
    \brief  Definition file for the ROL::Quadrature class.
    \author Created by D. Kouri and D. Ridzal.
*/

#include "ROL_QuadratureHelpers.hpp"

namespace ROL {

template <class Real> 
Quadrature<Real>::Quadrature(int dimension) : dimension_(dimension) {
  points_.clear(); weights_.clear(); accuracy_.clear();
  accuracy_.resize(dimension);
}

template <class Real>
void Quadrature<Real>::addPointAndWeight(std::vector<Real> point, Real weight, int loc) {
  this->points_.insert(std::pair<std::vector<Real>,int>(point,loc));
  this->weights_.push_back(weight);
}

template <class Real>
int Quadrature<Real>::getNumPoints() const {
  return this->weights_.size();
} 

template <class Real>
void Quadrature<Real>::getAccuracy(std::vector<int> & accuracy) const {
  accuracy = this->accuracy_;
}

template <class Real>
void Quadrature<Real>::getCubature(std::vector<std::vector<Real> >& points, 
                                   std::vector<Real>& weights) const {
  points.resize(this->weights_.size());
  weights.resize(this->weights_.size());
  typename std::map<std::vector<Real>,int>::const_iterator it;
  for (it=(this->points_).begin(); it!=(this->points_).end();it++) {
    points[it->second] = it->first;
    weights[it->second] = (this->weights_)[it->second];
  }
}

template <class Real>
int Quadrature<Real>::getDimension() const {
  return this->dimension_;
}

template <class Real> 
typename std::map<std::vector<Real>,int>::iterator Quadrature<Real>::begin() {
  return (this->points_).begin();
}

template <class Real> 
typename std::map<std::vector<Real>,int>::iterator Quadrature<Real>::end() {
  return (this->points_).end();
}

template <class Real> 
void Quadrature<Real>::insert(typename std::map<std::vector<Real>,int>::iterator it, 
                              std::vector<Real> point, Real weight) {
  (this->points_).insert(it,std::pair<std::vector<Real>,int>(point,(int)(this->points_).size()));
  (this->weights_).push_back(weight);
}

template <class Real> 
std::vector<Real> Quadrature<Real>::getNode(typename std::map<std::vector<Real>,int>::iterator it) {
  return it->first;
}

template <class Real> 
Real Quadrature<Real>::getWeight(int node) { 
  return (this->weights_)[node];
}

template <class Real> 
Real Quadrature<Real>::getWeight(std::vector<Real> point) {
  return (this->weights_)[(this->points_)[point]];
}

template <class Real>
void Quadrature<Real>::update(Real alpha, Quadrature<Real> &rule) {
  if ( rule.begin() != rule.end() ) {  
    // Initialize an iterator on std::map<std::vector<Real>,Real>
    typename std::map<std::vector<Real>,int>::iterator it;
    typename std::map<std::vector<Real>,int>::iterator it2;
  
    // Get location of last point and weight 
    int loc = points_.size();
 
    // Intersection of rule1 and rule2 and set difference rule2 \ rule1
    for ( it2=rule.begin(); it2!=rule.end(); ++it2 ) {
      it = points_.find(it2->first);
      if ( it != points_.end() ) {
        weights_[it->second] += alpha*rule.getWeight(it2->second);
      } 
      else {
        points_.insert(std::pair<std::vector<Real>, int>(it2->first,loc));
        weights_.push_back(alpha*rule.getWeight(it2->second));
        loc++;
      }
    }
  }
}

template <class Real>
void Quadrature<Real>::normalize(){
  Real sum = 0.0;
  typename std::vector<Real>::iterator it;
  for (it=weights_.begin(); it!=weights_.end(); it++) {
    sum += *it;
  }
  for (it=weights_.begin(); it!=weights_.end(); it++) {
    *it /= sum;
  }
}

// PRIVATE MEMBER FUNCTION DEFINITIONS
template<class Real>
void Quadrature<Real>::buildInitial(const int dimension,
                                    const int maxNumPoints,
                                    const std::vector<EQuadrature> &rule1D,
                                    const std::vector<EGrowth> &growth1D,
                                    const bool isNormalized) {
  accuracy_.clear();
  accuracy_.resize(dimension);
  std::vector<int> degree(1);
  Quadrature<Real> newRule(1); 
  for (int i=0; i<dimension; i++) {
    // Compute 1D rules   
    int numPoints = growthRule1D(maxNumPoints,growth1D[i],rule1D[i]);
    Quadrature<Real> rule1(rule1D[i],numPoints,isNormalized);
    rule1.getAccuracy(degree);
    accuracy_.push_back(degree[0]);
    // Build Tensor Rule
    newRule = kron_prod<Real>(newRule,rule1);
  }
  typename std::map<std::vector<Real>,int>::iterator it;
  int loc = 0;
  std::vector<Real> node;  
  for (it=newRule.begin(); it!=newRule.end(); it++) { 
    node = it->first;      
    addPointAndWeight(node,newRule.getWeight(node),loc);
    loc++;
  }
}

template<class Real>
void Quadrature<Real>::buildSGMGA(const int dimension,
                                  const int maxLevel,
                                  const std::vector<EQuadrature> &rule1D,
                                  const std::vector<EGrowth> &growth1D) {
  std::vector<Real> sparse_point;
  std::vector<Real> sparse_weight;
  int point_num = 0;

  Real tol = std::sqrt(ROL_EPSILON<Real>());

  GWPointer *gw_compute_points;
  gw_compute_points  = new GWPointer[dimension];
  GWPointer *gw_compute_weights;
  gw_compute_weights = new GWPointer[dimension];
  GWPointer2 *gw_compute_order;
  gw_compute_order   = new GWPointer2[dimension];

  int growth = 0;

  for (int i=0; i<dimension; ++i) {
    if(rule1D[i] == QUAD_CLENSHAWCURTIS) {
      gw_compute_points[i]  = &Quadrature<Real>::ClenshawCurtisPoints;
      gw_compute_weights[i] = &Quadrature<Real>::ClenshawCurtisWeights;
      gw_compute_order[i]   = &webbur::level_to_order_exp_cc;
    }
    else if (rule1D[i] == QUAD_FEJER2) {
      gw_compute_points[i]  = &Quadrature<Real>::Fejer2Points;
      gw_compute_weights[i] = &Quadrature<Real>::Fejer2Weights;
      gw_compute_order[i]   = &webbur::level_to_order_exp_f2;
    }
    else if (rule1D[i] == QUAD_LEGENDRE) {
      gw_compute_points[i]  = &Quadrature<Real>::LegendrePoints;
      gw_compute_weights[i] = &Quadrature<Real>::LegendreWeights;
      gw_compute_order[i]   = &webbur::level_to_order_exp_gauss;
    }
    else if (rule1D[i] == QUAD_PATTERSON) {
      gw_compute_points[i]  = &Quadrature<Real>::PattersonPoints;
      gw_compute_weights[i] = &Quadrature<Real>::PattersonWeights;
      gw_compute_order[i]   = &webbur::level_to_order_exp_gp;
    }
    else if (rule1D[i] == QUAD_HERMITE) {
      gw_compute_points[i]  = &Quadrature<Real>::HermitePoints;
      gw_compute_weights[i] = &Quadrature<Real>::HermiteWeights;
      gw_compute_order[i]   = &webbur::level_to_order_exp_gauss;
    }
    else if (rule1D[i] == QUAD_GENZKEISTER) {
      gw_compute_points[i]  = &Quadrature<Real>::GenzKeisterPoints;
      gw_compute_weights[i] = &Quadrature<Real>::GenzKeisterWeights;
      gw_compute_order[i]   = &webbur::level_to_order_exp_hgk;
    }
    else if (rule1D[i] == QUAD_LAGUERRE) {
      gw_compute_points[i]  = &Quadrature<Real>::LaguerrePoints;
      gw_compute_weights[i] = &Quadrature<Real>::LaguerreWeights;
      gw_compute_order[i]   = &webbur::level_to_order_exp_gauss;
    }

    if( growth1D[i] == GROWTH_DEFAULT ||
        growth1D[i] == GROWTH_FULLEXP    ) {
      growth = 2;
    }
    else if ( growth1D[i] == GROWTH_SLOWLIN    ||
              growth1D[i] == GROWTH_SLOWLINODD ||
              growth1D[i] == GROWTH_SLOWEXP       ) {
      growth = 0;
    }
    else if ( growth1D[i] == GROWTH_MODLIN ||
              growth1D[i] == GROWTH_MODEXP    ) {
      growth = 1;
    }
  }

  std::vector<Real> level_weight(dimension,1);
  std::vector<int> np(dimension,0);
  int np_sum = webbur::i4vec_sum(dimension,&np[0]);
  std::vector<Real> p(np_sum,0);

  //  Compute necessary data.
  int point_total_num
    = webbur::sandia_sgmga_size_total(dimension, &level_weight[0], maxLevel,
                                      growth, gw_compute_order);

  point_num
    = webbur::sandia_sgmga_size(dimension, &level_weight[0], maxLevel,
                                gw_compute_points, tol, growth,
                                gw_compute_order);

  std::vector<int> sparse_unique_index(point_total_num);
  webbur::sandia_sgmga_unique_index(dimension, &level_weight[0], maxLevel,
                                    gw_compute_points, tol, point_num,
                                    point_total_num, growth, gw_compute_order,
                                    &sparse_unique_index[0]);

  std::vector<int> sparse_order(dimension*point_num,0);
  std::vector<int> sparse_index(dimension*point_num,0);
  webbur::sandia_sgmga_index(dimension, &level_weight[0], maxLevel, point_num,
                             point_total_num, &sparse_unique_index[0], growth,
                             gw_compute_order, &sparse_order[0],
                             &sparse_index[0]);

  //  Compute points and weights.
  sparse_point.resize(dimension*point_num,0.0);
  webbur::sandia_sgmga_point(dimension, &level_weight[0], maxLevel,
                             gw_compute_points, point_num, &sparse_order[0],
                             &sparse_index[0], growth, gw_compute_order,
                             &sparse_point[0]);

  sparse_weight.resize(point_num,0.0);
  webbur::sandia_sgmga_weight(dimension, &level_weight[0], maxLevel,
                              gw_compute_weights, point_num, point_total_num,
                              &sparse_unique_index[0], growth,
                              gw_compute_order, &sparse_weight[0]);

  delete [] gw_compute_points;
  delete [] gw_compute_weights;
  delete [] gw_compute_order;
  /***********************************************************************/
  /******** END BUILD SPARSE GRID USING SANDIA_SGMGA *********************/
  /***********************************************************************/
  typename std::map<std::vector<Real>,int>::iterator it;
  Real weight(0);
  std::vector<Real> node(dimension);
  for (int i = 0; i < point_num; ++i) {
    weight = sparse_weight[i];
    for (int j = 0; j < dimension; ++j) {
      node[j] = sparse_point[j+i*dimension];
    }
    addPointAndWeight(node,weight,i);
  }
}

template<class Real>
void Quadrature<Real>::ClenshawCurtisPoints(int n, int dim, double x[]) {
  webbur::clenshaw_curtis_compute_points(n,x);
}

template<class Real>
void Quadrature<Real>::ClenshawCurtisWeights(int n, int dim, double w[]) {
  webbur::clenshaw_curtis_compute_weights(n,w);
}

template<class Real>
void Quadrature<Real>::Fejer2Points(int n, int dim, double x[]) {
  webbur::fejer2_compute_points(n,x);
}

template<class Real>
void Quadrature<Real>::Fejer2Weights(int n, int dim, double w[]) {
  webbur::fejer2_compute_weights(n,w);
}

template<class Real>
void Quadrature<Real>::LegendrePoints(int n, int dim, double x[]) {
  webbur::legendre_compute_points(n,x);
}

template<class Real>
void Quadrature<Real>::LegendreWeights(int n, int dim, double w[]) {
  webbur::legendre_compute_weights(n,w);
}

template<class Real>
void Quadrature<Real>::PattersonPoints(int n, int dim, double x[]) {
  webbur::patterson_lookup_points(n,x);
}

template<class Real>
void Quadrature<Real>::PattersonWeights(int n, int dim, double w[]) {
  webbur::patterson_lookup_weights(n,w);
}

template<class Real>
void Quadrature<Real>::HermitePoints(int n, int dim, double x[]) {
  webbur::hermite_compute_points(n,x);
}

template<class Real>
void Quadrature<Real>::HermiteWeights(int n, int dim, double w[]) {
  webbur::hermite_compute_weights(n,w);
}

template<class Real>
void Quadrature<Real>::GenzKeisterPoints(int n, int dim, double x[]) {
  webbur::hermite_genz_keister_lookup_points(n,x);
}

template<class Real>
void Quadrature<Real>::GenzKeisterWeights(int n, int dim, double w[]) {
  webbur::hermite_genz_keister_lookup_weights(n,w);
}

template<class Real>
void Quadrature<Real>::LaguerrePoints(int n, int dim, double x[]) {
  webbur::laguerre_compute_points(n,x);
}

template<class Real>
void Quadrature<Real>::LaguerreWeights(int n, int dim, double w[]) {
  webbur::laguerre_compute_weights(n,w);
}

} // end ROL namespace

#endif
