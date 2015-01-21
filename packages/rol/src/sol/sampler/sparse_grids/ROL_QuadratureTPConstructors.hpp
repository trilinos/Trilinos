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


/** \file   ROL_QuadratureDef.hpp
    \brief  Definition file for the ROL::Quadrature class.
    \author Created by D. Kouri and D. Ridzal.
*/

namespace ROL {

template <class Real> 
Quadrature<Real>::Quadrature(int dimension, std::vector<int> numPoints1D, std::vector<EROLBurkardt> rule1D, 
                             bool isNormalized) 
  : dimension_(dimension) {
  /*
    This constructor builds a tensor product rule according to quadInfo myRule.
  */  
  TEUCHOS_TEST_FOR_EXCEPTION((dimension!=(int)numPoints1D.size()||
		      dimension!=(int)rule1D.size()),std::out_of_range,
           ">>> ERROR (Quadrature): Dimension mismatch for inputs.");
  accuracy_.clear();
  std::vector<int> degree(1,0);
  Quadrature<Real> newRule(0,1);
  for (int i=0; i<dimension; i++) {
    // Compute 1D rules
    Quadrature<Real> rule1(rule1D[i],numPoints1D[i],isNormalized);
    rule1.getAccuracy(degree);
    accuracy_.push_back(degree[0]);
    // Build Tensor Rule
    newRule = kron_prod<Real>(newRule,rule1);
  }
  typename std::map<std::vector<Real>,int>::iterator it;
  int loc = 0;
  std::vector<Real> node(dimension,0.0);
  for (it=newRule.begin(); it!=newRule.end(); it++) {
    node = it->first;
    addPointAndWeight(node,newRule.getWeight(node),loc);
    loc++;
  }
  if (isNormalized) {
    normalize();
  }
} 

template <class Real> 
Quadrature<Real>::Quadrature(int dimension, std::vector<int> numPoints1D, std::vector<EROLBurkardt> rule1D, 
                             std::vector<EROLGrowth> growth1D, bool isNormalized) 
  : dimension_(dimension) {
  /*
    This constructor builds a tensor product rule according to quadInfo myRule.
  */  
  TEUCHOS_TEST_FOR_EXCEPTION((dimension!=(int)numPoints1D.size()||
		      dimension!=(int)rule1D.size()||
		      dimension!=(int)growth1D.size()),std::out_of_range,
           ">>> ERROR (Quadrature): Dimension mismatch for inputs.");
  accuracy_.clear();
  accuracy_.resize(dimension);
  std::vector<int> degree(1);
  Quadrature<Real> newRule(0,1);
  for (int i=0; i<dimension; i++) {
    // Compute 1D rules
    int numPoints = growthRule1D(numPoints1D[i],growth1D[i],rule1D[i]);
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
  if (isNormalized) {
    normalize();
  }
} 

template <class Real> 
Quadrature<Real>::Quadrature(int dimension, int maxNumPoints, std::vector<EROLBurkardt> rule1D, 
                             std::vector<EROLGrowth> growth1D, bool isNormalized)
  : dimension_(dimension) {
  /*
    This constructor builds a tensor product rule according to quadInfo myRule.
  */  
  TEUCHOS_TEST_FOR_EXCEPTION((dimension!=(int)rule1D.size()||
		      dimension!=(int)growth1D.size()),std::out_of_range,
            ">>> ERROR (Quadrature): Dimension mismatch for inputs.");
  accuracy_.clear();
  accuracy_.resize(dimension);
  std::vector<int> degree(1);
  Quadrature<Real> newRule(0,1);
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
  if (isNormalized) {
    normalize();
  }
} 

typedef void ( *GWPointer  ) ( int order, int np, double p[], double w[] );
typedef int  ( SandiaRules::*GWPointer2  ) ( int level, int growth );
typedef void ( SandiaRules2::*GWPointer1 ) ( int order, int dim, double w[] );

template <class Real>
Quadrature<Real>::Quadrature(int dim, int maxLevel, std::vector<EROLBurkardt> rule1D,
                                                 std::vector<EROLGrowth> growth1D, bool isNormalized,
                                                 bool useSandia) : dimension_(dim) {
  std::vector<Real> sparse_point;
  std::vector<Real> sparse_weight;
  int point_num = 0;

  if (!useSandia) {
    SandiaRules SR; sgmga webbur;

    Real tol = std::sqrt(SR.r8_epsilon());

    GWPointer *gw_compute_points;
    gw_compute_points  = new GWPointer[dim];
    GWPointer *gw_compute_weights;
    gw_compute_weights = new GWPointer[dim];

    std::vector<int> growth(dim,0);
    std::vector<int> rule(dim,0);
    for( unsigned i = 0; i < dim; i++ ) {
      switch( rule1D[i] ) {
        case BURK_CLENSHAWCURTIS: rule[i] = 1;  break;
        case BURK_FEJER2:         rule[i] = 2;  break;
        case BURK_PATTERSON:      rule[i] = 3;  break;
        case BURK_LEGENDRE:       rule[i] = 4;  break;
        case BURK_HERMITE:        rule[i] = 5;  break;
        case BURK_LAGUERRE:       rule[i] = 7;  break;
        case BURK_GENZKEISTER:    rule[i] = 10; break;
        default:                  rule[i] = 1;
      }
      switch( growth1D[i] ) {
        case GROWTH_DEFAULT:    growth[i] = 0; break;
        case GROWTH_SLOWLIN:    growth[i] = 1; break;
        case GROWTH_SLOWLINODD: growth[i] = 2; break;
        case GROWTH_MODLIN:     growth[i] = 3; break;
        case GROWTH_SLOWEXP:    growth[i] = 4; break;
        case GROWTH_MODEXP:     growth[i] = 5; break;
        case GROWTH_FULLEXP:    growth[i] = 6; break;
        default:                growth[i] = 0;
      }
    }

    std::vector<Real> level_weight(dim,1.0);
    std::vector<int> np(dim,0);
    int np_sum = SR.i4vec_sum(dim,&np[0]);
    std::vector<Real> p(np_sum,0.0);

    // Compute necessary data.
    int point_total_num = webbur.sgmga_size_total(dim,&level_weight[0],maxLevel,&rule[0],&growth[0]);

    point_num = webbur.sgmga_size(dim,&level_weight[0],maxLevel,&rule[0],&np[0],&p[0],
                                  gw_compute_points,tol,&growth[0]);

    std::vector<int> sparse_unique_index( point_total_num );
    webbur.sgmga_unique_index(dim,&level_weight[0],maxLevel,&rule[0],&np[0],&p[0],gw_compute_points,
                              tol,point_num,point_total_num,&growth[0],&sparse_unique_index[0]);

    std::vector<int> sparse_order(dim*point_num,0);
    std::vector<int> sparse_index(dim*point_num,0);
    webbur.sgmga_index(dim,&level_weight[0],maxLevel,&rule[0],point_num,point_total_num,
                       &sparse_unique_index[0],&growth[0],&sparse_order[0],&sparse_index[0]);

    // Compute points and weights.
    sparse_point.resize(dim*point_num,0.0);
    webbur.sgmga_point(dim,&level_weight[0],maxLevel,&rule[0],&np[0],&p[0],gw_compute_points,point_num, 
                       &sparse_order[0],&sparse_index[0],&growth[0],&sparse_point[0]);

    sparse_weight.resize(point_num,0.0);
    webbur.sgmga_weight(dim,&level_weight[0],maxLevel,&rule[0],&np[0],&p[0],gw_compute_weights,
                        point_num,point_total_num,&sparse_unique_index[0],&growth[0],&sparse_weight[0]);

    delete [] gw_compute_points;
    delete [] gw_compute_weights;

    /***********************************************************************/
    /******** END BUILD SPARSE GRID USING SGMGA ****************************/
    /***********************************************************************/
  }
  else {
    SandiaRules2 SR; SandiaSGMGA webbur;

    Real tol = std::sqrt(SR.r8_epsilon());

    GWPointer1 *gw_compute_points;
    gw_compute_points  = new GWPointer1[dim];
    GWPointer1 *gw_compute_weights;
    gw_compute_weights = new GWPointer1[dim];
    GWPointer2 *gw_compute_order;
    gw_compute_order   = new GWPointer2[dim];

    int growth = 0;

    for (int i=0; i<dim; i++) {
      if(rule1D[i] == BURK_CLENSHAWCURTIS) {
        gw_compute_points[i]  = &SandiaRules2::clenshaw_curtis_points;
        gw_compute_weights[i] = &SandiaRules2::clenshaw_curtis_weights;
        gw_compute_order[i]   = &SandiaRules::level_to_order_exp_cc;
      }
      else if (rule1D[i] == BURK_FEJER2) {
        gw_compute_points[i]  = &SandiaRules2::fejer2_points;
        gw_compute_weights[i] = &SandiaRules2::fejer2_weights;
        gw_compute_order[i]   = &SandiaRules::level_to_order_exp_f2;
      }
      else if (rule1D[i] == BURK_PATTERSON) {
        gw_compute_points[i]  = &SandiaRules2::patterson_points;
        gw_compute_weights[i] = &SandiaRules2::patterson_weights;
        gw_compute_order[i]   = &SandiaRules::level_to_order_exp_gp;
      }
      else if (rule1D[i] == BURK_LEGENDRE) {
        gw_compute_points[i]  = &SandiaRules2::legendre_points;
        gw_compute_weights[i] = &SandiaRules2::legendre_weights;
        gw_compute_order[i]   = &SandiaRules::level_to_order_exp_gauss;
      }
      else if (rule1D[i] == BURK_HERMITE) {
        gw_compute_points[i]  = &SandiaRules2::hermite_points;
        gw_compute_weights[i] = &SandiaRules2::hermite_weights;
        gw_compute_order[i]   = &SandiaRules::level_to_order_exp_gauss;
      }
      else if (rule1D[i] == BURK_LAGUERRE) {
        gw_compute_points[i]  = &SandiaRules2::laguerre_points;
        gw_compute_weights[i] = &SandiaRules2::laguerre_weights;
        gw_compute_order[i]   = &SandiaRules::level_to_order_exp_gauss;
      }
      else if (rule1D[i] == BURK_GENZKEISTER) {
        gw_compute_points[i]  = &SandiaRules2::hermite_genz_keister_points;
        gw_compute_weights[i] = &SandiaRules2::hermite_genz_keister_weights;
        gw_compute_order[i]   = &SandiaRules::level_to_order_exp_hgk;
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

    std::vector<Real> level_weight(dim,1.0);
    std::vector<int> np(dim,0);
    int np_sum = SR.i4vec_sum(dim,&np[0]);
    std::vector<Real> p(np_sum,0.0);

    //  Compute necessary data.
    int point_total_num = webbur.sgmga_size_total(dim,&level_weight[0],maxLevel,growth,
                                                  gw_compute_order);

    point_num = webbur.sgmga_size(dim,&level_weight[0],maxLevel,gw_compute_points,tol,growth,
                                  gw_compute_order);

    std::vector<int> sparse_unique_index(point_total_num);
    webbur.sgmga_unique_index(dim,&level_weight[0],maxLevel,gw_compute_points,tol,point_num,
                              point_total_num,growth,gw_compute_order,&sparse_unique_index[0]);

    std::vector<int> sparse_order(dim*point_num,0);
    std::vector<int> sparse_index(dim*point_num,0);
    webbur.sgmga_index(dim,&level_weight[0],maxLevel,point_num,point_total_num,&sparse_unique_index[0],
                       growth,gw_compute_order,&sparse_order[0],&sparse_index[0]);

    //  Compute points and weights.
    sparse_point.resize(dim*point_num,0.0);
    webbur.sgmga_point(dim,&level_weight[0],maxLevel,gw_compute_points,point_num,&sparse_order[0],
                       &sparse_index[0],growth,gw_compute_order,&sparse_point[0]);

    sparse_weight.resize(point_num,0.0);
    webbur.sgmga_weight(dim,&level_weight[0],maxLevel,gw_compute_weights,point_num,point_total_num,
                        &sparse_unique_index[0],growth,gw_compute_order,&sparse_weight[0]);

    delete [] gw_compute_points;
    delete [] gw_compute_weights;
    delete [] gw_compute_order;
    /***********************************************************************/
    /******** END BUILD SPARSE GRID USING SANDIA_SGMGA *********************/
    /***********************************************************************/
  }

  typename std::map<std::vector<Real>,int>::iterator it;
  Real weight = 0.0;
  std::vector<Real> node(dim,0.0);;
  for (int i = 0; i < point_num; i++) {
    weight = sparse_weight[i];
    for (int j = 0; j < dim; j++) {
      node[j] = sparse_point[j+i*dim];
    }
    addPointAndWeight(node,weight,i);
  }
  if (isNormalized) {
    normalize();
  }
}

template<class Real>
Quadrature<Real>::Quadrature(const char* SGinfo, const char* SGdata, bool isNormalized) {
  // Open SGdata
  std::fstream info; info.open(SGinfo);
  Real buf;

  // Get Dimension Number
  info >> buf;
  int dim_num = (int)buf;
  dimension_ = dim_num;

  // Get Number of Cubature Points/Weights
  info >> buf;
  int point_num = (int)buf;

  // Close SGinfo
  info.close();

  // Open SGdata File
  std::fstream file; file.open(SGdata);

  // Get Cubature Points/Weights
  Real weight = 0.0;
  std::vector<Real> node(dim_num, 0.0);
  int loc = 0;
  for ( int i = 0; i < point_num; i++ ) {
    // Get Cubature Point
    for ( int j = 0; j < dim_num; j++ ) {
      file >> buf;
      node[j] = buf;
    }
    // Get Cubature Weight
    file >> buf;
    weight = buf;
    // Insert Point/Weight into Sparse Grid Storage
    addPointAndWeight(node,weight,loc);
    loc++;
    // Move to Next Line in SGdata
    file.ignore(256,'\n');
  }
  file.close();
  if (isNormalized) {
    normalize();
  }
}
} // end ROL namespace
