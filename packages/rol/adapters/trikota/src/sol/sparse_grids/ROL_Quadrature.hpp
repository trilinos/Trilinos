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


/** \file   ROL_Quadrature.hpp
    \brief  Header file for the ROL::Quadrature class.
    \author Created by D. Kouri and D. Ridzal.
*/

#ifndef ROL_QUADRATURE_HPP
#define ROL_QUADRATURE_HPP

#include "TriKota_ConfigDefs.hpp"
#include "sandia_rules.hpp"
#include "sandia_sgmga.hpp"

#include "ROL_QuadratureTypes.hpp"
#include "ROL_QuadratureHelpers.hpp"

namespace ROL {

template<class Real>
class Quadrature {
private:  
  typename std::map<std::vector<Real>,int> points_; // keys = nodes, values = location of weights
  std::vector<Real> weights_;
  std::vector<int> accuracy_;
  int dimension_;

  typedef void ( *GWPointer  ) ( int order, int dim, double w[] );
  typedef int ( *GWPointer2  ) ( int level, int growth );

  void buildInitial(const int dimension,
                    const int maxNumPoints,
                    const std::vector<EQuadrature> &rule1D,
                    const std::vector<EGrowth> &growth1D,
                    const bool isNormalized);
  void buildSGMGA(const int dimension,
                  const int maxNumPoints,
                  const std::vector<EQuadrature> &rule1D,
                  const std::vector<EGrowth> &growth1D);

  static void ClenshawCurtisPoints(int n, int dim, double x[]);
  static void ClenshawCurtisWeights(int n, int dim, double w[]);
  static void Fejer2Points(int n, int dim, double x[]);
  static void Fejer2Weights(int n, int dim, double w[]);
  static void LegendrePoints(int n, int dim, double x[]);
  static void LegendreWeights(int n, int dim, double w[]);
  static void PattersonPoints(int n, int dim, double x[]);
  static void PattersonWeights(int n, int dim, double w[]);
  static void HermitePoints(int n, int dim, double x[]);
  static void HermiteWeights(int n, int dim, double w[]);
  static void GenzKeisterPoints(int n, int dim, double x[]);
  static void GenzKeisterWeights(int n, int dim, double w[]);
  static void LaguerrePoints(int n, int dim, double x[]);
  static void LaguerreWeights(int n, int dim, double w[]);

protected:
  void addPointAndWeight(std::vector<Real> point, Real weight, int loc);
 
public:
  virtual ~Quadrature() {}
  Quadrature(const int dimension = 1);

  // 1D Constructors
  Quadrature(const int degree,
             const EQuadrature rule,
             const bool isNormalized);
  Quadrature(const EQuadrature rule,
             const int numPoints,
             const bool isNormalized);
  Quadrature(const std::vector<Real>& points,
             const std::vector<Real>& weights);

  // Multi-D Constructors
  Quadrature(const int dimension,
             const std::vector<int> &numPoints1D,
             const std::vector<EQuadrature> &rule1D,
             const bool isNormalized);
  Quadrature(const int dimension,
             const std::vector<int> &numPoints1D,
             const std::vector<EQuadrature> &rule1D,
             const std::vector<EGrowth> &growth1D,
             const bool isNormalized);
  Quadrature(const int dimension,
             const int maxNumPoints,
             const std::vector<EQuadrature> &rule1D,
             const std::vector<EGrowth> &growth1D,
             const bool isNormalized,
             const bool adaptive);
  Quadrature(const QuadratureInfo &info);
  Quadrature(const char* SGinfo,
             const char* SGdata,
             const bool isNormalized);

  virtual void getCubature(std::vector<std::vector<Real> >& cubPoints,
                           std::vector<Real>& cubWeights) const;
  virtual int getNumPoints() const;
  virtual void getAccuracy(std::vector<int> & accuracy) const;
  virtual int getDimension() const;
  virtual typename std::map<std::vector<Real>,int>::iterator begin();
  virtual typename std::map<std::vector<Real>,int>::iterator end();
  virtual void insert(typename std::map<std::vector<Real>,int>::iterator it, 
                      std::vector<Real> point, Real weight);
  virtual std::vector<Real> getNode(typename std::map<std::vector<Real>,int>::iterator it);
  virtual Real getWeight(int node);
  virtual Real getWeight(std::vector<Real> point);
  virtual void update(Real alpha, Quadrature<Real> &rule); 
  virtual void normalize();

};

} // ROL Namespace

#include <ROL_QuadratureDef.hpp>
#include <ROL_Quadrature1dConstructors.hpp>
#include <ROL_QuadratureTPConstructors.hpp>

#endif
