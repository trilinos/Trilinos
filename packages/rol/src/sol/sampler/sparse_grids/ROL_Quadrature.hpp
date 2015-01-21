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

#include "ROL_BurkardtRules.hpp"
#include "ROL_SGMGA.hpp"
#include "ROL_SandiaSGMGA.hpp"
#include "Teuchos_Assert.hpp"

namespace ROL {

enum EROLGrowth {
  GROWTH_DEFAULT = 0,
  GROWTH_SLOWLIN,
  GROWTH_SLOWLINODD,
  GROWTH_MODLIN,
  GROWTH_SLOWEXP,
  GROWTH_MODEXP,
  GROWTH_FULLEXP,
  GROWTH_LAST
};

template<class Real>
class Quadrature {
private:  
  typename std::map<std::vector<Real>,int> points_; // keys = nodes, values = location of weights
  std::vector<Real> weights_;
  std::vector<int> accuracy_;
  int dimension_;

protected:
  void addPointAndWeight(std::vector<Real> point, Real weight, int loc);
 
public:
  virtual ~Quadrature() {}
  Quadrature(int dimension = 1);

  // 1D Constructors
  Quadrature(int degree, EROLBurkardt rule, bool isNormalized);
  Quadrature(EROLBurkardt rule, int numPoints, bool isNormalized);
  Quadrature( std::vector<Real>& points, std::vector<Real>& weights );

  // Multi-D Constructors
  Quadrature(int dimension, std::vector<int> numPoints1D, std::vector<EROLBurkardt> rule1D,
                       bool isNormalized);
  Quadrature(int dimension, std::vector<int> numPoints1D, std::vector<EROLBurkardt> rule1D,
                       std::vector<EROLGrowth> growth1D, bool isNormalized);
  Quadrature(int dimension, int maxNumPoints, std::vector<EROLBurkardt> rule1D,
                       std::vector<EROLGrowth> growth1D, bool isNormalized);
  Quadrature(int dimension, int maxNumPoints, std::vector<EROLBurkardt> rule1D,
                       std::vector<EROLGrowth> growth1D, bool isNormalized, bool useSandia);
  Quadrature(const char* SGinfo, const char* SGdata, bool isNormalized);

  virtual void getCubature(std::vector<std::vector<Real> >& cubPoints, std::vector<Real>& cubWeights) const;
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
#include <ROL_QuadratureHelpers.hpp>
#include <ROL_Quadrature1dConstructors.hpp>
#include <ROL_QuadratureTPConstructors.hpp>

#endif
