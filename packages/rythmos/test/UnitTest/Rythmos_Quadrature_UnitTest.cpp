//@HEADER
// ***********************************************************************
//
//                           Rythmos Package
//                 Copyright (2006) Sandia Corporation
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
// Questions? Contact Todd S. Coffey (tscoffe@sandia.gov)
//
// ***********************************************************************
//@HEADER

#include "Teuchos_UnitTestHarness.hpp"

#include "Rythmos_QuadratureBase.hpp"
#include "../PolynomialModel/PolynomialModel.hpp"
#include "Thyra_DetachedVectorView.hpp"

namespace Rythmos {

using Teuchos::rcp;

TEUCHOS_UNIT_TEST( Rythmos_Quadrature, createGauss ){
  GaussLegendreQuadrature1D<double> gq(2);
  TEST_EQUALITY_CONST( gq.getOrder(), 4 );
  RCP<const Array<double> > points = gq.getPoints();
  RCP<const Array<double> > weights = gq.getWeights();
  double tol = 1.0e-10;
  TEST_FLOATING_EQUALITY( (*points)[0], -0.5773502692, tol );
  TEST_FLOATING_EQUALITY( (*points)[1], +0.5773502692, tol );
  TEST_EQUALITY_CONST( (*weights)[0], 1.0 );
  TEST_EQUALITY_CONST( (*weights)[1], 1.0 );
  RCP<const TimeRange<double> > range = gq.getRange();
  TEST_EQUALITY_CONST( range->lower(), -1.0 );
  TEST_EQUALITY_CONST( range->upper(), +1.0 );
}

TEUCHOS_UNIT_TEST( Rythmos_Quadrature, translateTimeRange ) {
  TimeRange<double> dest(-1.0,1.0);
  TimeRange<double> source(1.0,2.0);
  double t_source = 1.0;
  double t_dest = translateTimeRange<double>(t_source,source,dest);
  TEST_EQUALITY_CONST( t_dest, -1.0 );
  t_source = 2.0;
  t_dest = translateTimeRange<double>(t_source,source,dest);
  TEST_EQUALITY_CONST( t_dest, 1.0 );
  t_source = 1.5;
  t_dest = translateTimeRange<double>(t_source,source,dest);
  TEST_EQUALITY_CONST( t_dest, 0.0 );
}

TEUCHOS_UNIT_TEST( Rythmos_Quadrature, computeLegendreArea) {
  RCP<Teuchos::Polynomial<double> > p;
  RCP<PolynomialModel> me;
  double T0 = -0.25;
  double T1 = +0.75;
  RCP<TimeRange<double> > tr = rcp(new TimeRange<double>(T0,T1) );
  RCP<GaussQuadrature1D<double> > gqPtr;
  unsigned int order = 0;
  double exact_area = 0.0;
  RCP<Thyra::VectorBase<double> > area;

  for (unsigned int numNodes = 2 ; numNodes <= 10 ; ++numNodes) {
    gqPtr = rcp( new GaussLegendreQuadrature1D<double>(numNodes) );
    order = gqPtr->getOrder();

    TEST_EQUALITY( order, order );
    exact_area = 0.0;
    for (unsigned int n=0 ; n < order ; ++n ) {
      TEST_EQUALITY_CONST( n, n );
      p = rcp(new Teuchos::Polynomial<double>(n,1.0) );
      me = polynomialModel(p);
      exact_area += (std::pow(double(T1),int(n+1))-std::pow(double(T0),int(n+1)))/(n+1.0);
      area = computeArea<double>(*me,*tr,*gqPtr);
      TEST_EQUALITY_CONST( Teuchos::is_null(area), false );
      { // scope to delete CosntDetachedVectorView
        Thyra::ConstDetachedVectorView<double> area_view(*area);
        double tol = 1.0e-13;
        TEST_FLOATING_EQUALITY( area_view[0], exact_area, tol );
      }
    }
    // Now verify that the next higher order polynomial is not exact
    unsigned int n = order;
    TEST_EQUALITY( n, n );
    p = rcp(new Teuchos::Polynomial<double>(n,1.0) );
    me = polynomialModel(p);
    exact_area += (std::pow(double(T1),int(n+1))-std::pow(double(T0),int(n+1)))/(n+1.0);
    area = computeArea<double>(*me,*tr,*gqPtr);
    TEST_EQUALITY_CONST( Teuchos::is_null(area), false );
    { // scope to delete CosntDetachedVectorView
      Thyra::ConstDetachedVectorView<double> area_view(*area);
      double tol = 1.0e-13;
      TEST_COMPARE( fabs((area_view[0]-exact_area)/exact_area), >, tol );
    }
  }
}

TEUCHOS_UNIT_TEST( Rythmos_Quadrature, computeLobattoArea) {
  RCP<Teuchos::Polynomial<double> > p;
  RCP<PolynomialModel> me;
  double T0 = -0.25;
  double T1 = +0.75;
  RCP<TimeRange<double> > tr = rcp(new TimeRange<double>(T0,T1) );
  RCP<GaussQuadrature1D<double> > gqPtr;
  unsigned int order = 0;
  double exact_area = 0.0;
  RCP<Thyra::VectorBase<double> > area;

  for (unsigned int numNodes = 3 ; numNodes <= 10 ; ++numNodes) {
    gqPtr = rcp( new GaussLobattoQuadrature1D<double>(numNodes) );
    order = gqPtr->getOrder();

    TEST_EQUALITY( order, order );
    exact_area = 0.0;
    for (unsigned int n=0 ; n < order ; ++n ) {
      TEST_EQUALITY_CONST( n, n );
      p = rcp(new Teuchos::Polynomial<double>(n,1.0) );
      me = polynomialModel(p);
      exact_area += (std::pow(double(T1),int(n+1))-std::pow(double(T0),int(n+1)))/(n+1.0);
      area = computeArea<double>(*me,*tr,*gqPtr);
      TEST_EQUALITY_CONST( Teuchos::is_null(area), false );
      { // scope to delete CosntDetachedVectorView
        Thyra::ConstDetachedVectorView<double> area_view(*area);
        double tol = 1.0e-13;
        TEST_FLOATING_EQUALITY( area_view[0], exact_area, tol );
      }
    }
    // Now verify that the next higher order polynomial is not exact
    unsigned int n = order;
    TEST_EQUALITY( n, n );
    p = rcp(new Teuchos::Polynomial<double>(n,1.0) );
    me = polynomialModel(p);
    exact_area += (std::pow(double(T1),int(n+1))-std::pow(double(T0),int(n+1)))/(n+1.0);
    area = computeArea<double>(*me,*tr,*gqPtr);
    TEST_EQUALITY_CONST( Teuchos::is_null(area), false );
    { // scope to delete CosntDetachedVectorView
      Thyra::ConstDetachedVectorView<double> area_view(*area);
      double tol = 1.0e-13;
      TEST_COMPARE( fabs((area_view[0]-exact_area)/exact_area), >, tol );
    }
  }
}

TEUCHOS_UNIT_TEST( Rythmos_Quadrature, eval_f_t ) {
  RCP<Teuchos::Polynomial<double> > p = rcp(new Teuchos::Polynomial<double>(2,0.0) );  
  p->setCoefficient(0,-1.0); // p = x^2-1
  p->setCoefficient(1,0.0);
  p->setCoefficient(2,1.0);
  RCP<PolynomialModel> me = polynomialModel(p);
  RCP<Thyra::VectorBase<double> > pval = eval_f_t<double>(*me,-1.0);
  { // scope to delete CosntDetachedVectorView
    Thyra::ConstDetachedVectorView<double> pval_view(*pval);
    double tol = 1.0e-10;
    TEST_FLOATING_EQUALITY( pval_view[0], 0.0, tol );
  }
  pval = eval_f_t<double>(*me,1.0);
  { // scope to delete CosntDetachedVectorView
    Thyra::ConstDetachedVectorView<double> pval_view(*pval);
    double tol = 1.0e-10;
    TEST_FLOATING_EQUALITY( pval_view[0], 0.0, tol );
  }
  pval = eval_f_t<double>(*me,0.0);
  { // scope to delete CosntDetachedVectorView
    Thyra::ConstDetachedVectorView<double> pval_view(*pval);
    double tol = 1.0e-10;
    TEST_FLOATING_EQUALITY( pval_view[0], -1.0, tol );
  }
}

} // namespace Rythmos

