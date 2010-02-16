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

#ifndef Rythmos_QUADRATURE_BASE_H
#define Rythmos_QUADRATURE_BASE_H

#include "Rythmos_TimeRange.hpp"
#include "Rythmos_Types.hpp"

#include "Thyra_ModelEvaluator.hpp"
#include "Thyra_VectorStdOps.hpp"


namespace Rythmos {

/** \brief Specific implementation of 1D Gaussian based quadrature formulas
 */
template<class Scalar> 
class GaussQuadrature1D : virtual public Teuchos::Describable
{
public:

  virtual RCP<const Array<Scalar> > getPoints() const =0;
  virtual RCP<const Array<Scalar> > getWeights() const =0;
  virtual RCP<const TimeRange<Scalar> > getRange() const =0;
  virtual int getOrder() const =0;
  
};

template<class Scalar>
class GaussLegendreQuadrature1D : virtual public GaussQuadrature1D<Scalar>
{
  public:
    GaussLegendreQuadrature1D(int numNodes);
    virtual ~GaussLegendreQuadrature1D() {}
    
    RCP<const Array<Scalar> > getPoints() const { return points_; }
    RCP<const Array<Scalar> > getWeights() const { return weights_; }
    RCP<const TimeRange<Scalar> > getRange() const { return range_; }
    int getOrder() const { return order_; }

  private:
    int numNodes_;
    void fixQuadrature_(int numNodes);
    RCP<Array<Scalar> > points_;
    RCP<Array<Scalar> > weights_;
    int order_;
    RCP<TimeRange<Scalar> > range_;
};

template<class Scalar>
GaussLegendreQuadrature1D<Scalar>::GaussLegendreQuadrature1D(int numNodes) {
  fixQuadrature_(numNodes);
}

template<class Scalar>
void GaussLegendreQuadrature1D<Scalar>::fixQuadrature_(int numNodes) {
  typedef Teuchos::ScalarTraits<Scalar> ST;
  TEST_FOR_EXCEPTION( numNodes < 2, std::out_of_range, "Error, numNodes < 2" );
  TEST_FOR_EXCEPTION( numNodes > 10, std::out_of_range, "Error, numNodes > 10" );
  numNodes_ = numNodes;
  range_ = Teuchos::rcp(new TimeRange<Scalar>(Scalar(-ST::one()),ST::one()));
  order_ = 2*numNodes_;
  points_ = rcp(new Array<Scalar>(numNodes_) );
  weights_ = rcp(new Array<Scalar>(numNodes_) );

  // These numbers are from David Day's matlab script
  if (numNodes_ == 2) {
    (*points_)[0] = Scalar( -0.57735026918963 );
    (*points_)[1] = Scalar( +0.57735026918963 );
    (*weights_)[0] = ST::one();
    (*weights_)[1] = ST::one();
  } else if (numNodes_ == 3) {
    (*points_)[0] = Scalar( -0.77459666924148 );
    (*points_)[1] = ST::zero();
    (*points_)[2] = Scalar( +0.77459666924148 );
    (*weights_)[0] = Scalar( 0.55555555555556 );
    (*weights_)[1] = Scalar( 0.88888888888889 );
    (*weights_)[2] = Scalar( 0.55555555555556 );
  } else if (numNodes_ == 4) {
    (*points_)[0] = Scalar( -0.86113631159405 );
    (*points_)[1] = Scalar( -0.33998104358486 );
    (*points_)[2] = Scalar( +0.33998104358486 );
    (*points_)[3] = Scalar( +0.86113631159405 );
    (*weights_)[0] = Scalar( 0.34785484513745 );
    (*weights_)[1] = Scalar( 0.65214515486255 );
    (*weights_)[2] = Scalar( 0.65214515486255 );
    (*weights_)[3] = Scalar( 0.34785484513745 );
  } else if (numNodes_ == 5) {
    (*points_)[0] = Scalar( -0.90617984593866 );
    (*points_)[1] = Scalar( -0.53846931010568 );
    (*points_)[2] = ST::zero();
    (*points_)[3] = Scalar( +0.53846931010568 );
    (*points_)[4] = Scalar( +0.90617984593866 );
    (*weights_)[0] = Scalar( 0.23692688505619 );
    (*weights_)[1] = Scalar( 0.47862867049937 );
    (*weights_)[2] = Scalar( 0.56888888888889 );
    (*weights_)[3] = Scalar( 0.47862867049937 );
    (*weights_)[4] = Scalar( 0.23692688505619 );
  } else if (numNodes_ == 6) {
    (*points_)[0] = Scalar( -0.93246951420315 );
    (*points_)[1] = Scalar( -0.66120938646626 );
    (*points_)[2] = Scalar( -0.23861918608320 );
    (*points_)[3] = Scalar( +0.23861918608320 );
    (*points_)[4] = Scalar( +0.66120938646626 );
    (*points_)[5] = Scalar( +0.93246951420315 );
    (*weights_)[0] = Scalar( 0.17132449237917 );
    (*weights_)[1] = Scalar( 0.36076157304814 );
    (*weights_)[2] = Scalar( 0.46791393457269 );
    (*weights_)[3] = Scalar( 0.46791393457269 );
    (*weights_)[4] = Scalar( 0.36076157304814 );
    (*weights_)[5] = Scalar( 0.17132449237917 );
  } else if (numNodes_ == 7) {
    (*points_)[0] = Scalar( -0.94910791234276 );
    (*points_)[1] = Scalar( -0.74153118559939 );
    (*points_)[2] = Scalar( -0.40584515137740 );
    (*points_)[3] = ST::zero();
    (*points_)[4] = Scalar( +0.40584515137740 );
    (*points_)[5] = Scalar( +0.74153118559939 );
    (*points_)[6] = Scalar( +0.94910791234276 );
    (*weights_)[0] = Scalar( 0.12948496616887 );
    (*weights_)[1] = Scalar( 0.27970539148928 );
    (*weights_)[2] = Scalar( 0.38183005050512 );
    (*weights_)[3] = Scalar( 0.41795918367347 );
    (*weights_)[4] = Scalar( 0.38183005050512 );
    (*weights_)[5] = Scalar( 0.27970539148928 );
    (*weights_)[6] = Scalar( 0.12948496616887 );
  } else if (numNodes_ == 8) {
    (*points_)[0] = Scalar( -0.96028985649754 );
    (*points_)[1] = Scalar( -0.79666647741363 );
    (*points_)[2] = Scalar( -0.52553240991633 );
    (*points_)[3] = Scalar( -0.18343464249565 );
    (*points_)[4] = Scalar( +0.18343464249565 );
    (*points_)[5] = Scalar( +0.52553240991633 );
    (*points_)[6] = Scalar( +0.79666647741363 );
    (*points_)[7] = Scalar( +0.96028985649754 );
    (*weights_)[0] = Scalar( 0.10122853629038 );
    (*weights_)[1] = Scalar( 0.22238103445337 );
    (*weights_)[2] = Scalar( 0.31370664587789 );
    (*weights_)[3] = Scalar( 0.36268378337836 );
    (*weights_)[4] = Scalar( 0.36268378337836 );
    (*weights_)[5] = Scalar( 0.31370664587789 );
    (*weights_)[6] = Scalar( 0.22238103445337 );
    (*weights_)[7] = Scalar( 0.10122853629038 );
  } else if (numNodes_ == 9) {
    (*points_)[0] = Scalar( -0.96816023950763 );
    (*points_)[1] = Scalar( -0.83603110732664 );
    (*points_)[2] = Scalar( -0.61337143270059 );
    (*points_)[3] = Scalar( -0.32425342340381 );
    (*points_)[4] = ST::zero();
    (*points_)[5] = Scalar( +0.32425342340381 );
    (*points_)[6] = Scalar( +0.61337143270059 );
    (*points_)[7] = Scalar( +0.83603110732664 );
    (*points_)[8] = Scalar( +0.96816023950763 );
    (*weights_)[0] = Scalar( 0.08127438836157 );
    (*weights_)[1] = Scalar( 0.18064816069486 );
    (*weights_)[2] = Scalar( 0.26061069640294 );
    (*weights_)[3] = Scalar( 0.31234707704000 );
    (*weights_)[4] = Scalar( 0.33023935500126 );
    (*weights_)[5] = Scalar( 0.31234707704000 );
    (*weights_)[6] = Scalar( 0.26061069640294 );
    (*weights_)[7] = Scalar( 0.18064816069486 );
    (*weights_)[8] = Scalar( 0.08127438836157 );
  } else if (numNodes_ == 10) {
    (*points_)[0] = Scalar( -0.97390652851717 );
    (*points_)[1] = Scalar( -0.86506336668898 );
    (*points_)[2] = Scalar( -0.67940956829902 );
    (*points_)[3] = Scalar( -0.43339539412925 );
    (*points_)[4] = Scalar( -0.14887433898163 );
    (*points_)[5] = Scalar( +0.14887433898163 );
    (*points_)[6] = Scalar( +0.43339539412925 );
    (*points_)[7] = Scalar( +0.67940956829902 );
    (*points_)[8] = Scalar( +0.86506336668898 );
    (*points_)[9] = Scalar( +0.97390652851717 );
    (*weights_)[0] = Scalar( 0.06667134430869 );
    (*weights_)[1] = Scalar( 0.14945134915058 );
    (*weights_)[2] = Scalar( 0.21908636251598 );
    (*weights_)[3] = Scalar( 0.26926671931000 );
    (*weights_)[4] = Scalar( 0.29552422471475 );
    (*weights_)[5] = Scalar( 0.29552422471475 );
    (*weights_)[6] = Scalar( 0.26926671931000 );
    (*weights_)[7] = Scalar( 0.21908636251598 );
    (*weights_)[8] = Scalar( 0.14945134915058 );
    (*weights_)[9] = Scalar( 0.06667134430869 );
  }
}

template<class Scalar>
class GaussLobattoQuadrature1D : virtual public GaussQuadrature1D<Scalar>
{
  public:
    GaussLobattoQuadrature1D(int numNodes);
    virtual ~GaussLobattoQuadrature1D() {}
    
    RCP<const Array<Scalar> > getPoints() const { return points_; }
    RCP<const Array<Scalar> > getWeights() const { return weights_; }
    RCP<const TimeRange<Scalar> > getRange() const { return range_; }
    int getOrder() const { return order_; }

  private:
    int numNodes_;
    void fixQuadrature_(int numNodes);
    RCP<Array<Scalar> > points_;
    RCP<Array<Scalar> > weights_;
    int order_;
    RCP<TimeRange<Scalar> > range_;
};

template<class Scalar>
GaussLobattoQuadrature1D<Scalar>::GaussLobattoQuadrature1D(int numNodes) {
  fixQuadrature_(numNodes);
}

template<class Scalar>
void GaussLobattoQuadrature1D<Scalar>::fixQuadrature_(int numNodes) {
  typedef Teuchos::ScalarTraits<Scalar> ST;
  TEST_FOR_EXCEPTION( numNodes < 3, std::out_of_range, "Error, numNodes < 3" );
  TEST_FOR_EXCEPTION( numNodes > 10, std::out_of_range, "Error, numNodes > 10" );
  numNodes_ = numNodes;
  range_ = Teuchos::rcp(new TimeRange<Scalar>(Scalar(-ST::one()),ST::one()));
  order_ = 2*numNodes_-2;
  points_ = rcp(new Array<Scalar>(numNodes_) );
  weights_ = rcp(new Array<Scalar>(numNodes_) );

  // These numbers are from David Day's matlab script
  if (numNodes_ == 3) {
    (*points_)[0] = Scalar(-ST::one());
    (*points_)[1] = ST::zero();
    (*points_)[2] = ST::one();
    (*weights_)[0] = Scalar( ST::one()/(3*ST::one()) ); 
    (*weights_)[1] = Scalar( 4*ST::one()/(3*ST::one()) );
    (*weights_)[2] = Scalar( ST::one()/(3*ST::one()) ); 
  } else if (numNodes_ == 4) {
    (*points_)[0] = Scalar(-ST::one());
    (*points_)[1] = Scalar( -0.44721359549996);
    (*points_)[2] = Scalar( +0.44721359549996);
    (*points_)[3] = ST::one();
    (*weights_)[0] = Scalar( ST::one()/(6*ST::one()) );
    (*weights_)[1] = Scalar( 5*ST::one()/(6*ST::one()) );
    (*weights_)[2] = Scalar( 5*ST::one()/(6*ST::one()) );
    (*weights_)[3] = Scalar( ST::one()/(6*ST::one()) );
  } else if (numNodes_ == 5) {
    (*points_)[0] = Scalar(-ST::one());
    (*points_)[1] = Scalar( -0.65465367070798 );
    (*points_)[2] = ST::zero();
    (*points_)[3] = Scalar( +0.65465367070798 );
    (*points_)[4] = ST::one();
    (*weights_)[0] = Scalar( ST::one()/(10*ST::one()) );
    (*weights_)[1] = Scalar( 49*ST::one()/(90*ST::one()) );
    (*weights_)[2] = Scalar( 32*ST::one()/(45*ST::one()) );
    (*weights_)[3] = Scalar( 49*ST::one()/(90*ST::one()) );
    (*weights_)[4] = Scalar( ST::one()/(10*ST::one()) );
  } else if (numNodes_ == 6) {
    (*points_)[0] = Scalar(-ST::one());
    (*points_)[1] = Scalar( -0.76505532392946 );
    (*points_)[2] = Scalar( -0.28523151648064 );
    (*points_)[3] = Scalar( +0.28523151648064 );
    (*points_)[4] = Scalar( +0.76505532392946 );
    (*points_)[5] = ST::one();
    (*weights_)[0] = Scalar( 0.06666666666667 );
    (*weights_)[1] = Scalar( 0.37847495629785 );
    (*weights_)[2] = Scalar( 0.55485837703549 );
    (*weights_)[3] = Scalar( 0.55485837703549 );
    (*weights_)[4] = Scalar( 0.37847495629785 );
    (*weights_)[5] = Scalar( 0.06666666666667 );
  } else if (numNodes_ == 7) {
    (*points_)[0] = Scalar(-ST::one());
    (*points_)[1] = Scalar( -0.83022389627857 );
    (*points_)[2] = Scalar( -0.46884879347071 );
    (*points_)[3] = ST::zero();
    (*points_)[4] = Scalar( +0.46884879347071 );
    (*points_)[5] = Scalar( +0.83022389627857 );
    (*points_)[6] = ST::one();
    (*weights_)[0] = Scalar( 0.04761904761905 );
    (*weights_)[1] = Scalar( 0.27682604736157 );
    (*weights_)[2] = Scalar( 0.43174538120986 );
    (*weights_)[3] = Scalar( 0.48761904761905 );
    (*weights_)[4] = Scalar( 0.43174538120986 );
    (*weights_)[5] = Scalar( 0.27682604736157 );
    (*weights_)[6] = Scalar( 0.04761904761905 );
  } else if (numNodes_ == 8) {
    (*points_)[0] = Scalar(-ST::one());
    (*points_)[1] = Scalar( -0.87174014850961 );
    (*points_)[2] = Scalar( -0.59170018143314 );
    (*points_)[3] = Scalar( -0.20929921790248 );
    (*points_)[4] = Scalar( +0.20929921790248 );
    (*points_)[5] = Scalar( +0.59170018143314 );
    (*points_)[6] = Scalar( +0.87174014850961 );
    (*points_)[7] = ST::one();
    (*weights_)[0] = Scalar( 0.03571428571429 );
    (*weights_)[1] = Scalar( 0.21070422714351 );
    (*weights_)[2] = Scalar( 0.34112269248350 );
    (*weights_)[3] = Scalar( 0.41245879465870 );
    (*weights_)[4] = Scalar( 0.41245879465870 );
    (*weights_)[5] = Scalar( 0.34112269248350 );
    (*weights_)[6] = Scalar( 0.21070422714351 );
    (*weights_)[7] = Scalar( 0.03571428571429 );
  } else if (numNodes_ == 9) {
    (*points_)[0] = Scalar(-ST::one());
    (*points_)[1] = Scalar( -0.89975799541146 );
    (*points_)[2] = Scalar( -0.67718627951074 );
    (*points_)[3] = Scalar( -0.36311746382618 );
    (*points_)[4] = ST::zero();
    (*points_)[5] = Scalar( +0.36311746382618 );
    (*points_)[6] = Scalar( +0.67718627951074 );
    (*points_)[7] = Scalar( +0.89975799541146 );
    (*points_)[8] = ST::one();
    (*weights_)[0] = Scalar( 0.02777777777778 );
    (*weights_)[1] = Scalar( 0.16549536156081 );
    (*weights_)[2] = Scalar( 0.27453871250016 );
    (*weights_)[3] = Scalar( 0.34642851097305 );
    (*weights_)[4] = Scalar( 0.37151927437642 );
    (*weights_)[5] = Scalar( 0.34642851097305 );
    (*weights_)[6] = Scalar( 0.27453871250016 );
    (*weights_)[7] = Scalar( 0.16549536156081 );
    (*weights_)[8] = Scalar( 0.02777777777778 );
  } else if (numNodes_ == 10) {
    (*points_)[0] = Scalar(-ST::one());
    (*points_)[1] = Scalar( -0.91953390816646 );
    (*points_)[2] = Scalar( -0.73877386510551 );
    (*points_)[3] = Scalar( -0.47792494981044 );
    (*points_)[4] = Scalar( -0.16527895766639 );
    (*points_)[5] = Scalar( +0.16527895766639 );
    (*points_)[6] = Scalar( +0.47792494981044 );
    (*points_)[7] = Scalar( +0.73877386510551 );
    (*points_)[8] = Scalar( +0.91953390816646 );
    (*points_)[9] = ST::one();
    (*weights_)[0] = Scalar( 0.02222222222222 );
    (*weights_)[1] = Scalar( 0.13330599085107 );
    (*weights_)[2] = Scalar( 0.22488934206313 );
    (*weights_)[3] = Scalar( 0.29204268367968 );
    (*weights_)[4] = Scalar( 0.32753976118390 );
    (*weights_)[5] = Scalar( 0.32753976118390 );
    (*weights_)[6] = Scalar( 0.29204268367968 );
    (*weights_)[7] = Scalar( 0.22488934206313 );
    (*weights_)[8] = Scalar( 0.13330599085107 );
    (*weights_)[9] = Scalar( 0.02222222222222 );
  }
}


template<class Scalar>
RCP<Thyra::VectorBase<Scalar> > eval_f_t(
    const Thyra::ModelEvaluator<Scalar>& me,
    Scalar t
    ) {
  typedef Teuchos::ScalarTraits<Scalar> ST;
  typedef Thyra::ModelEvaluatorBase MEB;
  MEB::InArgs<Scalar> inArgs = me.createInArgs();
  inArgs.set_t(t);
  MEB::OutArgs<Scalar> outArgs = me.createOutArgs();
  RCP<Thyra::VectorBase<Scalar> > f_out = Thyra::createMember(me.get_f_space());
  V_S(outArg(*f_out),ST::zero());
  outArgs.set_f(f_out);
  me.evalModel(inArgs,outArgs);
  return f_out;
}

template<class Scalar>
Scalar translateTimeRange(
    Scalar t,
    const TimeRange<Scalar>& sourceRange,
    const TimeRange<Scalar>& destinationRange
    ) {
  Scalar r = destinationRange.length()/sourceRange.length();
  return r*t+destinationRange.lower()-r*sourceRange.lower();
}

template<class Scalar>
RCP<Thyra::VectorBase<Scalar> > computeArea(
    const Thyra::ModelEvaluator<Scalar>& me, 
    const TimeRange<Scalar>& tr, 
    const GaussQuadrature1D<Scalar>& gq
    ) {
  typedef Teuchos::ScalarTraits<Scalar> ST;
  RCP<Thyra::VectorBase<Scalar> > area = Thyra::createMember(me.get_x_space());
  V_S(outArg(*area),ST::zero());
  RCP<const TimeRange<Scalar> > sourceRange = gq.getRange();
  RCP<const Array<Scalar> > sourcePts = gq.getPoints();
  RCP<const Array<Scalar> > sourceWts = gq.getWeights();
  Array<Scalar> destPts(*sourcePts);
  for (unsigned int i=0 ; i<sourcePts->size() ; ++i) {
    destPts[i] = translateTimeRange<Scalar>((*sourcePts)[i],*sourceRange,tr);
  }
  Scalar r = tr.length()/sourceRange->length();
  for (unsigned int i=0 ; i<destPts.size() ; ++i) {
    RCP<Thyra::VectorBase<Scalar> > tmpVec = eval_f_t<Scalar>(me,destPts[i]);
    Vp_StV(outArg(*area),r*(*sourceWts)[i],*tmpVec);
  }
  return area;
}


} // namespace Rythmos

#endif //Rythmos_QUADRATURE_BASE_H
