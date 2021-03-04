#pragma once
#ifndef ROL2_VALIDATEFUNCTION_DEF_HPP
#define ROL2_VALIDATEFUNCTION_DEF_HPP

namespace ROL2 {

template<typename Real>
ValidateFunction<Real>::ValidateFunction( int  order,
                                          int  numSteps,
                                          int  precision,
                                          bool printToStream,
                                          std::ostream& os ) 
  : os_(os), order_(order), numSteps_(numSteps), precision_(precision), 
  printToStream_(printToStream) {
  steps_.resize(numSteps_);
  float fact = 1.0;
  for( auto& e: steps_ ) {
    e = fact;
    fact *= 0.1;
  }
} 

template<typename Real>
ValidateFunction<Real>::ValidateFunction( ROL2::ParameterList& parlist, 
                                          std::ostream&        os ) 
  : os_(os) {
  auto& vflist = parlist.sublist("General").sublist("Function Validator");
  order_     = vflist.get("Finite Difference Order",1);
  numSteps_  = vflist.get("Number of Steps", 13 );
  precision_ = vflist.get("Display Precision",11);
  steps_.resize(numSteps_);
  float fact = 1.0;
  for( auto& e: steps_ ) {
    e = fact;
    fact *= 0.1;
  }    
}


  /* Check the directional derivative of a function that maps vectors to scalars */
template<typename Real>
template<typename ValueFunction,
         typename DerivativeFunction,
         typename UpdateFunction,
         typename PrimalVector,
         typename DualVector>
std::vector<std::vector<Real>>
ValidateFunction<Real>::
scalar_derivative_check( ValueFunction       f_value,
                         DerivativeFunction  f_derivative,
                         UpdateFunction      f_update,
                         const DualVector&   g,
                         const PrimalVector& v,
                         const PrimalVector& x,
                         const std::string&  label ) const {

  auto table = Table("",{"Step size",label,"FD approx","abs error"},
                     true, precision_ );

  int numVals = 4;
  std::vector<Real> tmp(numVals);
  std::vector<std::vector<Real>> vCheck(numSteps_,tmp);

  ROL2::NullStream oldFormatState;
  oldFormatState.copyfmt(os_);

  ROL2::NullStream no_output;

  Ptr<std::ostream> pStream;
  if( printToStream_ ) pStream = makePtrFromRef(os_);
  else                 pStream = makePtrFromRef(no_output);

  auto& os = *pStream;

  // Evaluate reference scalar at input
  f_update(x);

  auto gtmp_ptr = g.clone();
  auto& gtmp    = *gtmp_ptr;

  // Compute the derivative in the given direction
  f_derivative( gtmp, x );

  Real gv = gtmp.dot(v.dual());

  for (int i=0; i<numSteps_; i++) {
    vCheck[i][0] = steps_[i];
    vCheck[i][1] = gv;
    vCheck[i][2] = UFD::diff_scalar( f_value, f_update, v, x, steps_[i], order_ );
    vCheck[i][3] = std::abs(vCheck[i][2] - vCheck[i][1]);

    table.add_row(vCheck[i][0],vCheck[i][1],vCheck[i][2],vCheck[i][3]);
  }

   os << table;

   os_.copyfmt(oldFormatState);
   return vCheck;

} // scalar_derivative_check


template<typename Real>
template<typename ValueFunction,
         typename DerivativeFunction,
         typename UpdateFunction,
         typename DomainVector,
         typename RangeVector>
std::vector<std::vector<Real>>
ValidateFunction<Real>::
vector_derivative_check( ValueFunction       f_value,
                         DerivativeFunction  f_derivative,
                         UpdateFunction      f_update,
                         const RangeVector&  c,
                         const DomainVector& x,
                         const DomainVector& v,
                         const std::string&  label ) const {

  auto table = Table("",{"Step size",label,"FD approx","abs error"},
                     true, precision_ );
  int numVals = 4;
  std::vector<Real> tmp(numVals);
  std::vector<std::vector<Real>> vCheck(numSteps_,tmp);

  ROL2::NullStream oldFormatState;
  oldFormatState.copyfmt(os_);

  auto pStream = ROL2::makeStreamPtr( os_, printToStream_ );

  auto  dc_ptr = c.clone();
  auto  jv_ptr = c.clone();
  auto& dc     = *dc_ptr;
  auto& jv     = *jv_ptr;
  auto& os     = *pStream;

  // Evaluate reference scalar at input
  f_update(x);

  f_derivative(jv,v,x);
  auto norm_jv = jv.norm();
  for(int i=0; i<numSteps_; ++i) {
    vCheck[i][0] = steps_[i];
    vCheck[i][1] = norm_jv;

    UFD::diff_vector( f_value, f_update, dc, v, x, steps_[i], order_ );

    vCheck[i][2] = dc.norm();
    dc.axpy(-1.0,jv);
    vCheck[i][3] = dc.norm();

    table.add_row(vCheck[i][0],vCheck[i][1],vCheck[i][2],vCheck[i][3]);
  }
  os << table;

  os_.copyfmt(oldFormatState);
  return vCheck;
} // vector_derivative_check

} // namespace ROL2

#endif // ROL2_VALIDATEFUNCTION_DEF_HPP

