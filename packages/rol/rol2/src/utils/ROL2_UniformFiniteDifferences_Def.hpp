#pragma once
#ifndef ROL2_UNIFORMFINITEDIFFERENCES_DEF_HPP
#define ROL2_UNIFORMFINITEDIFFERENCES_DEF_HPP

namespace ROL2 {

template<typename ScalarValuedFunction,
         typename UpdateFunction,
         typename Real,
         typename VectorType>
Real UFD::diff_scalar( ScalarValuedFunction f_value,
                       UpdateFunction       f_update,
                       const VectorType&    v,
                       const VectorType&    x,
                       Real                 h,
                       int                  order ) {

  auto  xc_ptr = x.clone();
  int i = order - 1;
  auto s = [h,i]( int j ) { return h*static_cast<Real>(shift[i][j]); };
  auto w = [i]( int j )   { return static_cast<Real>(w_numer[i][j]); };

  auto& xc = *xc_ptr;

  auto phi = [&] ( Real alpha ) {
    if(alpha==0) xc.set(x);
    else xc.axpy(alpha,v);
    f_update(xc);
    return f_value(xc);
  };

  Real fx = 0;
  for( int j=0; j<=order; ++j ) fx += w(j) * phi( s(j) );

  return fx/(h*w_denom[i]);
} // UFD::diff_scalar



template<typename VectorValuedFunction,
         typename UpdateFunction,
         typename RangeVectorType,
         typename DomainVectorType,
         typename Real>
void UFD::diff_vector( VectorValuedFunction    f_value,
                       UpdateFunction          f_update,
                       RangeVectorType&        Jv,
                       const DomainVectorType& v,
                       const DomainVectorType& x,
                       Real                    h,
                       int                     order ) {

  auto xc_ptr  = x.clone();
  auto c_ptr   = Jv.clone();
  int  i = order - 1;

  auto s = [h,i]( int j ) { return h*static_cast<Real>(shift[i][j]); };
  auto w = [i]( int j )   { return static_cast<Real>(w_numer[i][j]); };

  auto& xc  = *xc_ptr;
  auto& c   = *c_ptr;

  Jv.zero();

  auto phi = [&xc,&x,&v,&f_value,&f_update] ( RangeVectorType& result, Real alpha ) {
    result.zero();
    if(alpha==0) xc.set(x);
    else xc.axpy(alpha,v);
    f_update(xc);
    f_value(result,xc);
  };

  for( int j=0; j<=order; ++j ) {
    c.zero();
    phi(c,s(j));
    Jv.axpy(w(j),c);
  }

  Jv.scale(1/(h*w_denom[i]));

} // UFD::diff_vector

} // namespace ROL2

#endif // ROL2_UNIFORMFINITEDIFFERENCES_DEF_HPP

