#pragma once
#ifndef ROL2_UNIFORMFINITEDIFFERENCES_DECL_HPP
#define ROL2_UNIFORMFINITEDIFFERENCES_DECL_HPP

namespace ROL2 {

struct UFD {

  /** Approximately compute the derivative of f(x) in the direction v
      using the step size h */
  template<typename ScalarValuedFunction,
           typename UpdateFunction,
           typename Real,
           typename VectorType>
  static Real diff_scalar( ScalarValuedFunction f_value,
                           UpdateFunction       f_update,
                           const VectorType&    v,
                           const VectorType&    x,
                           Real                 h,
                           int                  order = 1 );

  /** Approximately compute the Jacobian of f(x) applied to the direction v
      using the step size h */
  template<typename VectorValuedFunction,
           typename UpdateFunction,
           typename RangeVectorType,
           typename DomainVectorType,
           typename Real>
  static void diff_vector( VectorValuedFunction    f_value,
                           UpdateFunction          f_update,
                           RangeVectorType&        Jv,
                           const DomainVectorType& v,
                           const DomainVectorType& x,
                           Real                    h,
                           int                     order = 1 );

  static constexpr int shift[4][4] = { {  1,  0, 0, 0 },
                                       { -1,  2, 0, 0 },
                                       { -1,  2, 1, 0 },
                                       { -1, -1, 3, 1 } };

  static constexpr int w_numer[4][5] = { { -1,  1,  0,  0,  0 },
                                         {  0, -1,  1,  0,  0 },
                                         { -3, -2,  6, -1,  0 },
                                         {  0, -8,  1,  8, -1 } };

  static constexpr int w_denom[4] = { 1, 2, 6, 12 };

};  // UFD

//constexpr int UFD::shift[4][4];
//constexpr int UFD::w_numer[4][5];
//constexpr int UFD::w_denom[4];

} // namespace ROL2

#endif // ROL2_UNIFORMFINITEDIFFERENCES_DECL_HPP

