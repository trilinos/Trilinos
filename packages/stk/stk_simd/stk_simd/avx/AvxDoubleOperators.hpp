
namespace stk {
namespace simd {

namespace hidden {
static const simd::Bool TRUE_VAL(true);
}

STK_MATH_FORCE_INLINE simd::Double operator+ (const simd::Double& a, const simd::Double& b) {
  return simd::Double(_mm256_add_pd(a._data, b._data));
}

STK_MATH_FORCE_INLINE simd::Double operator- (const simd::Double& a, const simd::Double& b) {
  return simd::Double(_mm256_sub_pd(a._data, b._data));
}

STK_MATH_FORCE_INLINE simd::Double operator* (const simd::Double& a, const simd::Double& b) {
  return simd::Double(_mm256_mul_pd(a._data, b._data));
}
   
STK_MATH_FORCE_INLINE simd::Double operator/ (const simd::Double& a, const simd::Double& b) {
  return simd::Double(_mm256_div_pd(a._data, b._data));
}

STK_MATH_FORCE_INLINE simd::Double operator+ (const double a, const simd::Double& b) {
  return simd::Double(_mm256_add_pd(_mm256_set1_pd(a),b._data));
}

STK_MATH_FORCE_INLINE simd::Double operator+ (const simd::Double& a, const double b) {
  return simd::Double(_mm256_add_pd(a._data, _mm256_set1_pd(b)));
}

STK_MATH_FORCE_INLINE simd::Double operator- (const double a, const simd::Double& b) {
  return simd::Double(_mm256_sub_pd(_mm256_set1_pd(a),b._data));
}
  
STK_MATH_FORCE_INLINE simd::Double operator- (const simd::Double& a, const double b) {
  return simd::Double(_mm256_sub_pd(a._data, _mm256_set1_pd(b)));
}

STK_MATH_FORCE_INLINE simd::Double operator* (const double a, const simd::Double& b) {
  return simd::Double(_mm256_mul_pd(_mm256_set1_pd(a),b._data));
}

STK_MATH_FORCE_INLINE simd::Double operator* (const simd::Double& a, const double b) {
  return simd::Double(_mm256_mul_pd(a._data, _mm256_set1_pd(b)));
}

STK_MATH_FORCE_INLINE simd::Double operator/ (const double a, const simd::Double& b) {
  return simd::Double(_mm256_div_pd(_mm256_set1_pd(a),b._data));
}
  
STK_MATH_FORCE_INLINE simd::Double operator/ (const simd::Double& a, const double b) {
  return simd::Double(_mm256_div_pd(a._data, _mm256_set1_pd(b)));
}

STK_MATH_FORCE_INLINE simd::Bool operator< (const simd::Double& x, const simd::Double& y) {
  return simd::Bool(_mm256_cmp_pd(x._data, y._data, _CMP_LT_OQ));
}

STK_MATH_FORCE_INLINE simd::Bool operator<=(const simd::Double& x, const simd::Double& y) {
  return simd::Bool(_mm256_cmp_pd(x._data, y._data, _CMP_LE_OQ));
}

STK_MATH_FORCE_INLINE simd::Bool operator> (const simd::Double& x, const simd::Double& y) {
  return simd::Bool(_mm256_cmp_pd(y._data, x._data, _CMP_LT_OQ));
}

STK_MATH_FORCE_INLINE simd::Bool operator>=(const simd::Double& x, const simd::Double& y) {
  return simd::Bool(_mm256_cmp_pd(y._data, x._data, _CMP_LE_OQ));
}

STK_MATH_FORCE_INLINE simd::Bool operator==(const simd::Double& x, const simd::Double& y) {
  return simd::Bool(_mm256_cmp_pd(x._data, y._data, _CMP_EQ_OQ));
}
  
STK_MATH_FORCE_INLINE simd::Bool operator!=(const simd::Double& x, const simd::Double& y) {
  return simd::Bool(_mm256_cmp_pd(x._data, y._data, _CMP_NEQ_UQ));
}

STK_MATH_FORCE_INLINE simd::Bool operator&& (const simd::Bool& x, const simd::Bool& y) {
  return simd::Bool(_mm256_and_pd(x._data, y._data));
}

STK_MATH_FORCE_INLINE simd::Bool operator|| (const simd::Bool& x, const simd::Bool& y) {
  return simd::Bool(_mm256_or_pd(x._data, y._data));
}

STK_MATH_FORCE_INLINE simd::Bool operator! (const simd::Bool& x) {
  return simd::Bool( _mm256_xor_pd(hidden::TRUE_VAL._data, x._data) );
}

} // namespace simd
} // namespace stk
