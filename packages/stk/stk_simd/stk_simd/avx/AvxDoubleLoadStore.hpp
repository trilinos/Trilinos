
namespace stk {
namespace simd {

STK_MATH_FORCE_INLINE simd::Double load_aligned(const double* x) {
  return simd::Double(_mm256_load_pd(x));
}

STK_MATH_FORCE_INLINE simd::Double load(const double* x) {
  return simd::Double(_mm256_loadu_pd(x));
}
    
STK_MATH_FORCE_INLINE simd::Double load(const double* x, const int offset) {
  return simd::Double(_mm256_setr_pd(x[0],x[offset],x[2*offset],x[3*offset]));
}
  
STK_MATH_FORCE_INLINE void store_aligned(double* x, const simd::Double& z) {
  _mm256_store_pd(x,z._data);
}

STK_MATH_FORCE_INLINE void store(double* x, const simd::Double& z) {
  _mm256_storeu_pd(x,z._data);
}
  
STK_MATH_FORCE_INLINE void store(double* x, const simd::Double& z, const int offset) {
  for (int i=0; i < ndoubles; ++i) {
    x[offset*i] = z[i];
  }
}

} // namespace simd
} // namespace stk
