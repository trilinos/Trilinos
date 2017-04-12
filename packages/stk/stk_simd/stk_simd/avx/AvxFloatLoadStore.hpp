
namespace stk {
namespace simd {

STK_MATH_FORCE_INLINE simd::Float load_aligned(const float* x) {
  return simd::Float(_mm256_load_ps(x));
}

STK_MATH_FORCE_INLINE simd::Float load(const float* x) {
  return simd::Float(_mm256_loadu_ps(x));
}
    
STK_MATH_FORCE_INLINE simd::Float load(const float* x, const int offset) {
  return simd::Float(_mm256_setr_ps(x[0],       x[offset],  x[2*offset],x[3*offset],
                                     x[4*offset],x[5*offset],x[6*offset],x[7*offset]));
}

STK_MATH_FORCE_INLINE void store_aligned(float* x, const simd::Float& z) {
  _mm256_store_ps(x,z._data);
}

STK_MATH_FORCE_INLINE void store(float* x, const simd::Float& z) {
  _mm256_storeu_ps(x,z._data);
}
  
STK_MATH_FORCE_INLINE void store(float* x, const simd::Float& z, const int offset) {
  for (int i=0; i < nfloats; ++i) {
    x[offset*i] = z[i];
  }
}

} // namespace simd
} // namespace stk
