
namespace stk {
namespace simd {

inline simd::Float load_aligned(const float* x) {
  return simd::Float(_mm_load_ps(x));
}
    
inline simd::Float load(const float* x) {
  return simd::Float(_mm_loadu_ps(x));
}
    
inline simd::Float load(const float* x, const int offset) {
  return simd::Float(_mm_setr_ps(x[0],x[offset],x[2*offset],x[3*offset]));
}
    
inline void store_aligned(float* x, const simd::Float& z) {
  _mm_store_ps(x,z._data);
}
    
inline void store(float* x, const simd::Float& z) {
  _mm_storeu_ps(x,z._data);
}
    
inline void store(float* x, const simd::Float& z, const int offset) {
  for (int i=0; i < nfloats; ++i) {
    x[offset*i] = z[i];
  }
}

} // namespace simd
} // namespace stk
