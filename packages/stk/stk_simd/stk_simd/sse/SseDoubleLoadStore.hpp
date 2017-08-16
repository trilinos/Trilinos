
namespace stk {
namespace simd {

inline simd::Double load_aligned(const double* x) {
  return simd::Double(_mm_load_pd(x));
}

inline simd::Double load(const double* x) {
  return simd::Double(_mm_loadu_pd(x));
}
    
inline simd::Double load(const double* x, const int offset) {
  return simd::Double(_mm_setr_pd(x[0], x[offset]));
}

inline void store_aligned(double* x, const simd::Double& z) {
  _mm_store_pd(x,z._data);
}

inline void store(double* x, const simd::Double& z) {
  _mm_storeu_pd(x,z._data);
}
  
inline void store(double* x, const simd::Double& z, const int offset) {
  for (int i=0; i < simd::ndoubles; ++i) {
    x[offset*i] = z[i];
  }
}

} // namespace simd
} // namespace stk
