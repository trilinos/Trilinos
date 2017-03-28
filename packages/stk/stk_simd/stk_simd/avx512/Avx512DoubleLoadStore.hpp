
namespace stk {
namespace simd {

inline simd::Double load_aligned(const double* x) {
  return simd::Double(_mm512_load_pd(x));
}

inline simd::Double load(const double* x) {
  return simd::Double(_mm512_loadu_pd(x));
}
    
inline simd::Double load(const double* x, const int offset) {
  return simd::Double(_mm512_setr_pd(x[0],       x[offset],  x[2*offset],x[3*offset],
                                      x[4*offset],x[5*offset],x[6*offset],x[7*offset]));
}
  
inline void store_aligned(double* x, const simd::Double& z) {
  _mm512_store_pd(x, z._data);
}

inline void store(double* x, const simd::Double& z) {
  _mm512_storeu_pd(x, z._data);
}
  
inline void store(double* x, const simd::Double& z, const int offset) {
  for (int i=0; i < ndoubles; ++i) {
    x[offset*i] = z[i];
  }
}

} // namespace simd
} // namespace stk
