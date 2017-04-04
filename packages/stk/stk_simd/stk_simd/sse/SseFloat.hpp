// Copyright 2013 Sandia Corporation, Albuquerque, NM.

namespace stk {
namespace simd {

struct Float 
{

  inline Float() 
    : _data(_mm_setzero_ps()) {
  }
  
  inline Float(const float* x) 
    : _data(_mm_loadu_ps(x)) {
  }

  inline Float(const float* x, const int offset) 
    :_data(_mm_setr_ps(x[0],x[offset],x[2*offset],x[3*offset])) {
  }

  inline Float(const volatile float x)
    : _data(_mm_set1_ps(x)) {
  }

  inline Float(const __m128& x)
    : _data(x) {
  }

  inline Float(const Float& x)
    : _data(x._data) {
  }

  inline Float& operator= (const Float& x) {
    _data = x._data;
    return *this;
  }

  inline Float& operator= (const volatile float x) {
    _data = _mm_set1_ps(x);
    return *this;
  }

  inline Float& operator+= (const Float& a) {
    _data = _mm_add_ps(_data,a._data);
    return *this;
  }

  inline Float& operator-= (const Float& a) {
    _data = _mm_sub_ps(_data,a._data);
    return *this;
  }

  inline Float& operator*= (const Float& a) {
    _data = _mm_mul_ps(_data,a._data);
    return *this;
  }

  inline Float& operator/= (const Float& a) {
    _data = _mm_div_ps(_data,a._data);
    return *this;
  }

  inline Float& operator+= (const float a) {
    _data = _mm_add_ps(_data,_mm_set1_ps(a));
    return *this;
  }

  inline Float& operator-= (const float a) {
    _data = _mm_sub_ps(_data,_mm_set1_ps(a));
    return *this;
  }

  inline Float& operator*= (const float a) {
    _data = _mm_mul_ps(_data,_mm_set1_ps(a));
    return *this;
  }

  inline Float& operator/= (const float a) {
    _data = _mm_div_ps(_data,_mm_set1_ps(a));
    return *this;
  }

  inline Float operator-() const {
    return Float( _mm_sub_ps(Float(0.0)._data,_data) );
  }

  inline float& operator[](int i) {return (reinterpret_cast<float*>(&_data))[i];}
  inline const float& operator[](int i) const {return (reinterpret_cast<const float*>(&_data))[i];}
    
  inline int32_t& Int(int i) {return (reinterpret_cast<int32_t*>(&_data))[i];}
  inline const int32_t& Int(int i) const {return (reinterpret_cast<const int32_t*>(&_data))[i];}

  inline uint32_t& UInt(int i) {return (reinterpret_cast<uint32_t*>(&_data))[i];}
  inline const uint32_t& UInt(int i) const {return (reinterpret_cast<const uint32_t*>(&_data))[i];}

  __m128 _data; // the "_" means you should try not to use this directly
  // it is made public to avoid function call overhead 
  // and/or so the compiler doesn't have to use up one of
  // inlining depths (usually max inlining depth ~5)
  
};

} // namespace simd
} // namespace stk

