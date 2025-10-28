// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
// 
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 

#ifndef STK_MESH_BASE_FIELDBLAS_HPP
#define STK_MESH_BASE_FIELDBLAS_HPP

#include <stk_util/stk_config.h>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/GetNgpField.hpp>
#include <stk_mesh/base/GetNgpMesh.hpp>
#include <stk_mesh/base/DeviceField.hpp>
#include <stk_util/util/Blas.hpp>
#include <complex>
#include <string>
#include <algorithm>
#include <type_traits>

namespace stk {
namespace mesh {

template<typename T> struct is_complex_t : public std::false_type {};
template<typename U> struct is_complex_t<std::complex<U>> : public std::true_type {};
template<typename T> constexpr bool is_complex_v = is_complex_t<T>::value;

template <typename T1>
constexpr bool is_layout_right() {
  return (T1::layout == Layout::Right);
}

template <typename T1, typename T2>
constexpr bool is_layout_right() {
  return (T1::layout == Layout::Right) && (T2::layout == Layout::Right);
}

template <typename T1, typename T2, typename T3>
constexpr bool is_layout_right() {
  return (T1::layout == Layout::Right) && (T2::layout == Layout::Right) && (T3::layout == Layout::Right);
}

template <typename FieldX, typename FieldY>
inline bool is_compatible(const FieldX& x, const FieldY& y)
{
  if (&x.get_mesh() != &y.get_mesh()) return false;
  if (x.entity_rank() != y.entity_rank()) return false;
  if (x.data_traits().type_info != y.data_traits().type_info) return false;
  return true;
}

template <typename T, typename FieldX>
inline bool is_compatible(const FieldX& x)
{
  if (x.data_traits().type_info != typeid(T)) return false;
  return true;
}

template <typename XVALS, typename YVALS>
void check_matching_extents(const std::string& functionName, XVALS xValues, YVALS yValues)
{
  STK_ThrowAssertMsg(xValues.num_scalars() == yValues.num_scalars(),
                     "Cannot perform " << functionName << " operation on different-length Fields.  fieldX has " <<
                     xValues.num_scalars() << " scalars and fieldY has " << yValues.num_scalars() << " scalars.");
}

template <typename XVALS, typename YVALS, typename ZVALS>
void check_matching_extents(const std::string& functionName, XVALS xValues, YVALS yValues, ZVALS zValues)
{
  STK_ThrowAssertMsg((xValues.num_scalars() == yValues.num_scalars()) && (yValues.num_scalars() == zValues.num_scalars()),
                     "Cannot perform " << functionName << " operation on different-length Fields.  fieldX has " <<
                     xValues.num_scalars() << " scalars, fieldY has " << yValues.num_scalars() <<
                     " scalars, and fieldZ has " << zValues.num_scalars() << " scalars.");
}

struct MinMaxInfo {
  unsigned bucketId;
  int entityIndex;
  int scalarIndex;
};

template <typename T>
struct FieldBlasImpl
{
  template <typename XDATA, typename YDATA>
  inline static void axpy(const BucketVector& buckets, T alpha, const XDATA& xData, const YDATA& yData)
  {
#ifdef STK_USE_OPENMP
    #pragma omp parallel for schedule(static)
#endif
    for (unsigned bucketIndex = 0; bucketIndex < buckets.size(); ++bucketIndex) {
      const Bucket& bucket = *buckets[bucketIndex];
      auto xValues = xData.bucket_values(bucket);
      auto yValues = yData.bucket_values(bucket);
      check_matching_extents("field_axpy()", xValues, yValues);

      if constexpr (std::is_floating_point_v<T>) {
        if constexpr (is_layout_right<XDATA, YDATA>()) {
          const int numValues = yValues.num_entities()*yValues.num_scalars();
          const int stride = xValues.scalar_stride();

          if constexpr (std::is_same_v<T, double>) {
            SIERRA_FORTRAN(daxpy)(&numValues, &alpha, xValues.pointer(), &stride, yValues.pointer(), &stride);
          }
          else if constexpr (std::is_same_v<T, float>) {
            SIERRA_FORTRAN(saxpy)(&numValues, &alpha, xValues.pointer(), &stride, yValues.pointer(), &stride);
          }
        }
        else {  // Layout::Left or mixed-layout cases
          // Stride isn't uniform through the whole Bucket for one or the other Field, so we have to chop
          // it up into uniform segments
          for (ScalarIdx scalar : yValues.scalars()) {
            const int numValues = yValues.num_entities();
            const int xStride = xValues.entity_stride();
            const int yStride = yValues.entity_stride();
            const T* xPtr = xValues.pointer() + scalar*xValues.scalar_stride();
            T* yPtr = yValues.pointer() + scalar*yValues.scalar_stride();

            if constexpr (std::is_same_v<T, double>) {
              SIERRA_FORTRAN(daxpy)(&numValues, &alpha, xPtr, &xStride, yPtr, &yStride);
            }
            else if constexpr (std::is_same_v<T, float>) {
              SIERRA_FORTRAN(saxpy)(&numValues, &alpha, xPtr, &xStride, yPtr, &yStride);
            }
          }
        }
      }
      else {  // All non-floating-point types
        for (stk::mesh::EntityIdx entity : bucket.entities()) {
          for (stk::mesh::ScalarIdx scalar : yValues.scalars()) {
            yValues(entity, scalar) += alpha * xValues(entity, scalar);
          }
        }
      }
    }
  }

  template <typename XDATA, typename YDATA>
  inline static void axpby(const BucketVector& buckets, T alpha, const XDATA& xData, T beta, const YDATA& yData)
  {
    if constexpr (is_layout_right<XDATA, YDATA>()) {
#ifdef STK_USE_OPENMP
      #pragma omp parallel for schedule(static)
#endif
      for (unsigned bucketIndex = 0; bucketIndex < buckets.size(); ++bucketIndex) {
        const Bucket& bucket = *buckets[bucketIndex];
        auto xValues = xData.bucket_values(bucket);
        auto yValues = yData.bucket_values(bucket);
        check_matching_extents("field_axpby()", xValues, yValues);

        const int numValues = xValues.num_entities() * xValues.num_scalars();
        const T* xPtr = xValues.pointer();
        T* yPtr = yValues.pointer();

        for (int i = 0; i < numValues; ++i) {
          yPtr[i] = beta * yPtr[i] + alpha * xPtr[i];
        }
      }
    }
    else {  // Layout::Left or mixed-layout cases
#ifdef STK_USE_OPENMP
      #pragma omp parallel for schedule(static)
#endif
      for (unsigned bucketIndex = 0; bucketIndex < buckets.size(); ++bucketIndex) {
        const Bucket& bucket = *buckets[bucketIndex];
        auto xValues = xData.bucket_values(bucket);
        auto yValues = yData.bucket_values(bucket);
        check_matching_extents("field_axpby()", xValues, yValues);

        for (stk::mesh::EntityIdx entity : bucket.entities()) {
          for (stk::mesh::ScalarIdx scalar : yValues.scalars()) {
            yValues(entity, scalar) = beta * yValues(entity, scalar) + alpha * xValues(entity, scalar);
          }
        }
      }
    }
  }

  template <typename XDATA, typename YDATA, typename ZDATA>
  inline static void product(const BucketVector& buckets, const XDATA& xData, const YDATA& yData, const ZDATA& zData)
  {
    if constexpr (is_layout_right<XDATA, YDATA, ZDATA>()) {
#ifdef STK_USE_OPENMP
      #pragma omp parallel for schedule(static)
#endif
      for (unsigned bucketIndex = 0; bucketIndex < buckets.size(); ++bucketIndex) {
        const Bucket& bucket = *buckets[bucketIndex];
        auto xValues = xData.bucket_values(bucket);
        auto yValues = yData.bucket_values(bucket);
        auto zValues = zData.bucket_values(bucket);
        check_matching_extents("field_product()", xValues, yValues, zValues);

        const int numValues = xValues.num_entities() * xValues.num_scalars();
        const T* xPtr = xValues.pointer();
        const T* yPtr = yValues.pointer();
        T* zPtr = zValues.pointer();

        for (int i = 0; i < numValues; ++i) {
          zPtr[i] = xPtr[i] * yPtr[i];
        }
      }
    }
    else {  // Layout::Left or mixed-layout cases
#ifdef STK_USE_OPENMP
      #pragma omp parallel for schedule(static)
#endif
      for (unsigned bucketIndex = 0; bucketIndex < buckets.size(); ++bucketIndex) {
        const Bucket& bucket = *buckets[bucketIndex];
        auto xValues = xData.bucket_values(bucket);
        auto yValues = yData.bucket_values(bucket);
        auto zValues = zData.bucket_values(bucket);
        check_matching_extents("field_product()", xValues, yValues, zValues);

        for (stk::mesh::EntityIdx entity : bucket.entities()) {
          for (stk::mesh::ScalarIdx scalar : yValues.scalars()) {
            zValues(entity, scalar) = xValues(entity, scalar) * yValues(entity, scalar);
          }
        }
      }
    }
  }

  template <typename XDATA, typename YDATA>
  inline static void copy(const BucketVector& buckets, const XDATA& xData, const YDATA& yData)
  {
    if constexpr (is_layout_right<XDATA, YDATA>()) {
#ifdef STK_USE_OPENMP
      #pragma omp parallel for schedule(static)
#endif
      for (unsigned bucketIndex = 0; bucketIndex < buckets.size(); ++bucketIndex) {
        const Bucket& bucket = *buckets[bucketIndex];
        auto xValues = xData.bucket_values(bucket);
        auto yValues = yData.bucket_values(bucket);
        check_matching_extents("field_copy()", xValues, yValues);

        const int numValues = xValues.num_entities() * xValues.num_scalars();
        const T* xPtr = xValues.pointer();
        T* yPtr = yValues.pointer();

        for (int i = 0; i < numValues; ++i) {
          yPtr[i] = xPtr[i];
        }
      }
    }
    else {  // Layout::Left or mixed-layout cases
#ifdef STK_USE_OPENMP
      #pragma omp parallel for schedule(static)
#endif
      for (unsigned bucketIndex = 0; bucketIndex < buckets.size(); ++bucketIndex) {
        const Bucket& bucket = *buckets[bucketIndex];
        auto xValues = xData.bucket_values(bucket);
        auto yValues = yData.bucket_values(bucket);
        check_matching_extents("field_copy()", xValues, yValues);

        for (stk::mesh::EntityIdx entity : bucket.entities()) {
          for (stk::mesh::ScalarIdx scalar : yValues.scalars()) {
            yValues(entity, scalar) = xValues(entity, scalar);
          }
        }
      }
    }
  }

  template <typename XDATA, typename YDATA>
  inline static T dot(const BucketVector& buckets, const XDATA& xData, const YDATA& yData)
  {
    if constexpr (is_complex_v<T>) {
      using ComplexType = typename T::value_type;
      ComplexType localResultReal {};  // OpenMP can't do reductions on std::complex types
      ComplexType localResultImag {};

#ifdef STK_USE_OPENMP
      #pragma omp parallel for schedule(static) reduction(+:localResultReal,localResultImag)
#endif
      for (unsigned bucketIndex = 0; bucketIndex < buckets.size(); ++bucketIndex) {
        const Bucket& bucket = *buckets[bucketIndex];
        auto xValues = xData.bucket_values(bucket);
        auto yValues = yData.bucket_values(bucket);
        check_matching_extents("field_dot()", xValues, yValues);

        for (stk::mesh::EntityIdx entity : bucket.entities()) {
          for (stk::mesh::ScalarIdx scalar : yValues.scalars()) {
            T localResult = xValues(entity, scalar) * yValues(entity, scalar);
            localResultReal += localResult.real();
            localResultImag += localResult.imag();
          }
        }
      }

      return T(localResultReal, localResultImag);
    }
    else {  // All non-complex types
      T localResult {};

#ifdef STK_USE_OPENMP
      #pragma omp parallel for reduction(+:localResult) schedule(static)
#endif
      for (unsigned bucketIndex = 0; bucketIndex < buckets.size(); ++bucketIndex) {
        const Bucket& bucket = *buckets[bucketIndex];
        auto xValues = xData.bucket_values(bucket);
        auto yValues = yData.bucket_values(bucket);
        check_matching_extents("field_dot()", xValues, yValues);

        if constexpr (std::is_floating_point_v<T>) {
          if constexpr (is_layout_right<XDATA, YDATA>()) {
            const int numValues = yValues.num_entities()*yValues.num_scalars();
            const int stride = xValues.scalar_stride();

            if constexpr (std::is_same_v<T, double>) {
              localResult += SIERRA_FORTRAN(ddot)(&numValues, xValues.pointer(), &stride, yValues.pointer(), &stride);
            }
            else if constexpr (std::is_same_v<T, float>) {
              localResult += SIERRA_FORTRAN(sdot)(&numValues, xValues.pointer(), &stride, yValues.pointer(), &stride);
            }
          }
          else {  // Layout::Left or mixed-layout cases
            // Stride isn't uniform through the whole Bucket for one or the other Field, so we have to chop
            // it up into uniform segments
            for (ScalarIdx scalar : xValues.scalars()) {
              const int numValues = xValues.num_entities();
              const int xStride = xValues.entity_stride();
              const int yStride = yValues.entity_stride();
              const T* xPtr = xValues.pointer() + scalar*xValues.scalar_stride();
              const T* yPtr = yValues.pointer() + scalar*yValues.scalar_stride();

              if constexpr (std::is_same_v<T, double>) {
                localResult += SIERRA_FORTRAN(ddot)(&numValues, xPtr, &xStride, yPtr, &yStride);
              }
              else if constexpr (std::is_same_v<T, float>) {
                localResult += SIERRA_FORTRAN(sdot)(&numValues, xPtr, &xStride, yPtr, &yStride);
              }
            }
          }
        }
        else {  // All non-floating-point types
          for (stk::mesh::EntityIdx entity : bucket.entities()) {
            for (stk::mesh::ScalarIdx scalar : yValues.scalars()) {
              localResult += xValues(entity, scalar) * yValues(entity, scalar);
            }
          }
        }
      }

      return localResult;
    }
  }

  template <typename XDATA>
  inline static T nrm2(const BucketVector& buckets, const XDATA& xData)
  {
    if constexpr (is_complex_v<T>) {
      using ComplexType = typename T::value_type;
      ComplexType localResultReal {};

#ifdef STK_USE_OPENMP
      #pragma omp parallel for schedule(static) reduction(+:localResultReal)
#endif
      for (unsigned bucketIndex = 0; bucketIndex < buckets.size(); ++bucketIndex) {
        const Bucket& bucket = *buckets[bucketIndex];
        auto xValues = xData.bucket_values(bucket);

        for (stk::mesh::EntityIdx entity : bucket.entities()) {
          for (stk::mesh::ScalarIdx scalar : xValues.scalars()) {
            T localResult = std::pow(std::abs(xValues(entity, scalar)), 2);
            localResultReal += localResult.real();
          }
        }
      }

      return T(localResultReal, 0);
    }
    else {  // All non-complex types
      T localResult {};

#ifdef STK_USE_OPENMP
      #pragma omp parallel for reduction(+:localResult) schedule(static)
#endif
      for (unsigned bucketIndex = 0; bucketIndex < buckets.size(); ++bucketIndex) {
        const Bucket& bucket = *buckets[bucketIndex];
        auto xValues = xData.bucket_values(bucket);

        if constexpr (std::is_floating_point_v<T>) {
          if constexpr (is_layout_right<XDATA>()) {
            const int numValues = xValues.num_entities()*xValues.num_scalars();
            const int stride = xValues.scalar_stride();
            if constexpr (std::is_same_v<T, double>) {
              localResult += SIERRA_FORTRAN(ddot)(&numValues, xValues.pointer(), &stride, xValues.pointer(), &stride);
            }
            else if constexpr (std::is_same_v<T, float>) {
              localResult += SIERRA_FORTRAN(sdot)(&numValues, xValues.pointer(), &stride, xValues.pointer(), &stride);
            }
          }
          else {  // Layout::Left case
            // Stride isn't uniform through the whole Bucket, so we have to chop it up into uniform segments
            for (ScalarIdx scalar : xValues.scalars()) {
              const int numValues = xValues.num_entities();
              const int xStride = xValues.entity_stride();
              const T* xPtr = xValues.pointer() + scalar*xValues.scalar_stride();
              if constexpr (std::is_same_v<T, double>) {
                localResult += SIERRA_FORTRAN(ddot)(&numValues, xPtr, &xStride, xPtr, &xStride);
              }
              else if constexpr (std::is_same_v<T, float>) {
                localResult += SIERRA_FORTRAN(sdot)(&numValues, xPtr, &xStride, xPtr, &xStride);
              }
            }
          }
        }
        else {  // All non-floating-point types
          for (stk::mesh::EntityIdx entity : bucket.entities()) {
            for (stk::mesh::ScalarIdx scalar : xValues.scalars()) {
              localResult += xValues(entity, scalar) * xValues(entity, scalar);
            }
          }
        }
      }

      return localResult;
    }
  }

  template <typename XDATA>
  inline static void scale(const BucketVector& buckets, T alpha, const XDATA& xData)
  {
#ifdef STK_USE_OPENMP
    #pragma omp parallel for schedule(static)
#endif
    for (unsigned bucketIndex = 0; bucketIndex < buckets.size(); ++bucketIndex) {
      const Bucket& bucket = *buckets[bucketIndex];
      auto xValues = xData.bucket_values(bucket);

      if constexpr (std::is_floating_point_v<T>) {
        if constexpr (is_layout_right<XDATA>()) {
          const int numValues = xValues.num_entities()*xValues.num_scalars();
          const int stride = xValues.scalar_stride();

          if constexpr (std::is_same_v<T, double>) {
            SIERRA_FORTRAN(dscal)(&numValues, &alpha, xValues.pointer(), &stride);
          }
          else if constexpr (std::is_same_v<T, float>) {
            SIERRA_FORTRAN(sscal)(&numValues, &alpha, xValues.pointer(), &stride);
          }
        }
        else {  // Layout::Left case
          // Stride isn't uniform through the whole Bucket, so we have to chop it up into uniform segments
          for (ScalarIdx scalar : xValues.scalars()) {
            const int numValues = xValues.num_entities();
            const int xStride = xValues.entity_stride();
            T* xPtr = xValues.pointer() + scalar*xValues.scalar_stride();

            if constexpr (std::is_same_v<T, double>) {
              SIERRA_FORTRAN(dscal)(&numValues, &alpha, xPtr, &xStride);
            }
            else if constexpr (std::is_same_v<T, float>) {
              SIERRA_FORTRAN(sscal)(&numValues, &alpha, xPtr, &xStride);
            }
          }
        }
      }
      else {  // All non-floating-point types
        for (stk::mesh::EntityIdx entity : bucket.entities()) {
          for (stk::mesh::ScalarIdx scalar : xValues.scalars()) {
            xValues(entity, scalar) = alpha * xValues(entity, scalar);
          }
        }
      }
    }
  }

  template <typename XDATA>
  inline static void fill(const BucketVector& buckets, T alpha, const XDATA& xData)
  {
#ifdef STK_USE_OPENMP
    #pragma omp parallel for schedule(static)
#endif
    for (unsigned bucketIndex = 0; bucketIndex < buckets.size(); ++bucketIndex) {
      const Bucket& bucket = *buckets[bucketIndex];
      auto xValues = xData.bucket_values(bucket);

      if constexpr (is_layout_right<XDATA>()) {
        const int numValues = xValues.num_entities() * xValues.num_scalars();
        T* xPtr = xValues.pointer();
        std::fill(xPtr, xPtr+numValues, alpha);
      }
      else {  // Layout::Left case
        for (stk::mesh::EntityIdx entity : bucket.entities()) {
          for (stk::mesh::ScalarIdx scalar : xValues.scalars()) {
            xValues(entity, scalar) = alpha;
          }
        }
      }
    }
  }

  template <typename XDATA>
  inline static void fill_component(const BucketVector& buckets, const T* alpha, const XDATA& xData)
  {
#ifdef STK_USE_OPENMP
    #pragma omp parallel for schedule(static)
#endif
    for (unsigned bucketIndex = 0; bucketIndex < buckets.size(); ++bucketIndex) {
      const Bucket& bucket = *buckets[bucketIndex];
      auto xValues = xData.bucket_values(bucket);

      for (stk::mesh::EntityIdx entity : bucket.entities()) {
        for (stk::mesh::ScalarIdx scalar : xValues.scalars()) {
          xValues(entity, scalar) = alpha[scalar];
        }
      }
    }
  }

  template <typename XDATA, typename YDATA>
  inline static void swap(const BucketVector& buckets, const XDATA& xData, const YDATA& yData)
  {
#ifdef STK_USE_OPENMP
    #pragma omp parallel for schedule(static)
#endif
    for (unsigned bucketIndex = 0; bucketIndex < buckets.size(); ++bucketIndex) {
      const Bucket& bucket = *buckets[bucketIndex];
      auto xValues = xData.bucket_values(bucket);
      auto yValues = yData.bucket_values(bucket);
      check_matching_extents("field_swap()", xValues, yValues);

      if constexpr (std::is_floating_point_v<T>) {
        if constexpr (is_layout_right<XDATA, YDATA>()) {
          const int numValues = xValues.num_entities()*xValues.num_scalars();
          const int stride = xValues.scalar_stride();

          if constexpr (std::is_same_v<T, double>) {
            SIERRA_FORTRAN(dswap)(&numValues, xValues.pointer(), &stride, yValues.pointer(), &stride);
          }
          else if constexpr (std::is_same_v<T, float>) {
            SIERRA_FORTRAN(sswap)(&numValues, xValues.pointer(), &stride, yValues.pointer(), &stride);
          }
        }
        else {  // Layout::Left and mixed-layout cases
          // Stride isn't uniform through the whole Bucket for one or the other Field, so we have to chop
          // it up into uniform segments
          for (ScalarIdx scalar : yValues.scalars()) {
            const int numValues = yValues.num_entities();
            const int xStride = xValues.entity_stride();
            const int yStride = yValues.entity_stride();
            T* xPtr = xValues.pointer() + scalar*xValues.scalar_stride();
            T* yPtr = yValues.pointer() + scalar*yValues.scalar_stride();

            if constexpr (std::is_same_v<T, double>) {
              SIERRA_FORTRAN(dswap)(&numValues, xPtr, &xStride, yPtr, &yStride);
            }
            else if constexpr (std::is_same_v<T, float>) {
              SIERRA_FORTRAN(sswap)(&numValues, xPtr, &xStride, yPtr, &yStride);
            }
          }
        }
      }
      else {  // All non-floating-point types
        for (stk::mesh::EntityIdx entity : bucket.entities()) {
          for (stk::mesh::ScalarIdx scalar : yValues.scalars()) {
            T temp = yValues(entity, scalar);
            yValues(entity, scalar) = xValues(entity, scalar);
            xValues(entity, scalar) = temp;
          }
        }
      }
    }
  }

  template <typename XDATA>
  inline static T asum(const BucketVector& buckets, const XDATA& xData)
  {
    if constexpr (is_complex_v<T>) {
      using ComplexType = typename T::value_type;
      ComplexType localResultReal {};

#ifdef STK_USE_OPENMP
      #pragma omp parallel for schedule(static) reduction(+:localResultReal)
#endif
      for (unsigned bucketIndex = 0; bucketIndex < buckets.size(); ++bucketIndex) {
        const Bucket& bucket = *buckets[bucketIndex];
        auto xValues = xData.bucket_values(bucket);

        for (stk::mesh::EntityIdx entity : bucket.entities()) {
          for (stk::mesh::ScalarIdx scalar : xValues.scalars()) {
            T localResult = std::abs(xValues(entity, scalar));
            localResultReal += localResult.real();
          }
        }
      }

      return T(localResultReal, 0);
    }
    else {  // All non-complex types
      T localResult {};

#ifdef STK_USE_OPENMP
      #pragma omp parallel for schedule(static) reduction(+:localResult)
#endif
      for (unsigned bucketIndex = 0; bucketIndex < buckets.size(); ++bucketIndex) {
        const Bucket& bucket = *buckets[bucketIndex];
        auto xValues = xData.bucket_values(bucket);

        if constexpr (std::is_floating_point_v<T>) {
          if constexpr (is_layout_right<XDATA>()) {
            const int numValues = xValues.num_entities()*xValues.num_scalars();
            const int stride = xValues.scalar_stride();

            if constexpr (std::is_same_v<T, double>) {
              localResult += SIERRA_FORTRAN(dasum)(&numValues, xValues.pointer(), &stride);
            }
            else if constexpr (std::is_same_v<T, float>) {
              localResult += SIERRA_FORTRAN(sasum)(&numValues, xValues.pointer(), &stride);
            }
          }
          else {  // Layout::Left and mixed-layout cases
            // Stride isn't uniform through the whole Bucket, so we have to chop it up into uniform segments
            for (ScalarIdx scalar : xValues.scalars()) {
              const int numValues = xValues.num_entities();
              const int xStride = xValues.entity_stride();
              const T* xPtr = xValues.pointer() + scalar*xValues.scalar_stride();

              if constexpr (std::is_same_v<T, double>) {
                localResult += SIERRA_FORTRAN(dasum)(&numValues, xPtr, &xStride);
              }
              else if constexpr (std::is_same_v<T, float>) {
                localResult += SIERRA_FORTRAN(sasum)(&numValues, xPtr, &xStride);
              }
            }
          }
        }
        else {  // All non-floating-point types
          for (stk::mesh::EntityIdx entity : bucket.entities()) {
            for (stk::mesh::ScalarIdx scalar : xValues.scalars()) {
              localResult += std::abs(xValues(entity, scalar));
            }
          }
        }
      }

      return localResult;
    }
  }

  template <typename XDATA>
  inline static T amax(const BucketVector& buckets, const XDATA& xData)
  {
    if constexpr (is_complex_v<T>) {
      using ComplexType = typename T::value_type;
      ComplexType localMaxValueReal = {};

#ifdef STK_USE_OPENMP
      #pragma omp parallel for schedule(static) reduction(max:localMaxValueReal)
#endif
      for (unsigned bucketIndex = 0; bucketIndex < buckets.size(); ++bucketIndex) {
        const Bucket& bucket = *buckets[bucketIndex];
        auto xValues = xData.bucket_values(bucket);

        for (stk::mesh::EntityIdx entity : bucket.entities()) {
          for (stk::mesh::ScalarIdx scalar : xValues.scalars()) {
            ComplexType localValueReal = std::abs(xValues(entity, scalar));
            if (localValueReal > localMaxValueReal) {
              localMaxValueReal = localValueReal;
            }
          }
        }
      }

      return T(localMaxValueReal, 0);
    }
    else {  // All non-complex types
      T localMaxValue {};

#ifdef STK_USE_OPENMP
      #pragma omp parallel for schedule(static) reduction(max:localMaxValue)
#endif
      for (unsigned bucketIndex = 0; bucketIndex < buckets.size(); ++bucketIndex) {
        const Bucket& bucket = *buckets[bucketIndex];
        auto xValues = xData.bucket_values(bucket);
        if (xValues.is_field_defined()) {
          if constexpr (std::is_floating_point_v<T> && is_layout_right<XDATA>()) {
            const int numValues = xValues.num_entities()*xValues.num_scalars();
            const int stride = xValues.scalar_stride();
            int localIndex {};
            if constexpr (std::is_same_v<T, double>) {
              localIndex = SIERRA_FORTRAN(idamax)(&numValues, xValues.pointer(), &stride) - 1;
            }
            else if constexpr (std::is_same_v<T, float>) {
              localIndex = SIERRA_FORTRAN(isamax)(&numValues, xValues.pointer(), &stride) - 1;
            }
            EntityIdx localEntityIdx(localIndex / xValues.num_scalars());
            ScalarIdx localScalarIdx(localIndex % xValues.num_scalars());
            T localValue = std::abs(xValues(localEntityIdx, localScalarIdx));

            if (localValue > localMaxValue) {
              localMaxValue = localValue;
            }
          }
          else {  // All non-floating-point types and Layout::Left cases
            for (stk::mesh::EntityIdx entityIdx : bucket.entities()) {
              for (stk::mesh::ScalarIdx scalar : xValues.scalars()) {
                T localValue = std::abs(xValues(entityIdx, scalar));

                if (localValue > localMaxValue) {
                  localMaxValue = localValue;
                }
              }
            }
          }
        }
      }

      return localMaxValue;
    }
  }


  // Note: Can't include the maxValue information in the MinMaxInfo struct return type due to
  // un-suppressable gcc warnings about ABI changes for std::complex<float> in structs.  So,
  // return it separately as an out-parameter.

  template <typename XDATA>
  inline static MinMaxInfo iamax(const BucketVector& buckets, const XDATA& xData, T& maxValue)
  {
    if constexpr (is_complex_v<T>) {
      using ComplexType = typename T::value_type;
      ComplexType globalMaxValueReal {};
      MinMaxInfo globalMinMaxInfo {InvalidOrdinal, 0, 0};

#ifdef STK_USE_OPENMP
      #pragma omp parallel
#endif
      {
        ComplexType localMaxValueReal {};
        MinMaxInfo localMinMaxInfo {InvalidOrdinal, 0, 0};

#ifdef STK_USE_OPENMP
        #pragma omp for schedule(static)
#endif
        for (unsigned bucketIndex = 0; bucketIndex < buckets.size(); ++bucketIndex) {
          const Bucket& bucket = *buckets[bucketIndex];
          auto xValues = xData.bucket_values(bucket);

          for (stk::mesh::EntityIdx entityIdx : bucket.entities()) {
            for (stk::mesh::ScalarIdx scalar : xValues.scalars()) {
              ComplexType localValueReal = std::abs(xValues(entityIdx, scalar));
              if (localValueReal > localMaxValueReal) {
                localMaxValueReal = localValueReal;
                localMinMaxInfo.bucketId = bucket.bucket_id();
                localMinMaxInfo.entityIndex = entityIdx;
                localMinMaxInfo.scalarIndex = scalar;
              }
            }
          }
        }

#ifdef STK_USE_OPENMP
        #pragma omp critical
#endif
        if (localMaxValueReal > globalMaxValueReal) {
          globalMaxValueReal = localMaxValueReal;
          globalMinMaxInfo = localMinMaxInfo;
        }
      }

      maxValue = T(globalMaxValueReal, 0);
      return globalMinMaxInfo;
    }
    else {  // All non-complex types
      T globalMaxValue {};
      MinMaxInfo globalMinMaxInfo {InvalidOrdinal, 0, 0};

#ifdef STK_USE_OPENMP
      #pragma omp parallel
#endif
      {
        T localMaxValue {};
        MinMaxInfo localMinMaxInfo {InvalidOrdinal, 0, 0};

#ifdef STK_USE_OPENMP
        #pragma omp for schedule(static)
#endif
        for (unsigned bucketIndex = 0; bucketIndex < buckets.size(); ++bucketIndex) {
          const Bucket& bucket = *buckets[bucketIndex];
          auto xValues = xData.bucket_values(bucket);

          if (xValues.is_field_defined()) {
            if constexpr (std::is_floating_point_v<T> && is_layout_right<XDATA>()) {
              const int numValues = xValues.num_entities()*xValues.num_scalars();
              const int stride = xValues.scalar_stride();
              int localIndex {};
              if constexpr (std::is_same_v<T, double>) {
                localIndex = SIERRA_FORTRAN(idamax)(&numValues, xValues.pointer(), &stride) - 1;
              }
              else if constexpr (std::is_same_v<T, float>) {
                localIndex = SIERRA_FORTRAN(isamax)(&numValues, xValues.pointer(), &stride) - 1;
              }
              EntityIdx localEntityIdx(localIndex / xValues.num_scalars());
              ScalarIdx localScalarIdx(localIndex % xValues.num_scalars());
              T localValue = std::abs(xValues(localEntityIdx, localScalarIdx));

              if (localValue > localMaxValue) {
                localMaxValue = localValue;
                localMinMaxInfo.bucketId = bucket.bucket_id();
                localMinMaxInfo.entityIndex = localEntityIdx;
                localMinMaxInfo.scalarIndex = localScalarIdx;
              }
            }
            else {  // All non-floating-point types and Layout::Left cases
              for (stk::mesh::EntityIdx entityIdx : bucket.entities()) {
                for (stk::mesh::ScalarIdx scalar : xValues.scalars()) {
                  T localValue = std::abs(xValues(entityIdx, scalar));

                  if (localValue > localMaxValue) {
                    localMaxValue = localValue;
                    localMinMaxInfo.bucketId = bucket.bucket_id();
                    localMinMaxInfo.entityIndex = entityIdx;
                    localMinMaxInfo.scalarIndex = scalar;
                  }
                }
              }
            }
          }
        }

#ifdef STK_USE_OPENMP
        #pragma omp critical
#endif
        if (localMaxValue > globalMaxValue) {
          globalMaxValue = localMaxValue;
          globalMinMaxInfo = localMinMaxInfo;
        }
      }

      maxValue = globalMaxValue;
      return globalMinMaxInfo;
    }
  }

  template <typename XDATA>
  inline static T amin(const BucketVector& buckets, const XDATA& xData)
  {
    if constexpr (is_complex_v<T>) {
      using ComplexType = typename T::value_type;
      ComplexType localMinValueReal = std::numeric_limits<ComplexType>::max();

#ifdef STK_USE_OPENMP
      #pragma omp parallel for schedule(static) reduction(min:localMinValueReal)
#endif
      for (unsigned bucketIndex = 0; bucketIndex < buckets.size(); ++bucketIndex) {
        const Bucket& bucket = *buckets[bucketIndex];
        auto xValues = xData.bucket_values(bucket);

        for (stk::mesh::EntityIdx entity : bucket.entities()) {
          for (stk::mesh::ScalarIdx scalar : xValues.scalars()) {
            ComplexType localValueReal = std::abs(xValues(entity, scalar));
            if (localValueReal < localMinValueReal) {
              localMinValueReal = localValueReal;
            }
          }
        }
      }

      return T(localMinValueReal, 0);
    }
    else {  // All non-complex types
      T localMinValue = std::numeric_limits<T>::max();

#ifdef STK_USE_OPENMP
      #pragma omp parallel for schedule(static) reduction(min:localMinValue)
#endif
      for (unsigned bucketIndex = 0; bucketIndex < buckets.size(); ++bucketIndex) {
        const Bucket& bucket = *buckets[bucketIndex];
        auto xValues = xData.bucket_values(bucket);

        for (stk::mesh::EntityIdx entityIdx : bucket.entities()) {
          for (stk::mesh::ScalarIdx scalar : xValues.scalars()) {
            T localValue = std::abs(xValues(entityIdx, scalar));

            if (localValue < localMinValue) {
              localMinValue = localValue;
            }
          }
        }
      }

      return localMinValue;
    }
  }

  template <typename XDATA>
  inline static MinMaxInfo iamin(const BucketVector& buckets, const XDATA& xData, T& minValue)
  {
    if constexpr (is_complex_v<T>) {
      using ComplexType = typename T::value_type;
      ComplexType globalMinValueReal = std::numeric_limits<ComplexType>::max();
      MinMaxInfo globalMinMaxInfo {InvalidOrdinal, 0, 0};

#ifdef STK_USE_OPENMP
      #pragma omp parallel
#endif
      {
        ComplexType localMinValueReal = std::numeric_limits<ComplexType>::max();
        MinMaxInfo localMinMaxInfo {InvalidOrdinal, 0, 0};

#ifdef STK_USE_OPENMP
        #pragma omp for schedule(static)
#endif
        for (unsigned bucketIndex = 0; bucketIndex < buckets.size(); ++bucketIndex) {
          const Bucket& bucket = *buckets[bucketIndex];
          auto xValues = xData.bucket_values(bucket);

          for (stk::mesh::EntityIdx entityIdx : bucket.entities()) {
            for (stk::mesh::ScalarIdx scalar : xValues.scalars()) {
              ComplexType localValueReal = std::abs(xValues(entityIdx, scalar));

              if (localValueReal < localMinValueReal) {
                localMinValueReal = localValueReal;
                localMinMaxInfo.bucketId = bucket.bucket_id();
                localMinMaxInfo.entityIndex = entityIdx;
                localMinMaxInfo.scalarIndex = scalar;
              }
            }
          }
        }

#ifdef STK_USE_OPENMP
        #pragma omp critical
#endif
        if (localMinValueReal < globalMinValueReal) {
          globalMinValueReal = localMinValueReal;
          globalMinMaxInfo = localMinMaxInfo;
        }
      }

      minValue = T(globalMinValueReal, 0);
      return globalMinMaxInfo;
    }
    else {  // All non-complex types
      T globalMinValue = std::numeric_limits<T>::max();
      MinMaxInfo globalMinMaxInfo {InvalidOrdinal, 0, 0};

#ifdef STK_USE_OPENMP
      #pragma omp parallel
#endif
      {
        T localMinValue = std::numeric_limits<T>::max();
        MinMaxInfo localMinMaxInfo {InvalidOrdinal, 0, 0};

#ifdef STK_USE_OPENMP
        #pragma omp for schedule(static)
#endif
        for (unsigned bucketIndex = 0; bucketIndex < buckets.size(); ++bucketIndex) {
          const Bucket& bucket = *buckets[bucketIndex];
          auto xValues = xData.bucket_values(bucket);

          for (stk::mesh::EntityIdx entityIdx : bucket.entities()) {
            for (stk::mesh::ScalarIdx scalar : xValues.scalars()) {
              T localValue = std::abs(xValues(entityIdx, scalar));

              if (localValue < localMinValue) {
                localMinValue = localValue;
                localMinMaxInfo.bucketId = bucket.bucket_id();
                localMinMaxInfo.entityIndex = entityIdx;
                localMinMaxInfo.scalarIndex = scalar;
              }
            }
          }
        }

#ifdef STK_USE_OPENMP
        #pragma omp critical
#endif
        if (localMinValue < globalMinValue) {
          globalMinValue = localMinValue;
          globalMinMaxInfo = localMinMaxInfo;
        }
      }

      minValue = globalMinValue;
      return globalMinMaxInfo;
    }
  }

};


//==============================================================================
// axpy: y[i] = a*x[i] + y[i]
//
template <typename T>
inline
void field_axpy(T alpha, const FieldBase& xField, const FieldBase& yField, const Selector& selector)
{
  STK_ThrowRequire(is_compatible<T>(xField));
  STK_ThrowRequire(is_compatible(xField, yField));

  const BucketVector& buckets = xField.get_mesh().get_buckets(xField.entity_rank(), selector);

  if (xField.host_data_layout() == Layout::Right && yField.host_data_layout() == Layout::Right) {
    auto xData = xField.data<T, ConstUnsynchronized, stk::ngp::HostSpace, Layout::Right>();
    auto yData = yField.data<T, Unsynchronized, stk::ngp::HostSpace, Layout::Right>();
    FieldBlasImpl<T>::axpy(buckets, alpha, xData, yData);
  }
  else if (xField.host_data_layout() == Layout::Left && yField.host_data_layout() == Layout::Left) {
    auto xData = xField.data<T, ConstUnsynchronized, stk::ngp::HostSpace, Layout::Left>();
    auto yData = yField.data<T, Unsynchronized, stk::ngp::HostSpace, Layout::Left>();
    FieldBlasImpl<T>::axpy(buckets, alpha, xData, yData);
  }
  else if (xField.host_data_layout() == Layout::Left && yField.host_data_layout() == Layout::Right) {
    auto xData = xField.data<T, ConstUnsynchronized, stk::ngp::HostSpace, Layout::Left>();
    auto yData = yField.data<T, Unsynchronized, stk::ngp::HostSpace, Layout::Right>();
    FieldBlasImpl<T>::axpy(buckets, alpha, xData, yData);
  }
  else if (xField.host_data_layout() == Layout::Right && yField.host_data_layout() == Layout::Left) {
    auto xData = xField.data<T, ConstUnsynchronized, stk::ngp::HostSpace, Layout::Right>();
    auto yData = yField.data<T, Unsynchronized, stk::ngp::HostSpace, Layout::Left>();
    FieldBlasImpl<T>::axpy(buckets, alpha, xData, yData);
  }
  else {
    STK_ThrowErrorMsg("Unsupported Field data layout detected in field_axpy().  xField layout: "
                      << xField.host_data_layout() << " and yField layout: " << yField.host_data_layout());
  }
}

template <typename T>
inline
void field_axpy(T alpha, const FieldBase& xField, const FieldBase& yField)
{
  const Selector selector = selectField(xField) & selectField(yField);
  field_axpy(alpha, xField, yField, selector);
}


//==============================================================================
// axpby: y[i] = a*x[i] + b*y[i]
//
template <typename T>
inline
void field_axpby(T alpha, const FieldBase& xField, T beta, const FieldBase& yField, const Selector& selector)
{
  STK_ThrowRequire(is_compatible<T>(xField));
  STK_ThrowRequire(is_compatible(xField, yField));

  const BucketVector& buckets = xField.get_mesh().get_buckets(xField.entity_rank(), selector);

  if ( beta == T(1) ) {
    if (xField.host_data_layout() == Layout::Right && yField.host_data_layout() == Layout::Right) {
      auto xData = xField.data<T, ConstUnsynchronized, stk::ngp::HostSpace, Layout::Right>();
      auto yData = yField.data<T, Unsynchronized, stk::ngp::HostSpace, Layout::Right>();
      FieldBlasImpl<T>::axpy(buckets, alpha, xData, yData);
    }
    else if (xField.host_data_layout() == Layout::Left && yField.host_data_layout() == Layout::Left) {
      auto xData = xField.data<T, ConstUnsynchronized, stk::ngp::HostSpace, Layout::Left>();
      auto yData = yField.data<T, Unsynchronized, stk::ngp::HostSpace, Layout::Left>();
      FieldBlasImpl<T>::axpy(buckets, alpha, xData, yData);
    }
    else if (xField.host_data_layout() == Layout::Right && yField.host_data_layout() == Layout::Left) {
      auto xData = xField.data<T, ConstUnsynchronized, stk::ngp::HostSpace, Layout::Right>();
      auto yData = yField.data<T, Unsynchronized, stk::ngp::HostSpace, Layout::Left>();
      FieldBlasImpl<T>::axpy(buckets, alpha, xData, yData);
    }
    else if (xField.host_data_layout() == Layout::Left && yField.host_data_layout() == Layout::Right) {
      auto xData = xField.data<T, ConstUnsynchronized, stk::ngp::HostSpace, Layout::Left>();
      auto yData = yField.data<T, Unsynchronized, stk::ngp::HostSpace, Layout::Right>();
      FieldBlasImpl<T>::axpy(buckets, alpha, xData, yData);
    }
    else {
      STK_ThrowErrorMsg("Unsupported Field data layout detected in field_axpby().  xField layout: "
                        << xField.host_data_layout() << " and yField layout: " << yField.host_data_layout());
    }
  }
  else {
    if (xField.host_data_layout() == Layout::Right && yField.host_data_layout() == Layout::Right) {
      auto xData = xField.data<T, ConstUnsynchronized, stk::ngp::HostSpace, Layout::Right>();
      auto yData = yField.data<T, Unsynchronized, stk::ngp::HostSpace, Layout::Right>();
      FieldBlasImpl<T>::axpby(buckets, alpha, xData, beta, yData);
    }
    else if (xField.host_data_layout() == Layout::Left && yField.host_data_layout() == Layout::Left) {
      auto xData = xField.data<T, ConstUnsynchronized, stk::ngp::HostSpace, Layout::Left>();
      auto yData = yField.data<T, Unsynchronized, stk::ngp::HostSpace, Layout::Left>();
      FieldBlasImpl<T>::axpby(buckets, alpha, xData, beta, yData);
    }
    else if (xField.host_data_layout() == Layout::Right && yField.host_data_layout() == Layout::Left) {
      auto xData = xField.data<T, ConstUnsynchronized, stk::ngp::HostSpace, Layout::Right>();
      auto yData = yField.data<T, Unsynchronized, stk::ngp::HostSpace, Layout::Left>();
      FieldBlasImpl<T>::axpby(buckets, alpha, xData, beta, yData);
    }
    else if (xField.host_data_layout() == Layout::Left && yField.host_data_layout() == Layout::Right) {
      auto xData = xField.data<T, ConstUnsynchronized, stk::ngp::HostSpace, Layout::Left>();
      auto yData = yField.data<T, Unsynchronized, stk::ngp::HostSpace, Layout::Right>();
      FieldBlasImpl<T>::axpby(buckets, alpha, xData, beta, yData);
    }
    else {
      STK_ThrowErrorMsg("Unsupported Field data layout detected in field_axpby().  xField layout: "
                        << xField.host_data_layout() << " and yField layout: " << yField.host_data_layout());
    }
  }
}

template <typename T>
inline
void field_axpby(T alpha, const FieldBase& xField, T beta, const FieldBase& yField)
{
  const Selector selector = selectField(xField) & selectField(yField);
  field_axpby(alpha, xField, beta, yField, selector);
}


//==============================================================================
// product: z[i] = x[i] * y[i]
//
template <typename T>
inline
void field_product(const FieldBase& xField, const FieldBase& yField, const FieldBase& zField, const Selector& selector)
{
  STK_ThrowRequire(is_compatible<T>(xField));
  STK_ThrowRequire(is_compatible(xField, yField));
  STK_ThrowRequire(is_compatible(yField, zField));

  const BucketVector& buckets = xField.get_mesh().get_buckets(xField.entity_rank(), selector);

  const stk::mesh::Layout xLayout = xField.host_data_layout();
  const stk::mesh::Layout yLayout = yField.host_data_layout();
  const stk::mesh::Layout zLayout = zField.host_data_layout();

  if (xLayout == Layout::Right && yLayout == Layout::Right && zLayout == Layout::Right) {
    auto xData = xField.data<T, ConstUnsynchronized, stk::ngp::HostSpace, Layout::Right>();
    auto yData = yField.data<T, ConstUnsynchronized, stk::ngp::HostSpace, Layout::Right>();
    auto zData = zField.data<T, Unsynchronized,      stk::ngp::HostSpace, Layout::Right>();
    FieldBlasImpl<T>::product(buckets, xData, yData, zData);
  }
  else if (xLayout == Layout::Right && yLayout == Layout::Right && zLayout == Layout::Left) {
    auto xData = xField.data<T, ConstUnsynchronized, stk::ngp::HostSpace, Layout::Right>();
    auto yData = yField.data<T, ConstUnsynchronized, stk::ngp::HostSpace, Layout::Right>();
    auto zData = zField.data<T, Unsynchronized,      stk::ngp::HostSpace, Layout::Left>();
    FieldBlasImpl<T>::product(buckets, xData, yData, zData);
  }
  else if (xLayout == Layout::Right && yLayout == Layout::Left && zLayout == Layout::Right) {
    auto xData = xField.data<T, ConstUnsynchronized, stk::ngp::HostSpace, Layout::Right>();
    auto yData = yField.data<T, ConstUnsynchronized, stk::ngp::HostSpace, Layout::Left>();
    auto zData = zField.data<T, Unsynchronized,      stk::ngp::HostSpace, Layout::Right>();
    FieldBlasImpl<T>::product(buckets, xData, yData, zData);
  }
  else if (xLayout == Layout::Right && yLayout == Layout::Left && zLayout == Layout::Left) {
    auto xData = xField.data<T, ConstUnsynchronized, stk::ngp::HostSpace, Layout::Right>();
    auto yData = yField.data<T, ConstUnsynchronized, stk::ngp::HostSpace, Layout::Left>();
    auto zData = zField.data<T, Unsynchronized,      stk::ngp::HostSpace, Layout::Left>();
    FieldBlasImpl<T>::product(buckets, xData, yData, zData);
  }
  else if (xLayout == Layout::Left && yLayout == Layout::Right && zLayout == Layout::Right) {
    auto xData = xField.data<T, ConstUnsynchronized, stk::ngp::HostSpace, Layout::Left>();
    auto yData = yField.data<T, ConstUnsynchronized, stk::ngp::HostSpace, Layout::Right>();
    auto zData = zField.data<T, Unsynchronized,      stk::ngp::HostSpace, Layout::Right>();
    FieldBlasImpl<T>::product(buckets, xData, yData, zData);
  }
  else if (xLayout == Layout::Left && yLayout == Layout::Right && zLayout == Layout::Left) {
    auto xData = xField.data<T, ConstUnsynchronized, stk::ngp::HostSpace, Layout::Left>();
    auto yData = yField.data<T, ConstUnsynchronized, stk::ngp::HostSpace, Layout::Right>();
    auto zData = zField.data<T, Unsynchronized,      stk::ngp::HostSpace, Layout::Left>();
    FieldBlasImpl<T>::product(buckets, xData, yData, zData);
  }
  else if (xLayout == Layout::Left && yLayout == Layout::Left && zLayout == Layout::Right) {
    auto xData = xField.data<T, ConstUnsynchronized, stk::ngp::HostSpace, Layout::Left>();
    auto yData = yField.data<T, ConstUnsynchronized, stk::ngp::HostSpace, Layout::Left>();
    auto zData = zField.data<T, Unsynchronized,      stk::ngp::HostSpace, Layout::Right>();
    FieldBlasImpl<T>::product(buckets, xData, yData, zData);
  }
  else if (xLayout == Layout::Left && yLayout == Layout::Left && zLayout == Layout::Left) {
    auto xData = xField.data<T, ConstUnsynchronized, stk::ngp::HostSpace, Layout::Left>();
    auto yData = yField.data<T, ConstUnsynchronized, stk::ngp::HostSpace, Layout::Left>();
    auto zData = zField.data<T, Unsynchronized,      stk::ngp::HostSpace, Layout::Left>();
    FieldBlasImpl<T>::product(buckets, xData, yData, zData);
  }
  else {
    STK_ThrowErrorMsg("Unsupported Field data layout detected in field_product().  xField layout: " << xLayout <<
                      ", yField layout: " << yLayout << ", zField layout: " << zLayout);
  }
}

template <typename T>
inline
void field_product(const FieldBase& xField, const FieldBase& yField, const FieldBase& zField)
{
  const Selector selector = selectField(xField) & selectField(yField);
  field_product<T>(xField, yField, zField, selector);
}

inline
void field_product(const FieldBase& xField, const FieldBase& yField, const FieldBase& zField, const Selector& selector)
{
  if (xField.data_traits().type_info == typeid(double)) {
    field_product<double>(xField, yField, zField, selector);
  }
  else if (xField.data_traits().type_info == typeid(float)) {
    field_product<float>(xField, yField, zField, selector);
  }
  else if (xField.data_traits().type_info == typeid(std::complex<double>)) {
    field_product<std::complex<double>>(xField, yField, zField, selector);
  }
  else if (xField.data_traits().type_info == typeid(std::complex<float>)) {
    field_product<std::complex<float>>(xField, yField, zField, selector);
  }
  else if (xField.data_traits().type_info == typeid(int)) {
    field_product<int>(xField, yField, zField, selector);
  }
  else {
    STK_ThrowErrorMsg("Error in field_product(): Field is of type " << xField.data_traits().type_info.name() <<
                      ", which is not supported.");
  }
}

inline
void field_product(const FieldBase& xField, const FieldBase& yField, const FieldBase& zField)
{
  const Selector selector = selectField(xField) & selectField(yField);
  field_product(xField, yField, zField, selector);
}


//==============================================================================
// copy: y[i] = x[i]
//
template <typename XDATA, typename YDATA>
class LegacyDeviceFieldCopy {
public:
  LegacyDeviceFieldCopy(const XDATA& xData_, const YDATA& yData_)
    : xData(xData_),
      yData(yData_)
  {}

  KOKKOS_FUNCTION
  void operator()(const stk::mesh::FastMeshIndex& entityIndex) const
  {
    auto xValues = xData.entity_values(entityIndex);
    auto yValues = yData.entity_values(entityIndex);

    for (stk::mesh::ScalarIdx scalar : yValues.scalars()) {
      yValues(scalar) = xValues(scalar);
    }
  }

  XDATA xData;
  YDATA yData;
};

template <typename T>
inline
void field_copy(const FieldBase& xField, const FieldBase& yField, const Selector& selector)
{
  STK_ThrowRequire(is_compatible<T>(xField));
  STK_ThrowRequire(is_compatible(xField, yField));

  const BucketVector& buckets = xField.get_mesh().get_buckets(xField.entity_rank(), selector);

  // We need to remove the modify_on_host()/modify_on_device() behavior from this, as well as
  // the device-side copying.  Folks should use the NGP-flavored version of field_copy() if
  // they want to copy data in a different memory space, rather than relying on an odd
  // side-effect of this function.  None of the other functions in this file have this
  // behavior.
#ifdef STK_USE_DEVICE_MESH
  const bool upToDateOnHost = not xField.need_sync_to_host();
  if (upToDateOnHost) {
#endif
    const stk::mesh::Layout xLayout = xField.host_data_layout();
    const stk::mesh::Layout yLayout = yField.host_data_layout();

    if (xLayout == Layout::Right && yLayout == Layout::Right) {
      auto xData = xField.data<T, ConstUnsynchronized, stk::ngp::HostSpace, Layout::Right>();
      auto yData = yField.data<T, ReadWrite, stk::ngp::HostSpace, Layout::Right>();
      FieldBlasImpl<T>::copy(buckets, xData, yData);
    }
    else if (xLayout == Layout::Left && yLayout == Layout::Left) {
      auto xData = xField.data<T, ConstUnsynchronized, stk::ngp::HostSpace, Layout::Left>();
      auto yData = yField.data<T, ReadWrite, stk::ngp::HostSpace, Layout::Left>();
      FieldBlasImpl<T>::copy(buckets, xData, yData);
    }
    else if (xLayout == Layout::Right && yLayout == Layout::Left) {
      auto xData = xField.data<T, ConstUnsynchronized, stk::ngp::HostSpace, Layout::Right>();
      auto yData = yField.data<T, ReadWrite, stk::ngp::HostSpace, Layout::Left>();
      FieldBlasImpl<T>::copy(buckets, xData, yData);
    }
    else if (xLayout == Layout::Left && yLayout == Layout::Right) {
      auto xData = xField.data<T, ConstUnsynchronized, stk::ngp::HostSpace, Layout::Left>();
      auto yData = yField.data<T, ReadWrite, stk::ngp::HostSpace, Layout::Right>();
      FieldBlasImpl<T>::copy(buckets, xData, yData);
    }

#ifdef STK_USE_DEVICE_MESH
  }
  else {
    if constexpr (std::is_floating_point_v<T> || std::is_same_v<T, int>) {
      auto ngpMesh = stk::mesh::get_updated_ngp_mesh(xField.get_mesh());

      auto xData = xField.data<T, ConstUnsynchronized, stk::ngp::DeviceSpace>();
      auto yData = yField.data<T, ReadWrite, stk::ngp::DeviceSpace>();
      LegacyDeviceFieldCopy fieldCopy(xData, yData);
      stk::mesh::for_each_entity_run(ngpMesh, xField.entity_rank(), selector, fieldCopy);
    }
    else {
      STK_ThrowErrorMsg("Device Fields do not support std::complex datatypes!");
    }
  }
#endif
}

inline
void field_copy(const FieldBase& xField, const FieldBase& yField, const Selector& selector)
{
  if (xField.data_traits().type_info == typeid(double)) {
    field_copy<double>(xField, yField, selector);
  }
  else if (xField.data_traits().type_info == typeid(float)) {
    field_copy<float>(xField, yField, selector);
  }
  else if (xField.data_traits().type_info == typeid(std::complex<double>)) {
    field_copy<std::complex<double>>(xField, yField, selector);
  }
  else if (xField.data_traits().type_info == typeid(std::complex<float>)) {
    field_copy<std::complex<float>>(xField, yField, selector);
  }
  else if (xField.data_traits().type_info == typeid(int)) {
    field_copy<int>(xField, yField, selector);
  }
  else {
    STK_ThrowAssertMsg(false,"Error in field_copy(): Field is of type "<<xField.data_traits().type_info.name() <<
                       ", which is not supported");
  }
}

inline
void field_copy(const FieldBase& xField, const FieldBase& yField)
{
  const Selector selector = selectField(xField) & selectField(yField);
  field_copy(xField, yField, selector);
}


//==============================================================================
// dot: global_sum( sum_i( x[i]*y[i] ) )
//
template <typename T, Layout xLayout, Layout yLayout>
inline
T field_dot(const Field<T, xLayout>& xField, const Field<T, yLayout>& yField, const Selector& selector,
            const MPI_Comm comm)
{
  const BucketVector& buckets = xField.get_mesh().get_buckets(xField.entity_rank(),
                                                              selector & xField.mesh_meta_data().locally_owned_part());

  auto xData = xField.template data<ConstUnsynchronized>();
  auto yData = yField.template data<Unsynchronized>();
  T localResult = FieldBlasImpl<T>::dot(buckets, xData, yData);

  T globalResult = localResult;
  if constexpr (is_complex_v<T>) {
    using ComplexType = typename T::value_type;
    const ComplexType* localResultArray = reinterpret_cast<ComplexType*>(&localResult);
    ComplexType* globalResultArray = reinterpret_cast<ComplexType*>(&globalResult);
    stk::all_reduce_sum(comm, localResultArray, globalResultArray, 2u);
  }
  else {
    stk::all_reduce_sum(comm, &localResult, &globalResult, 1u);
  }

  return globalResult;
}

template <typename T, Layout xLayout, Layout yLayout>
inline
T field_dot(const Field<T, xLayout>& xField, const Field<T, yLayout>& yField, const Selector& selector)
{
  const MPI_Comm comm = xField.get_mesh().parallel();
  return field_dot(xField, yField, selector, comm);
}

template <typename T, Layout xLayout, Layout yLayout>
inline
T field_dot(const Field<T, xLayout>& xField, const Field<T, yLayout>& yField)
{
  const Selector selector = selectField(xField) & selectField(yField);
  return field_dot(xField, yField, selector);
}

template <typename T>
inline
void field_dot(T& globalResult, const FieldBase& xFieldBase, const FieldBase& yFieldBase,
               const Selector& selector, const MPI_Comm comm)
{
  STK_ThrowRequire(is_compatible<T>(xFieldBase));
  STK_ThrowRequire(is_compatible(xFieldBase, yFieldBase));

  if (xFieldBase.host_data_layout() == Layout::Right && yFieldBase.host_data_layout() == Layout::Right) {
    const auto& xField = dynamic_cast<const stk::mesh::Field<T, Layout::Right>&>(xFieldBase);
    const auto& yField = dynamic_cast<const stk::mesh::Field<T, Layout::Right>&>(yFieldBase);
    globalResult = field_dot(xField, yField, selector, comm);
  }
  else if (xFieldBase.host_data_layout() == Layout::Left && yFieldBase.host_data_layout() == Layout::Left) {
    const auto& xField = dynamic_cast<const stk::mesh::Field<T, Layout::Left>&>(xFieldBase);
    const auto& yField = dynamic_cast<const stk::mesh::Field<T, Layout::Left>&>(yFieldBase);
    globalResult = field_dot(xField, yField, selector, comm);
  }
  else if (xFieldBase.host_data_layout() == Layout::Left && yFieldBase.host_data_layout() == Layout::Right) {
    const auto& xField = dynamic_cast<const stk::mesh::Field<T, Layout::Left>&>(xFieldBase);
    const auto& yField = dynamic_cast<const stk::mesh::Field<T, Layout::Right>&>(yFieldBase);
    globalResult = field_dot(xField, yField, selector, comm);
  }
  else if (xFieldBase.host_data_layout() == Layout::Right && yFieldBase.host_data_layout() == Layout::Left) {
    const auto& xField = dynamic_cast<const stk::mesh::Field<T, Layout::Right>&>(xFieldBase);
    const auto& yField = dynamic_cast<const stk::mesh::Field<T, Layout::Left>&>(yFieldBase);
    globalResult = field_dot(xField, yField, selector, comm);
  }
  else {
    STK_ThrowErrorMsg("Unsupported Field data layout detected in field_dot().  xField layout: "
                      << xFieldBase.host_data_layout() << " and yField layout: " << yFieldBase.host_data_layout());
  }
}

template <typename T>
inline
void field_dot(T& result, const FieldBase& xField, const FieldBase& yField, const Selector& selector)
{
  const MPI_Comm comm = xField.get_mesh().parallel();
  field_dot(result, xField, yField, selector, comm);
}

template <typename T>
inline
void field_dot(T& result, const FieldBase& xField, const FieldBase& yField)
{
  const Selector selector = selectField(xField) & selectField(yField);
  field_dot(result, xField, yField, selector);
}


//==============================================================================
// nrm2: sqrt( global_sum( sum_i( x[i]*x[i] )))
//
template <typename T, Layout xLayout>
inline
T field_nrm2(const Field<T, xLayout>& xField, const Selector& selector, const MPI_Comm comm)
{
  const BucketVector& buckets = xField.get_mesh().get_buckets(xField.entity_rank(),
                                                              selector & xField.mesh_meta_data().locally_owned_part());

  auto xData = xField.template data<ConstUnsynchronized>();
  T localResult = FieldBlasImpl<T>::nrm2(buckets, xData);

  T globalResult = localResult;
  if constexpr (is_complex_v<T>) {
    using ComplexType = typename T::value_type;
    const ComplexType* localResultArray = reinterpret_cast<ComplexType*>(&localResult);
    ComplexType* globalResultArray = reinterpret_cast<ComplexType*>(&globalResult);
    stk::all_reduce_sum(comm, localResultArray, globalResultArray, 1u);
    globalResult = {std::sqrt(globalResult.real()), 0.0};  // Imaginary already zero
  }
  else {
    stk::all_reduce_sum(comm, &localResult, &globalResult, 1u);
    globalResult = std::sqrt(globalResult);
  }

  return globalResult;
}

template <typename T, Layout xLayout>
inline
T field_nrm2(const Field<T, xLayout>& xField, const Selector& selector)
{
  const MPI_Comm comm = xField.get_mesh().parallel();
  return field_nrm2(xField, selector, comm);
}

template <typename T, Layout xLayout>
inline
T field_nrm2(const Field<T, xLayout>& xField)
{
  const Selector selector = selectField(xField);
  return field_nrm2(xField, selector);
}

template <typename T>
inline
void field_nrm2(T& globalResult, const FieldBase& xFieldBase, const Selector& selector, const MPI_Comm comm)
{
  STK_ThrowRequire(is_compatible<T>(xFieldBase));

  if (xFieldBase.host_data_layout() == Layout::Right) {
    const auto& xField = dynamic_cast<const stk::mesh::Field<T, Layout::Right>&>(xFieldBase);
    globalResult = field_nrm2(xField, selector, comm);
  }
  else if (xFieldBase.host_data_layout() == Layout::Left) {
    const auto& xField = dynamic_cast<const stk::mesh::Field<T, Layout::Left>&>(xFieldBase);
    globalResult = field_nrm2(xField, selector, comm);
  }
  else {
    STK_ThrowErrorMsg("Unsupported Field data layout detected in field_nrm2().  xField layout: "
                      << xFieldBase.host_data_layout());
  }
}

template <typename T>
inline
void field_nrm2(T& result, const FieldBase& xField, const Selector& selector)
{
  const MPI_Comm comm = xField.get_mesh().parallel();
  field_nrm2(result, xField, selector, comm);
}

template <typename T>
inline
void field_nrm2(T& result, const FieldBase& xField)
{
  const Selector selector = selectField(xField);
  field_nrm2(result, xField, selector);
}


//==============================================================================
// scale: x[i] = a*x[i]
//
template <typename T>
inline
void field_scale(T alpha, const FieldBase& xField, const Selector& selector)
{
  STK_ThrowRequire(is_compatible<T>(xField));

  const BucketVector& buckets = xField.get_mesh().get_buckets(xField.entity_rank(), selector);

  if (xField.host_data_layout() == Layout::Right) {
    auto xData = xField.data<T, Unsynchronized, stk::ngp::HostSpace, Layout::Right>();
    FieldBlasImpl<T>::scale(buckets, alpha, xData);
  }
  else if (xField.host_data_layout() == Layout::Left) {
    auto xData = xField.data<T, Unsynchronized, stk::ngp::HostSpace, Layout::Left>();
    FieldBlasImpl<T>::scale(buckets, alpha, xData);
  }
  else {
    STK_ThrowErrorMsg("Unsupported Field data layout detected in field_scale().  xField layout: "
                      << xField.host_data_layout());
  }
}

template <typename T>
inline
void field_scale(T alpha, const FieldBase& xField)
{
  const Selector selector = selectField(xField);
  field_scale(alpha, xField, selector);
}


//==============================================================================
// fill: x[i] = a
//
template <typename T>
inline
void field_fill(T alpha, const FieldBase& xField, const Selector& selector)
{
  STK_ThrowRequire(is_compatible<T>(xField));

  const BucketVector& buckets = xField.get_mesh().get_buckets(xField.entity_rank(), selector);

  if (xField.host_data_layout() == Layout::Right) {
    auto xData = xField.data<T, Unsynchronized, stk::ngp::HostSpace, Layout::Right>();
    FieldBlasImpl<T>::fill(buckets, alpha, xData);
  }
  else if (xField.host_data_layout() == Layout::Left) {
    auto xData = xField.data<T, Unsynchronized, stk::ngp::HostSpace, Layout::Left>();
    FieldBlasImpl<T>::fill(buckets, alpha, xData);
  }
  else {
    STK_ThrowErrorMsg("Unsupported Field data layout detected in field_fill().  xField layout: "
                      << xField.host_data_layout());
  }
}

template <typename T>
inline
void field_fill(T alpha, const FieldBase& xField)
{
  const Selector selector = selectField(xField);
  field_fill(alpha, xField, selector);
}


template <typename T>
inline
void field_fill(T alpha, const std::vector<const FieldBase*>& xFields, const Selector& selector)
{
  for (const FieldBase* field : xFields) {
    STK_ThrowRequire(is_compatible<T>(*field));
    field_fill(alpha, *field, selector);
  }
}

template <typename T>
inline
void field_fill(T alpha, const std::vector<const FieldBase*>& xFields)
{
  for (const FieldBase* field : xFields) {
    const Selector selector = selectField(*field);
    field_fill(alpha, *field, selector);
  }
}


//==============================================================================
// fill_component: x[i, comp] = a[comp]
//
template <typename T>
inline
void field_fill_component(const T* alpha, const FieldBase& xField, const Selector& selector)
{
  STK_ThrowRequire(is_compatible<T>(xField));

  const BucketVector& buckets = xField.get_mesh().get_buckets(xField.entity_rank(), selector);

  if (xField.host_data_layout() == Layout::Right) {
    auto xData = xField.data<T, Unsynchronized, stk::ngp::HostSpace, Layout::Right>();
    FieldBlasImpl<T>::fill_component(buckets, alpha, xData);
  }
  else if (xField.host_data_layout() == Layout::Left) {
    auto xData = xField.data<T, Unsynchronized, stk::ngp::HostSpace, Layout::Left>();
    FieldBlasImpl<T>::fill_component(buckets, alpha, xData);
  }
  else {
    STK_ThrowErrorMsg("Unsupported Field data layout detected in field_fill_component().  xField layout: "
                      << xField.host_data_layout());
  }
}

template <typename T>
inline
void field_fill_component(const T* alpha, const FieldBase& xField)
{
  const Selector selector = selectField(xField);
  field_fill_component(alpha, xField, selector);
}


//==============================================================================
// swap: x[i] = y[i]
//       y[i] = x[i]
//
template <typename T>
inline
void field_swap(const FieldBase& xField, const FieldBase& yField, const Selector& selector)
{
  STK_ThrowRequire(is_compatible<T>(xField));
  STK_ThrowRequire(is_compatible(xField, yField));

  const BucketVector& buckets = xField.get_mesh().get_buckets(xField.entity_rank(), selector);

  const stk::mesh::Layout xLayout = xField.host_data_layout();
  const stk::mesh::Layout yLayout = yField.host_data_layout();

  if (xLayout == Layout::Right && yLayout == Layout::Right) {
    auto xData = xField.data<T, Unsynchronized, stk::ngp::HostSpace, Layout::Right>();
    auto yData = yField.data<T, Unsynchronized, stk::ngp::HostSpace, Layout::Right>();
    FieldBlasImpl<T>::swap(buckets, xData, yData);
  }
  else if (xLayout == Layout::Left && yLayout == Layout::Left) {
    auto xData = xField.data<T, Unsynchronized, stk::ngp::HostSpace, Layout::Left>();
    auto yData = yField.data<T, Unsynchronized, stk::ngp::HostSpace, Layout::Left>();
    FieldBlasImpl<T>::swap(buckets, xData, yData);
  }
  else if (xLayout == Layout::Right && yLayout == Layout::Left) {
    auto xData = xField.data<T, Unsynchronized, stk::ngp::HostSpace, Layout::Right>();
    auto yData = yField.data<T, Unsynchronized, stk::ngp::HostSpace, Layout::Left>();
    FieldBlasImpl<T>::swap(buckets, xData, yData);
  }
  else if (xLayout == Layout::Left && yLayout == Layout::Right) {
    auto xData = xField.data<T, Unsynchronized, stk::ngp::HostSpace, Layout::Left>();
    auto yData = yField.data<T, Unsynchronized, stk::ngp::HostSpace, Layout::Right>();
    FieldBlasImpl<T>::swap(buckets, xData, yData);
  }
}

inline
void field_swap(const FieldBase& xField, const FieldBase& yField, const Selector& selector)
{
  if (xField.data_traits().type_info == typeid(double)) {
    field_swap<double>(xField, yField, selector);
  }
  else if (xField.data_traits().type_info == typeid(float)) {
    field_swap<float>(xField, yField, selector);
  }
  else if (xField.data_traits().type_info == typeid(std::complex<double>)) {
    field_swap<std::complex<double>>(xField, yField, selector);
  }
  else if (xField.data_traits().type_info == typeid(std::complex<float>)) {
    field_swap<std::complex<float>>(xField, yField, selector);
  }
  else if (xField.data_traits().type_info == typeid(int)) {
    field_swap<int>(xField, yField, selector);
  }
  else {
    STK_ThrowErrorMsg("Error in field_swap(): Fields are of type " << xField.data_traits().type_info.name() <<
                       ", which is not supported.");
  }
}

inline
void field_swap(const FieldBase& xField, const FieldBase& yField)
{
  const Selector selector = selectField(xField) & selectField(yField);
  field_swap(xField, yField, selector);
}


//==============================================================================
// asum: global_sum( sum_i( abs(x[i]) ))
//
template <typename T, Layout xLayout>
inline
T field_asum(const Field<T, xLayout>& xField, const Selector& selector, const MPI_Comm comm)
{
  const BucketVector& buckets = xField.get_mesh().get_buckets(xField.entity_rank(),
                                                              selector & xField.mesh_meta_data().locally_owned_part());

  auto xData = xField.template data<ConstUnsynchronized>();
  T localResult = FieldBlasImpl<T>::asum(buckets, xData);

  T globalResult = localResult;
  if constexpr (is_complex_v<T>) {
    using ComplexType = typename T::value_type;
    const ComplexType* localResultArray = reinterpret_cast<ComplexType*>(&localResult);
    ComplexType* globalResultArray = reinterpret_cast<ComplexType*>(&globalResult);
    stk::all_reduce_sum(comm, localResultArray, globalResultArray, 2u);
  }
  else {
    stk::all_reduce_sum(comm, &localResult, &globalResult, 1u);
  }

  return globalResult;
}

template <typename T, Layout xLayout>
inline
T field_asum(const Field<T, xLayout>& xField, const Selector& selector)
{
  const MPI_Comm comm = xField.get_mesh().parallel();
  return field_asum(xField, selector, comm);
}

template <typename T, Layout xLayout>
inline
T field_asum(const Field<T, xLayout>& xField)
{
  const Selector selector = selectField(xField);
  return field_asum(xField,selector);
}


template <typename T>
inline
void field_asum(T& globalResult, const FieldBase& xFieldBase, const Selector& selector, const MPI_Comm comm)
{
  STK_ThrowRequire(is_compatible<T>(xFieldBase));

  if (xFieldBase.host_data_layout() == Layout::Right) {
    const auto& xField = dynamic_cast<const stk::mesh::Field<T, Layout::Right>&>(xFieldBase);
    globalResult = field_asum(xField, selector, comm);
  }
  else if (xFieldBase.host_data_layout() == Layout::Left) {
    const auto& xField = dynamic_cast<const stk::mesh::Field<T, Layout::Left>&>(xFieldBase);
    globalResult = field_asum(xField, selector, comm);
  }
  else {
    STK_ThrowErrorMsg("Unsupported Field data layout detected in field_asum().  xField layout: "
                      << xFieldBase.host_data_layout());
  }
}

template <typename T>
inline
void field_asum(T& result, const FieldBase& xFieldBase, const Selector& selector)
{
  const MPI_Comm comm = xFieldBase.get_mesh().parallel();
  field_asum(result, xFieldBase, selector, comm);
}

template <typename T>
inline
void field_asum(T& result, const FieldBase& xFieldBase)
{
  const Selector selector = selectField(xFieldBase);
  field_asum(result, xFieldBase, selector);
}


//==============================================================================
// amax: global_max( max_i( abs(x[i]) ))
//
template <typename T, Layout xLayout>
inline
T field_amax(const Field<T, xLayout>& xField, const Selector& selector, const MPI_Comm comm)
{
  const BucketVector& buckets = xField.get_mesh().get_buckets(xField.entity_rank(),
                                                              selector & xField.mesh_meta_data().locally_owned_part());

  auto xData = xField.template data<ConstUnsynchronized>();
  T localMaxValue = FieldBlasImpl<T>::amax(buckets, xData);

  T globalResult = localMaxValue;
  if constexpr (is_complex_v<T>) {
    using ComplexType = typename T::value_type;
    const ComplexType* localResultArray = reinterpret_cast<ComplexType*>(&localMaxValue);
    ComplexType* globalResultArray = reinterpret_cast<ComplexType*>(&globalResult);
    stk::all_reduce_max(comm, localResultArray, globalResultArray, 1u);  // Only the real part has a value
  }
  else {
    stk::all_reduce_max(comm, &localMaxValue, &globalResult, 1u);
  }

  return globalResult;
}

template <typename T, Layout xLayout>
inline
T field_amax(const Field<T, xLayout>& xField, const Selector& selector)
{
  const MPI_Comm comm = xField.get_mesh().parallel();
  return field_amax(xField, selector, comm);
}

template <typename T, Layout xLayout>
inline
T field_amax(const Field<T, xLayout>& xField)
{
  const Selector selector = selectField(xField);
  return field_amax(xField, selector);
}


template <typename T>
inline
void field_amax(T& result, const FieldBase& xFieldBase, const Selector& selector, const MPI_Comm comm)
{
  STK_ThrowRequire(is_compatible<T>(xFieldBase));

  if (xFieldBase.host_data_layout() == Layout::Right) {
    const auto& xField = dynamic_cast<const stk::mesh::Field<T, Layout::Right>&>(xFieldBase);
    result = field_amax(xField, selector, comm);
  }
  else if (xFieldBase.host_data_layout() == Layout::Left) {
    const auto& xField = dynamic_cast<const stk::mesh::Field<T, Layout::Left>&>(xFieldBase);
    result = field_amax(xField, selector, comm);
  }
  else {
    STK_ThrowErrorMsg("Unsupported Field data layout detected in field_amax().  xField layout: "
                      << xFieldBase.host_data_layout());
  }
}

template <typename T>
inline
void field_amax(T& result, const FieldBase& xFieldBase, const Selector& selector)
{
  const MPI_Comm comm = xFieldBase.get_mesh().parallel();
  field_amax(result, xFieldBase, selector, comm);
}

template <typename T>
inline
void field_amax(T& result, const FieldBase& xFieldBase)
{
  const Selector selector = selectField(xFieldBase);
  field_amax(result, xFieldBase, selector);
}


//==============================================================================
// eamax: Entity(loc:global_max( max_i( abs(x[i]) )))
//
template <typename T>
inline
Entity field_eamax(const FieldBase& xField, const Selector& selector)
{
  STK_ThrowRequire(is_compatible<T>(xField));

  const BucketVector& buckets = xField.get_mesh().get_buckets(xField.entity_rank(),
                                                              selector & xField.mesh_meta_data().locally_owned_part());

  const stk::mesh::Layout xLayout = xField.host_data_layout();

  MinMaxInfo localMinMaxInfo {};
  T localMaxValue {};

  if (xLayout == Layout::Right) {
    auto xData = xField.data<T, ConstUnsynchronized, stk::ngp::HostSpace, Layout::Right>();
    localMinMaxInfo = FieldBlasImpl<T>::iamax(buckets, xData, localMaxValue);
  }
  else if (xLayout == Layout::Left) {
    auto xData = xField.data<T, ConstUnsynchronized, stk::ngp::HostSpace, Layout::Left>();
    localMinMaxInfo = FieldBlasImpl<T>::iamax(buckets, xData, localMaxValue);
  }
  else {
    STK_ThrowErrorMsg("Unsupported Field data layout detected in field_eamax().  xField layout: " << xLayout);
  }

  const stk::mesh::BulkData& bulk = xField.get_mesh();
  T globalMaxValue;
  EntityId globalEntityId {};
  EntityId localEntityId {};
  STK_ThrowRequireMsg(localMinMaxInfo.bucketId != InvalidOrdinal, "No minimum value location found in field_eamax()");

  if constexpr (is_complex_v<T>) {
    using ComplexType = typename T::value_type;
    const ComplexType* localValueArray = reinterpret_cast<ComplexType*>(&localMaxValue);
    ComplexType* globalValueArray = reinterpret_cast<ComplexType*>(&globalMaxValue);
    const Bucket& localBucket = *bulk.buckets(xField.entity_rank())[localMinMaxInfo.bucketId];
    localEntityId = bulk.identifier(localBucket[localMinMaxInfo.entityIndex]);
    stk::all_reduce_maxloc(bulk.parallel(), localValueArray, &localEntityId, globalValueArray, &globalEntityId, 1u);
  }
  else {
    const Bucket& localBucket = *bulk.buckets(xField.entity_rank())[localMinMaxInfo.bucketId];
    localEntityId = bulk.identifier(localBucket[localMinMaxInfo.entityIndex]);
    stk::all_reduce_maxloc(bulk.parallel(), &localMaxValue, &localEntityId, &globalMaxValue, &globalEntityId, 1u);
  }

  if (globalEntityId == localEntityId) {
    return bulk.get_entity(xField.entity_rank(), globalEntityId);
  } else {
    return Entity();
  }
}

inline
Entity field_eamax(const FieldBase& xField, const Selector& selector)
{
  if (xField.data_traits().type_info == typeid(double)) {
    return field_eamax<double>(xField, selector);
  }
  else if (xField.data_traits().type_info == typeid(float)) {
    return field_eamax<float>(xField, selector);
  }
  else if (xField.data_traits().type_info == typeid(std::complex<double>)) {
    return field_eamax<std::complex<double>>(xField, selector);
  }
  else if (xField.data_traits().type_info == typeid(std::complex<float>)) {
    return field_eamax<std::complex<float>>(xField, selector);
  }
  else if (xField.data_traits().type_info == typeid(int)) {
    return field_eamax<int>(xField, selector);
  }
  else {
    STK_ThrowErrorMsg("Error in field_eamax(): Field is of type " << xField.data_traits().type_info.name() <<
                      ", which is not supported.");
  }
  return stk::mesh::Entity();
}

inline
Entity field_eamax(const FieldBase& xField)
{
  const Selector selector = selectField(xField);
  return field_eamax(xField, selector);
}

template <typename T>
inline
Entity field_eamax(const Field<T>& xField)
{
  const Selector selector = selectField(xField);
  return field_eamax(xField, selector);
}


//==============================================================================
// amin: global_min( min_i( abs(x[i]) ))
//
template <typename T, Layout xLayout>
inline
T field_amin(const Field<T, xLayout>& xField, const Selector& selector, const MPI_Comm comm)
{
  const BucketVector& buckets = xField.get_mesh().get_buckets(xField.entity_rank(),
                                                              selector & xField.mesh_meta_data().locally_owned_part());

  auto xData = xField.template data<ConstUnsynchronized>();
  T localMinValue = FieldBlasImpl<T>::amin(buckets, xData);

  T globalResult = localMinValue;
  if constexpr (is_complex_v<T>) {
    using ComplexType = typename T::value_type;
    const ComplexType* localResultArray = reinterpret_cast<ComplexType*>(&localMinValue);
    ComplexType* globalResultArray = reinterpret_cast<ComplexType*>(&globalResult);
    stk::all_reduce_min(comm, localResultArray, globalResultArray, 1u);  // Only the real part has a value
  }
  else {
    stk::all_reduce_min(comm, &localMinValue, &globalResult, 1u);
  }

  return globalResult;
}

template <typename T, Layout xLayout>
inline
T field_amin(const Field<T, xLayout>& xField, const Selector& selector)
{
  const MPI_Comm comm = xField.get_mesh().parallel();
  return field_amin(xField, selector, comm);
}

template <typename T, Layout xLayout>
inline
T field_amin(const Field<T, xLayout>& xField)
{
  const Selector selector = selectField(xField);
  return field_amin(xField, selector);
}


template <typename T>
inline
void field_amin(T& result, const FieldBase& xFieldBase, const Selector& selector, const MPI_Comm comm)
{
  STK_ThrowRequire(is_compatible<T>(xFieldBase));

  if (xFieldBase.host_data_layout() == Layout::Right) {
    const auto& xField = dynamic_cast<const stk::mesh::Field<T, Layout::Right>&>(xFieldBase);
    result = field_amin(xField, selector, comm);
  }
  else if (xFieldBase.host_data_layout() == Layout::Left) {
    const auto& xField = dynamic_cast<const stk::mesh::Field<T, Layout::Left>&>(xFieldBase);
    result = field_amin(xField, selector, comm);
  }
  else {
    STK_ThrowErrorMsg("Unsupported Field data layout detected in field_amin().  xField layout: "
                      << xFieldBase.host_data_layout());
  }
}

template <typename T>
inline
void field_amin(T& result, const FieldBase& xFieldBase, const Selector& selector)
{
  const MPI_Comm comm = xFieldBase.get_mesh().parallel();
  field_amin(result, xFieldBase, selector, comm);
}

template <typename T>
inline
void field_amin(T& result, const FieldBase& xFieldBase)
{
  const Selector selector = selectField(xFieldBase);
  field_amin(result, xFieldBase, selector);
}


//==============================================================================
// eamin: Entity(loc:global_min( min_i( abs(x[i]) )))
//
template <typename T>
inline
Entity field_eamin(const FieldBase& xField, const Selector& selector)
{
  STK_ThrowRequire(is_compatible<T>(xField));

  const BucketVector& buckets = xField.get_mesh().get_buckets(xField.entity_rank(),
                                                              selector & xField.mesh_meta_data().locally_owned_part());

  const stk::mesh::Layout xLayout = xField.host_data_layout();

  MinMaxInfo localMinMaxInfo {};
  T localMinValue {};

  if (xLayout == Layout::Right) {
    auto xData = xField.data<T, ConstUnsynchronized, stk::ngp::HostSpace, Layout::Right>();
    localMinMaxInfo = FieldBlasImpl<T>::iamin(buckets, xData, localMinValue);
  }
  else if (xLayout == Layout::Left) {
    auto xData = xField.data<T, ConstUnsynchronized, stk::ngp::HostSpace, Layout::Left>();
    localMinMaxInfo = FieldBlasImpl<T>::iamin(buckets, xData, localMinValue);
  }
  else {
    STK_ThrowErrorMsg("Unsupported Field data layout detected in field_eamin().  xField layout: " << xLayout);
  }

  const stk::mesh::BulkData& bulk = xField.get_mesh();
  T globalMinValue;
  EntityId globalEntityId {};
  EntityId localEntityId {};
  STK_ThrowRequireMsg(localMinMaxInfo.bucketId != InvalidOrdinal, "No minimum value location found in field_eamin()");

  if constexpr (is_complex_v<T>) {
    using ComplexType = typename T::value_type;
    const ComplexType* localValueArray = reinterpret_cast<ComplexType*>(&localMinValue);
    ComplexType* globalValueArray = reinterpret_cast<ComplexType*>(&globalMinValue);
    const Bucket& localBucket = *bulk.buckets(xField.entity_rank())[localMinMaxInfo.bucketId];
    localEntityId = bulk.identifier(localBucket[localMinMaxInfo.entityIndex]);
    stk::all_reduce_minloc(bulk.parallel(), localValueArray, &localEntityId, globalValueArray, &globalEntityId, 1u);
  }
  else {
    const Bucket& localBucket = *bulk.buckets(xField.entity_rank())[localMinMaxInfo.bucketId];
    localEntityId = bulk.identifier(localBucket[localMinMaxInfo.entityIndex]);
    stk::all_reduce_minloc(bulk.parallel(), &localMinValue, &localEntityId, &globalMinValue, &globalEntityId, 1u);
  }

  if (globalEntityId == localEntityId) {
    return bulk.get_entity(xField.entity_rank(), globalEntityId);
  } else {
    return Entity();
  }
}

inline
Entity field_eamin(const FieldBase& xField, const Selector& selector)
{
  if (xField.data_traits().type_info == typeid(double)) {
    return field_eamin<double>(xField, selector);
  }
  else if (xField.data_traits().type_info == typeid(float)) {
    return field_eamin<float>(xField, selector);
  }
  else if (xField.data_traits().type_info == typeid(std::complex<double>)) {
    return field_eamin<std::complex<double>>(xField, selector);
  }
  else if (xField.data_traits().type_info == typeid(std::complex<float>)) {
    return field_eamin<std::complex<float>>(xField, selector);
  }
  else if (xField.data_traits().type_info == typeid(int)) {
    return field_eamin<int>(xField, selector);
  }
  else {
    STK_ThrowErrorMsg("Error in field_eamin(): Field is of type " << xField.data_traits().type_info.name() <<
                      ", which is not supported.");
  }
  return stk::mesh::Entity();
}

inline
Entity field_eamin(const FieldBase& xField)
{
  const Selector selector = selectField(xField);
  return field_eamin(xField, selector);
}

template <typename T>
inline
Entity field_eamin(const Field<T>& xField)
{
  const Selector selector = selectField(xField);
  return field_eamin(xField, selector);
}


} // mesh
} // stk

#endif // STK_MESH_BASE_FIELDBLAS_HPP

