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

#ifndef STK_MESH_BASEIMPL_FIELDBLASIMPL_HPP
#define STK_MESH_BASEIMPL_FIELDBLASIMPL_HPP

#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_util/util/Blas.hpp>
#include <complex>
#include <algorithm>
#include <type_traits>

namespace stk::field_blas::impl {

template<typename T> struct is_complex_t : public std::false_type {};
template<typename U> struct is_complex_t<std::complex<U>> : public std::true_type {};
template<typename U> struct is_complex_t<Kokkos::complex<U>> : public std::true_type {};
template<typename T> constexpr bool is_complex_v = is_complex_t<T>::value;

template <typename T1>
constexpr bool is_layout_right() {
  return (T1::layout == mesh::Layout::Right);
}

template <typename T1, typename T2>
constexpr bool is_layout_right() {
  return (T1::layout == mesh::Layout::Right) &&
         (T2::layout == mesh::Layout::Right);
}

template <typename T1, typename T2, typename T3>
constexpr bool is_layout_right() {
  return (T1::layout == mesh::Layout::Right) &&
         (T2::layout == mesh::Layout::Right) &&
         (T3::layout == mesh::Layout::Right);
}

template <typename Scalar>
concept is_float_or_double = requires
{
  requires std::is_same_v<Scalar, float> || std::is_same_v<Scalar, double>;
};

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
void check_matching_extents([[maybe_unused]] const std::string& functionName,
    [[maybe_unused]] const XVALS& xValues,
    [[maybe_unused]] const YVALS& yValues)
{
  STK_ThrowAssertMsg(xValues.num_scalars() == yValues.num_scalars(),
                     "Cannot perform " << functionName << " operation on different-length Fields.  fieldX has " <<
                     xValues.num_scalars() << " scalars and fieldY has " << yValues.num_scalars() << " scalars.");
}

template <typename XVALS, typename YVALS, typename ZVALS>
void check_matching_extents([[maybe_unused]] const std::string& functionName,
    [[maybe_unused]] const XVALS& xValues,
    [[maybe_unused]] const YVALS& yValues,
    [[maybe_unused]] const ZVALS& zValues)
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
  // Note: Can't include the maxValue information in the MinMaxInfo struct return type due to
  // un-suppressable gcc warnings about ABI changes for std::complex<float> in structs.  So,
  // return it separately as an out-parameter.

  template <typename XDATA>
  inline static MinMaxInfo iamax(const mesh::BucketVector& buckets, const XDATA& xData, T& maxValue)
  {
    if constexpr (is_complex_v<T>) {
      using ComplexType = typename T::value_type;
      ComplexType globalMaxValueReal {};
      MinMaxInfo globalMinMaxInfo {mesh::InvalidOrdinal, 0, 0};

#ifdef STK_USE_OPENMP
      #pragma omp parallel
#endif
      {
        ComplexType localMaxValueReal {};
        MinMaxInfo localMinMaxInfo {mesh::InvalidOrdinal, 0, 0};

#ifdef STK_USE_OPENMP
        #pragma omp for schedule(static)
#endif
        for (unsigned bucketIndex = 0; bucketIndex < buckets.size(); ++bucketIndex) {
          const mesh::Bucket& bucket = *buckets[bucketIndex];
          auto xValues = xData.bucket_values(bucket);

          for (stk::mesh::EntityIdx entityIdx : bucket.entities()) {
            for (stk::mesh::ScalarIdx scalar : xValues.scalars()) {
              Kokkos::complex<ComplexType> kc = xValues(entityIdx, scalar);
              ComplexType localValueReal = Kokkos::abs(kc);
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
      MinMaxInfo globalMinMaxInfo {mesh::InvalidOrdinal, 0, 0};

#ifdef STK_USE_OPENMP
      #pragma omp parallel
#endif
      {
        T localMaxValue {};
        MinMaxInfo localMinMaxInfo {mesh::InvalidOrdinal, 0, 0};

#ifdef STK_USE_OPENMP
        #pragma omp for schedule(static)
#endif
        for (unsigned bucketIndex = 0; bucketIndex < buckets.size(); ++bucketIndex) {
          const mesh::Bucket& bucket = *buckets[bucketIndex];
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
              mesh::EntityIdx localEntityIdx(localIndex / xValues.num_scalars());
              mesh::ScalarIdx localScalarIdx(localIndex % xValues.num_scalars());
              T localValue = Kokkos::abs(xValues(localEntityIdx, localScalarIdx));

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
                  T localValue = Kokkos::abs(xValues(entityIdx, scalar));

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
  inline static MinMaxInfo iamin(const mesh::BucketVector& buckets, const XDATA& xData, T& minValue)
  {
    if constexpr (is_complex_v<T>) {
      using ComplexType = typename T::value_type;
      ComplexType globalMinValueReal = std::numeric_limits<ComplexType>::max();
      MinMaxInfo globalMinMaxInfo {mesh::InvalidOrdinal, 0, 0};

#ifdef STK_USE_OPENMP
      #pragma omp parallel
#endif
      {
        ComplexType localMinValueReal = std::numeric_limits<ComplexType>::max();
        MinMaxInfo localMinMaxInfo {mesh::InvalidOrdinal, 0, 0};

#ifdef STK_USE_OPENMP
        #pragma omp for schedule(static)
#endif
        for (unsigned bucketIndex = 0; bucketIndex < buckets.size(); ++bucketIndex) {
          const mesh::Bucket& bucket = *buckets[bucketIndex];
          auto xValues = xData.bucket_values(bucket);

          for (stk::mesh::EntityIdx entityIdx : bucket.entities()) {
            for (stk::mesh::ScalarIdx scalar : xValues.scalars()) {
              Kokkos::complex<ComplexType> kc = xValues(entityIdx, scalar);
              ComplexType localValueReal = Kokkos::abs(kc);

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
      MinMaxInfo globalMinMaxInfo {mesh::InvalidOrdinal, 0, 0};

#ifdef STK_USE_OPENMP
      #pragma omp parallel
#endif
      {
        T localMinValue = std::numeric_limits<T>::max();
        MinMaxInfo localMinMaxInfo {mesh::InvalidOrdinal, 0, 0};

#ifdef STK_USE_OPENMP
        #pragma omp for schedule(static)
#endif
        for (unsigned bucketIndex = 0; bucketIndex < buckets.size(); ++bucketIndex) {
          const mesh::Bucket& bucket = *buckets[bucketIndex];
          auto xValues = xData.bucket_values(bucket);

          for (stk::mesh::EntityIdx entityIdx : bucket.entities()) {
            for (stk::mesh::ScalarIdx scalar : xValues.scalars()) {
              T localValue = Kokkos::abs(xValues(entityIdx, scalar));

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

template <typename Scalar, typename XDATA, typename YDATA> requires (is_float_or_double<Scalar> &&
                                                                     is_layout_right<XDATA, YDATA>())
void field_axpy_on_host(const mesh::BucketVector& buckets, Scalar alpha, const XDATA& xData, const YDATA& yData)
{
#ifdef STK_USE_OPENMP
  #pragma omp parallel for schedule(static)
#endif
  for (unsigned bucketIndex = 0; bucketIndex < buckets.size(); ++bucketIndex) {
    const mesh::Bucket& bucket = *buckets[bucketIndex];
    auto xValues = xData.bucket_values(bucket);
    auto yValues = yData.bucket_values(bucket);
    check_matching_extents("field_axpy()", xValues, yValues);

    const int numValues = yValues.num_entities()*yValues.num_scalars();
    const int stride = xValues.scalar_stride();

    if constexpr (std::is_same_v<Scalar, float>) {
      SIERRA_FORTRAN(saxpy)(&numValues, &alpha, xValues.pointer(), &stride, yValues.pointer(), &stride);
    }
    else {
      SIERRA_FORTRAN(daxpy)(&numValues, &alpha, xValues.pointer(), &stride, yValues.pointer(), &stride);
    }
  }
}

template <typename Scalar, typename XDATA, typename YDATA> requires (is_float_or_double<Scalar> &&
                                                                     !is_layout_right<XDATA, YDATA>())
void field_axpy_on_host(const mesh::BucketVector& buckets, Scalar alpha, const XDATA& xData, const YDATA& yData)
{
#ifdef STK_USE_OPENMP
  #pragma omp parallel for schedule(static)
#endif
  for (unsigned bucketIndex = 0; bucketIndex < buckets.size(); ++bucketIndex) {
    const mesh::Bucket& bucket = *buckets[bucketIndex];
    auto xValues = xData.bucket_values(bucket);
    auto yValues = yData.bucket_values(bucket);
    check_matching_extents("field_axpy()", xValues, yValues);

    for (mesh::ScalarIdx scalar : yValues.scalars()) {
      const int numValues = yValues.num_entities();
      const int xStride = xValues.entity_stride();
      const int yStride = yValues.entity_stride();
      const Scalar* xPtr = xValues.pointer() + scalar*xValues.scalar_stride();
      Scalar* yPtr = yValues.pointer() + scalar*yValues.scalar_stride();

      if constexpr (std::is_same_v<Scalar, float>) {
        SIERRA_FORTRAN(saxpy)(&numValues, &alpha, xPtr, &xStride, yPtr, &yStride);
      }
      else {
        SIERRA_FORTRAN(daxpy)(&numValues, &alpha, xPtr, &xStride, yPtr, &yStride);
      }
    }
  }
}

template <typename Scalar, typename XDATA, typename YDATA>
void field_axpy_on_host(const mesh::BucketVector& buckets, Scalar alpha, const XDATA& xData, const YDATA& yData)
{
#ifdef STK_USE_OPENMP
  #pragma omp parallel for schedule(static)
#endif
  for (unsigned bucketIndex = 0; bucketIndex < buckets.size(); ++bucketIndex) {
    const mesh::Bucket& bucket = *buckets[bucketIndex];
    auto xValues = xData.bucket_values(bucket);
    auto yValues = yData.bucket_values(bucket);
    check_matching_extents("field_axpy()", xValues, yValues);

    for (mesh::EntityIdx entityIdx : bucket.entities()) {
      for (mesh::ScalarIdx scalarIdx : yValues.scalars()) {
        yValues(entityIdx, scalarIdx) += alpha * xValues(entityIdx, scalarIdx);
      }
    }
  }
}

template <typename T>
void field_axpy_impl(T alpha, const mesh::FieldBase& xField, const mesh::FieldBase& yField, const mesh::Selector& selector)
{
  const auto& buckets = xField.get_mesh().get_buckets(xField.entity_rank(), selector);

  field_data_execute<T, T, mesh::ReadOnly, mesh::ReadWrite>(xField,
                                                            yField,
                                                            [&buckets, alpha](auto& xData, auto& yData) {
                                                            field_axpy_on_host(buckets, alpha, xData, yData);
                                                            });
}

template <typename Scalar, typename XDATA, typename YDATA> requires (is_layout_right<XDATA, YDATA>())
void field_axpby_on_host(const mesh::BucketVector& buckets, Scalar alpha, const XDATA& xData, Scalar beta, const YDATA& yData)
{
#ifdef STK_USE_OPENMP
  #pragma omp parallel for schedule(static)
#endif
  for (unsigned bucketIndex = 0; bucketIndex < buckets.size(); ++bucketIndex) {
    const mesh::Bucket& bucket = *buckets[bucketIndex];
    auto xValues = xData.bucket_values(bucket);
    auto yValues = yData.bucket_values(bucket);
    check_matching_extents("field_axpby()", xValues, yValues);

    const int numValues = xValues.num_entities() * xValues.num_scalars();
    const Scalar* xPtr = xValues.pointer();
    Scalar* yPtr = yValues.pointer();

    for (int i = 0; i < numValues; ++i) {
      yPtr[i] = beta * yPtr[i] + alpha * xPtr[i];
    }
  }
}

template <typename Scalar, typename XDATA, typename YDATA>
void field_axpby_on_host(const mesh::BucketVector& buckets, Scalar alpha, const XDATA& xData, Scalar beta, const YDATA& yData)
{
#ifdef STK_USE_OPENMP
    #pragma omp parallel for schedule(static)
#endif
  for (unsigned bucketIndex = 0; bucketIndex < buckets.size(); ++bucketIndex) {
    const mesh::Bucket& bucket = *buckets[bucketIndex];
    auto xValues = xData.bucket_values(bucket);
    auto yValues = yData.bucket_values(bucket);
    check_matching_extents("field_axpby()", xValues, yValues);

    for (mesh::EntityIdx entityIdx : bucket.entities()) {
      for (mesh::ScalarIdx scalarIdx : yValues.scalars()) {
        yValues(entityIdx, scalarIdx) = beta * yValues(entityIdx, scalarIdx) + alpha * xValues(entityIdx, scalarIdx);
      }
    }
  }
}

template <typename T>
void field_axpby_impl(T alpha, const mesh::FieldBase& xField, T beta, const mesh::FieldBase& yField, const mesh::Selector& selector)
{
  if ( beta == T(1) ) {
    field_axpy_impl(alpha, xField, yField, selector);
  }
  else {
    const auto& buckets = xField.get_mesh().get_buckets(xField.entity_rank(), selector);
    field_data_execute<T, T, mesh::ReadOnly, mesh::ReadWrite>(xField,
                                                  yField,
                                                  [&buckets, alpha, beta](auto& xData, auto& yData) {
                                                    field_axpby_on_host(buckets, alpha, xData, beta, yData);
                                                  });
  }
}

template <typename Scalar, typename XDATA, typename YDATA, typename ZDATA> requires (is_layout_right<XDATA, YDATA, ZDATA>())
void field_axpbyz_on_host(const mesh::BucketVector& buckets, Scalar alpha, const XDATA& xData, Scalar beta, const YDATA& yData, const ZDATA& zData)
{
#ifdef STK_USE_OPENMP
  #pragma omp parallel for schedule(static)
#endif
  for (unsigned bucketIndex = 0; bucketIndex < buckets.size(); ++bucketIndex) {
    const mesh::Bucket& bucket = *buckets[bucketIndex];
    auto xValues = xData.bucket_values(bucket);
    auto yValues = yData.bucket_values(bucket);
    auto zValues = zData.bucket_values(bucket);
    check_matching_extents("field_axpbyz()", xValues, yValues);

    const int numValues = xValues.num_entities() * xValues.num_scalars();
    const Scalar* xPtr = xValues.pointer();
    const Scalar* yPtr = yValues.pointer();
    Scalar* zPtr = zValues.pointer();

    for (int i = 0; i < numValues; ++i) {
      zPtr[i] = beta * yPtr[i] + alpha * xPtr[i];
    }
  }
}

template <typename Scalar, typename XDATA, typename YDATA, typename ZDATA>
void field_axpbyz_on_host(const mesh::BucketVector& buckets, Scalar alpha, const XDATA& xData, Scalar beta, const YDATA& yData, const ZDATA& zData)
{
#ifdef STK_USE_OPENMP
  #pragma omp parallel for schedule(static)
#endif
  for (unsigned bucketIndex = 0; bucketIndex < buckets.size(); ++bucketIndex) {
    const mesh::Bucket& bucket = *buckets[bucketIndex];
    auto xValues = xData.bucket_values(bucket);
    auto yValues = yData.bucket_values(bucket);
    auto zValues = zData.bucket_values(bucket);
    check_matching_extents("field_axpbyz()", xValues, yValues);

    for (mesh::EntityIdx entityIdx : bucket.entities()) {
      for (mesh::ScalarIdx scalarIdx : yValues.scalars()) {
        zValues(entityIdx, scalarIdx) = beta * yValues(entityIdx, scalarIdx) + alpha * xValues(entityIdx, scalarIdx);
      }
    }
  }
}

template <typename T>
void field_axpbyz_impl(T alpha, const mesh::FieldBase& xField, T beta, const mesh::FieldBase& yField, const mesh::FieldBase& zField, const mesh::Selector& selector)
{
  const auto& buckets = xField.get_mesh().get_buckets(xField.entity_rank(), selector);

  field_data_execute<T, T, T, mesh::ReadOnly, mesh::ReadOnly, mesh::OverwriteAll>(xField,
                                                                                  yField,
                                                                                  zField,
                                                                                  [&buckets, alpha, beta](auto& xData, auto& yData, auto& zData) {
                                                                                    field_axpbyz_on_host(buckets, alpha, xData, beta, yData, zData);
                                                                                  });
}

template <typename XDATA, typename YDATA, typename ZDATA> requires (is_layout_right<XDATA, YDATA, ZDATA>())
void field_product_on_host(const mesh::BucketVector& buckets, const XDATA& xData, const YDATA& yData, const ZDATA& zData)
{
#ifdef STK_USE_OPENMP
  #pragma omp parallel for schedule(static)
#endif
  for (unsigned bucketIndex = 0; bucketIndex < buckets.size(); ++bucketIndex) {
    const mesh::Bucket& bucket = *buckets[bucketIndex];
    auto xValues = xData.bucket_values(bucket);
    auto yValues = yData.bucket_values(bucket);
    auto zValues = zData.bucket_values(bucket);
    check_matching_extents("field_product()", xValues, yValues, zValues);

    const int numValues = xValues.num_entities() * xValues.num_scalars();
    const auto* xPtr = xValues.pointer();
    const auto* yPtr = yValues.pointer();
    auto* zPtr = zValues.pointer();

    for (int i = 0; i < numValues; ++i) {
      zPtr[i] = xPtr[i] * yPtr[i];
    }
  }
}

template <typename XDATA, typename YDATA, typename ZDATA>
void field_product_on_host(const mesh::BucketVector& buckets, const XDATA& xData, const YDATA& yData, const ZDATA& zData)
{
#ifdef STK_USE_OPENMP
  #pragma omp parallel for schedule(static)
#endif
  for (unsigned bucketIndex = 0; bucketIndex < buckets.size(); ++bucketIndex) {
    const mesh::Bucket& bucket = *buckets[bucketIndex];
    auto xValues = xData.bucket_values(bucket);
    auto yValues = yData.bucket_values(bucket);
    auto zValues = zData.bucket_values(bucket);
    check_matching_extents("field_product()", xValues, yValues, zValues);

    for (mesh::EntityIdx entityIdx : bucket.entities()) {
      for (mesh::ScalarIdx scalarIdx : yValues.scalars()) {
        zValues(entityIdx, scalarIdx) = xValues(entityIdx, scalarIdx) * yValues(entityIdx, scalarIdx);
      }
    }
  }
}

template <typename T>
void apply_product_on_field(const mesh::FieldBase& xField,
                            const mesh::FieldBase& yField,
                            const mesh::FieldBase& zField,
                            const mesh::Selector& selector)
{
  const auto& buckets = xField.get_mesh().get_buckets(xField.entity_rank(), selector);

  field_data_execute<T, T, T, mesh::ReadOnly, mesh::ReadOnly, mesh::OverwriteAll>(xField,
                                                                                  yField,
                                                                                  zField,
                                                                                  [&buckets](auto& xData, auto& yData, auto& zData) {
                                                                                    field_product_on_host(buckets, xData, yData, zData);
                                                                                  });
}

inline
void field_product_impl(const mesh::FieldBase& xField, const mesh::FieldBase& yField, const mesh::FieldBase& zField,
                        const mesh::Selector& selector)
{
  field_datatype_execute(xField,
    [&]<typename T>(const stk::mesh::FieldBase&) {
      apply_product_on_field<T>(xField, yField, zField, selector);
    }
  );
}

template <typename XDATA, typename YDATA> requires (is_layout_right<XDATA, YDATA>())
void field_copy_on_host(const mesh::BucketVector& buckets, const XDATA& xData, const YDATA& yData)
{
#ifdef STK_USE_OPENMP
  #pragma omp parallel for schedule(static)
#endif
  for (unsigned bucketIndex = 0; bucketIndex < buckets.size(); ++bucketIndex) {
    const mesh::Bucket& bucket = *buckets[bucketIndex];
    auto xValues = xData.bucket_values(bucket);
    auto yValues = yData.bucket_values(bucket);
    check_matching_extents("field_copy()", xValues, yValues);

    const int numValues = xValues.num_entities() * xValues.num_scalars();
    const auto* xPtr = xValues.pointer();
    auto* yPtr = yValues.pointer();

    for (int i = 0; i < numValues; ++i) {
      yPtr[i] = xPtr[i];
    }
  }
}

template <typename XDATA, typename YDATA>
void field_copy_on_host(const mesh::BucketVector& buckets, const XDATA& xData, const YDATA& yData)
{
#ifdef STK_USE_OPENMP
  #pragma omp parallel for schedule(static)
#endif
  for (unsigned bucketIndex = 0; bucketIndex < buckets.size(); ++bucketIndex) {
    const mesh::Bucket& bucket = *buckets[bucketIndex];
    auto xValues = xData.bucket_values(bucket);
    auto yValues = yData.bucket_values(bucket);
    check_matching_extents("field_copy()", xValues, yValues);

    for (mesh::EntityIdx entityIdx : bucket.entities()) {
      for (mesh::ScalarIdx scalarIdx : yValues.scalars()) {
        yValues(entityIdx, scalarIdx) = xValues(entityIdx, scalarIdx);
      }
    }
  }
}

template <typename T>
void apply_copy_on_field(const mesh::FieldBase& xField, const mesh::FieldBase& yField, const mesh::Selector& selector)
{
  const auto& buckets = xField.get_mesh().get_buckets(xField.entity_rank(), selector);

  field_data_execute<T, T, mesh::ReadOnly, mesh::OverwriteAll>(xField,
                                                               yField,
                                                               [&buckets](auto xData, auto yData) {
                                                                 field_copy_on_host(buckets, xData, yData);
                                                               });
}

inline
void field_copy_impl(const mesh::FieldBase& xField, const mesh::FieldBase& yField, const mesh::Selector& selector)
{
  field_datatype_execute(xField,
    [&]<typename T>(const stk::mesh::FieldBase&) {
      apply_copy_on_field<T>(xField, yField, selector);
    }
  );
}

template <typename Scalar, typename XDATA, typename YDATA> requires is_complex_v<Scalar>
void field_dot_on_host(Scalar& localResult, const mesh::BucketVector& buckets, const XDATA& xData, const YDATA& yData)
{
  using ComplexType = typename Scalar::value_type;
  ComplexType localResultReal {};  // OpenMP can't do reductions on std::complex types
  ComplexType localResultImag {};

#ifdef STK_USE_OPENMP
  #pragma omp parallel for schedule(static) reduction(+:localResultReal,localResultImag)
#endif
  for (unsigned bucketIndex = 0; bucketIndex < buckets.size(); ++bucketIndex) {
    const mesh::Bucket& bucket = *buckets[bucketIndex];
    auto xValues = xData.bucket_values(bucket);
    auto yValues = yData.bucket_values(bucket);
    check_matching_extents("field_dot()", xValues, yValues);

    for (mesh::EntityIdx entityIdx : bucket.entities()) {
      for (mesh::ScalarIdx scalarIdx : yValues.scalars()) {
        Scalar value = xValues(entityIdx, scalarIdx) * yValues(entityIdx, scalarIdx);
        localResultReal += value.real();
        localResultImag += value.imag();
      }
    }
  }

  localResult = Scalar(localResultReal, localResultImag);
}

template <typename Scalar, typename XDATA, typename YDATA> requires (is_float_or_double<Scalar> && is_layout_right<XDATA, YDATA>())
void field_dot_on_host(Scalar& localResult, const mesh::BucketVector& buckets, const XDATA& xData, const YDATA& yData)
{
  localResult = Scalar{0};
#ifdef STK_USE_OPENMP
  #pragma omp parallel for reduction(+:localResult) schedule(static)
#endif
  for (unsigned bucketIndex = 0; bucketIndex < buckets.size(); ++bucketIndex) {
    const mesh::Bucket& bucket = *buckets[bucketIndex];
    auto xValues = xData.bucket_values(bucket);
    auto yValues = yData.bucket_values(bucket);
    check_matching_extents("field_dot()", xValues, yValues);

    const int numValues = yValues.num_entities()*yValues.num_scalars();
    const int stride = xValues.scalar_stride();

    if constexpr (std::is_same_v<Scalar, float>) {
      localResult += SIERRA_FORTRAN(sdot)(&numValues, xValues.pointer(), &stride, yValues.pointer(), &stride);
    }
    else {
      localResult += SIERRA_FORTRAN(ddot)(&numValues, xValues.pointer(), &stride, yValues.pointer(), &stride);
    }
  }
}

template <typename Scalar, typename XDATA, typename YDATA> requires (is_float_or_double<Scalar> && !is_layout_right<XDATA, YDATA>())
void field_dot_on_host(Scalar& localResult, const mesh::BucketVector& buckets, const XDATA& xData, const YDATA& yData)
{
  localResult = Scalar{0};
#ifdef STK_USE_OPENMP
  #pragma omp parallel for reduction(+:localResult) schedule(static)
#endif
  for (unsigned bucketIndex = 0; bucketIndex < buckets.size(); ++bucketIndex) {
    const mesh::Bucket& bucket = *buckets[bucketIndex];
    auto xValues = xData.bucket_values(bucket);
    auto yValues = yData.bucket_values(bucket);
    check_matching_extents("field_dot()", xValues, yValues);

    // Stride isn't uniform through the whole Bucket for one or the other Field, so we have to chop
    // it up into uniform segments
    for (mesh::ScalarIdx scalar : xValues.scalars()) {
      const int numValues = xValues.num_entities();
      const int xStride = xValues.entity_stride();
      const int yStride = yValues.entity_stride();
      const Scalar* xPtr = xValues.pointer() + scalar*xValues.scalar_stride();
      const Scalar* yPtr = yValues.pointer() + scalar*yValues.scalar_stride();

      if constexpr (std::is_same_v<Scalar, float>) {
        localResult += SIERRA_FORTRAN(sdot)(&numValues, xPtr, &xStride, yPtr, &yStride);
      }
      else {
        localResult += SIERRA_FORTRAN(ddot)(&numValues, xPtr, &xStride, yPtr, &yStride);
      }
    }
  }
}

template <typename Scalar, typename XDATA, typename YDATA>
void field_dot_on_host(Scalar& localResult, const mesh::BucketVector& buckets, const XDATA& xData, const YDATA& yData)
{
  localResult = Scalar{0};
#ifdef STK_USE_OPENMP
  #pragma omp parallel for reduction(+:localResult) schedule(static)
#endif
  for (unsigned bucketIndex = 0; bucketIndex < buckets.size(); ++bucketIndex) {
    const mesh::Bucket& bucket = *buckets[bucketIndex];
    auto xValues = xData.bucket_values(bucket);
    auto yValues = yData.bucket_values(bucket);
    check_matching_extents("field_dot()", xValues, yValues);

    for (mesh::EntityIdx entityIdx : bucket.entities()) {
      for (mesh::ScalarIdx scalarIdx : yValues.scalars()) {
        localResult += xValues(entityIdx, scalarIdx) * yValues(entityIdx, scalarIdx);
      }
    }
  }
}

template <typename T>
void field_dot_impl(T& result, const mesh::FieldBase& xField, const mesh::FieldBase& yField, const mesh::Selector& selector)
{
  const auto& buckets = xField.get_mesh().get_buckets(xField.entity_rank(),
                                                              selector & xField.mesh_meta_data().locally_owned_part());
  const MPI_Comm comm = xField.get_mesh().parallel();

  field_data_execute<T, T, mesh::ReadOnly, mesh::ReadOnly>(xField, yField, [&](auto& xData, auto& yData) {
    T localResult;
    field_dot_on_host<T>(localResult, buckets, xData, yData);

    result = localResult;
    if constexpr (is_complex_v<T>) {
      using ComplexType = typename T::value_type;
      const ComplexType* localResultArray = reinterpret_cast<ComplexType*>(&localResult);
      ComplexType* globalResultArray = reinterpret_cast<ComplexType*>(&result);
      stk::all_reduce_sum(comm, localResultArray, globalResultArray, 2u);
    }
    else {
      stk::all_reduce_sum(comm, &localResult, &result, 1u);
    }
  });
}

template <typename Scalar, typename XDATA> requires is_complex_v<Scalar>
void field_nrm2_on_host(Scalar& localResult, const mesh::BucketVector& buckets, const XDATA& xData)
{
  using ComplexType = typename Scalar::value_type;
  ComplexType localResultReal {};

#ifdef STK_USE_OPENMP
  #pragma omp parallel for schedule(static) reduction(+:localResultReal)
#endif
  for (unsigned bucketIndex = 0; bucketIndex < buckets.size(); ++bucketIndex) {
    const mesh::Bucket& bucket = *buckets[bucketIndex];
    auto xValues = xData.bucket_values(bucket);

    for (mesh::EntityIdx entityIdx : bucket.entities()) {
      for (mesh::ScalarIdx scalarIdx : xValues.scalars()) {
        Scalar value = std::pow(Kokkos::abs(xValues(entityIdx, scalarIdx)), 2);
        localResultReal += value.real();
      }
    }
  }

  localResult = Scalar(localResultReal, 0);
}

template <typename Scalar, typename XDATA> requires (is_float_or_double<Scalar> &&
                                                     is_layout_right<XDATA>())
void field_nrm2_on_host(Scalar& localResult, const mesh::BucketVector& buckets, const XDATA& xData)
{
  localResult = Scalar{0};

#ifdef STK_USE_OPENMP
  #pragma omp parallel for reduction(+:localResult) schedule(static)
#endif
  for (unsigned bucketIndex = 0; bucketIndex < buckets.size(); ++bucketIndex) {
    const mesh::Bucket& bucket = *buckets[bucketIndex];
    auto xValues = xData.bucket_values(bucket);

    const int numValues = xValues.num_entities()*xValues.num_scalars();
    const int stride = xValues.scalar_stride();

    if constexpr (std::is_same_v<Scalar, float>) {
      localResult += SIERRA_FORTRAN(sdot)(&numValues, xValues.pointer(), &stride, xValues.pointer(), &stride);
    }
    else {
      localResult += SIERRA_FORTRAN(ddot)(&numValues, xValues.pointer(), &stride, xValues.pointer(), &stride);
    }
  }
}

template <typename Scalar, typename XDATA> requires (is_float_or_double<Scalar> &&
                                                     !is_layout_right<XDATA>())
void field_nrm2_on_host(Scalar& localResult, const mesh::BucketVector& buckets, const XDATA& xData)
{
  localResult = Scalar{0};

#ifdef STK_USE_OPENMP
  #pragma omp parallel for reduction(+:localResult) schedule(static)
#endif
  for (unsigned bucketIndex = 0; bucketIndex < buckets.size(); ++bucketIndex) {
    const mesh::Bucket& bucket = *buckets[bucketIndex];
    auto xValues = xData.bucket_values(bucket);

    // Stride isn't uniform through the whole Bucket, so we have to chop it up into uniform segments
    for (mesh::ScalarIdx scalar : xValues.scalars()) {
      const int numValues = xValues.num_entities();
      const int xStride = xValues.entity_stride();
      const Scalar* xPtr = xValues.pointer() + scalar*xValues.scalar_stride();

      if constexpr (std::is_same_v<Scalar, float>) {
        localResult += SIERRA_FORTRAN(sdot)(&numValues, xPtr, &xStride, xPtr, &xStride);
      }
      else {
        localResult += SIERRA_FORTRAN(ddot)(&numValues, xPtr, &xStride, xPtr, &xStride);
      }
    }
  }
}

template <typename Scalar, typename XDATA>
void field_nrm2_on_host(Scalar& localResult, const mesh::BucketVector& buckets, const XDATA& xData)
{
  localResult = Scalar{0};

#ifdef STK_USE_OPENMP
  #pragma omp parallel for reduction(+:localResult) schedule(static)
#endif
  for (unsigned bucketIndex = 0; bucketIndex < buckets.size(); ++bucketIndex) {
    const mesh::Bucket& bucket = *buckets[bucketIndex];
    auto xValues = xData.bucket_values(bucket);

    for (mesh::EntityIdx entityIdx : bucket.entities()) {
      for (mesh::ScalarIdx scalarIdx : xValues.scalars()) {
        localResult += xValues(entityIdx, scalarIdx) * xValues(entityIdx, scalarIdx);
      }
    }
  }
}

template <typename T>
void field_nrm2_impl(T& result, const mesh::FieldBase& xField, const mesh::Selector& selector)
{
  const auto& buckets = xField.get_mesh().get_buckets(xField.entity_rank(),
                                                              selector & xField.mesh_meta_data().locally_owned_part());
  const MPI_Comm comm = xField.get_mesh().parallel();

  field_data_execute<T, mesh::ReadOnly>(xField, [&](auto& xData) {
    T localResult;
    field_nrm2_on_host(localResult, buckets, xData);

    result = localResult;
    if constexpr (is_complex_v<T>) {
      using ComplexType = typename T::value_type;
      const ComplexType* localResultArray = reinterpret_cast<ComplexType*>(&localResult);
      ComplexType* globalResultArray = reinterpret_cast<ComplexType*>(&result);
      stk::all_reduce_sum(comm, localResultArray, globalResultArray, 1u);
      result = T{Kokkos::sqrt(result.real()), 0.0};  // Imaginary already zero
    }
    else {
      stk::all_reduce_sum(comm, &localResult, &result, 1u);
      result = Kokkos::sqrt(result);
    }
  });
}

template <typename Scalar, typename XDATA> requires (is_float_or_double<Scalar> &&
                                                     is_layout_right<XDATA>())
void field_scale_on_host(const mesh::BucketVector& buckets, Scalar alpha, const XDATA& xData)
{
#ifdef STK_USE_OPENMP
  #pragma omp parallel for schedule(static)
#endif
  for (unsigned bucketIndex = 0; bucketIndex < buckets.size(); ++bucketIndex) {
    const mesh::Bucket& bucket = *buckets[bucketIndex];
    auto xValues = xData.bucket_values(bucket);

    const int numValues = xValues.num_entities()*xValues.num_scalars();
    const int stride = xValues.scalar_stride();

    if constexpr (std::is_same_v<Scalar, float>) {
      SIERRA_FORTRAN(sscal)(&numValues, &alpha, xValues.pointer(), &stride);
    }
    else {
      SIERRA_FORTRAN(dscal)(&numValues, &alpha, xValues.pointer(), &stride);
    }
  }
}

template <typename Scalar, typename XDATA> requires (is_float_or_double<Scalar> &&
                                                     !is_layout_right<XDATA>())
void field_scale_on_host(const mesh::BucketVector& buckets, Scalar alpha, const XDATA& xData)
{
#ifdef STK_USE_OPENMP
  #pragma omp parallel for schedule(static)
#endif
  for (unsigned bucketIndex = 0; bucketIndex < buckets.size(); ++bucketIndex) {
    const mesh::Bucket& bucket = *buckets[bucketIndex];
    auto xValues = xData.bucket_values(bucket);

    // Stride isn't uniform through the whole Bucket, so we have to chop it up into uniform segments
    for (mesh::ScalarIdx scalar : xValues.scalars()) {
      const int numValues = xValues.num_entities();
      const int xStride = xValues.entity_stride();
      Scalar* xPtr = xValues.pointer() + scalar*xValues.scalar_stride();

      if constexpr (std::is_same_v<Scalar, float>) {
        SIERRA_FORTRAN(sscal)(&numValues, &alpha, xPtr, &xStride);
      }
      else {
        SIERRA_FORTRAN(dscal)(&numValues, &alpha, xPtr, &xStride);
      }
    }
  }
}

template <typename Scalar, typename XDATA>
void field_scale_on_host(const mesh::BucketVector& buckets, Scalar alpha, const XDATA& xData)
{
#ifdef STK_USE_OPENMP
  #pragma omp parallel for schedule(static)
#endif
  for (unsigned bucketIndex = 0; bucketIndex < buckets.size(); ++bucketIndex) {
    const mesh::Bucket& bucket = *buckets[bucketIndex];
    auto xValues = xData.bucket_values(bucket);

    for (mesh::EntityIdx entityIdx : bucket.entities()) {
      for (mesh::ScalarIdx scalarIdx : xValues.scalars()) {
        xValues(entityIdx, scalarIdx) = alpha * xValues(entityIdx, scalarIdx);
      }
    }
  }
}

template <typename T>
void field_scale_impl(T alpha, const mesh::FieldBase& xField, const mesh::Selector& selector)
{
  const auto& buckets = xField.get_mesh().get_buckets(xField.entity_rank(), selector);

  field_data_execute<T, mesh::ReadWrite>(xField,
                                   [&buckets, alpha](auto& xData) {
                                     field_scale_on_host(buckets, alpha, xData);
                                   });
}

template <typename Scalar, typename XDATA> requires (is_layout_right<XDATA>())
void field_fill_on_host(const mesh::BucketVector& buckets, Scalar alpha, const XDATA& xData)
{
#ifdef STK_USE_OPENMP
  #pragma omp parallel for schedule(static)
#endif
  for (unsigned bucketIndex = 0; bucketIndex < buckets.size(); ++bucketIndex) {
    const mesh::Bucket& bucket = *buckets[bucketIndex];
    auto xValues = xData.bucket_values(bucket);

    const int numValues = xValues.num_entities() * xValues.num_scalars();
    Scalar* xPtr = xValues.pointer();
    std::fill(xPtr, xPtr+numValues, alpha);
  }
}

template <typename Scalar, typename XDATA>
void field_fill_on_host(const mesh::BucketVector& buckets, Scalar alpha, const XDATA& xData)
{
#ifdef STK_USE_OPENMP
  #pragma omp parallel for schedule(static)
#endif
  for (unsigned bucketIndex = 0; bucketIndex < buckets.size(); ++bucketIndex) {
    const mesh::Bucket& bucket = *buckets[bucketIndex];
    auto xValues = xData.bucket_values(bucket);

    for (mesh::EntityIdx entityIdx : bucket.entities()) {
      for (mesh::ScalarIdx scalarIdx : xValues.scalars()) {
        xValues(entityIdx, scalarIdx) = alpha;
      }
    }
  }
}

template <typename T>
void field_fill_impl(T alpha, const mesh::FieldBase& xField, const mesh::Selector& selector)
{
  const auto& buckets = xField.get_mesh().get_buckets(xField.entity_rank(), selector);

  field_data_execute<T, mesh::OverwriteAll>(xField,
                                            [&](auto& xData) {
                                              field_fill_on_host(buckets, alpha, xData);
                                            });
}

template <typename Scalar, typename XDATA>
void field_fill_on_host(const mesh::BucketVector& buckets, Scalar alpha, const XDATA& xData, int component)
{
#ifdef STK_USE_OPENMP
  #pragma omp parallel for schedule(static)
#endif
  for (unsigned bucketIndex = 0; bucketIndex < buckets.size(); ++bucketIndex) {
    const mesh::Bucket& bucket = *buckets[bucketIndex];
    auto xValues = xData.bucket_values(bucket);

    for (mesh::EntityIdx entityIdx : bucket.entities()) {
      xValues(entityIdx, mesh::ScalarIdx(component)) = alpha;
    }
  }
}

template <typename T>
void field_fill_impl(T alpha, const mesh::FieldBase& xField, int component, const mesh::Selector& selector)
{
  const auto& buckets = xField.get_mesh().get_buckets(xField.entity_rank(), selector);

  field_data_execute<T, mesh::ReadWrite>(xField,
                                         [&](auto& xData) {
                                           field_fill_on_host(buckets, alpha, xData, component);
                                         });
}

template <typename Scalar, typename XDATA>
void field_fill_component_on_host(const mesh::BucketVector& buckets, const Scalar* alpha, const XDATA& xData)
{
#ifdef STK_USE_OPENMP
  #pragma omp parallel for schedule(static)
#endif
  for (unsigned bucketIndex = 0; bucketIndex < buckets.size(); ++bucketIndex) {
    const mesh::Bucket& bucket = *buckets[bucketIndex];
    auto xValues = xData.bucket_values(bucket);

    for (mesh::EntityIdx entityIdx : bucket.entities()) {
      for (mesh::ScalarIdx scalarIdx : xValues.scalars()) {
        xValues(entityIdx, scalarIdx) = alpha[scalarIdx];
      }
    }
  }
}

template <typename T>
void field_fill_component_impl(const T* alpha, const mesh::FieldBase& xField, const mesh::Selector& selector)
{
  const auto& buckets = xField.get_mesh().get_buckets(xField.entity_rank(), selector);

  field_data_execute<T, mesh::OverwriteAll>(xField,
                                            [&](auto& xData) {
                                              field_fill_component_on_host(buckets, alpha, xData);
                                            });
}

template <typename Scalar, typename XDATA, typename YDATA> requires (is_float_or_double<Scalar> &&
                                                                     is_layout_right<XDATA, YDATA>())
void field_swap_on_host(const mesh::BucketVector& buckets, const XDATA& xData, const YDATA& yData)
{
#ifdef STK_USE_OPENMP
  #pragma omp parallel for schedule(static)
#endif
  for (unsigned bucketIndex = 0; bucketIndex < buckets.size(); ++bucketIndex) {
    const mesh::Bucket& bucket = *buckets[bucketIndex];
    auto xValues = xData.bucket_values(bucket);
    auto yValues = yData.bucket_values(bucket);
    check_matching_extents("field_swap()", xValues, yValues);

    const int numValues = xValues.num_entities()*xValues.num_scalars();
    const int stride = xValues.scalar_stride();

    if constexpr (std::is_same_v<Scalar, float>) {
      SIERRA_FORTRAN(sswap)(&numValues, xValues.pointer(), &stride, yValues.pointer(), &stride);
    }
    else {
      SIERRA_FORTRAN(dswap)(&numValues, xValues.pointer(), &stride, yValues.pointer(), &stride);
    }
  }
}

template <typename Scalar, typename XDATA, typename YDATA> requires (is_float_or_double<Scalar> &&
                                                                     !is_layout_right<XDATA, YDATA>())
void field_swap_on_host(const mesh::BucketVector& buckets, const XDATA& xData, const YDATA& yData)
{
#ifdef STK_USE_OPENMP
  #pragma omp parallel for schedule(static)
#endif
  for (unsigned bucketIndex = 0; bucketIndex < buckets.size(); ++bucketIndex) {
    const mesh::Bucket& bucket = *buckets[bucketIndex];
    auto xValues = xData.bucket_values(bucket);
    auto yValues = yData.bucket_values(bucket);
    check_matching_extents("field_swap()", xValues, yValues);

    // Stride isn't uniform through the whole Bucket for one or the other Field, so we have to chop
    // it up into uniform segments
    for (mesh::ScalarIdx scalar : yValues.scalars()) {
      const int numValues = yValues.num_entities();
      const int xStride = xValues.entity_stride();
      const int yStride = yValues.entity_stride();
      Scalar* xPtr = xValues.pointer() + scalar*xValues.scalar_stride();
      Scalar* yPtr = yValues.pointer() + scalar*yValues.scalar_stride();

      if constexpr (std::is_same_v<Scalar, float>) {
        SIERRA_FORTRAN(sswap)(&numValues, xPtr, &xStride, yPtr, &yStride);
      }
      else {
        SIERRA_FORTRAN(dswap)(&numValues, xPtr, &xStride, yPtr, &yStride);
      }
    }
  }
}

template <typename Scalar, typename XDATA, typename YDATA>
void field_swap_on_host(const mesh::BucketVector& buckets, const XDATA& xData, const YDATA& yData)
{
#ifdef STK_USE_OPENMP
  #pragma omp parallel for schedule(static)
#endif
  for (unsigned bucketIndex = 0; bucketIndex < buckets.size(); ++bucketIndex) {
    const mesh::Bucket& bucket = *buckets[bucketIndex];
    auto xValues = xData.bucket_values(bucket);
    auto yValues = yData.bucket_values(bucket);
    check_matching_extents("field_swap()", xValues, yValues);

    for (mesh::EntityIdx entityIdx : bucket.entities()) {
      for (mesh::ScalarIdx scalarIdx : yValues.scalars()) {
        Scalar temp = yValues(entityIdx, scalarIdx);
        yValues(entityIdx, scalarIdx) = xValues(entityIdx, scalarIdx);
        xValues(entityIdx, scalarIdx) = temp;
      }
    }
  }
}

template <typename T>
void apply_swap_on_field(const mesh::FieldBase& xField, const mesh::FieldBase& yField, const mesh::Selector& selector)
{
  const auto& buckets = xField.get_mesh().get_buckets(xField.entity_rank(), selector);

  field_data_execute<T, T, mesh::ReadWrite, mesh::ReadWrite>(xField,
                                                 yField,
                                                 [&buckets](auto& xData, auto& yData) {
                                                   field_swap_on_host<T>(buckets, xData, yData);
                                                 });
}

inline
void field_swap_impl(const mesh::FieldBase& xField, const mesh::FieldBase& yField, const mesh::Selector& selector)
{
  field_datatype_execute(xField,
    [&]<typename T>(const stk::mesh::FieldBase&) {
      apply_swap_on_field<T>(xField, yField, selector);
    }
  );
}

template <typename Scalar, typename XDATA> requires is_complex_v<Scalar>
void field_asum_on_host(Scalar& localResult, const mesh::BucketVector& buckets, const XDATA& xData)
{
  using ComplexType = typename Scalar::value_type;
  ComplexType localResultReal {};

#ifdef STK_USE_OPENMP
  #pragma omp parallel for schedule(static) reduction(+:localResultReal)
#endif
  for (unsigned bucketIndex = 0; bucketIndex < buckets.size(); ++bucketIndex) {
    const mesh::Bucket& bucket = *buckets[bucketIndex];
    auto xValues = xData.bucket_values(bucket);

    for (mesh::EntityIdx entityIdx : bucket.entities()) {
      for (mesh::ScalarIdx scalarIdx : xValues.scalars()) {
        Scalar value = Kokkos::abs(xValues(entityIdx, scalarIdx));
        localResultReal += value.real();
      }
    }
  }

  localResult = Scalar(localResultReal, 0);
}

template <typename Scalar, typename XDATA> requires (is_float_or_double<Scalar> &&
                                                     is_layout_right<XDATA>())
void field_asum_on_host(Scalar& localResult, const mesh::BucketVector& buckets, const XDATA& xData)
{
  localResult = Scalar{};

#ifdef STK_USE_OPENMP
  #pragma omp parallel for schedule(static) reduction(+:localResult)
#endif
  for (unsigned bucketIndex = 0; bucketIndex < buckets.size(); ++bucketIndex) {
    const mesh::Bucket& bucket = *buckets[bucketIndex];
    auto xValues = xData.bucket_values(bucket);

    const int numValues = xValues.num_entities()*xValues.num_scalars();
    const int stride = xValues.scalar_stride();

    if constexpr (std::is_same_v<Scalar, double>) {
      localResult += SIERRA_FORTRAN(dasum)(&numValues, xValues.pointer(), &stride);
    }
    else if constexpr (std::is_same_v<Scalar, float>) {
      localResult += SIERRA_FORTRAN(sasum)(&numValues, xValues.pointer(), &stride);
    }
  }
}

template <typename Scalar, typename XDATA> requires (is_float_or_double<Scalar> &&
                                                     !is_layout_right<XDATA>())
void field_asum_on_host(Scalar& localResult, const mesh::BucketVector& buckets, const XDATA& xData)
{
  localResult = Scalar{};

#ifdef STK_USE_OPENMP
  #pragma omp parallel for schedule(static) reduction(+:localResult)
#endif
  for (unsigned bucketIndex = 0; bucketIndex < buckets.size(); ++bucketIndex) {
    const mesh::Bucket& bucket = *buckets[bucketIndex];
    auto xValues = xData.bucket_values(bucket);

    // Stride isn't uniform through the whole Bucket, so we have to chop it up into uniform segments
    for (mesh::ScalarIdx scalar : xValues.scalars()) {
      const int numValues = xValues.num_entities();
      const int xStride = xValues.entity_stride();
      const Scalar* xPtr = xValues.pointer() + scalar*xValues.scalar_stride();

      if constexpr (std::is_same_v<Scalar, float>) {
        localResult += SIERRA_FORTRAN(sasum)(&numValues, xPtr, &xStride);
      }
      else {
        localResult += SIERRA_FORTRAN(dasum)(&numValues, xPtr, &xStride);
      }
    }
  }
}

template <typename Scalar, typename XDATA>
void field_asum_on_host(Scalar& localResult, const mesh::BucketVector& buckets, const XDATA& xData)
{
  localResult = Scalar{};

#ifdef STK_USE_OPENMP
  #pragma omp parallel for schedule(static) reduction(+:localResult)
#endif
  for (unsigned bucketIndex = 0; bucketIndex < buckets.size(); ++bucketIndex) {
    const mesh::Bucket& bucket = *buckets[bucketIndex];
    auto xValues = xData.bucket_values(bucket);

    for (mesh::EntityIdx entityIdx : bucket.entities()) {
      for (mesh::ScalarIdx scalarIdx : xValues.scalars()) {
        localResult += Kokkos::abs(xValues(entityIdx, scalarIdx));
      }
    }
  }
}

template <typename T>
void field_asum_impl(T& result, const mesh::FieldBase& xField, const mesh::Selector& selector)
{
  const auto& buckets = xField.get_mesh().get_buckets(xField.entity_rank(),
                                                      selector & xField.mesh_meta_data().locally_owned_part());
  const MPI_Comm comm = xField.get_mesh().parallel();

  field_data_execute<T, mesh::ReadOnly>(xField, [&](auto& xData) {
    T localResult;
    field_asum_on_host(localResult, buckets, xData);

    result = localResult;
    if constexpr (is_complex_v<T>) {
      using ComplexType = typename T::value_type;
      const ComplexType* localResultArray = reinterpret_cast<ComplexType*>(&localResult);
      ComplexType* globalResultArray = reinterpret_cast<ComplexType*>(&result);
      stk::all_reduce_sum(comm, localResultArray, globalResultArray, 2u);
    }
    else {
      stk::all_reduce_sum(comm, &localResult, &result, 1u);
    }
  });
}

template <typename Scalar, typename XDATA> requires is_complex_v<Scalar>
void field_amax_on_host(Scalar& localMaxValue, const mesh::BucketVector& buckets, const XDATA& xData)
{
  using ComplexType = typename Scalar::value_type;
  ComplexType localMaxValueReal = {};

#ifdef STK_USE_OPENMP
#pragma omp parallel for schedule(static) reduction(max:localMaxValueReal)
#endif
  for (unsigned bucketIndex = 0; bucketIndex < buckets.size(); ++bucketIndex) {
    const mesh::Bucket& bucket = *buckets[bucketIndex];
    auto xValues = xData.bucket_values(bucket);

    for (mesh::EntityIdx entityIdx : bucket.entities()) {
      for (mesh::ScalarIdx scalarIdx : xValues.scalars()) {
        ComplexType localValueReal = Kokkos::abs(xValues(entityIdx, scalarIdx));
        if (localValueReal > localMaxValueReal) {
          localMaxValueReal = localValueReal;
        }
      }
    }
  }
  localMaxValue = Scalar(localMaxValueReal, 0);
}

template <typename Scalar, typename XDATA> requires (is_float_or_double<Scalar> &&
                                                     is_layout_right<XDATA>())
void field_amax_on_host(Scalar& localMaxValue, const mesh::BucketVector& buckets, const XDATA& xData)
{
  localMaxValue = Scalar{};

#ifdef STK_USE_OPENMP
#pragma omp parallel for schedule(static) reduction(max:localMaxValue)
#endif
  for (unsigned bucketIndex = 0; bucketIndex < buckets.size(); ++bucketIndex) {
    const mesh::Bucket& bucket = *buckets[bucketIndex];
    auto xValues = xData.bucket_values(bucket);
    if (xValues.is_field_defined()) {
      const int numValues = xValues.num_entities()*xValues.num_scalars();
      const int stride = xValues.scalar_stride();
      const auto* xPtr = xValues.pointer();
      int localIndex {};
      if constexpr (std::is_same_v<Scalar, float>) {
        localIndex = SIERRA_FORTRAN(isamax)(&numValues, xPtr, &stride) - 1;
      }
      else {
        localIndex = SIERRA_FORTRAN(idamax)(&numValues, xPtr, &stride) - 1;
      }
      Scalar localValue = Kokkos::abs(xPtr[localIndex]);

      if (localValue > localMaxValue) {
        localMaxValue = localValue;
      }
    }
  }
}

template <typename Scalar, typename XDATA>
void field_amax_on_host(Scalar& localMaxValue, const mesh::BucketVector& buckets, const XDATA& xData)
{
  localMaxValue = Scalar{};

#ifdef STK_USE_OPENMP
#pragma omp parallel for schedule(static) reduction(max:localMaxValue)
#endif
  for (unsigned bucketIndex = 0; bucketIndex < buckets.size(); ++bucketIndex) {
    const mesh::Bucket& bucket = *buckets[bucketIndex];
    auto xValues = xData.bucket_values(bucket);
    if (xValues.is_field_defined()) {
      for (mesh::EntityIdx entityIdx : bucket.entities()) {
        for (mesh::ScalarIdx scalarIdx : xValues.scalars()) {
          Scalar localValue = Kokkos::abs(xValues(entityIdx, scalarIdx));

          if (localValue > localMaxValue) {
            localMaxValue = localValue;
          }
        }
      }
    }
  }
}

template <typename T>
void field_amax_impl(T& result, const mesh::FieldBase& xField, const mesh::Selector& selector)
{
  const auto& buckets = xField.get_mesh().get_buckets(xField.entity_rank(),
                                                      selector & xField.mesh_meta_data().locally_owned_part());
  const MPI_Comm comm = xField.get_mesh().parallel();

  field_data_execute<T, mesh::ReadOnly>(xField, [&](auto& xData) {
    T localMaxValue;
    field_amax_on_host(localMaxValue, buckets, xData);

    result = localMaxValue;
    if constexpr (is_complex_v<T>) {
      using ComplexType = typename T::value_type;
      const ComplexType* localResultArray = reinterpret_cast<ComplexType*>(&localMaxValue);
      ComplexType* globalResultArray = reinterpret_cast<ComplexType*>(&result);
      stk::all_reduce_max(comm, localResultArray, globalResultArray, 1u);  // Only the real part has a value
    }
    else {
      stk::all_reduce_max(comm, &localMaxValue, &result, 1u);
    }
  });
}

template <typename Scalar, typename XDATA> requires is_complex_v<Scalar>
void field_amin_on_host(Scalar& localMinValue, const mesh::BucketVector& buckets, const XDATA& xData)
{
  using ComplexType = typename Scalar::value_type;
  ComplexType localMinValueReal = std::numeric_limits<ComplexType>::max();

#ifdef STK_USE_OPENMP
  #pragma omp parallel for schedule(static) reduction(min:localMinValueReal)
#endif
  for (unsigned bucketIndex = 0; bucketIndex < buckets.size(); ++bucketIndex) {
    const mesh::Bucket& bucket = *buckets[bucketIndex];
    auto xValues = xData.bucket_values(bucket);

    for (mesh::EntityIdx entityIdx : bucket.entities()) {
      for (mesh::ScalarIdx scalarIdx : xValues.scalars()) {
        ComplexType localValueReal = Kokkos::abs(xValues(entityIdx, scalarIdx));
        if (localValueReal < localMinValueReal) {
          localMinValueReal = localValueReal;
        }
      }
    }
  }

  localMinValue = Scalar(localMinValueReal, 0);
}

template <typename Scalar, typename XDATA>
void field_amin_on_host(Scalar& localMinValue, const mesh::BucketVector& buckets, const XDATA& xData)
{
  localMinValue = std::numeric_limits<Scalar>::max();

#ifdef STK_USE_OPENMP
    #pragma omp parallel for schedule(static) reduction(min:localMinValue)
#endif
  for (unsigned bucketIndex = 0; bucketIndex < buckets.size(); ++bucketIndex) {
    const mesh::Bucket& bucket = *buckets[bucketIndex];
    auto xValues = xData.bucket_values(bucket);

    for (mesh::EntityIdx entityIdx : bucket.entities()) {
      for (mesh::ScalarIdx scalarIdx : xValues.scalars()) {
        Scalar localValue = Kokkos::abs(xValues(entityIdx, scalarIdx));

        if (localValue < localMinValue) {
          localMinValue = localValue;
        }
      }
    }
  }
}

template <typename T>
void field_amin_impl(T& result, const mesh::FieldBase& xField, const mesh::Selector& selector)
{
  const auto& buckets = xField.get_mesh().get_buckets(xField.entity_rank(),
                                                      selector & xField.mesh_meta_data().locally_owned_part());
  const MPI_Comm comm = xField.get_mesh().parallel();

  field_data_execute<T, mesh::ReadOnly>(xField, [&](auto& xData) {
    T localMinValue;
    field_amin_on_host(localMinValue, buckets, xData);

    result = localMinValue;
    if constexpr (is_complex_v<T>) {
      using ComplexType = typename T::value_type;
      const ComplexType* localResultArray = reinterpret_cast<ComplexType*>(&localMinValue);
      ComplexType* globalResultArray = reinterpret_cast<ComplexType*>(&result);
      stk::all_reduce_min(comm, localResultArray, globalResultArray, 1u);  // Only the real part has a value
    }
    else {
      stk::all_reduce_min(comm, &localMinValue, &result, 1u);
    }
  });
}

} // namespace stk::field_blas::impl

#endif // STK_MESH_BASEIMPL_FIELDBLASIMPL_HPP
