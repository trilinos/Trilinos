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

#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_util/util/Fortran.hpp> // For SIERRA_FORTRAN
#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_mesh/base/MetaData.hpp>

#include <complex>
#include <string>
#include <iostream>
#include <algorithm>

#if defined(_OPENMP) && !defined(__INTEL_COMPILER)
#define OPEN_MP_ACTIVE_FIELDBLAS_HPP
// there seems to be an issue with OpenMP combined with GoogleTest macros with Intel compilers
// example of error:
//    openMP.C(206): internal error: assertion failed at: "shared/cfe/edgcpfe/checkdir.c", line 5531
#include <omp.h>

#endif

extern "C"
{
void SIERRA_FORTRAN(daxpy)(const int *n, const double *dscale, const double x[], const int *incx, double y[],const int *incy); // y=y+dscale*x
void SIERRA_FORTRAN(dscal)(const int *n, const double *dscale, double *vect, const int *inc); //vect = dscale * vect
double SIERRA_FORTRAN(ddot)(const int * n, const double* x, const int * incx, const double* y, const int * incy); // < x , y >
void SIERRA_FORTRAN(sscal)(const int *n, const float *sscale, float *vect, const int *inc); //vect = sscale * vect
void SIERRA_FORTRAN(dswap)(const int* n, double* d, const int* inc, double* d1, const int* inc1); // switch d1 , d // D.N.E.
double SIERRA_FORTRAN(dasum)(const int * n,const double * x,const int * incx);
int SIERRA_FORTRAN(idamax)(const int *n, const double *vect, const int *inc);
void SIERRA_FORTRAN(saxpy)(const int *n, const float *xscale, const float x[], const int *incx, float y[],const int *incy); // y=y+sscale*x
float SIERRA_FORTRAN(sdot)(const int * n, const float* x, const int * incx, const float* y, const int * incy); // < x , y >
float SIERRA_FORTRAN(sasum)(const int * n,const float * x,const int * incx);
void SIERRA_FORTRAN(sswap)(const int* n, float* s, const int* inc, float* s1, const int* inc1); // switch s1 , s // D.N.E.
int SIERRA_FORTRAN(isamax)(const int *n, const float *vect, const int *inc);
}

namespace stk {
namespace mesh {

template<typename Scalar>
struct BucketSpan {
  Scalar* ptr     = nullptr;
  unsigned length = 0;

  unsigned size() const { return length; }

  const Scalar* data() const { return ptr; }
        Scalar* data()       { return ptr; }

        Scalar& operator[](unsigned i)       noexcept { return ptr[i]; }
  const Scalar& operator[](unsigned i) const noexcept { return ptr[i]; }

  template<typename Field_t>
  BucketSpan(const Field_t & f, const Bucket& b)
  {
    const FieldMetaData& field_meta_data = f.get_meta_data_for_field()[b.bucket_id()];
    ptr    = static_cast<Scalar*>(reinterpret_cast<typename FieldTraits<Field_t>::data_type*>(field_meta_data.m_data));
    length = b.size() * field_meta_data.m_bytes_per_entity/f.data_traits().size_of;
  }

  BucketSpan& operator=(const BucketSpan& B) noexcept
  {
    for(size_t i = 0; i < length; ++i) {
      ptr[i] = B.ptr[i];
    }
    return *this;
  }
};

template<typename Field_t>
inline bool is_compatible(const Field_t& x, const Field_t& y)
{
  if(&x.get_mesh() != &y.get_mesh()) return false;
  if(x.entity_rank() != y.entity_rank()) return false;
  if(x.data_traits().type_info != y.data_traits().type_info) return false;
  return true;
}

template<typename Scalar, typename Field_t>
inline bool is_compatible(const Field_t& x)
{
    if( x.data_traits().type_info != typeid(Scalar) ) return false;
    return true;
}

template<class Scalar>
struct FortranBLAS
{
    inline static void axpy(int kmax, const Scalar alpha, const Scalar x[], Scalar y[])
    {
        for(int k = 0; k < kmax; ++k) {
            y[k] += alpha * x[k];
        }
    }

    inline static void axpby(int kmax, const Scalar alpha, const Scalar x[], const Scalar beta, Scalar y[])
    {
        for(int k = 0; k < kmax; ++k) {
          y[k] *= beta;
        }
        for(int k = 0; k < kmax; ++k) {
            y[k] += alpha * x[k];
        }
    }

    inline static void copy(int kmax, const Scalar x[], Scalar y[])
    {
        for(int k = 0; k < kmax; ++k) {
            y[k] = x[k];
        }
    }

    inline static void product(int kmax, const Scalar x[], const Scalar y[], Scalar z[])
    {
        for (int k = 0; k < kmax; ++k) {
            z[k] = x[k]*y[k];
        }
    }

    inline static Scalar dot(int kmax, const Scalar x[], const Scalar y[])
    {
        Scalar result(0.0);
        for(int k = 0; k < kmax; ++k) {
            result += y[k] * x[k];
        }
        return result;
    }

    inline static Scalar nrm2(int kmax, const Scalar x[])
    {
        Scalar result(0.0);
        for(int k = 0; k < kmax; ++k) {
            result += std::pow(std::abs(x[k]),2);
        }
        return Scalar(std::sqrt(result));
    }

    inline static void scal(int kmax, const Scalar alpha, Scalar x[])
    {
        for(int k = 0; k < kmax; ++k) {
            x[k] = alpha * x[k];
        }
    }

    inline static void fill(int kmax, const Scalar alpha, Scalar x[], const int inc=1)
    {
        auto ke = kmax*inc;
        for(int k = 0; k < ke ; k += inc) {
            x[k] = alpha;
        }
    }

    inline static void swap(int kmax, Scalar x[], Scalar y[])
    {
        Scalar temp;
        for(int k = 0; k < kmax; ++k) {
            temp = y[k];
            y[k] = x[k];
            x[k] = temp;
        }
    }

    inline static Scalar asum(int kmax, const Scalar x[])
    {
        Scalar result(0.0);
        for(int k = 0; k < kmax; ++k) {
            result += std::abs(x[k]);
        }
        return Scalar(result);
    }

    inline static int iamax(int kmax, const Scalar x[])
    {
        double amax = 0.0;
        int result = 0;
        for(int k = 0; k < kmax; ++k) {
            if (amax < std::abs(x[k])) {
                result = k;
                amax = std::abs(x[k]);
            }
        }
        return result;
    }

    inline static int iamin(int kmax, const Scalar x[])
    {
        int result = 0;
        double amin = std::abs(x[0]);
        for(int k = 0; k < kmax; ++k) {
            if (std::abs(x[k]) < amin) {
                result = k;
                amin = std::abs(x[k]);
            }
        }
        return result;
    }
};

template<class Scalar>
struct FortranBLAS<std::complex<Scalar> >
{
    inline static void axpy(int kmax, const std::complex<Scalar> alpha, const std::complex<Scalar> x[], std::complex<Scalar> y[])
    {
        for(int k = 0; k < kmax; ++k) {
            y[k] += alpha * x[k];
        }
    }

    inline static void axpby(int kmax, const std::complex<Scalar> alpha, const std::complex<Scalar> x[], const std::complex<Scalar>& beta, std::complex<Scalar> y[])
    {
        for(int k = 0; k < kmax; ++k) {
          y[k] *= beta;
        }
        for(int k = 0; k < kmax; ++k) {
            y[k] += alpha * x[k];
        }
    }

    inline static void product(int kmax, const std::complex<Scalar> x[], const std::complex<Scalar> y[], std::complex<Scalar> z[])
    {
        for (int k = 0; k < kmax; ++k)
        {
            z[k] = x[k]*y[k];
        }
    }

    inline static void copy(int kmax, const std::complex<Scalar> x[], std::complex<Scalar> y[])
    {
        for(int k = 0; k < kmax; ++k) {
            y[k] = x[k];
        }
    }

    inline static std::complex<Scalar> dot(int kmax, const std::complex<Scalar> x[], const std::complex<Scalar> y[]) 
    {
        std::complex<Scalar> result = std::complex<Scalar>(0.0);
        for(int k = 0; k < kmax; ++k) {
            result += y[k] * x[k];
        }
        return result;
    }

    inline static std::complex<Scalar> nrm2(int kmax, const std::complex<Scalar> x[]) 
    {
        Scalar result(0.0);
        for(int k = 0; k < kmax; ++k) {
            result += std::pow(std::abs(x[k]),2);
        }
        return std::complex<Scalar>(std::sqrt(result));
    }

    inline static void scal(int kmax, const std::complex<Scalar> alpha, std::complex<Scalar> x[])
    {
        for(int k = 0; k < kmax; ++k) {
            x[k] = alpha * x[k];
        }
    }

    inline static void fill(int kmax, const std::complex<Scalar> alpha, std::complex<Scalar> x[], const int inc=1)
    {
        auto ke = kmax*inc;
        for(int k = 0; k < ke ; k += inc) {
            x[k] = alpha;
        }
    }

    inline static void swap(int kmax, std::complex<Scalar> x[], std::complex<Scalar> y[])
    {
        std::complex<Scalar> temp;
        for(int k = 0; k < kmax; ++k) {
            temp = y[k];
            y[k] = x[k];
            x[k] = temp;
        }
    }

    inline static std::complex<Scalar> asum(int kmax, const std::complex<Scalar> x[])
    {
        Scalar result(0.0);
        for(int k = 0; k < kmax; ++k) {
            result += std::abs(x[k]);
        }
        return std::complex<Scalar>(result,0.0);
    }

    inline static int iamax(int kmax, const std::complex<Scalar> x[]) 
    {
        Scalar amax(0.0);
        int result = 0;
        for(int k = 0; k < kmax; ++k) {
            if (amax < std::norm(x[k])) {
                result = k;
                amax = std::norm(x[k]);
            }
        }
        return result;
    }

    inline static int iamin(int kmax, const std::complex<Scalar> x[]) 
    {
        int result = 0;
        Scalar amin = std::norm(x[0]);
        for(int k = 0; k < kmax; ++k) {
            if (std::norm(x[k]) < amin) {
                result = k;
                amin = std::norm(x[k]);
            }
        }
        return result;
    }

};

template<>
struct FortranBLAS<double>
{
    inline static void axpy(int kmax, const double alpha, const double x[], double y[])
    {
        const int one = 1;
        SIERRA_FORTRAN(daxpy)(&kmax,&alpha,x,&one,y,&one);
    }

    inline static void axpby(int kmax, const double alpha, const double x[], const double beta, double y[])
    {
        for(int k = 0; k < kmax; ++k) {
          y[k] *= beta;
        }
        for(int k = 0; k < kmax; ++k) {
            y[k] += alpha * x[k];
        }
    }

    inline static void product(int kmax, const double x[], const double y[], double z[])
    {
        for (int k = 0; k < kmax; ++k)
        {
            z[k] = x[k]*y[k];
        }
    }

    inline static void copy(int kmax, const double x[], double y[])
    {
        for (int k = 0; k < kmax; ++k) {
            y[k] = x[k];
        }
    }

    inline static double dot(int kmax, const double x[], const double y[])
    {
        const int one = 1;
        return SIERRA_FORTRAN(ddot)(&kmax,x,&one,y,&one);
    }

    inline static double nrm2(int kmax, const double x[])
    {
        return std::sqrt(dot(kmax, x, x));
    }

    inline static void scal(int kmax, const double alpha, double x[])
    {
        const int one = 1;
        SIERRA_FORTRAN(dscal)(&kmax,&alpha,x,&one);
    }

    inline static void fill(int kmax, const double alpha, double x[], const int inc=1)
    {
        auto ke = kmax*inc;
        for(int k = 0; k < ke ; k += inc) {
            x[k] = alpha;
        }
    }

    inline static void swap(int kmax, double x[], double y[])
    {
        const int one = 1;
        SIERRA_FORTRAN(dswap)(&kmax,x,&one,y,&one);
    }

    inline static double asum(int kmax, const double x[])
    {
        const int one = 1;
        return SIERRA_FORTRAN(dasum)(&kmax,x,&one);
    }

    inline static int iamax(int kmax, const double x[])
    {
        const int one = 1;
        return (SIERRA_FORTRAN(idamax)(&kmax, x, &one) - 1);
    }

    inline static int iamin(int kmax, const double x[])
    {
        int result = 0;
        double amin = std::abs(x[0]);
        for(int k = 0; k < kmax; ++k) {
            if (std::abs(x[k]) < amin) {
                result = k;
                amin = std::abs(x[k]);
            }
        }
        return result;
    }

};

template<>
struct FortranBLAS<float>
{
    inline static void axpy(int kmax, const float alpha, const float x[], float y[])
    {
        const int one = 1;
        SIERRA_FORTRAN(saxpy)(&kmax,&alpha,x,&one,y,&one);
    }

    inline static void axpby(int kmax, const float alpha, const float x[], const float& beta, float y[])
    {
        for(int k = 0; k < kmax; ++k) {
          y[k] *= beta;
        }
        for(int k = 0; k < kmax; ++k) {
            y[k] += alpha * x[k];
        }
    }

    inline static void product(int kmax, const float x[], const float y[], float z[])
    {
        for (int k = 0; k < kmax; ++k) {
            z[k] = x[k]*y[k];
        }
    }

    inline static void copy(int kmax, const float x[], float y[])
    {
        for (int k = 0; k < kmax; ++k) {
            y[k] = x[k];
        }
    }

    inline static float dot(int kmax, const float x[], const float y[])
    {
        const int one = 1;
        return static_cast<float>(SIERRA_FORTRAN(sdot)(&kmax,x,&one,y,&one));
    }

    inline static float nrm2(int kmax, const float x[])
    {
        return std::sqrt(dot(kmax, x, x));
    }

    inline static void scal(int kmax, const float alpha, float x[])
    {
        const int one = 1;
        SIERRA_FORTRAN(sscal)(&kmax,&alpha,x,&one);
    }

    inline static void fill(int kmax, const float alpha, float x[], const int inc=1)
    {
        auto ke = kmax*inc;
        for(int k = 0; k < ke ; k += inc) {
            x[k] = alpha;
        }
    }

    inline static void swap(int kmax, float x[], float y[])
    {
        const int one = 1;
        SIERRA_FORTRAN(sswap)(&kmax,x,&one,y,&one);
    }

    inline static float asum(int kmax, const float x[])
    {
        const int one = 1;
        return static_cast<float>(SIERRA_FORTRAN(sasum)(&kmax,x,&one));
    }

    inline static int iamax(int kmax, const float x[])
    {
        const int one = 1;
        return (SIERRA_FORTRAN(isamax)(&kmax, x, &one) - 1);
    }

    inline static int iamin(int kmax, const float x[])
    {
        int result = 0;
        float amin = std::abs(x[0]);
        for(int k = 0; k < kmax; ++k) {
            if (std::abs(x[k]) < amin) {
                result = k;
                amin = std::abs(x[k]);
            }
        }
        return result;
    }

};

inline int fix_omp_threads()
{
  int orig_thread_count = 0;
#ifdef OPEN_MP_ACTIVE_FIELDBLAS_HPP
    orig_thread_count = omp_get_max_threads();
    if(omp_get_max_threads() >= omp_get_num_procs()) {omp_set_num_threads(omp_get_num_procs());}
#endif
    return orig_thread_count;
}

inline void unfix_omp_threads(int orig_thread_count)
{
#ifdef OPEN_MP_ACTIVE_FIELDBLAS_HPP
    omp_set_num_threads(orig_thread_count);
#endif
}

template<class Scalar>
inline
void field_axpy(const Scalar alpha, const FieldBase& xField, const FieldBase& yField, const Selector& selector)
{
    ThrowAssert(is_compatible(xField, yField));

    BucketVector const& buckets = xField.get_mesh().get_buckets( xField.entity_rank(), selector );

    int orig_thread_count = fix_omp_threads();
#ifdef OPEN_MP_ACTIVE_FIELDBLAS_HPP
#pragma omp parallel for schedule(static)
#endif
    for(size_t i=0; i < buckets.size(); i++) {
        Bucket & b = *buckets[i];
        BucketSpan<Scalar> x(xField, b);
        BucketSpan<Scalar> y(yField, b);
        ThrowAssert(x.size() == y.size());
        FortranBLAS<Scalar>::axpy(x.size(),alpha,x.data(),y.data());
    }
    unfix_omp_threads(orig_thread_count);
}

template<class Scalar>
inline
void field_axpy(const Scalar alpha, const FieldBase& xField, const FieldBase& yField)
{
    const Selector selector = selectField(xField) & selectField(yField);
    field_axpy(alpha,xField,yField,selector);
}

template<class Scalar>
inline
void field_axpby(const Scalar alpha, const FieldBase& xField, const Scalar beta, const FieldBase& yField, const Selector& selector)
{
    ThrowAssert(is_compatible(xField, yField));

    BucketVector const& buckets = xField.get_mesh().get_buckets( xField.entity_rank(), selector );

    int orig_thread_count = fix_omp_threads();
#ifdef OPEN_MP_ACTIVE_FIELDBLAS_HPP
#pragma omp parallel for schedule(static)
#endif
    for(size_t i=0; i < buckets.size(); i++) {
        Bucket & b = *buckets[i];
        BucketSpan<Scalar> x(xField, b);
        BucketSpan<Scalar> y(yField, b);
        ThrowAssert(x.size() == y.size());
        FortranBLAS<Scalar>::axpby(x.size(),alpha,x.data(),beta,y.data());
    }
    unfix_omp_threads(orig_thread_count);
}

template<class Scalar>
inline
void field_axpby(const Scalar alpha, const FieldBase& xField, const Scalar beta, const FieldBase& yField)
{
    const Selector selector = selectField(xField) & selectField(yField);
    field_axpby(alpha,xField,beta,yField,selector);
}

template<class Scalar>
inline
void INTERNAL_field_product(const FieldBase& xField, const FieldBase& yField, const FieldBase& zField, const Selector& selector)
{
    BucketVector const& buckets = xField.get_mesh().get_buckets( xField.entity_rank(), selector );

    int orig_thread_count = fix_omp_threads();
#ifdef OPEN_MP_ACTIVE_FIELDBLAS_HPP
#pragma omp parallel for schedule(static)
#endif
    for(size_t i=0; i < buckets.size(); i++) {
        Bucket & b = *buckets[i];
        BucketSpan<Scalar> x(xField, b);
        BucketSpan<Scalar> y(yField, b);
        BucketSpan<Scalar> z(zField, b);
        FortranBLAS<Scalar>::product(x.size(),x.data(),y.data(),z.data());
    }
    unfix_omp_threads(orig_thread_count);
}

inline
void field_product(const FieldBase& xField, const FieldBase& yField, const FieldBase& zField, const Selector& selector)
{
    ThrowAssert(is_compatible(xField, yField));
    ThrowAssert(is_compatible(yField, zField));

    if (xField.data_traits().type_info == typeid(double)) {
        INTERNAL_field_product<double>(xField,yField,zField,selector);
    } else if (xField.data_traits().type_info == typeid(float)) {
        INTERNAL_field_product<float>(xField,yField,zField,selector);
    } else if (xField.data_traits().type_info == typeid(std::complex<double>)) {
        INTERNAL_field_product<std::complex<double> >(xField,yField,zField,selector);
    } else if (xField.data_traits().type_info == typeid(std::complex<float>)) {
        INTERNAL_field_product<std::complex<float> >(xField,yField,zField,selector);
    } else if (xField.data_traits().type_info == typeid(int)) {
        INTERNAL_field_product<int>(xField,yField,zField,selector);
    } else {
        ThrowAssertMsg(false,"Error in field_product; field is of type "<<xField.data_traits().type_info.name()<<" which is not supported");
    }
}

inline
void field_product(const FieldBase& xField, const FieldBase& yField, const FieldBase& zField)
{
    const Selector selector = selectField(xField) & selectField(yField);
    field_product(xField,yField,zField,selector);
}

template<class Scalar>
inline
void INTERNAL_field_copy(const FieldBase& xField, const FieldBase& yField, const Selector& selector)
{
    BucketVector const& buckets = xField.get_mesh().get_buckets( xField.entity_rank(), selector );

    int orig_thread_count = fix_omp_threads();
#ifdef OPEN_MP_ACTIVE_FIELDBLAS_HPP
#pragma omp parallel for schedule(static)
#endif
    for(size_t i=0; i < buckets.size(); i++) {
        Bucket & b = *buckets[i];
        BucketSpan<Scalar> x(xField, b);
        BucketSpan<Scalar> y(yField, b);
        y = x;
    }
    unfix_omp_threads(orig_thread_count);
}

inline
void field_copy(const FieldBase& xField, const FieldBase& yField, const Selector& selector)
{
    ThrowAssert(is_compatible(xField, yField));

    if (xField.data_traits().type_info == typeid(double)) {
        INTERNAL_field_copy<double>(xField,yField,selector);
    } else if (xField.data_traits().type_info == typeid(float)) {
        INTERNAL_field_copy<float>(xField,yField,selector);
    } else if (xField.data_traits().type_info == typeid(std::complex<double>)) {
        INTERNAL_field_copy<std::complex<double> >(xField,yField,selector);
    } else if (xField.data_traits().type_info == typeid(std::complex<float>)) {
        INTERNAL_field_copy<std::complex<float> >(xField,yField,selector);
    } else if (xField.data_traits().type_info == typeid(int)) {
        INTERNAL_field_copy<int>(xField,yField,selector);
    } else {
        ThrowAssertMsg(false,"Error in field_copy; field is of type "<<xField.data_traits().type_info.name()<<" which is not supported");
    }
}

inline
void field_copy(const FieldBase& xField, const FieldBase& yField)
{
    const Selector selector = selectField(xField) & selectField(yField);
    field_copy(xField,yField,selector);
}

template<class Scalar,class T1,class T2,class T3,class T4,class T5,class T6,class T7>
inline
Scalar field_dot(const Field<Scalar,T1,T2,T3,T4,T5,T6,T7> & xField, const Field<Scalar,T1,T2,T3,T4,T5,T6,T7> & yField, const Selector& selector, const MPI_Comm comm)
{
    ThrowAssert(is_compatible(xField, yField));

    BucketVector const& buckets = xField.get_mesh().get_buckets(xField.entity_rank(), selector & xField.mesh_meta_data().locally_owned_part());

    Scalar local_result(0.0);

    int orig_thread_count = fix_omp_threads();
#ifdef OPEN_MP_ACTIVE_FIELDBLAS_HPP
#pragma omp parallel for schedule(static) reduction(+:local_result)
#endif
    for(size_t i=0; i < buckets.size(); i++) {
        Bucket & b = *buckets[i];
        BucketSpan<Scalar> x(xField, b);
        BucketSpan<Scalar> y(yField, b);
        ThrowAssert(x.size() == y.size());
        local_result += FortranBLAS<Scalar>::dot(x.size(),x.data(),y.data());
    }

    Scalar glob_result = local_result;
    stk::all_reduce_sum(comm,&local_result,&glob_result,1u);
    unfix_omp_threads(orig_thread_count);
    return glob_result;
}

template<class Scalar,class T1,class T2,class T3,class T4,class T5,class T6,class T7>
inline
std::complex<Scalar> field_dot(const Field<std::complex<Scalar>,T1,T2,T3,T4,T5,T6,T7>& xField, const Field<std::complex<Scalar>,T1,T2,T3,T4,T5,T6,T7>& yField, const Selector& selector, const MPI_Comm comm) 
{
    ThrowAssert(is_compatible(xField, yField));

    BucketVector const& buckets = xField.get_mesh().get_buckets(xField.entity_rank(), selector & xField.mesh_meta_data().locally_owned_part());

    Scalar local_result_r (0.0);
    Scalar local_result_i (0.0);
    std::complex<Scalar> priv_tmp;

    int orig_thread_count = fix_omp_threads();
#ifdef OPEN_MP_ACTIVE_FIELDBLAS_HPP
#pragma omp parallel for reduction(+:local_result_r,local_result_i) schedule(static) private(priv_tmp)
#endif
    for(size_t i=0; i < buckets.size(); i++) {
        Bucket & b = *buckets[i];
        BucketSpan<std::complex<Scalar>> x(xField, b);
        BucketSpan<std::complex<Scalar>> y(yField, b);
        ThrowAssert(x.size() == y.size());
        priv_tmp = FortranBLAS<std::complex<Scalar> >::dot(x.size(),x.data(),y.data());

        local_result_r += priv_tmp.real();
        local_result_i += priv_tmp.imag();
    }

    Scalar local_result_ri [2] = { local_result_r     , local_result_i     };
    Scalar  glob_result_ri [2] = { local_result_ri[0] , local_result_ri[1] };
    stk::all_reduce_sum(comm,local_result_ri,glob_result_ri,2u);
    unfix_omp_threads(orig_thread_count);
    return std::complex<Scalar> (glob_result_ri[0],glob_result_ri[1]);
}

template<class Scalar,class T1,class T2,class T3,class T4,class T5,class T6,class T7>
inline
Scalar field_dot(const Field<Scalar,T1,T2,T3,T4,T5,T6,T7> & xField, const Field<Scalar,T1,T2,T3,T4,T5,T6,T7> & yField, const Selector& selector)
{
    const MPI_Comm comm = xField.get_mesh().parallel();
    return field_dot(xField,yField,selector,comm);
}

template<class Scalar,class T1,class T2,class T3,class T4,class T5,class T6,class T7>
inline
Scalar field_dot(const Field<Scalar,T1,T2,T3,T4,T5,T6,T7> & xField, const Field<Scalar,T1,T2,T3,T4,T5,T6,T7> & yField)
{
    const Selector selector = selectField(xField) & selectField(yField);
    return field_dot(xField,yField,selector);
}

template<class Scalar>
inline
void field_dot(std::complex<Scalar>& global_result, const FieldBase& xField, const FieldBase& yField, const Selector& selector, const MPI_Comm comm)
{
    ThrowAssert(is_compatible(xField, yField));

    BucketVector const& buckets = xField.get_mesh().get_buckets(xField.entity_rank(), selector & xField.mesh_meta_data().locally_owned_part());

    Scalar local_result_r (0.0);
    Scalar local_result_i (0.0);
    std::complex<Scalar> priv_tmp;

    int orig_thread_count = fix_omp_threads();
#ifdef OPEN_MP_ACTIVE_FIELDBLAS_HPP
#pragma omp parallel for reduction(+:local_result_r,local_result_i) schedule(static) private(priv_tmp)
#endif
    for(size_t i=0; i < buckets.size(); i++) {
        Bucket & b = *buckets[i];
        BucketSpan<std::complex<Scalar>> x(xField, b);
        BucketSpan<std::complex<Scalar>> y(yField, b);
        ThrowAssert(x.size() == y.size());
        priv_tmp = FortranBLAS<std::complex<Scalar> >::dot(x.size(),x.data(),y.data());

        local_result_r += priv_tmp.real();
        local_result_i += priv_tmp.imag();
    }

    Scalar local_result_ri [2] = { local_result_r     , local_result_i     };
    Scalar  glob_result_ri [2] = { local_result_ri[0] , local_result_ri[1] };
    stk::all_reduce_sum(comm,local_result_ri,glob_result_ri,2u);
    global_result = std::complex<Scalar> (glob_result_ri[0],glob_result_ri[1]);
    unfix_omp_threads(orig_thread_count);
}

template<class Scalar>
inline
void field_dot(Scalar& glob_result, const FieldBase& xField, const FieldBase& yField, const Selector& selector, const MPI_Comm comm)
{
    ThrowAssert(is_compatible(xField, yField));

    BucketVector const& buckets = xField.get_mesh().get_buckets(xField.entity_rank(), selector & xField.mesh_meta_data().locally_owned_part());

    Scalar local_result(0.0);

    int orig_thread_count = fix_omp_threads();
#ifdef OPEN_MP_ACTIVE_FIELDBLAS_HPP
#pragma omp parallel for reduction(+:local_result) schedule(static)
#endif
    for(size_t i=0; i < buckets.size(); i++) {
        Bucket & b = *buckets[i];
        BucketSpan<Scalar> x(xField, b);
        BucketSpan<Scalar> y(yField, b);
        ThrowAssert(x.size() == y.size());
        local_result += FortranBLAS<Scalar>::dot(x.size(),x.data(),y.data());
    }

    glob_result = local_result;
    stk::all_reduce_sum(comm,&local_result,&glob_result,1u);
    unfix_omp_threads(orig_thread_count);
}

template<class Scalar>
inline
void field_dot(Scalar& result, const FieldBase& xField, const FieldBase& yField, const Selector& selector)
{
    const MPI_Comm comm = xField.get_mesh().parallel();
    field_dot(result,xField,yField,selector,comm);
}

template<class Scalar>
inline
void field_dot(Scalar& result, const FieldBase& xField, const FieldBase& yField)
{
    const Selector selector = selectField(xField) & selectField(yField);
    field_dot(result,xField,yField,selector);
}

template<class Scalar>
inline
void field_scale(const Scalar alpha, const FieldBase& xField, const Selector& selector)
{
    ThrowAssert( is_compatible<Scalar>(xField) );

    BucketVector const& buckets = xField.get_mesh().get_buckets(xField.entity_rank(),selector);

    int orig_thread_count = fix_omp_threads();
#ifdef OPEN_MP_ACTIVE_FIELDBLAS_HPP
#pragma omp parallel for schedule(static)
#endif
    for(size_t i=0; i < buckets.size(); i++) {
        BucketSpan<Scalar> x(xField, *buckets[i]);
        FortranBLAS<Scalar>::scal(x.size(),alpha,x.data());
    }
    unfix_omp_threads(orig_thread_count);
}

template<class Scalar>
inline
void field_scale(const Scalar alpha, const FieldBase& xField)
{
    const Selector selector = selectField(xField);
    field_scale(alpha,xField,selector);
}

template<class Scalar>
inline
void field_fill_component(const Scalar* alpha, const FieldBase& xField, const Selector& selector)
{
    ThrowAssert( is_compatible<Scalar>(xField) );

    BucketVector const& buckets = xField.get_mesh().get_buckets(xField.entity_rank(),selector);

    int orig_thread_count = fix_omp_threads();
#ifdef OPEN_MP_ACTIVE_FIELDBLAS_HPP
#pragma omp parallel for schedule(static)
#endif
    for(size_t i = 0; i < buckets.size(); i++) {
        const Bucket & b = *buckets[i];
        const int length = b.size();
        const unsigned int fieldSize = field_scalars_per_entity(xField, b);
        Scalar * x = static_cast<Scalar*>(field_data(xField, b));

        for (unsigned int j = 0; j < fieldSize; j++) {
            FortranBLAS<Scalar>::fill(length,alpha[j],x+j,fieldSize);
        }
    }
    unfix_omp_threads(orig_thread_count);
}

template<class Scalar>
inline
void field_fill_component(const Scalar* alpha, const FieldBase& xField)
{
    const Selector selector = selectField(xField);
    field_fill_component(alpha,xField,selector);
}

template<class Scalar>
inline
void field_fill(const Scalar alpha, const FieldBase& xField, const Selector& selector)
{
    ThrowAssert( is_compatible<Scalar>(xField) );

    BucketVector const& buckets = xField.get_mesh().get_buckets(xField.entity_rank(),selector);

    int orig_thread_count = fix_omp_threads();
#ifdef OPEN_MP_ACTIVE_FIELDBLAS_HPP
#pragma omp parallel for schedule(static)
#endif
    for(size_t i = 0; i < buckets.size(); i++) {
        BucketSpan<Scalar> x(xField, *buckets[i]);
        FortranBLAS<Scalar>::fill(x.size(),alpha,x.data());
    }
    unfix_omp_threads(orig_thread_count);
}

template<class Scalar>
inline
void field_fill(const Scalar alpha, const FieldBase& xField)
{
    const Selector selector = selectField(xField);
    field_fill(alpha,xField,selector);
}

template<class Scalar>
inline
void INTERNAL_field_swap(const FieldBase& xField, const FieldBase& yField, const Selector& selector)
{
    BucketVector const& buckets = xField.get_mesh().get_buckets( xField.entity_rank(), selector );

    int orig_thread_count = fix_omp_threads();
#ifdef OPEN_MP_ACTIVE_FIELDBLAS_HPP
#pragma omp parallel for schedule(static)
#endif
    for(size_t i=0; i < buckets.size(); i++) {
        Bucket & b = *buckets[i];
        BucketSpan<Scalar> x(xField, b);
        BucketSpan<Scalar> y(yField, b);
        ThrowAssert(x.size() == y.size());
        FortranBLAS<Scalar>::swap(x.size(),x.data(),y.data());
    }
    unfix_omp_threads(orig_thread_count);
}

inline
void field_swap(const FieldBase& xField, const FieldBase& yField, const Selector& selector)
{
    ThrowAssert(is_compatible(xField, yField));

    if (xField.data_traits().type_info == typeid(double)) {
        INTERNAL_field_swap<double>(xField,yField,selector);
    } else if (xField.data_traits().type_info == typeid(float)) {
        INTERNAL_field_swap<float>(xField,yField,selector);
    } else if (xField.data_traits().type_info == typeid(std::complex<double>)) {
        INTERNAL_field_swap<std::complex<double> >(xField,yField,selector);
    } else if (xField.data_traits().type_info == typeid(std::complex<float>)) {
        INTERNAL_field_swap<std::complex<float> >(xField,yField,selector);
    } else if (xField.data_traits().type_info == typeid(int)) {
        INTERNAL_field_swap<int>(xField,yField,selector);
    } else {
        ThrowAssertMsg(false,"Error in field_swap; field is of type "<<xField.data_traits().type_info.name()<<" which is not supported");
    }
}

inline
void field_swap(const FieldBase& xField, const FieldBase& yField)
{
    const Selector selector = selectField(xField) & selectField(yField);
    field_swap(xField,yField,selector);
}

template<class Scalar,class T1,class T2,class T3,class T4,class T5,class T6,class T7>
inline
Scalar field_nrm2(const Field<Scalar,T1,T2,T3,T4,T5,T6,T7> & xField, const Selector& selector, const MPI_Comm comm)
{
    BucketVector const& buckets = xField.get_mesh().get_buckets(xField.entity_rank(), selector & xField.mesh_meta_data().locally_owned_part());

    Scalar local_result(0.0);

    int orig_thread_count = fix_omp_threads();
#ifdef OPEN_MP_ACTIVE_FIELDBLAS_HPP
#pragma omp parallel for reduction(+:local_result) schedule(static)
#endif
    for(size_t i=0; i < buckets.size(); i++) {
        BucketSpan<Scalar> x(xField, *buckets[i]);
        local_result += FortranBLAS<Scalar>::dot(x.size(),x.data(),x.data());
    }

    Scalar glob_result = local_result;
    stk::all_reduce_sum(comm,&local_result,&glob_result,1u);
    unfix_omp_threads(orig_thread_count);
    return std::sqrt(glob_result);
}

template<class Scalar,class T1,class T2,class T3,class T4,class T5,class T6,class T7>
inline
std::complex<Scalar> field_nrm2(const Field< std::complex<Scalar>,T1,T2,T3,T4,T5,T6,T7> & xField, const Selector& selector, const MPI_Comm comm) 
{
    BucketVector const& buckets = xField.get_mesh().get_buckets(xField.entity_rank(), selector & xField.mesh_meta_data().locally_owned_part());

    Scalar local_result_r(0.0);

    int orig_thread_count = fix_omp_threads();
#ifdef OPEN_MP_ACTIVE_FIELDBLAS_HPP
#pragma omp parallel for reduction(+:local_result_r) schedule(static)
#endif
    for(size_t i=0; i < buckets.size(); i++) {
        BucketSpan<std::complex<Scalar>> x(xField, *buckets[i]);
        local_result_r += std::pow(FortranBLAS<std::complex<Scalar> >::nrm2(x.size(),x.data()).real(),2);
    }

    Scalar glob_result = local_result_r;
    stk::all_reduce_sum(comm,&local_result_r,&glob_result,1u);
    unfix_omp_threads(orig_thread_count);
    return std::complex<Scalar>(std::sqrt(glob_result),0.0);
}

template<class Scalar,class T1,class T2,class T3,class T4,class T5,class T6,class T7>
inline
Scalar field_nrm2(const Field<Scalar,T1,T2,T3,T4,T5,T6,T7> & xField, const Selector& selector)
{
    const MPI_Comm comm = xField.get_mesh().parallel();
    return field_nrm2(xField,selector,comm);
}

template<class Scalar,class T1,class T2,class T3,class T4,class T5,class T6,class T7>
inline
Scalar field_nrm2(const Field<Scalar,T1,T2,T3,T4,T5,T6,T7> & xField)
{
    const Selector selector = selectField(xField);
    return field_nrm2(xField,selector);
}

template<class Scalar>
inline
void field_nrm2(Scalar& glob_result, const FieldBase& xField, const Selector& selector, const MPI_Comm comm)
{
    ThrowAssert( is_compatible<Scalar>(xField) );

    BucketVector const& buckets = xField.get_mesh().get_buckets(xField.entity_rank(), selector & xField.mesh_meta_data().locally_owned_part());

    Scalar local_result(0.0);

    int orig_thread_count = fix_omp_threads();
#ifdef OPEN_MP_ACTIVE_FIELDBLAS_HPP
#pragma omp parallel for reduction(+:local_result) schedule(static)
#endif
    for(size_t i=0; i < buckets.size(); i++) {
        BucketSpan<Scalar> x(xField, *buckets[i]);
        local_result += FortranBLAS<Scalar>::dot(x.size(),x.data(),x.data());
    }

    glob_result = local_result;
    stk::all_reduce_sum(comm,&local_result,&glob_result,1u);
    glob_result = std::sqrt(glob_result);
    unfix_omp_threads(orig_thread_count);
}

template<class Scalar>
inline
void field_nrm2(std::complex<Scalar>& result, const FieldBase& xField, const Selector& selector, const MPI_Comm comm)
{
    ThrowAssert( is_compatible<std::complex<Scalar>>(xField) );

    BucketVector const& buckets = xField.get_mesh().get_buckets(xField.entity_rank(), selector & xField.mesh_meta_data().locally_owned_part());

    Scalar local_result_r(0.0);

    int orig_thread_count = fix_omp_threads();
#ifdef OPEN_MP_ACTIVE_FIELDBLAS_HPP
#pragma omp parallel for reduction(+:local_result_r) schedule(static)
#endif
    for(size_t i=0; i < buckets.size(); i++) {
        BucketSpan<std::complex<Scalar>> x(xField, *buckets[i]);
        local_result_r += std::pow(FortranBLAS<std::complex<Scalar> >::nrm2(x.size(),x.data()).real(),2);
    }

    Scalar glob_result = local_result_r;
    stk::all_reduce_sum(comm,&local_result_r,&glob_result,1u);
    result = std::complex<Scalar>(std::sqrt(glob_result),0.0);
    unfix_omp_threads(orig_thread_count);
}

template<class Scalar>
inline
void field_nrm2(Scalar& result, const FieldBase& xField, const Selector& selector)
{
    const MPI_Comm comm = xField.get_mesh().parallel();
    field_nrm2(result,xField,selector,comm);
}

template<class Scalar>
inline
void field_nrm2(Scalar& result, const FieldBase& xField)
{
    const Selector selector = selectField(xField);
    field_nrm2(result,xField,selector);
}

template<class Scalar,class T1,class T2,class T3,class T4,class T5,class T6,class T7>
inline
Scalar field_asum(const Field<Scalar,T1,T2,T3,T4,T5,T6,T7> & xField, const Selector& selector, const MPI_Comm comm)
{
    BucketVector const& buckets = xField.get_mesh().get_buckets(xField.entity_rank(), selector & xField.mesh_meta_data().locally_owned_part());

    Scalar local_result(0.0);

    int orig_thread_count = fix_omp_threads();
#ifdef OPEN_MP_ACTIVE_FIELDBLAS_HPP
#pragma omp parallel for reduction(+:local_result) schedule(static)
#endif
    for(size_t i=0; i < buckets.size(); i++) {
        BucketSpan<Scalar> x(xField, *buckets[i]);
        local_result += FortranBLAS<Scalar>::asum(x.size(),x.data());
    }

    Scalar glob_result = local_result;
    stk::all_reduce_sum(comm,&local_result,&glob_result,1u);
    unfix_omp_threads(orig_thread_count);
    return glob_result;
}

template<class Scalar,class T1,class T2,class T3,class T4,class T5,class T6,class T7>
inline
std::complex<Scalar> field_asum(const Field<std::complex<Scalar>,T1,T2,T3,T4,T5,T6,T7> & xField, const Selector& selector, const MPI_Comm comm) 
{
    BucketVector const& buckets = xField.get_mesh().get_buckets(xField.entity_rank(), selector & xField.mesh_meta_data().locally_owned_part());

    Scalar local_result(0.0);

    int orig_thread_count = fix_omp_threads();
#ifdef OPEN_MP_ACTIVE_FIELDBLAS_HPP
#pragma omp parallel for reduction(+:local_result) schedule(static)
#endif
    for(size_t i=0; i < buckets.size(); i++) {
        BucketSpan<std::complex<Scalar>> x(xField, *buckets[i]);
        local_result += FortranBLAS<std::complex<Scalar> >::asum(x.size(),x.data()).real();
    }

    Scalar glob_result = local_result;
    stk::all_reduce_sum(comm,&local_result,&glob_result,1u);
    unfix_omp_threads(orig_thread_count);
    return std::complex<Scalar>(glob_result,0.0);
}

template<class Scalar,class T1,class T2,class T3,class T4,class T5,class T6,class T7>
inline
Scalar field_asum(const Field<Scalar,T1,T2,T3,T4,T5,T6,T7> & xField, const Selector& selector)
{
    const MPI_Comm comm = xField.get_mesh().parallel();
    return field_asum(xField,selector,comm);
}

template<class Scalar,class T1,class T2,class T3,class T4,class T5,class T6,class T7>
inline
Scalar field_asum(const Field<Scalar,T1,T2,T3,T4,T5,T6,T7> & xField)
{
    const Selector selector = selectField(xField);
    return field_asum(xField,selector);
}

template<class Scalar>
inline
void field_asum(Scalar& glob_result, const FieldBase& xField, const Selector& selector, const MPI_Comm comm)
{
    ThrowAssert( is_compatible<Scalar>(xField) );

    BucketVector const& buckets = xField.get_mesh().get_buckets(xField.entity_rank(), selector & xField.mesh_meta_data().locally_owned_part());

    Scalar local_result(0.0);

    int orig_thread_count = fix_omp_threads();
#ifdef OPEN_MP_ACTIVE_FIELDBLAS_HPP
#pragma omp parallel for reduction(+:local_result) schedule(static)
#endif
    for(size_t i=0; i < buckets.size(); i++) {
        BucketSpan<Scalar> x(xField, *buckets[i]);
        local_result += FortranBLAS<Scalar>::asum(x.size(),x.data());
    }

    glob_result = local_result;
    stk::all_reduce_sum(comm,&local_result,&glob_result,1u);
    unfix_omp_threads(orig_thread_count);
}

template<class Scalar>
inline
void field_asum(std::complex<Scalar>& result, const FieldBase& xField, const Selector& selector, const MPI_Comm comm)
{
    ThrowAssert( is_compatible<std::complex<Scalar>>(xField) );

    BucketVector const& buckets = xField.get_mesh().get_buckets(xField.entity_rank(), selector & xField.mesh_meta_data().locally_owned_part());

    Scalar local_result(0.0);

    int orig_thread_count = fix_omp_threads();
#ifdef OPEN_MP_ACTIVE_FIELDBLAS_HPP
#pragma omp parallel for reduction(+:local_result) schedule(static)
#endif
    for(size_t i=0; i < buckets.size(); i++) {
        BucketSpan<std::complex<Scalar>> x(xField, *buckets[i]);
        local_result += FortranBLAS<std::complex<Scalar> >::asum(x.size(),x.data()).real();
    }

    Scalar glob_result = local_result;
    stk::all_reduce_sum(comm,&local_result,&glob_result,1u);
    result = std::complex<Scalar>(glob_result,0.0);
    unfix_omp_threads(orig_thread_count);
}

template<class Scalar>
inline
void field_asum(Scalar& result, const FieldBase& xField, const Selector& selector)
{
    const MPI_Comm comm = xField.get_mesh().parallel();
    field_asum(result,xField,selector,comm);
}

template<class Scalar>
inline
void field_asum(Scalar& result, const FieldBase& xField)
{
    const Selector selector = selectField(xField);
    field_asum(result,xField,selector);
}

template<class Scalar,class T1,class T2,class T3,class T4,class T5,class T6,class T7>
inline
Scalar field_amax(const Field<Scalar,T1,T2,T3,T4,T5,T6,T7> & xField, const Selector& selector, const MPI_Comm comm)
{
    BucketVector const& buckets = xField.get_mesh().get_buckets(xField.entity_rank(), selector & xField.mesh_meta_data().locally_owned_part());

    Scalar priv_tmp;
    Scalar local_amax(-1.0);

    int orig_thread_count = fix_omp_threads();
#ifdef OPEN_MP_ACTIVE_FIELDBLAS_HPP
#pragma omp parallel for schedule(static) reduction(max:local_amax) private(priv_tmp)
#endif
    for(size_t i=0; i < buckets.size(); i++) {
        BucketSpan<Scalar> x(xField, *buckets[i]);
        priv_tmp = std::abs(x[FortranBLAS<Scalar>::iamax(x.size(),x.data())]);

        if (local_amax < priv_tmp) {
            local_amax = priv_tmp;
        }
    }

    Scalar global_amax = local_amax;
    stk::all_reduce_max(comm,&local_amax,&global_amax,1u);
    unfix_omp_threads(orig_thread_count);
    return global_amax;
}

template<class Scalar,class T1,class T2,class T3,class T4,class T5,class T6,class T7>
inline
std::complex<Scalar> field_amax(const Field<std::complex<Scalar>,T1,T2,T3,T4,T5,T6,T7> & xField, const Selector& selector, const MPI_Comm comm) 
{
    BucketVector const& buckets = xField.get_mesh().get_buckets(xField.entity_rank(), selector & xField.mesh_meta_data().locally_owned_part());

    Scalar priv_tmp;
    Scalar local_amax(0.0);

    int orig_thread_count = fix_omp_threads();
#ifdef OPEN_MP_ACTIVE_FIELDBLAS_HPP
#pragma omp parallel for schedule(static) reduction(max:local_amax) private(priv_tmp)
#endif
    for(size_t i=0; i < buckets.size(); i++) {
        BucketSpan<std::complex<Scalar>> x(xField, *buckets[i]);
        priv_tmp = std::abs(x[FortranBLAS<std::complex<Scalar> >::iamax(x.size(),x.data())]);

        if (local_amax < priv_tmp) {
            local_amax = priv_tmp;
        }
    }

    Scalar glob_amax = local_amax;
    stk::all_reduce_max(comm,&local_amax,&glob_amax,1u);
    unfix_omp_threads(orig_thread_count);
    return std::complex<Scalar>(glob_amax,0.0);
}

template<class Scalar,class T1,class T2,class T3,class T4,class T5,class T6,class T7>
inline
Scalar field_amax(const Field<Scalar,T1,T2,T3,T4,T5,T6,T7> & xField, const Selector& selector)
{
    const MPI_Comm comm = xField.get_mesh().parallel();
    return field_amax(xField,selector,comm);
}

template<class Scalar,class T1,class T2,class T3,class T4,class T5,class T6,class T7>
inline
Scalar field_amax(const Field<Scalar,T1,T2,T3,T4,T5,T6,T7> & xField)
{
    const Selector selector = selectField(xField);
    return field_amax(xField,selector);
}

template<class Scalar>
inline
Entity INTERNAL_field_eamax_complex(const FieldBase& xField, const Selector& selector)
{
    BucketVector const& buckets = xField.get_mesh().get_buckets(xField.entity_rank(), selector & xField.mesh_meta_data().locally_owned_part());

    int priv_iamax;
    Scalar priv_amax;
    Entity priv_result;
    Scalar local_amax(-2.0);
    Entity local_result = Entity();

    int orig_thread_count = fix_omp_threads();
#ifdef OPEN_MP_ACTIVE_FIELDBLAS_HPP
#pragma omp parallel private(priv_iamax,priv_amax,priv_result)
    {
#endif
        priv_amax = Scalar(-1.0);
#ifdef OPEN_MP_ACTIVE_FIELDBLAS_HPP
#pragma omp for schedule(static) reduction(max:local_amax)
#endif
        for(size_t i=0; i < buckets.size(); i++) {
            Bucket & b = *buckets[i];
            BucketSpan<std::complex<Scalar> > x(xField, b);
            priv_iamax = FortranBLAS<std::complex<Scalar> >::iamax(x.size(),x.data());

            if (priv_amax < std::norm(x[priv_iamax])) {
                priv_result = b[priv_iamax];
                priv_amax = std::norm(x[priv_iamax]);
                local_amax = priv_amax;
            }
        }
#ifdef OPEN_MP_ACTIVE_FIELDBLAS_HPP
        if (local_amax == priv_amax) {
#endif
            local_result = priv_result;
#ifdef OPEN_MP_ACTIVE_FIELDBLAS_HPP
        }
    }
#endif

    EntityId local_EntityId = xField.get_mesh().identifier(local_result);
    EntityId  glob_EntityId = local_EntityId;
    Scalar glob_amax = local_amax;
    stk::all_reduce_maxloc(xField.get_mesh().parallel(),&local_amax,&local_EntityId,&glob_amax,&glob_EntityId,1u);
    if (glob_EntityId == local_EntityId) {
        local_result = xField.get_mesh().get_entity(xField.entity_rank(),glob_EntityId);
    } else {
        local_result = Entity();
    }
    unfix_omp_threads(orig_thread_count);
    return local_result;
}

template<class Scalar>
inline
Entity INTERNAL_field_eamax(const FieldBase& xField, const Selector& selector)
{
    BucketVector const& buckets = xField.get_mesh().get_buckets(xField.entity_rank(), selector & xField.mesh_meta_data().locally_owned_part());

    int priv_iamax;
    Scalar priv_amax;
    Entity priv_result;
    Scalar local_amax(-2.0);
    Entity local_result = Entity();

    int orig_thread_count = fix_omp_threads();
#ifdef OPEN_MP_ACTIVE_FIELDBLAS_HPP
#pragma omp parallel private(priv_iamax,priv_amax,priv_result)
    {
#endif
        priv_amax = Scalar(-1.0);
#ifdef OPEN_MP_ACTIVE_FIELDBLAS_HPP
#pragma omp for schedule(static) reduction(max:local_amax)
#endif
        for(size_t i=0; i < buckets.size(); i++) {
            Bucket & b = *buckets[i];
            BucketSpan<Scalar> x(xField, b);
            priv_iamax = FortranBLAS<Scalar>::iamax(x.size(),x.data());

            if (priv_amax < std::abs(x[priv_iamax])) {
                priv_result = b[priv_iamax];
                priv_amax = std::abs(x[priv_iamax]);
                local_amax = priv_amax;
            }
        }
#ifdef OPEN_MP_ACTIVE_FIELDBLAS_HPP
        if (local_amax==priv_amax) {
#endif
            local_result = priv_result;
#ifdef OPEN_MP_ACTIVE_FIELDBLAS_HPP
        }
    }
#endif

    EntityId local_EntityId = xField.get_mesh().identifier(local_result);
    EntityId glob_EntityId = local_EntityId;
    Scalar glob_amax = local_amax;
    stk::all_reduce_maxloc(xField.get_mesh().parallel(),&local_amax,&local_EntityId,&glob_amax,&glob_EntityId,1u);
    if (glob_EntityId==local_EntityId) {
        local_result = xField.get_mesh().get_entity(xField.entity_rank(),glob_EntityId);
    } else {
        local_result = Entity();
    }
    unfix_omp_threads(orig_thread_count);
    return local_result;
}

inline
Entity field_eamax(const FieldBase& xField, const Selector& selector)
{
    if (xField.data_traits().type_info == typeid(double)) {
        return INTERNAL_field_eamax<double>(xField,selector);
    } else if (xField.data_traits().type_info == typeid(float)) {
        return INTERNAL_field_eamax<float>(xField,selector);
    } else if (xField.data_traits().type_info == typeid(std::complex<double>)) {
        return INTERNAL_field_eamax_complex<double>(xField,selector);
    } else if (xField.data_traits().type_info == typeid(std::complex<float>)) {
        return INTERNAL_field_eamax_complex<float>(xField,selector);
    } else if (xField.data_traits().type_info == typeid(int)) {
        return INTERNAL_field_eamax<int>(xField,selector);
    } else {
        ThrowAssertMsg(false,"Error in field_eamax; field is of type "<<xField.data_traits().type_info.name()<<" which is not supported");
    }
    return stk::mesh::Entity();
}

inline
Entity field_eamax(const FieldBase& xField)
{
    const Selector selector = selectField(xField);
    return field_eamax(xField,selector);
}

template<class Scalar,class T1,class T2,class T3,class T4,class T5,class T6,class T7>
inline
Entity field_eamax(const Field<Scalar,T1,T2,T3,T4,T5,T6,T7> & xField)
{
    const Selector selector = selectField(xField);
    return field_eamax(xField,selector);
}

template<class Scalar>
inline
void field_amax(std::complex<Scalar>& result, const FieldBase& xField, const Selector& selector, const MPI_Comm comm)
{
    ThrowAssert( is_compatible<std::complex<Scalar>>(xField) );

    BucketVector const& buckets = xField.get_mesh().get_buckets(xField.entity_rank(), selector & xField.mesh_meta_data().locally_owned_part());

    Scalar priv_tmp;
    Scalar local_amax(-1.0);

    int orig_thread_count = fix_omp_threads();
#ifdef OPEN_MP_ACTIVE_FIELDBLAS_HPP
#pragma omp parallel for schedule(static) reduction(max:local_amax) private(priv_tmp)
#endif
    for(size_t i=0; i < buckets.size(); i++) {
        BucketSpan<std::complex<Scalar>> x(xField, *buckets[i]);
        priv_tmp = std::abs(x[FortranBLAS<std::complex<Scalar> >::iamax(x.size(),x.data())]);

        if (local_amax < priv_tmp) {
            local_amax = priv_tmp;
        }
    }

    Scalar glob_amax = local_amax;
    stk::all_reduce_max(comm,&local_amax,&glob_amax,1u);
    result = std::complex<Scalar>(glob_amax,0.0);
    unfix_omp_threads(orig_thread_count);
}

template<class Scalar>
inline
void field_amax(Scalar& result, const FieldBase& xField, const Selector& selector, const MPI_Comm comm)
{
    ThrowAssert( is_compatible<Scalar>(xField) );

    BucketVector const& buckets = xField.get_mesh().get_buckets(xField.entity_rank(), selector & xField.mesh_meta_data().locally_owned_part());

    Scalar priv_tmp;
    Scalar local_amax(-1.0);

    int orig_thread_count = fix_omp_threads();
#ifdef OPEN_MP_ACTIVE_FIELDBLAS_HPP
#pragma omp parallel for schedule(static) reduction(max:local_amax) private(priv_tmp)
#endif
    for(size_t i=0; i < buckets.size(); i++) {
        BucketSpan<Scalar> x(xField, *buckets[i]);
        priv_tmp = std::abs(x[FortranBLAS<Scalar>::iamax(x.size(),x.data())]);
        if (local_amax < priv_tmp) {
            local_amax = priv_tmp;
        }
    }

    Scalar glob_amax = local_amax;
    stk::all_reduce_max(comm,&local_amax,&glob_amax,1u);
    result = glob_amax;
    unfix_omp_threads(orig_thread_count);
}

template<class Scalar>
inline
void field_amax(Scalar& result, const FieldBase& xField, const Selector& selector)
{
    const MPI_Comm comm = xField.get_mesh().parallel();
    field_amax(result,xField,selector,comm);
}

template<class Scalar>
inline
void field_amax(Scalar& result, const FieldBase& xField)
{
    const Selector selector = selectField(xField);
    field_amax(result,xField,selector);
}

template<class Scalar,class T1,class T2,class T3,class T4,class T5,class T6,class T7>
inline
Entity field_eamin(const Field<std::complex<Scalar>,T1,T2,T3,T4,T5,T6,T7> & xField, const Selector& selector)
{
    BucketVector const& buckets = xField.get_mesh().get_buckets(xField.entity_rank(), selector & xField.mesh_meta_data().locally_owned_part());

    int priv_iamin;
    Scalar priv_amin;
    Entity priv_result;
    Scalar local_amin = std::numeric_limits<Scalar>::max();
    Entity local_result;
    Entity glob_result;

    int orig_thread_count = fix_omp_threads();
#ifdef OPEN_MP_ACTIVE_FIELDBLAS_HPP
#pragma omp parallel private(priv_iamin,priv_amin,priv_result)
    {
#endif
        priv_amin = std::numeric_limits<Scalar>::max();
#ifdef OPEN_MP_ACTIVE_FIELDBLAS_HPP
#pragma omp for schedule(static) reduction(min:local_amin)
#endif
        for(size_t i=0; i < buckets.size(); i++) {
            Bucket & b = *buckets[i];
            BucketSpan<std::complex<Scalar>> x(xField, b);
            priv_iamin = FortranBLAS<std::complex<Scalar> >::iamin(x.size(),x.data());

            if (std::abs(x[priv_iamin]) < priv_amin) {
                priv_result = b[priv_iamin];
                priv_amin = std::abs(x[priv_iamin]);
                local_amin = priv_amin;
            }
        }
#ifdef OPEN_MP_ACTIVE_FIELDBLAS_HPP
        if (priv_amin==local_amin) {
#endif
            local_result = priv_result;
#ifdef OPEN_MP_ACTIVE_FIELDBLAS_HPP
        }
    }
#endif

    EntityId local_EntityId = xField.get_mesh().identifier(local_result);
    EntityId glob_EntityId = local_EntityId;
    Scalar glob_amin = local_amin;
    stk::all_reduce_minloc(xField.get_mesh().parallel(),&local_amin,&local_EntityId,&glob_amin,&glob_EntityId,1u);
    if (glob_EntityId == local_EntityId) {
        glob_result = xField.get_mesh().get_entity(xField.entity_rank(),glob_EntityId);
    } else {
        glob_result = Entity();
    }
    unfix_omp_threads(orig_thread_count);
    return glob_result;
}

template<class Scalar,class T1,class T2,class T3,class T4,class T5,class T6,class T7>
inline
Entity field_eamin(const Field<Scalar,T1,T2,T3,T4,T5,T6,T7> & xField, const Selector& selector)
{
    BucketVector const& buckets = xField.get_mesh().get_buckets(xField.entity_rank(), selector & xField.mesh_meta_data().locally_owned_part());

    int priv_iamin;
    Scalar priv_amin;
    Entity priv_result;
    Scalar local_amin = std::numeric_limits<Scalar>::max();
    Entity local_result;
    Entity glob_result;

    int orig_thread_count = fix_omp_threads();
#ifdef OPEN_MP_ACTIVE_FIELDBLAS_HPP
#pragma omp parallel private(priv_iamin,priv_amin,priv_result)
    {
#endif
        priv_amin = std::numeric_limits<Scalar>::max();
#ifdef OPEN_MP_ACTIVE_FIELDBLAS_HPP
#pragma omp for schedule(static) reduction(min:local_amin)
#endif
        for(size_t i=0; i < buckets.size(); i++) {
            Bucket & b = *buckets[i];
            BucketSpan<Scalar> x(xField, b);
            priv_iamin = FortranBLAS<Scalar>::iamin(x.size(),x.data());

            if (std::abs(x[priv_iamin]) < priv_amin) {
                priv_result = b[priv_iamin];
                priv_amin = std::abs(x[priv_iamin]);
                local_amin = priv_amin;
            }
        }
#ifdef OPEN_MP_ACTIVE_FIELDBLAS_HPP
        if (priv_amin == local_amin) {
#endif
            local_result = priv_result;
#ifdef OPEN_MP_ACTIVE_FIELDBLAS_HPP
        }
    }
#endif

    EntityId local_EntityId = xField.get_mesh().identifier(local_result);
    EntityId glob_EntityId = local_EntityId;
    Scalar glob_amin = local_amin;
    stk::all_reduce_minloc(xField.get_mesh().parallel(),&local_amin,&local_EntityId,&glob_amin,&glob_EntityId,1u);
    if (glob_EntityId == local_EntityId) {
        glob_result = xField.get_mesh().get_entity(xField.entity_rank(),glob_EntityId);
    } else {
        glob_result = Entity();
    }
    unfix_omp_threads(orig_thread_count);
    return glob_result;
}

template<class Scalar,class T1,class T2,class T3,class T4,class T5,class T6,class T7>
inline
Entity field_eamin(const Field<Scalar,T1,T2,T3,T4,T5,T6,T7> & xField)
{
    const Selector selector = selectField(xField);
    return field_eamin(xField,selector);
}

template<class Scalar,class T1,class T2,class T3,class T4,class T5,class T6,class T7>
inline
std::complex<Scalar> field_amin(const Field<std::complex<Scalar>,T1,T2,T3,T4,T5,T6,T7> & xField, const Selector& selector, const MPI_Comm comm) 
{
    BucketVector const& buckets = xField.get_mesh().get_buckets(xField.entity_rank(), selector & xField.mesh_meta_data().locally_owned_part());

    Scalar priv_tmp;
    Scalar local_amin = std::numeric_limits<Scalar>::max();

    int orig_thread_count = fix_omp_threads();
#ifdef OPEN_MP_ACTIVE_FIELDBLAS_HPP
#pragma omp parallel for schedule(static) reduction(min:local_amin) private(priv_tmp)
#endif
    for(size_t i=0; i < buckets.size(); i++) {
        BucketSpan<std::complex<Scalar>> x(xField, *buckets[i]);
        priv_tmp = std::norm(x[FortranBLAS<std::complex<Scalar> >::iamin(x.size(),x.data())]);

        if (priv_tmp < local_amin) {
            local_amin = priv_tmp;
        }
    }

    Scalar glob_amin = local_amin;
    stk::all_reduce_min(comm,&local_amin,&glob_amin,1u);
    unfix_omp_threads(orig_thread_count);
    return std::sqrt(glob_amin);
}

template<class Scalar,class T1,class T2,class T3,class T4,class T5,class T6,class T7>
inline
Scalar field_amin(const Field<Scalar,T1,T2,T3,T4,T5,T6,T7> & xField, const Selector& selector, const MPI_Comm comm)
{
    BucketVector const& buckets = xField.get_mesh().get_buckets(xField.entity_rank(), selector & xField.mesh_meta_data().locally_owned_part());

    Scalar priv_tmp;
    Scalar local_amin = std::numeric_limits<Scalar>::max();

    int orig_thread_count = fix_omp_threads();
#ifdef OPEN_MP_ACTIVE_FIELDBLAS_HPP
#pragma omp parallel for schedule(static) reduction(min:local_amin) private(priv_tmp)
#endif
    for(size_t i=0; i < buckets.size(); i++) {
        BucketSpan<Scalar> x(xField, *buckets[i]);
        priv_tmp = std::abs(x[FortranBLAS<Scalar>::iamin(x.size(),x.data())]);

        if (priv_tmp < local_amin) {
            local_amin = priv_tmp;
        }
    }

    Scalar glob_amin = local_amin;
    stk::all_reduce_min(comm,&local_amin,&glob_amin,1u);
    unfix_omp_threads(orig_thread_count);
    return glob_amin;
}

template<class Scalar,class T1,class T2,class T3,class T4,class T5,class T6,class T7>
inline
Scalar field_amin(const Field<Scalar,T1,T2,T3,T4,T5,T6,T7> & xField, const Selector& selector)
{
    const MPI_Comm comm = xField.get_mesh().parallel();
    return field_amin(xField,selector,comm);
}

template<class Scalar,class T1,class T2,class T3,class T4,class T5,class T6,class T7>
inline
Scalar field_amin(const Field<Scalar,T1,T2,T3,T4,T5,T6,T7> & xField)
{
    const Selector selector = selectField(xField);
    return field_amin(xField,selector);
}

template<class Scalar>
inline
Entity INTERNAL_field_eamin_complex(const FieldBase& xField, const Selector& selector) 
{
    BucketVector const& buckets = xField.get_mesh().get_buckets(xField.entity_rank(), selector & xField.mesh_meta_data().locally_owned_part());

    int priv_iamin;
    double priv_amin;
    Entity priv_result;
    double local_amin = std::numeric_limits<double>::max();
    Entity local_result;
    Entity glob_result;

    int orig_thread_count = fix_omp_threads();
#ifdef OPEN_MP_ACTIVE_FIELDBLAS_HPP
#pragma omp parallel private(priv_iamin,priv_amin,priv_result)
    {
#endif
        priv_amin = std::numeric_limits<double>::max();
#ifdef OPEN_MP_ACTIVE_FIELDBLAS_HPP
#pragma omp for schedule(static) reduction(min:local_amin)
#endif
        for(size_t i=0; i < buckets.size(); i++) {
            Bucket & b = *buckets[i];
            BucketSpan<std::complex<Scalar>> x(xField, b);
            priv_iamin = FortranBLAS<std::complex<Scalar> >::iamin(x.size(),x.data());

            if (std::norm(x[priv_iamin]) < priv_amin) {
                priv_result = b[priv_iamin];
                priv_amin = std::norm(x[priv_iamin]);
                local_amin = priv_amin;
            }
        }
#ifdef OPEN_MP_ACTIVE_FIELDBLAS_HPP
        if (priv_amin == local_amin) {
#endif
            local_result = priv_result;
#ifdef OPEN_MP_ACTIVE_FIELDBLAS_HPP
        }
    }
#endif

    EntityId local_EntityId = xField.get_mesh().identifier(local_result);
    EntityId glob_EntityId = local_EntityId;
    double glob_amin = local_amin;
    stk::all_reduce_minloc(xField.get_mesh().parallel(),&local_amin,&local_EntityId,&glob_amin,&glob_EntityId,1u);
    if (glob_EntityId == local_EntityId) {
        glob_result = xField.get_mesh().get_entity(xField.entity_rank(),glob_EntityId);
    } else {
        glob_result = Entity();
    }
    unfix_omp_threads(orig_thread_count);
    return glob_result;
}

template<class Scalar>
inline
Entity INTERNAL_field_eamin(const FieldBase& xField, const Selector& selector) 
{
    BucketVector const& buckets = xField.get_mesh().get_buckets(xField.entity_rank(), selector & xField.mesh_meta_data().locally_owned_part());

    int priv_iamin;
    double priv_amin;
    Entity priv_result;
    double local_amin = std::numeric_limits<double>::max();
    Entity local_result;
    Entity glob_result;

    int orig_thread_count = fix_omp_threads();
#ifdef OPEN_MP_ACTIVE_FIELDBLAS_HPP
#pragma omp parallel private(priv_iamin,priv_amin,priv_result)
    {
#endif
        priv_amin = std::numeric_limits<double>::max();
#ifdef OPEN_MP_ACTIVE_FIELDBLAS_HPP
#pragma omp for schedule(static) reduction(min:local_amin)
#endif
        for(size_t i=0; i < buckets.size(); i++) {
            Bucket & b = *buckets[i];
            BucketSpan<Scalar> x(xField, b);
            priv_iamin = FortranBLAS<Scalar>::iamin(x.size(),x.data());

            if ( std::abs(x[priv_iamin]) < priv_amin) {
                priv_result = b[priv_iamin];
                priv_amin = std::abs(x[priv_iamin]);
                local_amin = priv_amin;
            }
        }
#ifdef OPEN_MP_ACTIVE_FIELDBLAS_HPP
        if (priv_amin==local_amin) {
#endif
            local_result = priv_result;
#ifdef OPEN_MP_ACTIVE_FIELDBLAS_HPP
        }
    }
#endif

    EntityId local_EntityId = xField.get_mesh().identifier(local_result);
    EntityId glob_EntityId = local_EntityId;
    double glob_amin = local_amin;
    stk::all_reduce_minloc(xField.get_mesh().parallel(),&local_amin,&local_EntityId,&glob_amin,&glob_EntityId,1u);
   if (glob_EntityId == local_EntityId) {
        glob_result = xField.get_mesh().get_entity(xField.entity_rank(),glob_EntityId);
    } else {
        glob_result = Entity();
    }
    unfix_omp_threads(orig_thread_count);
    return glob_result;
}

inline
Entity field_eamin(const FieldBase& xField, const Selector& selector)
{
    if (xField.data_traits().type_info == typeid(double)) {
        return INTERNAL_field_eamin<double>(xField,selector);
    } else if (xField.data_traits().type_info == typeid(float)) {
        return INTERNAL_field_eamin<float>(xField,selector);
    } else if (xField.data_traits().type_info == typeid(std::complex<double>)) {
        return INTERNAL_field_eamin_complex<double>(xField,selector);
    } else if (xField.data_traits().type_info == typeid(std::complex<float>)) {
        return INTERNAL_field_eamin_complex<float>(xField,selector);
    } else if (xField.data_traits().type_info == typeid(int)) {
        return INTERNAL_field_eamin<int>(xField,selector);
    } else {
        ThrowAssertMsg(false,"Error in field_eamin; field is of type "<<xField.data_traits().type_info.name()<<" which is not supported");
    }
    return stk::mesh::Entity();
}

inline
Entity field_eamin(const FieldBase& xField)
{
    const Selector selector = selectField(xField);
    return field_eamin(xField,selector);
}

template<class Scalar>
inline
void field_amin(std::complex<Scalar>& result, const FieldBase& xField, const Selector& selector, const MPI_Comm comm)
{
    ThrowAssert( is_compatible<std::complex<Scalar>>(xField) );

    BucketVector const& buckets = xField.get_mesh().get_buckets(xField.entity_rank(), selector & xField.mesh_meta_data().locally_owned_part());

    Scalar priv_tmp;
    Scalar local_amin = std::numeric_limits<Scalar>::max();

    int orig_thread_count = fix_omp_threads();
#ifdef OPEN_MP_ACTIVE_FIELDBLAS_HPP
#pragma omp parallel for schedule(static) reduction(min:local_amin) private(priv_tmp)
#endif
    for(size_t i=0; i < buckets.size(); i++) {
        BucketSpan<std::complex<Scalar>> x(xField, *buckets[i]);
        priv_tmp = std::abs(x[FortranBLAS<std::complex<Scalar> >::iamin(x.size(),x.data())]);

        if (priv_tmp < local_amin) {
            local_amin = priv_tmp;
        }
    }

    Scalar glob_amin = local_amin;
    stk::all_reduce_min(comm,&local_amin,&glob_amin,1u);
    result = std::complex<Scalar>(glob_amin,0.0);
    unfix_omp_threads(orig_thread_count);
}

template<class Scalar>
inline
void field_amin(Scalar& result, const FieldBase& xField, const Selector& selector, const MPI_Comm comm)
{
    ThrowAssert( is_compatible<Scalar>(xField) );

    BucketVector const& buckets = xField.get_mesh().get_buckets(xField.entity_rank(), selector & xField.mesh_meta_data().locally_owned_part());

    Scalar priv_tmp;
    Scalar local_amin = std::numeric_limits<Scalar>::max();

    int orig_thread_count = fix_omp_threads();
#ifdef OPEN_MP_ACTIVE_FIELDBLAS_HPP
#pragma omp parallel for schedule(static) reduction(min:local_amin) private(priv_tmp)
#endif
    for(size_t i=0; i < buckets.size(); i++) {
        Bucket & b = *buckets[i];
        BucketSpan<Scalar> x(xField, b);
        priv_tmp = std::abs(x[FortranBLAS<Scalar>::iamin(x.size(),x.data())]);

        if (priv_tmp < local_amin) {
            local_amin = priv_tmp;
        }
    }
    Scalar glob_amin = local_amin;
    stk::all_reduce_min(comm,&local_amin,&glob_amin,1u);
    result = glob_amin;
    unfix_omp_threads(orig_thread_count);
}

template<class Scalar>
inline
void field_amin(Scalar& result, const FieldBase& xField, const Selector& selector)
{
    const MPI_Comm comm = xField.get_mesh().parallel();
    field_amin(result,xField,selector,comm);
}

template<class Scalar>
inline
void field_amin(Scalar& result, const FieldBase& xField)
{
    const Selector selector = selectField(xField);
    field_amin(result,xField,selector);
}

} // mesh
} // stk

#endif // STK_MESH_BASE_FIELDBLAS_HPP

