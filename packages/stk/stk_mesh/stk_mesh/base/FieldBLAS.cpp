#include "FieldBLAS.hpp"
#include <stk_util/util/BlasLapack.hpp>

namespace stk {
namespace mesh {

//-------------------------------------------------------------------------
// double instantiation
//
void FortranBLAS<double>::axpy(int kmax, const double alpha, const double x[], double y[])
{
    const int one = 1;
    SIERRA_FORTRAN(daxpy)(&kmax,&alpha,x,&one,y,&one);
}

void FortranBLAS<double>::axpby(int kmax, const double alpha, const double x[], const double beta, double y[])
{
    for(int k = 0; k < kmax; ++k) {
      y[k] *= beta;
    }
    for(int k = 0; k < kmax; ++k) {
        y[k] += alpha * x[k];
    }
}

void FortranBLAS<double>::product(int kmax, const double x[], const double y[], double z[])
{
    for (int k = 0; k < kmax; ++k)
    {
        z[k] = x[k]*y[k];
    }
}

void FortranBLAS<double>::copy(int kmax, const double x[], double y[])
{
    for (int k = 0; k < kmax; ++k) {
        y[k] = x[k];
    }
}

double FortranBLAS<double>::dot(int kmax, const double x[], const double y[])
{
    const int one = 1;
    return SIERRA_FORTRAN(ddot)(&kmax,x,&one,y,&one);
}

double FortranBLAS<double>::nrm2(int kmax, const double x[])
{
    return std::sqrt(dot(kmax, x, x));
}

void FortranBLAS<double>::scal(int kmax, const double alpha, double x[])
{
    const int one = 1;
    SIERRA_FORTRAN(dscal)(&kmax,&alpha,x,&one);
}

void FortranBLAS<double>::fill(int numVals, double alpha, double x[], int stride)
{
    if (stride == 1) {
        std::fill(x, x+numVals, alpha);
    }
    else {
        for (double * end = x+(numVals*stride); x < end; x+=stride) {
            *x = alpha;
        }
    }
}

void FortranBLAS<double>::swap(int kmax, double x[], double y[])
{
    const int one = 1;
    SIERRA_FORTRAN(dswap)(&kmax,x,&one,y,&one);
}

double FortranBLAS<double>::asum(int kmax, const double x[])
{
    const int one = 1;
    return SIERRA_FORTRAN(dasum)(&kmax,x,&one);
}

int FortranBLAS<double>::iamax(int kmax, const double x[])
{
    const int one = 1;
    return (SIERRA_FORTRAN(idamax)(&kmax, x, &one) - 1);
}

int FortranBLAS<double>::iamin(int kmax, const double x[])
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

//------------------------------------------------------------------------
// float instantiation

void FortranBLAS<float>::axpy(int kmax, const float alpha, const float x[], float y[])
{
    const int one = 1;
    SIERRA_FORTRAN(saxpy)(&kmax,&alpha,x,&one,y,&one);
}

void FortranBLAS<float>::axpby(int kmax, const float alpha, const float x[], const float& beta, float y[])
{
    for(int k = 0; k < kmax; ++k) {
      y[k] *= beta;
    }
    for(int k = 0; k < kmax; ++k) {
        y[k] += alpha * x[k];
    }
}

void FortranBLAS<float>::product(int kmax, const float x[], const float y[], float z[])
{
    for (int k = 0; k < kmax; ++k) {
        z[k] = x[k]*y[k];
    }
}

void FortranBLAS<float>::copy(int kmax, const float x[], float y[])
{
    for (int k = 0; k < kmax; ++k) {
        y[k] = x[k];
    }
}

float FortranBLAS<float>::dot(int kmax, const float x[], const float y[])
{
    const int one = 1;
    return static_cast<float>(SIERRA_FORTRAN(sdot)(&kmax,x,&one,y,&one));
}

float FortranBLAS<float>::nrm2(int kmax, const float x[])
{
    return std::sqrt(dot(kmax, x, x));
}

void FortranBLAS<float>::scal(int kmax, const float alpha, float x[])
{
    const int one = 1;
    SIERRA_FORTRAN(sscal)(&kmax,&alpha,x,&one);
}

void FortranBLAS<float>::fill(int numVals, float alpha, float x[], int stride)
{
    if (stride == 1) {
        std::fill(x, x+numVals, alpha);
    }
    else {
        for (float * end = x+(numVals*stride); x < end; x+=stride) {
            *x = alpha;
        }
    }
}

void FortranBLAS<float>::swap(int kmax, float x[], float y[])
{
    const int one = 1;
    SIERRA_FORTRAN(sswap)(&kmax,x,&one,y,&one);
}

float FortranBLAS<float>::asum(int kmax, const float x[])
{
    const int one = 1;
    return static_cast<float>(SIERRA_FORTRAN(sasum)(&kmax,x,&one));
}

int FortranBLAS<float>::iamax(int kmax, const float x[])
{
    const int one = 1;
    return (SIERRA_FORTRAN(isamax)(&kmax, x, &one) - 1);
}

int FortranBLAS<float>::iamin(int kmax, const float x[])
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


}
}

