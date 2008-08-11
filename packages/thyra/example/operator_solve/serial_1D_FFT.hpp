// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef THYRA_SERIAL_1D_FFT_HPP
#define THYRA_SERIAL_1D_FFT_HPP

#include "Thyra_LinearOpWithSolveBaseDecl.hpp"
#include "Thyra_SingleRhsLinearOpWithSolveBaseDecl.hpp"

/** \brief Simple templated function that implements a serial 1D
 * complex-to-complex FFT.
 *
 * Replaces data[1..2*nn] by its discrete Fourier transform, if isign is input
 * as 1; or replaces data[1..2*nn] by nn times its inverse discrete Fourier
 * transform, if isign is input as -1.  data is a complex array of length nn
 * or, equivalently, a real array of length 2*nn. nn MUST be an integer power
 * of 2 (this is not checked for!).
 */
template<class RealScalar>
void serial_1D_FFT( RealScalar data[], unsigned long nn, int isign )
{
  unsigned long n,mmax,m,j,istep,i;
  RealScalar wtemp,wr,wpr,wpi,wi,theta;
  RealScalar tempr,tempi;
  n=nn << 1;
  j=1;
  for (i=1;i<n;i+=2) { // This is the bit-reversal section of the routine.
    if (j > i) {
      std::swap(data[j],data[i]);
      std::swap(data[j+1],data[i+1]);
    }
    m=nn;
    while (m >= 2 && j > m) {
      j -= m;
      m >>= 1;
    }
    j += m;
  }
  // Here begins the Danielson-Lanczos section of the routine.
  mmax=2;
  while (n > mmax) { // Outer loop executed log2(nn) times.
    istep=mmax << 1;
    theta=isign*(6.28318530717959/mmax); // Initialize the trigonometric recurrence.
    wtemp=std::sin(0.5*theta);
    wpr = -2.0*wtemp*wtemp;
    wpi=std::sin(theta);
    wr=1.0;
    wi=0.0;
    for (m=1;m<mmax;m+=2) {
      for (i=m;i<=n;i+=istep) {
        j=i+mmax;
        tempr=wr*data[j]-wi*data[j+1];
        tempi=wr*data[j+1]+wi*data[j];
        data[j]=data[i]-tempr;
        data[j+1]=data[i+1]-tempi;
        data[i] += tempr;
        data[i+1] += tempi;
      }
      wr=(wtemp=wr)*wpr-wi*wpi+wr; // Trigonometric recurrence.
      wi=wi*wpr+wtemp*wpi+wi;
    }
    mmax=istep;
  }
}

#endif	// THYRA_SERIAL_1D_FFT_HPP
