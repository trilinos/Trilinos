/*
// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER
*/

// -*- c-file-style: "thyra" -*-

/** Same as sample.thyra-default.cpp except that the "thyra" style is used
 * explicitly!
 */

namespace NamespaceA {

void func1( int a, int b, int c,
  int d, int e, int f,
  int g, int h, int i
  );


void func2(
  int a, int b, int c,
  int d, int e, int f,
  int g, int h, int i
  );


} // namespace NamespaceA

void NamespaceA::func1( int a, int b, int c,
  int d, int e, int f,
  int g, int h, int i
  )
{
  
  double aa, bb, cc,
    dd;
  
  {
    std::vector<double> va(a);

    for ( int i = 0; i < a; ++i ) {
      if ( i*a < b ) {
        va[i] = 2.0;
      }
      else if ( i*b < c ) {
        va[i] = 2.5;
      }
      else {
        va[i] = 3.0;
      }
    }

    for ( int i = 0; i < a; ++i )
    {
      if ( i*a < b )
      {
        va[i] = 2.0;
      }
      else if ( i*b < c )
      {
        va[i] = 2.5;
      }
      else
      {
        va[i] = 3.0;
      }
    }

    switch(d) {
      case 0:
        aa = 4.0;
        break;
      case 1:
        aa = 5.0;
        break;
      case 2:
        aa = 6.0;
        break;
      default:
        TEST_FOR_EXCEPT("Should never get here!");
    }

    if(
      a < b
      && c > d
      && f < g
      )
    {
      bb = 8.0;
    }
    else if( h < i ) {
      bb = 9.0;
    }
    else
    {
      cc = 10.0;
    }
    
  }
  
}


void NamespaceA::func2(
  int a, int b, int c,
  int d, int e, int f,
  int g, int h, int i
  )
{

  func1( a, b, c, d, e,
    f, g, h, i );

  func2(
    a, b, c, d, e,
    f, g, h, i
    );

}
