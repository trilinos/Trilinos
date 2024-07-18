// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// -*- c-file-style: "krlong" -*-

/** Same as sample.thyra-default.cpp except that the "krlong" style is used
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
      TEUCHOS_TEST_FOR_EXCEPT("Should never get here!");
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
