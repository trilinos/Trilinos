// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/**
 * This file shows examples of indentatation of the "thyra" emacs style that
 *  is ment to follow the style guidelines in "Code Complete", 2nd edition.
 */


//
// Prototypes
//


//
// Don't indent for namespaces enclosures
//

namespace NamespaceA {


// Indent arguments on continuation lines one offset from from the beginning
// of the function type+name line.  This style does deviate from the
// recommended style in "Code Complete" in that the paren lines up with the
// argument list, not the type+name line.
void func1( int a, int b, int c,
  int d, int e, int f,
  int g, int h, int i
  );


// Same as above, except the the argument list starts on the line below the
// opening paren.  This case can be handled differently in emacs.
void func2(
  int a, int b, int c,
  int d, int e, int f,
  int g, int h, int i
  );


} // namespace NamespaceA


//
// Defintions:
//


// The following function definitions shows a few things:
// 1) The defintions are indented two spaces from other entities
// 2) The function begin '{' and end '}' are both indented from
//    the rest of the code by one space.  This sets off the
//    boundaries for the function.

void NamespaceA::func1( int a, int b, int c,
  int d, int e, int f,
  int g, int h, int i
  )
{
  
  // Indent continuation lines on variable declarations one offset.
  double aa, bb, cc,
    dd;
  
  {
    std::vector<double> va(a);

    // Use "pure block emulation" for one-line control statements

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

    // Uses "unindented begin-end pairs" (not recommended, but see below).

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

    // Above, not that (x)emacs shows the match for the opening '{' plus the
    // line above it when the '{' is not in the screen!
    
    // Indent case labels within switch statements

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

    // For control statements that extend over one line, use "unindented
    // begin-end pairs".  This breaks with the advise of "Code Complete", 2nd
    // edition but this is more standard within the C++ community than the
    // recommended indented "begin/end pairs".  To be consistent with the
    // initial 'if' statement block within the same if/else if/else structure,
    // I put the '{' on the next line from the 'else' statement.

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


// Indented two spaces from above end of function '}'.
void NamespaceA::func2(
  int a, int b, int c,
  int d, int e, int f,
  int g, int h, int i
  )
{

  // The function arguments on continuation lines in a function call are
  // indented one offset instead of aligning them with the opening '('.
  func1( a, b, c, d, e,
    f, g, h, i );

  // Same as above, except that the arguments start on the next line from the
  // opening '('.  This is handled differently in emacs.  Also, note that the
  // closing ')' is aligned with the arguments and not the 'func2(' beginning
  // 
  func2(
    a, b, c, d, e,
    f, g, h, i
    );

}

