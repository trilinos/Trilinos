
/* Copyright (2001) Sandia Corportation. Under the terms of Contract 
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this 
 * work by or on behalf of the U.S. Government.  Export of this program
 * may require a license from the United States Government. */


/* NOTICE:  The United States Government is granted for itself and others
 * acting on its behalf a paid-up, nonexclusive, irrevocable worldwide
 * license in ths data to reproduce, prepare derivative works, and
 * perform publicly and display publicly.  Beginning five (5) years from
 * July 25, 2001, the United States Government is granted for itself and
 * others acting on its behalf a paid-up, nonexclusive, irrevocable
 * worldwide license in this data to reproduce, prepare derivative works,
 * distribute copies to the public, perform publicly and display
 * publicly, and to permit others to do so.
 * 
 * NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED STATES DEPARTMENT
 * OF ENERGY, NOR SANDIA CORPORATION, NOR ANY OF THEIR EMPLOYEES, MAKES
 * ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR
 * RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY
 * INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS
 * THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS. */

#ifndef _EPETRA_OBJECT_H_
#define _EPETRA_OBJECT_H_


#ifndef __cplusplus
#define __cplusplus
#endif

#if defined(SGI) || defined(SGI64) || defined(SGI32) || defined(CPLANT)

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <iostream>
#include <math.h>
#include <string>
using namespace std;

#elif defined(TFLOP)

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string>
using std::string;
#include <iostream>
#include <iomanip>
using std::istream;
using std::ostream;
using std::cerr;
using std::cout;
using std::endl;

#else

#include <cstdlib>
#include <cstdio>
#include <cassert>
#include <iostream>
#include <cmath>
#include <string>
using namespace std;

#endif



#ifdef EPETRA_SIMULATE_BOOL
#ifdef bool
#undef bool
#endif
#ifdef true
#undef true
#endif
#ifdef false
#undef false
#endif

#define bool int
#define true 1
#define false 0

#endif

#define EPETRA_MAX(x,y) (( (x) > (y) ) ? x : y)     /* max function  */
#define EPETRA_MIN(x,y) (( (x) < (y) ) ? x : y)     /* min function  */
#define EPETRA_SGN(x) (((x) < 0.0) ? -1.0 : 1.0)  /* sign function */


const double Epetra_MinDouble = 1.0E-100;
const double Epetra_MaxDouble = 1.0E+100;
const double Epetra_Overflow = 1.79E308; // Used to test if equilibration should be done.
const double Epetra_Underflow = 2.23E-308;

// RAB: 2002/1/25: Define EPETRA_ANSI_CPP as an argument to the compiler!
//#undef EPETRA_ANSI_CPP // Do not use ANSI/ISO C++ (curently just checked for I/O functions)
#ifdef EPETRA_ANSI_CPP
typedef std::ios_base::fmtflags   Epetra_fmtflags;
#else
typedef long int                  Epetra_fmtflags;
#endif

const bool Epetra_FormatStdout = true; // Set true if the ostream << operator should format output

// Delete any previous definition of EPETRA_NO_ERROR_REPORTS

#ifdef EPETRA_CHK_ERR
#undef EPETRA_CHK_ERR
#endif
#ifdef EPETRA_CHK_PTR
#undef EPETRA_CHK_PTR
#endif
#ifdef EPETRA_CHK_REF
#undef EPETRA_CHK_REF
#endif

// Make error report silent by defining EPETRA_NO_ERROR_REPORTS

#ifdef EPETRA_NO_ERROR_REPORTS
#define EPETRA_CHK_ERR(a) { int epetra_err = a; if (epetra_err != 0)  return(epetra_err);}
#define EPETRA_CHK_PTR(a) { return(a);}
#define EPETRA_CHK_REF(a) { return(a);}
#else


// Great little macro obtained from Alan Williams (modified for dynamic switching on/off)

const int DefaultTracebackMode = 1; // Default value for traceback behavior

#define EPETRA_CHK_ERR(a) { { int epetra_err = a; \
                              if ((epetra_err < 0 && Epetra_Object::GetTracebackMode() > 0) || \
                                  (epetra_err > 0 && Epetra_Object::GetTracebackMode() > 1)) { \
                      cerr << "Epetra ERROR " << epetra_err << ", " \
                           << __FILE__ << ", line " << __LINE__ << endl; }\
                      if (epetra_err !=0) return(epetra_err);  }\
                   }

// Extension of same macro for pointer, returns zero if bad

#define EPETRA_CHK_PTR(a) { if (a == 0 && Epetra_Object::GetTracebackMode() > 0) { \
                      cerr << "Epetra returning zero pointer " << ", " \
                           << __FILE__ << ", line " << __LINE__ << endl; } \
                      return(a); \
                   }
// Extension of same macro for reference, returns a default reference

#define EPETRA_CHK_REF(a) { if (Epetra_Object::GetTracebackMode() > 0) {\
                      cerr << "Epetra returning default reference " << ", " \
                           << __FILE__ << ", line " << __LINE__ << endl; } \
                      return(a); \
                   }
#endif

#include "Epetra_CombineMode.h"
#include "Epetra_DataAccess.h"


//! Epetra_Object:  The base Epetra class.
/*! The Epetra_Object class provides capabilities common to all Epetra objects,
    such as a label that identifies an object instance, constant definitions,
    enum types.
  
*/
class Epetra_Object {
    
  public:
  //@{ \name Constructors/destructor.
  //! Epetra_Object Constructor.
  /*! Epetra_Object is the primary base class in Epetra.  All Epetra class
      are derived from it, directly or indirectly.  This class is seldom
      used explictly.
  */
  Epetra_Object(int TracebackModeIn = -1);

  //! Epetra_Object Constructor.
  /*! Creates a Epetra_Object with the given label.
  */
  Epetra_Object(const char * const Label, int TracebackModeIn = -1);

  //! Epetra_Object Copy Constructor.
  /*! Makes an exact copy of an existing Epetra_Object instance.
  */
  Epetra_Object(const Epetra_Object& Object);

  //! Epetra_Object Destructor.
  /*! Completely deletes a Epetra_Object object.  
  */
  virtual ~Epetra_Object();
  //@}
  
  //@{ \name Attribute set/get methods.

  //! Epetra_Object Label definition using char *.
  /*! Defines the label used to describe the \e this object.  
  */
  virtual void SetLabel(const char * const Label);

  //! Epetra_Object Label access funtion.
  /*! Returns the string used to define this object.  
  */
  virtual char * Label() const;

  //! Set the value of the Epetra_Object error traceback report mode.
  /*! Sets the integer error traceback behavior.  
      TracebackMode controls whether or not traceback information is printed when run time 
      integer errors are detected:

      <= 0 - No information report

       = 1 - Fatal (negative) values are reported

      >= 2 - All values (except zero) reported.

      Default is set to 1.
  */
  static void SetTracebackMode(int TracebackModeValue);

  //! Get the value of the Epetra_Object error report mode.
  static int GetTracebackMode();
  //@}

  //@{ \name Miscellaneous

  //! Print object to an output stream
  //! Print method
  virtual void Print(ostream & os) const;

  //! Error reporting method
  virtual int ReportError(const string Message, int ErrorCode) const;
  //@}

  
// TracebackMode controls how much traceback information is printed when run time 
// integer errors are detected:
// = 0 - No information report
// = 1 - Fatal (negative) values are reported
// = 2 - All values (except zero) reported.

// Default is set to 2.  Can be set to different value using SetTracebackMode() method in
// Epetra_Object class
  static int TracebackMode;


 protected:
  string toString(const int& x) const {
     char s[100];
     std::sprintf(s, "%d", x);
     return string(s);
}

  string toString(const double& x) const {
     char s[100];
     std::sprintf(s, "%g", x);
     return string(s);
}
  

 private:

  char * Label_;

};
inline ostream& operator<<(ostream& os, const Epetra_Object& obj)
{
  if (Epetra_FormatStdout) {
    const Epetra_fmtflags  olda = os.setf(ios::right,ios::adjustfield);
    const Epetra_fmtflags  oldf = os.setf(ios::scientific,ios::floatfield);
    const int              oldp = os.precision(12);

    os << obj.Label() << endl;
    obj.Print(os);

    os.setf(olda,ios::adjustfield);
    os.setf(oldf,ios::floatfield);
    os.precision(oldp);
  }
  else {

    os << obj.Label();
    obj.Print(os);
  }
  
  return os;
}


#endif /* _EPETRA_OBJECT_H_ */
