
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

#if defined(SGI) || defined(SGI64) || defined(SGI32)

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <iostream>
#include <math.h>
#include <string>
using namespace std;
#define MPI_NO_CPPBIND

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


// Great little macro obtained from Alan Williams

#define EPETRA_CHK_ERR(a) { { int epetra_err = a; if (epetra_err != 0) { \
                      cerr << "Epetra ERROR " << epetra_err << ", " \
                           << __FILE__ << ", line " << __LINE__ << endl; \
                      return(epetra_err);  } }\
                   }

// Extension of same macro for pointer, returns zero if bad

#define EPETRA_CHK_PTR(a) { if (a == 0) { \
                      cerr << "Epetra returning zero pointer " << ", " \
                           << __FILE__ << ", line " << __LINE__ << endl; } \
                      return(a); \
                   }
// Extension of same macro for reference, returns a default reference

#define EPETRA_CHK_REF(a) { \
                      cerr << "Epetra returning default reference " << ", " \
                           << __FILE__ << ", line " << __LINE__ << endl; \
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
  //! Epetra_Object Constructor.
  /*! Creates a Epetra_Object instance.
  */
  Epetra_Object();

  //! Epetra_Object Constructor.
  /*! Creates a Epetra_Object with the given label.
  */
  Epetra_Object(const char * const Label);

  //! Epetra_Object Copy Constructor.
  /*! Makes an exact copy of an existing Epetra_Object instance.
  */
  Epetra_Object(const Epetra_Object& Object);


  //! Epetra_Object Label definition using char *.
  /*! Defines the label used to describe the \e this object.  
  */
  virtual void SetLabel(const char * const Label);

  //! Epetra_Object Label access funtion.
  /*! Returns the string used to define this object.  
  */
  virtual char * Label() const;

  //! Print method
  virtual void Print(ostream & os) const;

  //! Error reporting method
  virtual int ReportError(const string Message, int ErrorCode) const;


  //! Epetra_Object Destructor.
  /*! Completely deletes a Epetra_Object object.  
  */
  virtual ~Epetra_Object();

 protected:
  string toString(const int& x) const {
     char s[100];
     sprintf(s, "%d", x);
     return string(s);
}

  string toString(const double& x) const {
     char s[100];
     sprintf(s, "%g", x);
     return string(s);
}
  

 private:

  char * Label_;
  
};

inline ostream& operator<<(ostream& os, const Epetra_Object& obj)
{
  os << obj.Label();
  obj.Print(os);
  return os;
}


#endif /* _EPETRA_OBJECT_H_ */
