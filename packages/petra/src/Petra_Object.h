#ifndef _PETRA_OBJECT_H_
#define _PETRA_OBJECT_H_

//! Petra_Object:  The base Petra class.
/*! The Petra_Object class provides capabilities common to all Petra objects,
    such as a label that identifies an object instance, constant definitions,
    enum types.
  
*/

#undef PETRA_LEVELSCHEDULING

#ifndef __cplusplus
#define __cplusplus
#endif

#undef NDEBUG             // make sure asserts are enabled

#if defined(SGI) || defined(SGI64) || defined(SGI32) || defined(SOLARIS) || defined(TFLOP)
#include <stdlib.h>
#include <assert.h>
#include <iostream.h>
#include <strstream.h>
#include <string.h>
#include <math.h>
#include <iomanip.h>

#else

#include <cstdlib>
#include <cassert>
#include <iostream>
#include <strstream>
#include <cstring>
#include <string>
#include <cmath>
#include <iomanip>
using namespace std;

#endif



#ifdef PETRA_MPI
#include <mpi.h>
#endif

#ifdef PETRA_SIMULATE_BOOL
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

#define PETRA_MAX(x,y) (( x > y ) ? x : y)     /* max function  */
#define PETRA_MIN(x,y) (( x < y ) ? x : y)     /* min function  */
#define PETRA_SGN(x) ((x < 0.0) ? -1.0 : 1.0)  /* sign function */

// Make error report silent by defining PETRA_NO_ERROR_REPORTS
#ifdef PETRA_NO_ERROR_REPORTS
#undef PETRA_CHK_ERR
#undef PETRA_CHK_PTR
#undef PETRA_CHK_REF

#else

#ifdef PETRA_CHK_ERR
#undef PETRA_CHK_ERR
#endif
#ifdef PETRA_CHK_PTR
#undef PETRA_CHK_PTR
#endif
#ifdef PETRA_CHK_REF
#undef PETRA_CHK_REF
#endif

// Great little macro obtained from Alan Williams

#define PETRA_CHK_ERR(a) { int petra_err = a; if (petra_err != 0) { \
                      cerr << "Petra ERROR " << petra_err << ", " \
                           << __FILE__ << ", line " << __LINE__ << endl; } \
                      return(petra_err); \
                   }

// Extension of same macro for pointer, returns zero if bad

#define PETRA_CHK_PTR(a) { if (a == 0) { \
                      cerr << "Petra returning zero pointer " << ", " \
                           << __FILE__ << ", line " << __LINE__ << endl; } \
                      return(a); \
                   }
// Extension of same macro for reference, returns a default reference

#define PETRA_CHK_REF(a) { \
                      cerr << "Petra returning default reference " << ", " \
                           << __FILE__ << ", line " << __LINE__ << endl; \
                      return(a); \
                   }
#endif



class Petra_Object {
    
  public:
  //! Petra_Object Constructor.
  /*! Creates a Petra_Object instance.
  */
  Petra_Object(void);

  //! Petra_Object Constructor.
  /*! Creates a Petra_Object with the given label.
  */
  Petra_Object(const char * const Label);

  //! Petra_Object Copy Constructor.
  /*! Makes an exact copy of an existing Petra_Object instance.
  */
  Petra_Object(const Petra_Object& Object);


  //! Petra_Object Label definition using char *.
  /*! Defines the label used to describe the \e this object.  
  */
  virtual void SetLabel(const char * const Label);

  //! Petra_Object Label access funtion.
  /*! Returns the string used to define this object.  
  */
  virtual char * Label(void) const;

  //! Print method
  virtual void Print(ostream & os) const;


  //! Petra_Object Destructor.
  /*! Completely deletes a Petra_Object object.  
  */
  virtual ~Petra_Object(void);

 protected:
  

 private:

  char * Label_;
  
};
inline ostream& operator<<(ostream& os, const Petra_Object& obj)
{
  os << obj.Label();
  obj.Print(os);
  return os;
}


#endif /* _PETRA_OBJECT_H_ */
