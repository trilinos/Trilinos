/*Paul
27-May-2002 General cleanup. Changed method names to fit namingConvention. Commented out reportError method (see Tpetra_Object.cpp for details).
*/

#ifndef _TPETRA_OBJECT_H_
#define _TPETRA_OBJECT_H_

#include "Tpetra_Common.h"
#include "Tpetra_CombineMode.h"
#include "Tpetra_DataAccess.h"


//! Tpetra::Object:  The base Tpetra class.
/*! The Tpetra::Object class provides capabilities common to all Tpetra objects,
    such as a label that identifies an object instance, constant definitions,
    enum types.
  
*/

namespace Tpetra
{

class Object
{
  public:
  //@{ \name Constructors/destructor.
  //! Tpetra::Object Constructor.
  /*! Tpetra::Object is the primary base class in Tpetra.  All Tpetra class
      are derived from it, directly or indirectly.  This class is seldom
      used explictly.
  */
  Object(int TracebackModeIn = -1);

  //! Tpetra::Object Constructor.
  /*! Creates a Tpetra::Object with the given label.
  */
  Object(const char* const Label, int TracebackModeIn = -1);

  //! Tpetra::Object Copy Constructor.
  /*! Makes an exact copy of an existing Tpetra::Object instance.
  */
  Object(const Object& Obj);

  //! Tpetra::Object Destructor.
  /*! Completely deletes a Tpetra::Object object.  
  */
  virtual ~Object();
  //@}
  
  //@{ \name Attribute set/get methods.

  //! Tpetra::Object Label definition using char *.
  /*! Defines the label used to describe the \e this object.  
  */
  virtual void setLabel(const char* const Label);

  //! Tpetra::Object Label access funtion.
  /*! Returns the string used to define this object.  
  */
  virtual char* label() const;

  //! Set the value of the Tpetra::Object error traceback report mode.
  /*! Sets the integer error traceback behavior.  
      TracebackMode controls whether or not traceback information is printed when run time 
      integer errors are detected:

      <= 0 - No information report

       = 1 - Fatal (negative) values are reported

      >= 2 - All values (except zero) reported.

      Default is set to 1.
  */
  static void setTracebackMode(int TracebackModeValue);

  //! Get the value of the Tpetra::Object error report mode.
  static int getTracebackMode();
  //@}

  //@{ \name Miscellaneous

  //! Print object to an output stream
  //! Print method
  virtual void print(ostream& os) const;

  //! Error reporting method.
  virtual int reportError(const string Message, int ErrorCode) const;
  //@}

  
// TracebackMode controls how much traceback information is printed when run time 
// integer errors are detected:
// = 0 - No information report
// = 1 - Fatal (negative) values are reported
// = 2 - All values (except zero) reported.

// Default is set to 2.  Can be set to different value using setTracebackMode() method in
// Tpetra::Object class
  static int TracebackMode;


 protected:
  string toString(const int& x) const
{
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

  char* Label_;

}; // class Object

} // namespace Tpetra

#include "Tpetra_Object.cpp"

inline ostream& operator<<(ostream& os, const Tpetra::Object& obj)
{
  if (Tpetra_FormatStdout)
  {
    const Tpetra_fmtflags  olda = os.setf(ios::right,ios::adjustfield);
    const Tpetra_fmtflags  oldf = os.setf(ios::scientific,ios::floatfield);
    const int              oldp = os.precision(12);

    os << obj.label() << endl;
    obj.print(os);

    os.setf(olda,ios::adjustfield);
    os.setf(oldf,ios::floatfield);
    os.precision(oldp);
  }
  else
  {

    os << obj.label();
    obj.print(os);
  }
  
  return os;
}


#endif /* _TPETRA_OBJECT_H_ */
