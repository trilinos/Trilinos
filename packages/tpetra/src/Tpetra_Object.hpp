/*Paul
27-May-2002 General cleanup. Changed method names to fit namingConvention. Commented out reportError method (see Tpetra_Object.cpp for details).
06-August-2002 Changed to images (nothing changed). Also touched up the documentation a bit, and updated some naming conventions.
12-Oct-2002 Updated for Common->Compiler_Directives renaming.
30-Oct-2002 Updated for Compiler_Directives -> ConfigDefs renaming.
04-Feb-2003 Moved reportError method into seperate file.
*/

#ifndef _TPETRA_OBJECT_HPP_
#define _TPETRA_OBJECT_HPP_

#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_CombineMode.hpp"
#include "Tpetra_DataAccess.hpp"
//#include "Tpetra_ReportError.hpp"

namespace Tpetra
{

//! Tpetra::Object:  The base Tpetra class.
/*! The Object class provides capabilities common to all Tpetra objects,
    such as a label that identifies an object instance, constant definitions,
    enum types.
*/

class Object
{
  public:
  //@{ \name Constructors/destructor.
  //! Object Constructor.
  /*! Object is the primary base class in Tpetra.  All Tpetra class
      are derived from it, directly or indirectly.  This class is seldom
      used explictly.
  */
  Object(int tracebackModeIn = -1);

  //! Object Constructor.
  /*! Creates an Object with the given label.
  */
  Object(const char* const label, int tracebackModeIn = -1);

  //! Object Copy Constructor.
  /*! Makes an exact copy of an existing Object instance.
  */
  Object(const Object& obj);

  //! Object Destructor.
  /*! Completely deletes an Object object.  
  */
  virtual ~Object();
  //@}
  
  //@{ \name Attribute set/get methods.

  //! Object Label definition using char *.
  /*! Defines the label used to describe the \e this object.  
  */
  virtual void setLabel(const char* const label);

  //! Object Label access funtion.
  /*! Returns the string used to define this object.  
  */
  virtual char* label() const;

  //! Set the value of the Object error traceback report mode.
  /*! Sets the integer error traceback behavior.  
      TracebackMode controls whether or not traceback information is printed when run time 
      integer errors are detected:

      <= 0 - No information report

       = 1 - Fatal (negative) values are reported

      >= 2 - All values (except zero) reported.

      Default is set to 1.
  */
  static void setTracebackMode(int tracebackModeValue);

  //! Get the value of the Object error report mode.
  static int getTracebackMode();
  //@}

  //@{ \name Miscellaneous

  //! Print object to an output stream
  //! Print method
  virtual void print(ostream& os) const;

	//! Error reporting method.
	virtual int reportError(const string message, int errorCode) const {
  // NOTE:  We are extracting a C-style string from Message because 
  //        the SGI compiler does not have a real string class with 
  //        the << operator.  Some day we should get rid of ".c_str()"
#ifndef TPETRA_NO_ERROR_REPORTS
		cerr << endl << "Error in Tpetra Object with label: " << label_ << endl 
				 << "Tpetra Error:  " << message.c_str() << "  Error Code:  " << errorCode << endl;
#endif
		return(errorCode);
	}
  
// tracebackMode controls how much traceback information is printed when run time 
// integer errors are detected:
// = 0 - No information report
// = 1 - Fatal (negative) values are reported
// = 2 - All values (except zero) reported.

// Default is set to 2.  Can be set to different value using setTracebackMode() method in
// Object class
  static int tracebackMode;


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

  char* label_;

}; // class Object

} // namespace Tpetra

#include "Tpetra_Object.cpp"

inline ostream& operator<<(ostream& os, const Tpetra::Object& Obj)
{
  os << Obj.label() << endl;
  Obj.print(os);
 
  return os;
}


#endif /* _TPETRA_OBJECT_HPP_ */
