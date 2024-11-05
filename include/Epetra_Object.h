/*
//@HEADER
// ************************************************************************
//
//               Epetra: Linear Algebra Services Package
//                 Copyright 2011 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// ************************************************************************
//@HEADER
*/

#ifndef EPETRA_OBJECT_H
#define EPETRA_OBJECT_H

#if defined(Epetra_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Epetra package is deprecated"
#endif
#endif



#include "Epetra_CombineMode.h"
#include "Epetra_DataAccess.h"
#include "Epetra_ConfigDefs.h"

//! Epetra_Object:  The base Epetra class.
/*! The Epetra_Object class provides capabilities common to all Epetra objects,
    such as a label that identifies an object instance, constant definitions,
    enum types.

*/
class EPETRA_LIB_DLL_EXPORT Epetra_Object {

  public:
    //! @name Constructors/destructor
  //@{
  //! Epetra_Object Constructor.
  /*! Epetra_Object is the primary base class in Epetra.  All Epetra class
      are derived from it, directly or indirectly.  This class is seldom
      used explictly.
  */
  Epetra_Object(int TracebackModeIn = -1, bool set_label = true);

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

  //! @name Attribute set/get methods
  //@{

  //! Epetra_Object Label definition using char *.
  /*! Defines the label used to describe the \e this object.
  */
  virtual void SetLabel(const char * const Label);

  //! Epetra_Object Label access funtion.
  /*! Returns the string used to define this object.
  */
  virtual const char * Label() const;

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

  //! Get the output stream for error reporting
  static std::ostream& GetTracebackStream();

  //@}

  //! @name Miscellaneous
  //@{

  //! Print object to an output stream
  //! Print method
  virtual void Print(std::ostream & os) const;

  //! Error reporting method
  virtual int ReportError(const std::string Message, int ErrorCode) const;
  //@}


// TracebackMode controls how much traceback information is printed when run time
// integer errors are detected:
// = 0 - No information report
// = 1 - Fatal (negative) values are reported
// = 2 - All values (except zero) reported.

// Default is set to 1.  Can be set to different value using SetTracebackMode() method in
// Epetra_Object class
  static int TracebackMode;


 protected:
  std::string toString(const int& x) const {
     char s[100];
     snprintf(s, sizeof(s), "%d", x);
     return std::string(s);
}
  std::string toString(const long long& x) const {
     char s[100];
     snprintf(s, sizeof(s), "%lld", x);
     return std::string(s);
}

  std::string toString(const double& x) const {
     char s[100];
     snprintf(s, sizeof(s), "%g", x);
     return std::string(s);
}


 private:
  Epetra_Object& operator=(const Epetra_Object& src) {
    SetLabel(src.Label());
    return *this;
  }

  char * Label_;

};

inline std::ostream& operator<<(std::ostream& os, const Epetra_Object& obj)
{
  if (Epetra_FormatStdout) {
/*    const Epetra_fmtflags  olda = os.setf(ios::right,ios::adjustfield);
    const Epetra_fmtflags  oldf = os.setf(ios::scientific,ios::floatfield);
    const int              oldp = os.precision(12); */

    os << obj.Label() << std::endl;
    obj.Print(os);

/*    os.setf(olda,ios::adjustfield);
    os.setf(oldf,ios::floatfield);
    os.precision(oldp); */
  }
  else {

    os << obj.Label();
    obj.Print(os);
  }

  return os;
}

/** \brief Macro for testing for and throwing and int exception for objects
 * derived from Epetra_Object.
 *
 * This macro adds the file name and line number to teh
 */
#define EPETRA_TEST_FOR_EXCEPTION(throw_exception_test,errCode,msg) \
{ \
    const bool throw_exception = (throw_exception_test); \
    if(throw_exception) { \
        std::ostringstream omsg; \
	    omsg \
        << __FILE__ << ":" << __LINE__ << ":" \
        << " Throw test that evaluated to true: "#throw_exception_test << ":" \
        << "Error message : " << msg; \
	    throw ReportError(omsg.str(),errCode); \
    } \
}

#endif /* EPETRA_OBJECT_H */
