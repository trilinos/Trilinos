// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
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

// Kris
// 07.08.03 -- Move into Teuchos package/namespace

#ifndef _TEUCHOS_OBJECT_HPP_
#define _TEUCHOS_OBJECT_HPP_

/*! \file Teuchos_Object.hpp
    \brief The base Teuchos object.
*/

#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_DataAccess.hpp"

// 2007/11/26: rabartl: This class has to change from using 'char*' to
// std::string!

/*! \class Teuchos::Object
    \brief The base Teuchos class.

    The Object class provides capabilities common to all Teuchos objects,
    such as a label that identifies an object instance, constant definitions,
    enum types.
*/

namespace Teuchos {

class TEUCHOSNUMERICS_LIB_DLL_EXPORT Object {
public:
  //! @name Constructors/Destructor.
  //@{
  //! Default Constructor.
  /*! Object is the primary base class in Teuchos.  All Teuchos class [sic]
      are derived from it, directly or indirectly.  This class is seldom
      used explictly.

     \warning (mfh 23 Nov 2014) Pretty much just ignore the above
       description of this class.  Very few classes in Teuchos inherit
       from Object.  Such inheritance should be treated as deprecated
       legacy behavior.  It's not 1990 and C++ is not Java 1.0; we
       don't need a common base class for everything.
  */
  Object (int tracebackModeIn = -1);

  //! Labeling Constructor.
  /*! Creates an Object with the given label.
   *
   * LEGACY; DEPRECATE.
  */
  Object (const char* label, int tracebackModeIn = -1);

  /// \brief Create an Object with the given label, and optionally,
  ///   with the given traceback mode.
  Object (const std::string& label, int tracebackModeIn = -1);

  //! Destructor (virtual, for safety of derived classes).
  virtual ~Object () {}

  //@}
  //! @name Set methods.
  //@{

  // LEGACY; REPLACE "const char*" with std::string.
  virtual void setLabel (const char* theLabel);

  /// \brief Set the value of the Object error traceback report mode.
  ///
  /// TracebackMode controls whether or not traceback information is
  /// printed when run time integer errors are detected:
  ///
  /// <= 0 - No information report
  ///
  /// = 1 - Fatal (negative) values are reported
  ///
  /// >= 2 - All values (except zero) reported.
  ///
  /// \note Default is set to -1 when object is constructed.
  static void setTracebackMode (int tracebackModeValue);

  //@}
  //! @name Accessor methods.
  //@{

  //! Access the object's label (LEGACY; return std::string instead).
  virtual const char* label () const;

  //! Get the value of the Object error traceback report mode.
  static int getTracebackMode();

  //@}
  //! @name I/O method.
  //@{

  //! Print the object to the given output stream.
  virtual void print (std::ostream& os) const;

  //@}
  //! @name Error reporting method.
  //@{

  //! Report an error with this Object.
  virtual int reportError (const std::string message, int errorCode) const;

  //@}

  static int tracebackMode;

private:
  //! The Object's current label.
  std::string label_;
};

/// \brief Print the given Object to the given output stream.
/// \relates Object
std::ostream& operator<< (std::ostream& os, const Teuchos::Object& obj);

} // namespace Teuchos

#endif /* _TEUCHOS_OBJECT_HPP_ */
