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


#ifndef TEUCHOS_C_WRAPPER_SUPPORT_CPP_HPP
#define TEUCHOS_C_WRAPPER_SUPPORT_CPP_HPP


#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_FancyOStream.hpp"


/** \defgroup Teuchos_CWrapperSupport_grp .*/


namespace Teuchos {


/** \brief Static C Wrapper Error Handling Policy Class.
 *
 */
class TEUCHOS_LIB_DLL_EXPORT CWrapperErrorHandling {
public:

  /** \brief Set the ostream that will be printed to when errors occur.
   *
   * \ingroup Teuchos_CWrapperSupport_grp
   */
  static void setPrintErrorOStream(const RCP<FancyOStream> &errorOStream);
  
  /** \brief Get the ostream that will be printed when errors occur.
   *
   * \ingroup Teuchos_CWrapperSupport_grp
   */
  static RCP<FancyOStream> getPrintErrorOStream();

  /** \brief Set if the stacktrace should be shown on every caught
   * exception.
   *
   * This only is meaningful if stacktracing is configured and enabled.
   */
  static void setShowStackTraceOnException(const bool showStrackTraceOnException);

  /** \brief Get if the stacktrace should be shown on every caught exception.
   */
  static bool getShowStackTraceOnException();
  

}; // class CWrapperErrorHandling


} // Teuchos



/** \brief Define a try block.
 *
 * \ingroup Teuchos_CWrapperSupport_grp
 */
#define TEUCHOS_CWRAPPER_TRY(IERR) \
  TEUCHOS_TEST_FOR_EXCEPT((IERR) == 0); \
  (*(IERR)) = 0; \
  bool cwrapper_try_success = true; \
  try
// Above: We must set *ierr = 0 in case a return statement will return the
// value.


/** \brief Define the catch blocks and set the error code.
 *
 * \ingroup Teuchos_CWrapperSupport_grp
 */
#define TEUCHOS_CWRAPPER_CATCH_ERROR_CODE(IERR) \
  TEUCHOS_STANDARD_CATCH_STATEMENTS_IMPL( \
    true, *Teuchos::CWrapperErrorHandling::getPrintErrorOStream(), \
    Teuchos::CWrapperErrorHandling::getShowStackTraceOnException(), \
    cwrapper_try_success ); \
  if (!cwrapper_try_success) { (*(IERR)) = -1; }
// Above: We make sure and set the error code in case there is a failure.



/** \brief Set the error code.
 *
 * This function is to be used inside of the try/catch blocks to set the error
 * code.
 *
 * \ingroup Teuchos_CWrapperSupport_grp
 */
#define TEUCHOS_CWRAPPER_SET_ERROR_CODE(IERR, IERR_VALUE) \
  if ((IERR)) { \
    (*(IERR)) = (IERR_VALUE); \
  }


#endif // TEUCHOS_C_WRAPPER_SUPPORT_CPP_HPP
