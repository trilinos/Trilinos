// Kris
// 07.08.03 -- Move into Teuchos package/namespace

#ifndef _TEUCHOS_REPORTERROR_HPP_
#define _TEUCHOS_REPORTERROR_HPP_

/*! \file Teuchos_ReportError.hpp 
    \brief Teuchos::reportError error reporting method.
*/

namespace Teuchos {
	
	/*! \enum reportError
      A static method called by Teuchos modules to throw an exception.
			reportError is called inside of the throw call, and errorCode
			is the integer value that the exception actually throws.

			The preprocessor directive TEUCHOS_NO_ERROR_REPORTS can be used
			to prevent reportError from outputing anything.
			This will not prevent exceptions from being thrown.
	*/

  // NOTE:  We are extracting a C-style string from Message because 
  //        the SGI compiler does not have a real string class with 
  //        the << operator.  Some day we should get rid of ".c_str()"

	static inline int reportError(const string message, int errorCode) {
#ifndef TEUCHOS_NO_ERROR_REPORTS
		cerr << endl << "Error in Teuchos Module:  " << __FILE__ << " on line " << __LINE__ << endl 
				 << "Teuchos Error:  " << message.c_str() << "  Error Code:  " << errorCode << endl;
#endif
		return(errorCode);
	}
	
} // namespace Teuchos
#endif /* _TEUCHOS_REPORTERROR_HPP_ */
