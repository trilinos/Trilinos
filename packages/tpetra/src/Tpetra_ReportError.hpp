/*Paul
Modified 04-Feb-2003
*/

#ifndef _TPETRA_REPORTERROR_HPP_
#define _TPETRA_REPORTERROR_HPP_

/*! \file Tpetra_ReportError.hpp 
    \brief Tpetra::reportError error reporting method.
*/

namespace Tpetra {
	
	/*! \enum reportError
      A static method called by Tpetra modules to throw an exception.
			reportError is called inside of the throw call, and errorCode
			is the integer value that the exception actually throws.

			The preprocessor directive TPETRA_NO_ERROR_REPORTS can be used
			to prevent reportError from outputing anything.
			This will not prevent exceptions from being thrown.
	*/

  // NOTE:  We are extracting a C-style string from Message because 
  //        the SGI compiler does not have a real string class with 
  //        the << operator.  Some day we should get rid of ".c_str()"

	static inline int reportError(const string message, int errorCode) {
#ifndef TPETRA_NO_ERROR_REPORTS
		cerr << endl << "Error in Tpetra Module:  " << __FILE__ << " on line " << __LINE__ << endl 
				 << "Tpetra Error:  " << message.c_str() << "  Error Code:  " << errorCode << endl;
#endif
		return(errorCode);
	}
	
} // namespace Tpetra
#endif /* _TPETRA_REPORTERROR_HPP_ */
