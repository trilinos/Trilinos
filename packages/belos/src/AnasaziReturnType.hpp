// file AnasaziReturnType.hpp
#ifndef ANASAZI_RETURN_TYPE_HPP
#define ANASAZI_RETURN_TYPE_HPP

/*!
	\brief The apply methods in the AnasaziEigenproblem or Eigensolver class
	may fail or not be defined, this information needs to be passed back to 
	the algorithm.  This will be used in the algorithm to decide what should 
	be done and determine if the information provided by the AnasaziEigenproblem 
	is sufficient to solve the eigenvalue problem.  The enumerated list is being
	put in the Anasazi namespace to distinguish it from any other enumerated list
	that may be included in Anasazi or its interfaces.
*/
	enum Anasazi_ReturnType {	Ok, 		//! Computation completed sucessfully
					Undefined, 	//! This operation is not defined
					Unconverged,    //! This operation returned unconverged
					Failed		//! Any other numerical failure in thie computation 
	};

#endif
// end of file AnasaziReturnType.hpp
