// /////////////////////////////////////////////
// Teuchos_BLAS_types.hpp

#ifndef TEUCHOS_BLAS_TYPES_HPP
#define TEUCHOS_BLAS_TYPES_HPP

/*! 	\file Teuchos_BLAS_types.hpp
	\brief Enumerated types for BLAS input characters.
	
	These enumerated lists are used in compile time checking of the input characters
	for BLAS methods.  

	\note Any other input other than those specified here will result
	in an error at compile time and are not supported by the templated BLAS interface.

	<ul>
	<li><b>Teuchos::ESide</b> : Enumerated list for BLAS character input "SIDE".
		<ul>
		<li>LEFT_SIDE : The matrix/vector is on the left side of the equation
		<li>RIGHT_SIDE : The matrix/vector is on the right side of the equation
		</ul><br>
	<li><b>Teuchos::ETransp</b> : Enumerated list for BLAS character input "TRANS".
		<ul>
		<li>NO_TRANS : The matrix/vector is not transposed
		<li>TRANS : The matrix/vector is transposed
		<li>CONJ_TRANS : The matrix/vector is conjugate transposed
		</ul><br>
	<li><b>Teuchos::EUplo</b> : Enumerated list for BLAS character input "UPLO".
		<ul>
		<li>UPPER_TRI : The matrix is upper triangular
		<li>LOWER_TRI : The matrix is lower triangular
		</ul><br>
	<li><b>Teuchos::EDiag</b> : Enumerated list for BLAS character input "DIAG".
		<ul>
		<li>UNIT_DIAG : The matrix has all ones on its diagonal
		<li>NON_UNIT_DIAG : The matrix does not have all ones on its diagonal
		</ul>
	</ul>
*/

namespace Teuchos {
	enum ESide { 	
		 LEFT_SIDE,	/*!< Left side */ 
		 RIGHT_SIDE 	/*!< Right side */
		 };

	enum ETransp { 	
		NO_TRANS,	/*!< Not transposed */ 
		TRANS, 		/*!< Transposed */
		CONJ_TRANS 	/*!< Conjugate transposed */
		};

	enum EUplo { 	
		UPPER_TRI,	/*!< Upper triangular */ 
		LOWER_TRI 	/*!< Lower triangular */
		};

	enum EDiag { 	
		UNIT_DIAG,	/*!< Unit diagaonal */ 
		NON_UNIT_DIAG	/*!< Not unit diagonal */ 
		};
}

#endif // TEUCHOS_BLAS_TYPES_HPP
