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

// Author: Ian Karlin ikarlin@sandia.gov 05-28-2008

#ifndef EPETRA_OSKIMATRIX_H
#define EPETRA_OSKIMATRIX_H

#if defined(Epetra_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Epetra package is deprecated"
#endif
#endif



#include "Epetra_OskiMultiVector.h"
#include "Epetra_OskiVector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_OskiPermutation.h"
#include "Teuchos_ParameterList.hpp"
#include "Epetra_Comm.h"
extern "C" {
  #include "oski/oski.h"
}

class Epetra_OskiVector;
class Epetra_OskiMultiVector;
class Teuchos_ParameterList;
class Epetra_OskiPermutation;

//! Epetra_OskiMatrix: A class for constructing and using OSKI Matrices within Epetra.  For information on known issues with OSKI see the detailed description.

/*! OSKI is a high-performance sparse matrix kernel package written by the UC Berkeley Benchmarking and
    Optimization Group.  The Epetra_OskiMatrix class is a lightweight interface to allow Epetra users access to OSKI functionality.
    This interface includes lightweight conversion of Epetra matrices to OSKI matrices, runtime tuning based on parameter list hints,
    and access to high-performance computational kernels to perform matrix-vector/matrix multi-vector calculations.

    The Epetra_OskiMatrix class provides access to the entire OSKI interface.  However, the following features are not fully implemented in OSKI
    version oski-1.0.1h, and therefore are unsupported in the Epetra_OskiMatrix class:

    - OSKI does not provide stock composed kernels. Hence, the tune function must be called in order to see performance gains when using the MatTransMatMultiply and MultiplyAndMatTransMultiply kernels.
    - The MatPowMultiply kernel does not work.
    - Optimized multivector kernels are not created by default when installing OSKI.
    - The tune function cannot transform a (nearly) symmetric matrix to be stored as such.
    - In order to use the \f$A^TA\f$ OSKI kernel (MatTransMatMultiply), in oski/src/MBCSR/ata.c you must replace the lines
    \code const char *kernel_name = (opA == OP_AT_A)
        ? "SubmatRperTransSubmatRperMult" : "SubmatRperHermSubmatRperMult";
    \endcode
    with
    \code const char *kernel_name = (opA == OP_AT_A)
        ? "SubmatReprTransSubmatReprMult" : "SubmatReprHermSubmatReprMult";
    \endcode
    - OSKI does not convert between CSR and CSC when it could be profitable, such as when performing
    \f$AA^T\f$ on a CSR matrix.
    - OSKI may be incompatible with the following architectures:
        - Barcelona (quad-core Opteron): errors during "make install" (confirmed with OSKI developers)
        - single core Xeon: OSKI installs, but never transforms matrices.  This includes cases where other machines will transform the same
          matrices, and where one would expect the matrix to be transformed, based on OSKI tuning data.
*/

class Epetra_OskiMatrix: public Epetra_CrsMatrix{
 public:
	//! @name Constructors/Destructor
	//@{
        //! Copy constructor.
        Epetra_OskiMatrix(const Epetra_OskiMatrix& Source); //not in use for now

	//! Constructor creates an Epetra_OskiMatrix from an Epetra_CrsMatrix.
	/*! \param Source (In) An Epetra_CrsMatrix that is to be wrapped as an Epetra_OskiMatrix.
	    \param List (In) Any options or data wanted or needed for the conversion.
            \return Pointer to an Epetra_OskiMatrix.

            Options that can be passed to the List are presented below.  They are: "<type> <option name> <default value>: <description of purpose>"

            - bool autotune false: If true, Epetra tries to set as many hints as possible based on its knowledge of the matrix.
	    - string matrixtype general: Other types that can be taken are: uppertri, lowertri, uppersymm, lowersymm, fullsymm, upperherm, lowerherm and fullherm.
            - bool diagstored false: If true, the diagonal entries are not stored in the matrix and are all assumed to be 1.
	    - bool zerobased false: If true, the array is zero based, as in C. Otherwise, it is 1 based, as in Fortran.
            - bool sorted false: If true, all elements in the passed in array are sorted.
            - bool unique false: If true, a value in a column only appears once in each row.
            - bool deepcopy false: If true, when the OSKI matrix is created it will be a deepcopy of the data in the function.
	*/
	Epetra_OskiMatrix(const Epetra_CrsMatrix& Source, const Teuchos::ParameterList& List);
	
	//! Destructor	
	virtual ~Epetra_OskiMatrix();
	//@}

	//! @name Extract/Replace Values
	//@{
	//! Replace current values with this list of entries for a given local row of the matrix.  Warning this could be expensive.
	/*! The reason this function could be expensive is its underlying implementation.
	    Both the OSKI and Epetra versions of the matrix must be changed when the matrix
	    has been permuted.  When this is the case, a call must be made to the Epetra
	    ReplaceMyValues, and NumEntries calls must be made to a function that changes the
	    OSKI matrix's values one at a time.
	    \param MyRow (In) Row number (in local coordinates) to put elements.
	    \param NumEntries (In) Number of entries.
	    \param Values (In) Values to enter.
	    \param Indices (In) Local column indices corresponding to the values.
       	    \return Integer error code, set to 0 if successful. Note that if the
    		    allocated length of the row has to be expanded, Oski will fail.
                    A positive warning code may be returned, but this should be treated
		    as a fatal error; part of the data will be changed, and OSKI cannot
                    support adding in new data values.
    	    \pre IndicesAreLocal()==true
    	    \post The given Values at the given Indices have been summed into the
    	   	  entries of MyRow.

	*/
	int ReplaceMyValues(int MyRow,
			    int NumEntries,
			    double* Values,
			    int* Indices);

   	//! Add this list of entries to existing values for a given local row of the matrix.  WARNING: this could be expensive.
  	/*! The reason this function could be expensive is its underlying implementation.
            Both the OSKI and Epetra versions of the Matrix must be changed when the matrix
            has been permuted.  When this is the case, a call must be made to the Epetra
            SumIntoMyValues, and NumEntries calls must be made to a function that changes the
            OSKI matrix's values one at a time.
	    \param MyRow - (In) Row number (in local coordinates) to put elements.
    	    \param NumEntries - (In) Number of entries.
    	    \param Values - (In) Values to enter.
    	    \param Indices - (In) Local column indices corresponding to values.
    	    \return Integer error code, set to 0 if successful. Note that if the
    		    allocated length of the row has to be expanded, a positive warning code
    		    may be returned.  This should be treated as a fatal error, as part of
		    the data will be changed, and OSKI cannot support adding in new data values.
    	    \pre IndicesAreLocal()==true
    	    \post The given Values at the given Indices have been summed into the
    		  entries of MyRow.
	*/
	int SumIntoMyValues(int MyRow,
			    int NumEntries,
			    double* Values,
			    int* Indices);


	//! Returns a copy of the main diagonal in a user-provided vector.
  	/*! \param Diagonal - (Out) Extracted main diagonal.
    	    \return Integer error code, set to 0 if successful and non-zero if not.
    	    \pre Filled()==true
    	    \post Unchanged.
  	*/
	int ExtractDiagonalCopy(Epetra_Vector& Diagonal) const;
	
	//! Replaces diagonal values of the matrix with those in the user-provided vector.
  	/*! This routine is meant to allow replacement of {\bf existing} diagonal values.
     	    If a diagonal value does not exist for a given row, the corresponding value in
            the input Epetra_OskiVector will be ignored, and the return code will be set to 1.
    	    \param Diagonal - (In) New values to be placed in the main diagonal.
    	    \return Integer error code, set to 0 if successful, set to 1 on the calling processor
	       	    if one or more diagonal entries not present in matrix.  Other error codes can
		    be returned as well, indicating improperly constructed matrices or vectors.
    	    \pre Filled()==true
    	    \post Diagonal values have been replaced with the values of Diagonal.
	*/
	int ReplaceDiagonalValues(const Epetra_OskiVector& Diagonal);
	//@}	

	//! @name Computational methods
	//@{
	//! Performs a matrix vector multiply of y = this^TransA*x
	/*! The vectors x and y can be either Epetra_Vectors or Epetra_OskiVectors.
	    \param TransA (In) If TransA = TRUE then use the transpose of the matrix in
	           computing the product.
	    \param x (In) The vector the matrix is multiplied by.
	    \param y (Out) The vector where the calculation result is stored.
	    \return Integer error code, set to 0 if successful.
            \pre Filled()==true
            \post Unchanged
	*/
  	int Multiply(bool TransA,
	    	     const Epetra_Vector& x,
		     Epetra_Vector& y) const;
	
	//! Performs a matrix vector multiply of y = Alpha*this^TransA*x + Beta*y
	/*! The vectors x and y can be either Epetra_Vectors or Epetra_OskiVectors.
	    \param TransA (In) If TransA = TRUE then use the transpose of the matrix in
	           computing the product.
	    \param x (In) The vector the matrix is multiplied by.
	    \param y (In/Out) The vector where the calculation result is stored.
	    \param Alpha (In) A scalar constant used to scale x.
	    \param Beta  (In) A scalar constant used to scale y.
	    \return Integer error code, set to 0 if successful.
            \pre Filled()==true
            \post Unchanged
	*/
  	int Multiply(bool TransA,
	    	     const Epetra_Vector& x,
		     Epetra_Vector& y,
		     double Alpha,
		     double Beta = 0.0) const;
	
	//! Performs a matrix multi-vector multiply of Y = this^TransA*X
	/*! The multi-vectors X and Y can be either Epetra_MultiVectors or Epetra_OskiMultiVectors.
	    \param TransA (In) If Trans = TRUE then use the transpose of the matrix in
	           computing the product.
	    \param X (In) The multi-vector the matrix is multiplied by.
	    \param Y (Out) The multi-vector where the calculation result is stored.
	    \return Integer error code, set to 0 if successful.
	    \pre Filled()==true
	    \post Unchanged
	*/
  	int Multiply(bool TransA,
   	    	     const Epetra_MultiVector& X,
	    	     Epetra_MultiVector& Y) const;

	//! Performs a matrix multi-vector multiply of Y = Alpha*this^TransA*X + Beta*Y
	/*! The multi-vectors X and Y can be either Epetra_MultiVectors or Epetra_OskiMultiVectors.
	    \param TransA (In) If Trans = TRUE then use the transpose of the matrix in
	           computing the product.
	    \param X (In) The multi-vector the matrix is multiplied by.
	    \param Y (In/Out) The multi-vector where the calculation result is stored.
	    \param Alpha (In) A scalar constant used to scale X.
	    \param Beta  (In) A scalar constant used to scale Y.
	    \return Integer error code, set to 0 if successful.
	    \pre Filled()==true
	    \post Unchanged
	*/
  	int Multiply(bool TransA,
   	    	     const Epetra_MultiVector& X,
	    	     Epetra_MultiVector& Y,
		     double Alpha,
	    	     double Beta = 0.0) const;

	//! Performs a triangular solve of y = (this^TransA)^-1*x where this is a triangular matrix.
	/*! The vector x can be either be an Epetra_Vector or Epetra_OskiVector.  The OskiMatrix must already be triangular, and the UnitDiagonal setting associated with it will be used.
	    \param Upper (In) This parameter is ignored, and is here only to match the Epetra_CrsMatrix::Solve syntax.
	    \param TransA (In) If TransA = TRUE then use the transpose of the matrix in
	           solving the equations.
	    \param UnitDiagonal (In) This parameter is ignored only here to match the Epetra_CrsMatrix::Solve syntax.
	    \param x (In) The vector solved against.
	    \param y (Out) The solution vector.
	    \return Integer error code, set to 0 if successful.
            \pre Filled()==true
            \post Unchanged
	*/
        int Solve(bool Upper, bool TransA, bool UnitDiagonal, const Epetra_Vector& x, Epetra_Vector &y) const;

	//! Performs a triangular solve of y = Alpha*(this^TransA)^-1*x where this is a triangular matrix.
	/*! The vector x can be either be an Epetra_Vector or Epetra_OskiVector.
	    \param TransA (In) If TransA = TRUE then use the transpose of the matrix in
	           solving the equations.
	    \param x (In) The vector solved against.
	    \param y (Out) The solution vector.
	    \param Alpha (In) A scalar constant used to scale x.
	    \return Integer error code, set to 0 if successful.
            \pre Filled()==true
            \post Unchanged
	*/
        int Solve(bool TransA, const Epetra_Vector& x, Epetra_Vector& y, double Alpha = 1.0) const;

	//! Performs a triangular solve of Y = (this^TransA)^-1*X where this is a triangular matrix.
	/*! The vector X can be either be an Epetra_MultiVector or Epetra_OskiMultiVector.  The OskiMatrix must already be triangular, and the UnitDiagonal setting associated with it will be used.
	    \param Upper (In) This parameter is ignored only here to match the Epetra_CrsMatrix::Solve syntax.
	    \param TransA (In) If TransA = TRUE then use the transpose of the matrix in
	           solving the equations.
	    \param UnitDiagonal (In) This parameter is ignored only here to match the Epetra_CrsMatrix::Solve syntax.
	    \param X (In) The multi-vector solved against.
	    \param Y (Out) The solution multi-vector.
	    \return Integer error code, set to 0 if successful.
            \pre Filled()==true
            \post Unchanged
	*/
        int Solve(bool Upper, bool TransA, bool UnitDiagonal, const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

	//! Performs a triangular solve of Y = Alpha*(this^TransA)^-1*X where this is a triangular matrix.
	/*! The vector X can be either be an Epetra_MultiVector or Epetra_OskiMultiVector.
	    \param TransA (In) If TransA = TRUE then use the transpose of the matrix in
	           solving the equations.
	    \param X (In) The multi-vector solved against.
	    \param Y (Out) The solution multi-vector.
	    \param Alpha (In) A scalar constant used to scale X.
	    \return Integer error code, set to 0 if successful.
            \pre Filled()==true
            \post Unchanged
	*/
        int Solve(bool TransA, const Epetra_MultiVector& X, Epetra_MultiVector& Y, double Alpha = 1.0) const;


	//! Performs two matrix vector multiplies of y = Alpha*this^TransA*this*x + Beta*y or y = Alpha*this*this^TransA*x + Beta*y.
	/*! The vectors x, y and t can be either Epetra_Vectors or Epetra_OskiVectors.
	    This composed routine is most commonly used in linear least squares and
	    bidiagonalization methods.  The parallel version of y = Alpha*this*this^TransA*x + Beta*y
            uses calls to the Multiply routine under the hood, as it is not possible to perform
            both multiplies automatically.
	    \param ATA (In) If TransA = TRUE then compute this^T*this*x otherwise compute
		   this*this^T*x.
	    \param x (In) The vector the matrix is multiplied by.
	    \param y (In/Out) The vector where the calculation result is stored.
	    \param t (Out) The vector where the result of the this*x is stored if
		   TransA = true and this^T*x is stored otherwise.
	    \param Alpha (In) A scalar constant used to scale x.
	    \param Beta  (In) A scalar constant used to scale y.
	    \return Integer error code, set to 0 if successful.
            \pre Filled()==true
            \post Unchanged
	*/
	int MatTransMatMultiply(bool ATA,
			 	const Epetra_Vector& x,
				Epetra_Vector& y,
				Epetra_Vector* t,
				double Alpha = 1.0,
				double Beta = 0.0) const;

	//! Performs two matrix multi-vector multiplies of Y = Alpha*this^TransA*this*X + Beta*Y or Y = Alpha*this*this^TransA*X + Beta*Y.
	/*! The multi-vectors X, Y and T can be either Epetra_MultiVectors or Epetra_OskiMultiVectors.
	    This composed routine is most commonly used in linear least squares and
	    bidiagonalization methods.  The parallel version of Y = Alpha*this*this^TransA*X + Beta*Y
            uses calls to the Multiply routine under the hood, as it is not possible to perform
            both multiplies automatically.
	    \param ATA (In) If TransA = TRUE then compute this^T*this*X otherwise compute
		   this*this^T*X.
	    \param X (In) The vector the matrix is multiplied by.
	    \param Y (In/Out) The vector where the calculation result is stored.
	    \param T (Out) The multi-vector where the result of the this*X is stored if
		   TransA = true and this^T*X is stored otherwise.
	    \param Alpha (In) A scalar constant used to scale X.
	    \param Beta  (In) A scalar constant used to scale Y.
	    \return Integer error code, set to 0 if successful.
            \pre Filled()==true
            \post Unchanged
	*/
	int MatTransMatMultiply(bool ATA,
				const Epetra_MultiVector& X,
				Epetra_MultiVector& Y,
				Epetra_MultiVector* T,
				double Alpha = 1.0,
				double Beta = 0.0) const;
	
	//! Performs the two matrix vector multiplies of y = Alpha*this*x + Beta*y and z = Omega*this^TransA*w + Zeta*z.
	/*! The vectors x, y, w and z can be either Epetra_Vectors or Epetra_OskiVectors.
	    This composed routine is most commonly used in bi-conjugate gradient calculations.
	    \param TransA (In) If TransA = TRUE then use the transpose of the matrix in
	           computing the second product.
	    \param x (In) A vector the matrix is multiplied by.
	    \param y (In/Out) A vector where the calculation result of the first multiply is stored.
	    \param w (In) A vector the matrix is multiplied by.
	    \param z (In/Out) A vector where the calculation result of the second multiply is stored.
	    \param Alpha (In) A scalar constant used to scale x.
	    \param Beta  (In) A scalar constant used to scale y.
	    \param Omega  (In) A scalar constant used to scale w.
	    \param Zeta  (In) A scalar constant used to scale z.
	    \return Integer error code, set to 0 if successful.
            \pre Filled()==true
            \post Unchanged
	*/
	int MultiplyAndMatTransMultiply(bool TransA,
				 	const Epetra_Vector& x,
					Epetra_Vector& y,
					const Epetra_Vector& w,
					Epetra_Vector& z,
					double Alpha = 1.0,
					double Beta = 0.0,
					double Omega = 1.0,
					double Zeta = 0.0) const;

	//! Performs the two matrix multi-vector multiplies of Y = Alpha*this*X + Beta*Y and Z = Omega*this^TransA*W + Zeta*Z.
	/*! The multi-vectors X, Y, W and Z can be either Epetra_MultiVectors or Epetra_OskiMultiVectors.
	    This composed routine is most commonly used in bi-conjugate gradient calculations.
	    \param TransA (In) If TransA = TRUE then use the transpose of the matrix in
	           computing the second product.
	    \param X (In) A multi-vector the matrix is multiplied by.
	    \param Y (In/Out) A multi-vector where the calculation result of the first multiply is stored.
	    \param W (In) A multi-vector the matrix is multiplied by.
	    \param Z (In/Out) A multi-vector where the calculation result of the second multiply is stored.
	    \param Alpha (In) A scalar constant used to scale X.
	    \param Beta  (In) A scalar constant used to scale Y.
	    \param Omega  (In) A scalar constant used to scale W.
	    \param Zeta  (In) A scalar constant used to scale Z.
	    \return Integer error code, set to 0 if successful.
            \pre Filled()==true
            \post Unchanged
	*/
	int MultiplyAndMatTransMultiply(bool TransA,
					const Epetra_MultiVector& X,
					Epetra_MultiVector& Y,
					const Epetra_MultiVector& W,
					Epetra_MultiVector& Z,
					double Alpha = 1.0,
					double Beta = 0.0,
					double Omega = 1.0,
					double Zeta = 0.0) const;

	//! Performs a matrix vector multiply of y = Alpha*(this^TransA)^Power*x + Beta*y.  This is not implemented as described in the detailed description.
	/*! The vectors x and y can be either Epetra_Vectors or Epetra_OskiVectors.
	    The vector T can be either an Epetra_MultiVector or Epetra_OskiMultiVector.
	    This composed routine is used in power and S-step methods.  This routine is
            not implemented due a bug in the oski-1.01h kernel that makes testing of correctness
            impossible.
	    \param TransA (In) If TransA = TRUE then use the transpose of the matrix in
	           computing the product.
	    \param x (In) The vector the matrix is multiplied by.
	    \param y (In/Out) The vector where the calculation result is stored.
	    \param T (Out) The multi-vector where the result of each subsequent multiplication
		   this*x ... this^(Power-1)*x is stored.
	    \param Power (In) The power to raise the matrix to in the calculation.
	    \param Alpha (In) A scalar constant used to scale x.
	    \param Beta  (In) A scalar constant used to scale y.
	    \return Integer error code, set to 0 if successful.
            \pre Filled()==true
            \post Unchanged
	*/
	int MatPowMultiply(bool TransA,
			   const Epetra_Vector& x,
			   Epetra_Vector& y,
			   Epetra_MultiVector& T,
 			   int Power = 2,
			   double Alpha = 1.0,
			   double Beta = 0.0) const;

	//! Performs a matrix vector multiply of y = Alpha*(this^TransA)^Power*x + Beta*y.  This is not implemented as described in the detailed description.
	/*! The vectors x and y can be either Epetra_Vectors or Epetra_OskiVectors.
	    This composed routine is used in power and S-step methods.  This routine is
            not implemented due a bug in the oski-1.01h kernel that makes testing of correctness
            impossible.
	    \param TransA (In) If TransA = TRUE then use the transpose of the matrix in
	           computing the product.
	    \param x (In) The vector the matrix is multiplied by.
	    \param y (In/Out) The vector where the calculation result is stored.
	    \param Power (In) The power to raise the matrix to in the calculation.
	    \param Alpha (In) A scalar constant used to scale x.
	    \param Beta  (In) A scalar constant used to scale y.
	    \return Integer error code, set to 0 if successful.
            \pre Filled()==true
            \post Unchanged
	*/
	int MatPowMultiply(bool TransA,
			   const Epetra_Vector& x,
			   Epetra_Vector& y,
 			   int Power = 2,
			   double Alpha = 1.0,
			   double Beta = 0.0) const;
	//@}

	//! @name Tuning
	//@{
	//! Stores the hints in List in the matrix structure.
	/*! \param List (In) A list of hints and options to register along with the matrix
		   used for tuning purposes.  The full list is given below.  It may
		   be moved to the user guide in the future.
	    \return On successful storage of the hint 0 is returned.  On failure an error code
		    is returned.

            Options that can be passed to the List are presented below.  They are: "<type> <option name>: <description of purpose>".  The
available hints are grouped by section, and only one hint from each section can be true for a given matrix.  For options where multiple arguments can be passed in at once the interface only supports up to 5.  This means only 5 block sizes or 5 diaganols can be passed in at once.  If you have more changing the code to support your needs should be simple, but every case you use must be enumerated.  Of course you can just say there are diagonals and blocks and not pass in specific sizes as wells.

	    - bool noblocks: If true, the matrix has no block structure
	    - bool singleblocksize: If true, the matrix structure is dominated by blocks of the size of the next two parameters.
              - int row: The number of rows in each block.
              - int col: The number of columns in each block.
	    - bool multipleblocksize: If true, the matrix consists of multiple block sizes.  The next 3 parameters describe these and are optional.
	      - int blocks: The number of block sizes in the matrix.
	      - int row<x>: Where x is the block number, and x goes from 1 to blocks.  This is the number of rows in block x.
	      - int col<x>: Where x is the block number, and x goes from 1 to blocks.  This is the number of cols in block x.

	    - bool alignedblocks: If true, all blocks are aligned to a grid.
	    - bool unalignedblocks: If true, blocks are not aligned to a grid.
	
	    - bool symmetricpattern: If true, the matrix is either symmetric or nearly symmetric.
	    - bool nonsymmetricpattern: If true, the matrix has a very unsymmetric pattern.
	
	    - bool randompattern: If true, the matrix's non-zeros are distributed in a random pattern.
	    - bool correlatedpattern: If true, the row and column indices for
	      non-zeros are highly correlated.

	    - bool nodiags : If true, the matrix has little or no diagonal structure.
	    - bool diags: If true, the matrix consists of diagonal structure described the next two optional parameters.
              - int numdiags: The number of diagonal sizes known to be present others not listed could be present.
	      - int diag<x>: Where x is the diagonal number, and x goes from 1 to numdiags.  This is the size of the diagonal.
 	*/
	int SetHint(const Teuchos::ParameterList& List);

	
	//! Workload hints for computing a matrix-vector multiply used by OskiTuneMat to optimize the data structure storage, and the routine to compute the calculation.
	/*! In parallel the routine uses symbolic vectors.  This is done for two reasons.  Doing this
            saves on data allocation and potentially communication overhead.  For a matrix-vector
            routine there should be no advantage to having the actual vector, as its size must be the same
            as a matrix dimension.  For a matrix-multivector routine there could be gains from knowing the
            number of vectors in the multi-vector. However, OSKI does not perform multi-vector optimizations,
            so there is no need to add the overhead.
            \param Trans (In) If Trans = true then the transpose of the matrix will be used in
	           computing the product.
	    \param Alpha (In) A scalar constant used to scale InVec.
	    \param InVec (In) The vector the matrix is multiplied by or whether it is a single vector or multi-vector.
	    \param Beta  (In) A scalar constant used to scale OutVec.
	    \param OutVec (In) The vector where the calculation result is stored or whether it is a single vector or multi-vector.
	    \param NumCalls (In) The number of times the operation is called or the tuning level wanted.
	    \param List (In) Used for denoting the use of symbolic vectors for both InVec
		   and OutVec, as well as for level of aggressive tuning if either NumCalls not
		   known or to be overridden.  Options are shown below. It should be noted that
		   by using these options the associated vector or NumCalls becomes invalid.
	    \return Stores the workload hint in the matrix if the operation is valid.  If the
	       	    operation is not valid an error code is returned.
	
            Options that can be passed to the List are presented below.  They are: "<type> <option name>: <description of purpose>".  The
available hints are grouped by section, and only one hint from each section can be true for a given matrix.

	    These replace InVec.
	    - bool symminvec: If true, use a symbolic vector rather than the vector passed in for tuning purposes.
	    - bool symminmultivec: If true, use a symbolic multi-vector rather than the multi-vector passed in for tuning purposes.

            These replace OutVec.
	    - bool symmoutvec: If true, use a symbolic vector rather than the vector passed in for tuning purposes.
	    - bool symmoutmultivec: If true, use a symbolic multi-vector rather than the multi-vector passed in for tuning purposes.

	    - bool tune: If true, have OSKI tune moderately rather than using the number of calls passed in.
	    - bool tuneaggressive: If true, have OSKI tune aggressively rather than using the number of calls passed in.
	*/
	int SetHintMultiply(bool TransA,
			    double Alpha,
			    const Epetra_OskiMultiVector& InVec,
			    double Beta,
			    const Epetra_OskiMultiVector& OutVec,
			    int NumCalls,
			    const Teuchos::ParameterList& List);
	
	//! Workload hints for computing a triangular solve used by OskiTuneMat to optimize the data structure storage, and the routine to compute the calculation.
	/*! In parallel the routine uses symbolic vectors.  This is done for two reasons.  Doing this
            saves on data allocation and potentially communication overhead.  For a matrix-vector
            routine there should be no advantage to having the actual vector, as its size must be the same
            as a matrix dimension.  For a matrix-multivector routine there could be gains from knowing the
            number of vectors in the multi-vector. However, OSKI does not perform multi-vector optimizations,
            so there is no need to add the overhead.

            \param Trans (In) If Trans = true then the transpose of the matrix will be used in
	           computing the product.
	    \param Alpha (In) A scalar constant used to scale InVec.
	    \param Vector (In) The vector being used in the solve and to store the solution.
	    \param NumCalls (In) The number of times the operation is called or the tuning level wanted.
	    \param List (In) Used for denoting the use of a symbolic vectors, as well as for
		   level of aggressive tuning if either NumCalls not
		   known or to be overridden.  Options are shown below. It should be noted that
		   by using these options the associated vector or NumCalls becomes invalid.
	    \return Stores the workload hint in the matrix if the operation is valid.  If the
	       	    operation is not valid an error code is returned.
	
            Options that can be passed to the List are presented below.  They are: "<type> <option name>: <description of purpose>".  The
available hints are grouped by section, and only one hint from each section can be true for a given matrix.

	    These replace Vector.
	    - bool symmvec: If true, use a symbolic vector rather than the vector passed in for tuning purposes.
	    - bool symmmultivec: If true, use a symbolic multi-vector rather than the multi-vector passed in for tuning purposes.

	    - bool tune: If true, have OSKI tune moderately rather than using the number of calls passed in.
	    - bool tuneaggressive: If true, have OSKI tune aggressively rather than using the number of calls passed in.
	*/
	int SetHintSolve(bool TransA,
		         double Alpha,
		         const Epetra_OskiMultiVector& Vector,
		   	 int NumCalls,
		   	 const Teuchos::ParameterList& List);
	
	//! Workload hints for computing a two matrix-vector multiplies that are composed used by OskiTuneMat to optimize the data structure storage, and the routine to compute the calculation.
	/*! In parallel the routine uses symbolic vectors.  This is done for two reasons.  Doing this
            saves on data allocation and potentially communication overhead.  For a matrix-vector
            routine there should be no advantage to having the actual vector, as its size must be the same
            as a matrix dimension.  For a matrix-multivector routine there could be gains from knowing the
            number of vectors in the multi-vector. However, OSKI does not perform multi-vector optimizations,
            so there is no need to add the overhead.
            \param ATA (In) If ATA = true then this^T*this*x will be computed otherwise this*this^T*x will be.
	    \param Alpha (In) A scalar constant used to scale InVec.
	    \param InVec (In) The vector the matrix is multiplied by or whether it is a single vector or multi-vector.
	    \param Beta  (In) A scalar constant used to scale OutVec.
	    \param OutVec (In) The vector where the calculation result is stored or whether it is a single vector or multi-vector.
	    \param Intermediate (In) The vector where result of the first product can be stored
		   or whether it is a single vector or multi-vector.  If this quantity is NULL
		   then the intermediate product is not stored.
	    \param NumCalls (In) The number of times the operation is called or the tuning level wanted.
	    \param List (In) Used for denoting the use of symbolic vectors for InVec,
		   OutVec and Intermediate, along with the level of aggressive tuning if either NumCalls not
		   known or to be overridden.  Options are shown below. It should be noted that
		   by using these options the associated vector or NumCalls becomes invalid.
	    \return Stores the workload hint in the matrix if the operation is valid.  If the
	       	    operation is not valid an error code is returned.
	
            Options that can be passed to the List are presented below.  They are: "<type> <option name>: <description of purpose>".  The
available hints are grouped by section, and only one hint from each section can be true for a given matrix.

	    These replace InVec.
	    - bool symminvec: If true, use a symbolic vector rather than the vector passed in for tuning purposes.
	    - bool symminmultivec: If true, use a symbolic multi-vector rather than the multi-vector passed in for tuning purposes.

            These replace OutVec.
	    - bool symmoutvec: If true, use a symbolic vector rather than the vector passed in for tuning purposes.
	    - bool symmoutmultivec: If true, use a symbolic multi-vector rather than the multi-vector passed in for tuning purposes.

            These replace Intermediate.
	    - bool symmintervec: If true, use a symbolic vector rather than the vector passed in for tuning purposes.
	    - bool symmintermultivec: If true, use a symbolic multi-vector rather than the multi-vector passed in for tuning purposes.

	    - bool tune: If true, have OSKI tune moderately rather than using the number of calls passed in.
	    - bool tuneaggressive: If true, have OSKI tune aggressively rather than using the number of calls passed in.
	*/
 	int SetHintMatTransMatMultiply(bool ATA,
			               double Alpha,
			   	       const Epetra_OskiMultiVector& InVec,
			   	       double Beta,
				       const Epetra_OskiMultiVector& OutVec,
				       const Epetra_OskiMultiVector& Intermediate,
			   	       int NumCalls,
			   	       const Teuchos::ParameterList& List);

	//! Workload hints for computing two matrix-vector multiplies used by OskiTuneMat to optimize the data structure storage, and the routine to compute the calculation.
	/*! In parallel the routine uses symbolic vectors.  This is done for two reasons.  Doing this
            saves on data allocation and potentially communication overhead.  For a matrix-vector
            routine there should be no advantage to having the actual vector, as its size must be the same
            as a matrix dimension.  For a matrix-multivector routine there could be gains from knowing the
            number of vectors in the multi-vector. However, OSKI does not perform multi-vector optimizations,
            so there is no need to add the overhead.

            \param Trans (In) If Trans = true then the transpose of the matrix will be used in
	           computing the product.
	    \param Alpha (In) A scalar constant used to scale InVec.
	    \param InVec (In) The vector the matrix is multiplied by or whether it is a single vector or multi-vector.
	    \param Beta  (In) A scalar constant used to scale OutVec.
	    \param OutVec (In) The vector where the calculation result is stored or whether it is a single vector or multi-vector.
	    \param Omega (In) A scalar constant used to scale InVec2.
	    \param InVec2 (In) The vector the matrix is multiplied by or whether it is a single vector or multi-vector.
	    \param Zeta  (In) A scalar constant used to scale OutVec2.
	    \param OutVec2 (In) The vector where the calculation result is stored or whether it is a single vector or multi-vector.
	    \param NumCalls (In) The number of times the operation is called or the tuning level wanted.
	    \param List (In) Used for denoting the use of symbolic vectors for both InVec
		   and OutVec, as well as for level of aggressive tuning if either NumCalls not
		   known or to be overridden.  Options are shown below. It should be noted that
		   by using these options the associated vector or NumCalls becomes invalid.
	    \return Stores the workload hint in the matrix if the operation is valid.  If the
	       	    operation is not valid an error code is returned.
	
            Options that can be passed to the List are presented below.  They are: "<type> <option name>: <description of purpose>".  The
available hints are grouped by section, and only one hint from each section can be true for a given matrix.

	    These replace InVec.
	    - bool symminvec: If true, use a symbolic vector rather than the vector passed in for tuning purposes.
	    - bool symminmultivec: If true, use a symbolic multi-vector rather than the multi-vector passed in for tuning purposes.

            These replace OutVec.
	    - bool symmoutvec: If true, use a symbolic vector rather than the vector passed in for tuning purposes.
	    - bool symmoutmultivec: If true, use a symbolic multi-vector rather than the multi-vector passed in for tuning purposes.

	    These replace InVec2.
	    - bool symminvec2: If true, use a symbolic vector rather than the vector passed in for tuning purposes.
	    - bool symminmultivec2: If true, use a symbolic multi-vector rather than the multi-vector passed in for tuning purposes.

            These replace OutVec2.
	    - bool symmoutvec2: If true, use a symbolic vector rather than the vector passed in for tuning purposes.
	    - bool symmoutmultivec2: If true, use a symbolic multi-vector rather than the multi-vector passed in for tuning purposes.

	    - bool tune: If true, have OSKI tune moderately rather than using the number of calls passed in.
	    - bool tuneaggressive: If true, have OSKI tune aggressively rather than using the number of calls passed in.
	*/
	int SetHintMultiplyAndMatTransMultiply(bool TransA,
			          	       double Alpha,
			   	  	       const Epetra_OskiMultiVector& InVec,
			   	  	       double Beta,
			   	  	       const Epetra_OskiMultiVector& OutVec,
			          	       double Omega,
			   	  	       const Epetra_OskiMultiVector& InVec2,
			   	  	       double Zeta,
			   	  	       const Epetra_OskiMultiVector& OutVec2,
			   	  	       int NumCalls,
			   	  	       const Teuchos::ParameterList& List);

	//! Workload hints for computing a matrix-vector multiply performed Power times used by OskiTuneMat to optimize the data structure storage and the routine to compute the calculation.
	/*! In parallel the routine uses symbolic vectors.  This is done for two reasons.  Doing this
            saves on data allocation and potentially communication overhead.  For a matrix-vector
            routine there should be no advantage to having the actual vector, as its size must be the same
            as a matrix dimension.  For a matrix-multivector routine there could be gains from knowing the
            number of vectors in the multi-vector. However, OSKI does not perform multi-vector optimizations,
            so there is no need to add the overhead.
            \param Trans (In) If Trans = true then the transpose of the matrix will be used in
	           computing the product.
	    \param Alpha (In) A scalar constant used to scale InVec.
	    \param InVec (In) The vector the matrix is multiplied by or whether it is a single vector or multi-vector.
	    \param Beta  (In) A scalar constant used to scale OutVec.
	    \param OutVec (In) The vector where the calculation result is stored or whether it is a single vector or multi-vector.
	    \param Intermediate (In) The multi-vector where result of the first product can be stored
		   or whether it is a single vector or multi-vector.  If this quantity is NULL
		   then the intermediate product is not stored.
	    \param Power (In) The power to raise the matrix to in the calculation.
	    \param NumCalls (In) The number of times the operation is called or the tuning level wanted.
	    \param List (In) Used for denoting the use of symbolic vectors for both InVec
		   and OutVec, as well as for level of aggressive tuning if either NumCalls not
		   known or to be overridden.  Options are shown below. It should be noted that
		   by using these options the associated vector or NumCalls becomes invalid.
	    \return Stores the workload hint in the matrix if the operation is valid.  If the
	       	    operation is not valid an error code is returned.
	
            Options that can be passed to the List are presented below.  They are: "<type> <option name>: <description of purpose>".  The
available hints are grouped by section, and only one hint from each section can be true for a given matrix.

	    These replace InVec.
	    - bool symminvec: If true, use a symbolic vector rather than the vector passed in for tuning purposes.
	    - bool symminmultivec: If true, use a symbolic multi-vector rather than the multi-vector passed in for tuning purposes.

            These replace OutVec.
	    - bool symmoutvec: If true, use a symbolic vector rather than the vector passed in for tuning purposes.
	    - bool symmoutmultivec: If true, use a symbolic multi-vector rather than the multi-vector passed in for tuning purposes.

            This replaces Intermediate.
	    - bool symmintermultivec: If true, use a symbolic multi-vector rather than the multi-vector passed in for tuning purposes.

	    - bool tune: If true, have OSKI tune moderately rather than using the number of calls passed in.
	    - bool tuneaggressive: If true, have OSKI tune aggressively rather than using the number of calls passed in.
	*/
	int SetHintPowMultiply(bool TransA,
			       double Alpha,
			       const Epetra_OskiMultiVector& InVec,
			       double Beta,
			       const Epetra_OskiMultiVector& OutVec,
			       const Epetra_OskiMultiVector& Intermediate,
			       int Power,
			       int NumCalls,
			       const Teuchos::ParameterList& List);

	//! Tunes the matrix multiply if its deemed profitable.
	/*! The routine tunes based upon user provided hints if given.  If hints are not given the
	    tuning is performed based on expected future workload for the calculation.
	    \return On success returns a non-negative status code of the transformations
		    performed.  On failure an error code is returned.
	*/
	int TuneMatrix();
	//@}	

	//! @name Data Structure Transformation Methods
	//@{
	//! Returns 1 if the matrix has been reordered by tuning, and 0 if it has not been.
	int IsMatrixTransformed() const;

	//! Returns the transformed version of InMat if InMat has been transformed.  If InMat has not been transformed then the return will equal InMat.
	const Epetra_OskiMatrix& ViewTransformedMat() const;

	//! Returns a read only row/left permutation of the Matrix.
	const Epetra_OskiPermutation& ViewRowPermutation() const;
	
	//! Returns a read only column/right permutation of the Matrix.
	const Epetra_OskiPermutation& ViewColumnPermutation() const;
	
	//! Returns a string holding the transformations performed on the matrix when it was tuned.
	/*! \return Upon success returns a newly-allocated string that stores the
		    transformations applied to the matrix during tuning.  NULL is returned
		    upon an error.  It is the users responsibility to deallocate the returned
                    string.
	*/
	char* GetMatrixTransforms() const;

	//! Replaces the current data structure of the matrix with the one specified in Transforms.
	/*! If a previously tuned copy of the Matrix existed it is now replaced by one
	    specified in Transforms.
	    \param Transforms (In) A string that holds the transformations to be applied
		   to the matrix.  If Transforms is NULL or the empty string then no
	           changes to the data structure are made.
	    \return If the transformation was successful 0 is returned.  Otherwise an error
		    code is returned.
	*/
	int ApplyMatrixTransforms(const char* Transforms);
	//@}

 protected:

 private:
	const Epetra_CrsMatrix* Epetra_View_;
	oski_matrix_t A_tunable_;
	bool Copy_Created_;
};
#endif /* EPETRA_OSKIMATRIX_H */
