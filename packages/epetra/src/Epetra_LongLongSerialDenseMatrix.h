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

#ifndef EPETRA_LONGLONGSERIALDENSEMATRIX_H
#define EPETRA_LONGLONGSERIALDENSEMATRIX_H

#include "Epetra_Object.h" 

//! Epetra_LongLongSerialDenseMatrix: A class for constructing and using general dense integer matrices.

/*! The Epetra_LongLongSerialDenseMatrix class enables the construction and use of integer-valued, general
    dense matrices. 

The Epetra_LongLongSerialDenseMatrix class is intended to provide very basic support for dense rectangular matrices.


<b>Constructing Epetra_LongLongSerialDenseMatrix Objects</b>

There are four Epetra_LongLongSerialDenseMatrix constructors.  The first constructs a zero-sized object which should be made
to appropriate length using the Shape() or Reshape() functions and then filled with the [] or () operators. 
The second constructs an object sized to the dimensions specified, which should be filled with the [] or () operators.
The third is a constructor that accepts user
data as a 2D array, and the fourth is a copy constructor. The third constructor has
two data access modes (specified by the Epetra_DataAccess argument):
<ol>
  <li> Copy mode - Allocates memory and makes a copy of the user-provided data. In this case, the
       user data is not needed after construction.
  <li> View mode - Creates a "view" of the user data. In this case, the
       user data is required to remain intact for the life of the object.
</ol>

\warning View mode is \e extremely dangerous from a data hiding perspective.
Therefore, we strongly encourage users to develop code using Copy mode first and 
only use the View mode in a secondary optimization phase.

Epetra_LongLongSerialDenseMatrix constructors will throw an exception if an error occurrs.  
These exceptions will alway be negative integer values as follows:
<ol>
  <li> -1  Invalid dimension specified.
  <li> -2  Shape returned non-zero.
  <li> -3  Null pointer specified for user's data.
  <li> -99 Internal Epetra_LongLongSerialDenseMatrix error.  Contact developer.
</ol>

Other Epetra_LongLongSerialDenseMatrix functions that do not return an integer error code
(such as operators () and [] ) will throw an exception if an error occurrs. 
These exceptions will be integer values as follows:
<ol>
  <li> -1  Invalid row specified.
  <li> -2  Invalid column specified.
	<li> -5  Invalid assignment (type mismatch).
  <li> -99 Internal Epetra_LongLongSerialDenseMatrix error.  Contact developer.
</ol>


b<b>Extracting Data from Epetra_LongLongSerialDenseMatrix Objects</b>

Once a Epetra_LongLongSerialDenseMatrix is constructed, it is possible to view the data via access functions.

\warning Use of these access functions cam be \e extremely dangerous from a data hiding perspective.


<b>Vector and Utility Functions</b>

Once a Epetra_LongLongSerialDenseMatrix is constructed, several mathematical functions can be applied to
the object.  Specifically:
<ul>
  <li> Multiplication.
  <li> Norms.
</ul>


*/


//=========================================================================
class EPETRA_LIB_DLL_EXPORT Epetra_LongLongSerialDenseMatrix : public Epetra_Object {

  public:
  
    //! @name Constructor/Destructor Methods
  //@{ 
  //! Default constructor; defines a zero size object.
  /*!
    Epetra_LongLongSerialDenseMatrix objects defined by the default constructor should be sized with the 
    Shape() or Reshape functions.  
    Values should be defined by using the [] or () operators.
   */
  Epetra_LongLongSerialDenseMatrix();
  
  //! Shaped constructor; defines a variable-sized object
  /*!
    \param In 
           NumRows - Number of rows in object.
    \param In 
           NumCols - Number of columns in object.

    Epetra_SerialDenseMatrix objects defined by the shaped constructor are already shaped to the
		dimensions given as a parameters. All values are initialized to 0. Calling this constructor 
		is equivalent to using the default constructor, and then calling the Shape function on it.
    Values should be defined by using the [] or () operators.
   */
  Epetra_LongLongSerialDenseMatrix(int NumRows, int NumCols);

  //! Set object values from two-dimensional array.
  /*!
    \param In 
           Epetra_DataAccess - Enumerated type set to Copy or View.
    \param In
           A - Pointer to an array of integer numbers.  The first vector starts at A.
	   The second vector starts at A+LDA, the third at A+2*LDA, and so on.
    \param In
           LDA - The "Leading Dimension", or stride between vectors in memory.
    \param In 
           NumRows - Number of rows in object.
    \param In 
           NumCols - Number of columns in object.

	   See Detailed Description section for further discussion.
  */
  Epetra_LongLongSerialDenseMatrix(Epetra_DataAccess CV, long long* A, int LDA, int NumRows, int NumCols);
  
  //! Epetra_LongLongSerialDenseMatrix copy constructor.
	/*!
		This matrix will take on the data access mode of the Source matrix.
	*/
  Epetra_LongLongSerialDenseMatrix(const Epetra_LongLongSerialDenseMatrix& Source);

  //! Epetra_LongLongSerialDenseMatrix destructor.  
  virtual ~Epetra_LongLongSerialDenseMatrix ();
  //@}

  //! @name Shaping/sizing Methods
  //@{ 
  //! Set dimensions of a Epetra_LongLongSerialDenseMatrix object; init values to zero.
  /*!
    \param In 
           NumRows - Number of rows in object.
    \param In 
           NumCols - Number of columns in object.

	   Allows user to define the dimensions of a Epetra_LongLongSerialDenseMatrix at any point. This function can
	   be called at any point after construction.  Any values that were previously in this object are
	   destroyed and the resized matrix starts off with all zero values.

    \return Integer error code, set to 0 if successful.
  */
  int Shape(int NumRows, int NumCols);
  
  //! Reshape a Epetra_LongLongSerialDenseMatrix object.
  /*!
    \param In 
           NumRows - Number of rows in object.
    \param In 
           NumCols - Number of columns in object.

	   Allows user to define the dimensions of a Epetra_LongLongSerialDenseMatrix at any point. This function can
	   be called at any point after construction.  Any values that were previously in this object are
	   copied into the new shape.  If the new shape is smaller than the original, the upper left portion
	   of the original matrix (the principal submatrix) is copied to the new matrix.

    \return Integer error code, set to 0 if successful.
  */
  int Reshape(int NumRows, int NumCols);
  //@}
  
  //! @name Data Accessor methods
  //@{ 

  //! Computes the 1-Norm of the \e this matrix.
  /*!
    \return Integer error code, set to 0 if successful.
  */
  virtual long long OneNorm();

  //! Computes the Infinity-Norm of the \e this matrix.
  virtual long long InfNorm();

  //! Copy from one matrix to another.
  /*!
    The operator= allows one to copy the values from one existing LongLongSerialDenseMatrix to another.
		The left hand side matrix will take on the data access mode of the right hand side matrix. 

    \return Values of the left hand side matrix are modified by the values of the right hand side matrix.
  */
    Epetra_LongLongSerialDenseMatrix& operator = (const Epetra_LongLongSerialDenseMatrix& Source);

    //! Comparison operator.
    /*! operator== compares two Epetra_LongLongSerialDenseMatrix objects, returns false if sizes are different,
      or if any coefficients differ.
    */
    bool operator==(const Epetra_LongLongSerialDenseMatrix& rhs) const;

    //! Inequality operator
    /*! operator!= simply returns the negation of operator==.
     */
    bool operator!=(const Epetra_LongLongSerialDenseMatrix& rhs) const
    { return !(*this == rhs); }

  //! Element access function.
  /*!
    The parentheses operator returns the element in the ith row and jth column if A(i,j) is
    specified, the expression A[j][i] (note that i and j are reversed) will return the same element.
    Thus, A(i,j) = A[j][i] for all valid i and j.

    \return Element from the specified row and column.

		\warning No bounds checking is done unless Epetra is compiled with HAVE_EPETRA_ARRAY_BOUNDS_CHECK.
  */
    long long& operator () (int RowIndex, int ColIndex);

  //! Element access function.
  /*!
    The parentheses operator returns the element in the ith row and jth column if A(i,j) is
    specified, the expression A[j][i] (note that i and j are reversed) will return the same element.
    Thus, A(i,j) = A[j][i] for all valid i and j.

    \return Element from the specified row and column.

		\warning No bounds checking is done unless Epetra is compiled with HAVE_EPETRA_ARRAY_BOUNDS_CHECK.
  */
    const long long& operator () (int RowIndex, int ColIndex) const;

  //! Column access function.
  /*!
    The parentheses operator returns the element in the ith row and jth column if A(i,j) is
    specified, the expression A[j][i] (note that i and j are reversed) will return the same element.
    Thus, A(i,j) = A[j][i] for all valid i and j.

    \return Pointer to address of specified column.

    \warning No bounds checking can be done for the index i in the expression A[j][i].
		\warning No bounds checking is done unless Epetra is compiled with HAVE_EPETRA_ARRAY_BOUNDS_CHECK.
  */
    long long* operator [] (int ColIndex);

  //! Column access function.
  /*!
    The parentheses operator returns the element in the ith row and jth column if A(i,j) is
    specified, the expression A[j][i] (note that i and j are reversed) will return the same element.
    Thus, A(i,j) = A[j][i] for all valid i and j.

    \return Pointer to address of specified column.

    \warning No bounds checking can be done for the index i in the expression A[j][i].
		\warning No bounds checking is done unless Epetra is compiled with HAVE_EPETRA_ARRAY_BOUNDS_CHECK.
  */
    const long long* operator [] (int ColIndex) const;

  //! Set matrix values to random numbers.
  /*! 
		LongLongSerialDenseMatrix uses the random number generator provided by Epetra_Util.
		The matrix values will be set to random values on the interval (0, 2^31 - 1).

    \return Integer error code, set to 0 if successful.
  */
  int Random();
    
  //! Returns row dimension of system.
  int M() const {return(M_);};

  //! Returns column dimension of system.
  int N() const {return(N_);};

  //! Returns const pointer to the \e this matrix.
  const long long* A() const {return(A_);};

  //! Returns pointer to the \e this matrix.
  long long* A() {return(A_);};

  //! Returns the leading dimension of the \e this matrix.
  int LDA() const {return(LDA_);};

	//! Returns the data access mode of the \e this matrix.
	Epetra_DataAccess CV() const {return(CV_);};
  //@}
  
  //! @name I/O methods
  //@{ 
  //! Print service methods; defines behavior of ostream << operator.
  virtual void Print(ostream& os) const;
  //@}

  //! @name Expert-only unsupported methods
  //@{ 

  //! Reset an existing LongLongSerialDenseMatrix to point to another Matrix.
	/*! Allows an existing LongLongSerialDenseMatrix to become a View of another
		matrix's data, regardless of the DataAccess mode of the Source matrix.
		It is assumed that the Source matrix is an independent matrix, and 
		no checking is done to verify this.

		This is used by Epetra_CrsGraph in the OptimizeStorage method. It is used so that
		an existing (Copy) matrix can be converted to a View. This frees up
		memory that CrsGraph no longer needs.
		
		@param Source The LongLongSerialDenseMatrix this will become a view of.
		
		\return Integer error code, set to 0 if successful, and set to -1 
		if a type mismatch occured.
		
		\warning This method is extremely dangerous and should only be used by experts.
	*/
	
	int MakeViewOf(const Epetra_LongLongSerialDenseMatrix& Source);
	//@}

 protected:

	void CopyMat(long long* Source, int Source_LDA, int NumRows, int NumCols, long long* Target, int Target_LDA);
  void CleanupData();

	Epetra_DataAccess CV_;
	bool A_Copied_;
  int M_;
  int N_;
  int LDA_;
  long long* A_;

};

// inlined definitions of op() and op[]
//=========================================================================
inline long long& Epetra_LongLongSerialDenseMatrix::operator () (int RowIndex, int ColIndex) {
#ifdef HAVE_EPETRA_ARRAY_BOUNDS_CHECK
  if(RowIndex >= M_ || RowIndex < 0) 
		throw ReportError("Row index = " + toString(RowIndex) + 
											" Out of Range 0 - " + toString(M_-1),-1);
  if(ColIndex >= N_ || ColIndex < 0) 
		throw ReportError("Column index = " + toString(ColIndex) + 
											" Out of Range 0 - " + toString(N_-1),-2);
#endif
  return(A_[ColIndex*LDA_ + RowIndex]);
}
//=========================================================================
inline const long long& Epetra_LongLongSerialDenseMatrix::operator () (int RowIndex, int ColIndex) const {
#ifdef HAVE_EPETRA_ARRAY_BOUNDS_CHECK
  if(RowIndex >= M_ || RowIndex < 0) 
		throw ReportError("Row index = " + toString(RowIndex) + 
											" Out of Range 0 - " + toString(M_-1),-1);
  if(ColIndex >= N_ || ColIndex < 0) 
		throw ReportError("Column index = " + toString(ColIndex) + 
											" Out of Range 0 - " + toString(N_-1),-2);
#endif
	return(A_[ColIndex * LDA_ + RowIndex]);
}
//=========================================================================
inline long long* Epetra_LongLongSerialDenseMatrix::operator [] (int ColIndex) {
#ifdef HAVE_EPETRA_ARRAY_BOUNDS_CHECK
  if(ColIndex >= N_ || ColIndex < 0) 
		throw ReportError("Column index = " + toString(ColIndex) + 
											" Out of Range 0 - " + toString(N_-1),-2);
#endif
  return(A_+ ColIndex * LDA_);
}
//=========================================================================
inline const long long* Epetra_LongLongSerialDenseMatrix::operator [] (int ColIndex) const {
#ifdef HAVE_EPETRA_ARRAY_BOUNDS_CHECK
  if(ColIndex >= N_ || ColIndex < 0) 
		throw ReportError("Column index = " + toString(ColIndex) + 
											" Out of Range 0 - " + toString(N_-1),-2);
#endif
  return(A_ + ColIndex * LDA_);
}
//=========================================================================

#endif /* EPETRA_LONGLONGSERIALDENSEMATRIX_H */
