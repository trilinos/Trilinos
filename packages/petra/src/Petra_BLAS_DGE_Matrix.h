#ifndef _PETRA_BLAS_DGE_MATRIX_H_
#define _PETRA_BLAS_DGE_MATRIX_H_

//! Petra_BLAS_DGE_Matrix: A class for constructing and using real double precision general dense matrices.

/*! The Petra_BLAS_DGE_Matrix class enables the construction and use of real-valued, general, 
    double-precision dense matrices.  It is built on the BLAS, and derives from the Petra_BLAS. 

The Petra_BLAS_DGE_Matrix class is intended to provide very basic support for dense rectangular matrices.


<b>Constructing Petra_BLAS_DGE_Matrix Objects</b>

There are three Petra_BLAS_DGE_Matrix constructors.  The first constructs a zero-sized object which should be made
to appropriate length using the Shape() or Reshape() functions and then filled with the [] or () operators. 
The second is a constructor that accepts user
data as a 2D array, the third is a copy constructor. The second constructor has
two data access modes (specified by the Petra_DataAccess argument):
<ol>
  <li> Copy mode - Allocates memory and makes a copy of the user-provided data. In this case, the
       user data is not needed after construction.
  <li> View mode - Creates a "view" of the user data. In this case, the
       user data is required to remain intact for the life of the object.
</ol>

\warning View mode is \e extremely dangerous from a data hiding perspective.
Therefore, we strongly encourage users to develop code using Copy mode first and 
only use the View mode in a secondary optimization phase.

<b>Extracting Data from Petra_BLAS_DGE_Matrix Objects</b>

Once a Petra_BLAS_DGE_Matrix is constructed, it is possible to view the data via access functions.

\warning Use of these access functions cam be \e extremely dangerous from a data hiding perspective.


<b>Vector and Utility Functions</b>

Once a Petra_BLAS_DGE_Matrix is constructed, several mathematical functions can be applied to
the object.  Specifically:
<ul>
  <li> Multiplication.
  <li> Norms.
</ul>

The final useful function is Flops().  Each Petra_BLAS_DGE_Matrix object keep track of the number
of \e serial floating point operations performed using the specified object as the \e this argument
to the function.  The Flops() function returns this number as a double precision number.  Using this 
information, in conjunction with the Petra_Time class, one can get accurate parallel performance
numbers.


*/
#include "Petra_Object.h" 
#include "Petra_Flops.h"
#include "Petra_BLAS.h"


//=========================================================================
class Petra_BLAS_DGE_Matrix : public Petra_Flops, public Petra_Object, public Petra_BLAS {

  public:
  
  //@{ \name Constructor/Destructor Methods
  //! Default constructor; defines a zero size object.
  /*!
    Petra_BLAS_DGE_Matrix objects defined by the default constructor should be sized with the 
    Shape() or Reshape functions.  
    Values should be defined by using the [] or () operators.
   */
  Petra_BLAS_DGE_Matrix(void);
  
  //! Set object values from two-dimensional array.
  /*!
    \param In 
           Petra_DataAccess - Enumerated type set to Copy or View.
    \param In
           A - Pointer to an array of double precision numbers.  The first vector starts at A.
	   The second vector starts at A+LDA, the third at A+2*LDA, and so on.
    \param In
           LDA - The "Leading Dimension", or stride between vectors in memory.
    \param In 
           NumRows - Number of rows in object.
    \param In 
           NumCols - Number of columns in object.

	   See Detailed Description section for further discussion.
  */
  Petra_BLAS_DGE_Matrix(Petra_DataAccess CV, double *A, int LDA, int NumRows, int NumCols);
  
  //! Petra_BLAS_DGE_Matrix copy constructor.
  
  Petra_BLAS_DGE_Matrix(const Petra_BLAS_DGE_Matrix& Source);

  //! Petra_BLAS_DGE_Matrix destructor.  
  virtual ~Petra_BLAS_DGE_Matrix ();
  //@}

  //@{ \name Shaping/sizing Methods
  //! Set dimensions of a Petra_BLAS_DGE_Matrix object; init values to zero.
  /*!
    \param In 
           NumRows - Number of rows in object.
    \param In 
           NumCols - Number of columns in object.

	   Allows user to define the dimensions of a Petra_BLAS_DGE_Matrix at any point. This function can
	   be called at any point after construction.  Any values that were previously in this object are
	   destroyed and the resized matrix starts off with all zero values.

    \return Integer error code, set to 0 if successful.
  */
  int Shape(int NumRows, int NumCols);
  
  //! Reshape a Petra_BLAS_DGE_Matrix object.
  /*!
    \param In 
           NumRows - Number of rows in object.
    \param In 
           NumCols - Number of columns in object.

	   Allows user to define the dimensions of a Petra_BLAS_DGE_Matrix at any point. This function can
	   be called at any point after construction.  Any values that were previously in this object are
	   copied into the new shape.  If the new shape is smaller than the original, the upper left portion
	   of the original matrix (the principal submatrix) is copied to the new matrix.

    \return Integer error code, set to 0 if successful.
  */
  int Reshape(int NumRows, int NumCols);
  //@}
  
  //@{ \name Mathematical methods

  //! Computes the 1-Norm of the \e this matrix.
  /*!
    \return Integer error code, set to 0 if successful.
  */
  double OneNorm();

  //! Computes the Infinity-Norm of the \e this matrix.
  /*!
    \return Integer error code, set to 0 if successful.
  */
  double InfNorm();

  //! Matrix-Matrix multiplication, \e this = Scalar*\e this + ScalarAB*A*B.
  /*! This function performs a variety of matrix-matrix multiply operations.

  \param In
         TransA - Operate with the transpose of A if = 'T', else no transpose if = 'N'.
  \param In
         TransB - Operate with the transpose of B if = 'T', else no transpose if = 'N'.

  \param In
         ScalarAB - Scalar to multiply with A*B.
  \param In
         A - Dense Matrix.
  \param In
         B - Dense Matrix.
  \param In
         Scalar - Scalar to multiply with \e this.

    \return Integer error code, set to 0 if successful.
	 
  */
  int  Multiply (char TransA, char TransB, double ScalarAB, 
                 const Petra_BLAS_DGE_Matrix& A, 
                 const Petra_BLAS_DGE_Matrix& B,
                 double Scalar );
  //@}

  //@{ \name Data Accessor methods
  //! Element access function.
  /*!
    The parentheses operator returns the element in the ith row and jth column if A(i,j) is
    specified, the expression A[j][i] (note that i and j are reversed) will return the same element.
    Thus, A(i,j) = A[j][i] for all valid i and j.

    \return Element from the specified row and column.
  */
    double& operator () (int RowIndex, int ColIndex);

  //! Element access function.
  /*!
    The parentheses operator returns the element in the ith row and jth column if A(i,j) is
    specified, the expression A[j][i] (note that i and j are reversed) will return the same element.
    Thus, A(i,j) = A[j][i] for all valid i and j.

    \return Element from the specified row and column.
  */
    const double& operator () (int RowIndex, int ColIndex) const;

  //! Column access function.
  /*!
    The parentheses operator returns the element in the ith row and jth column if A(i,j) is
    specified, the expression A[j][i] (note that i and j are reversed) will return the same element.
    Thus, A(i,j) = A[j][i] for all valid i and j.

    \return Pointer to address of specified column.

    \warning No bounds checking can be done for the index i in the expression A[j][i].
  */
    double* operator [] (int ColIndex);

  //! Column access function.
  /*!
    The parentheses operator returns the element in the ith row and jth column if A(i,j) is
    specified, the expression A[j][i] (note that i and j are reversed) will return the same element.
    Thus, A(i,j) = A[j][i] for all valid i and j.

    \return Pointer to address of specified column.

    \warning No bounds checking can be done for the index i in the expression A[j][i].
  */
    const double* operator [] (int ColIndex) const;
    
  //! Returns row dimension of system.
  int M()  const {return(M_);};

  //! Returns column dimension of system.
  int N()  const {return(N_);};

  //! Returns pointer to the \e this matrix.
  double * A()  const {return(A_);};

  //! Returns pointer to the \e this matrix.
  double * A() {return(A_);};

  //! Returns the leading dimension of the \e this matrix.
  int LDA()  const {return(LDA_);};
  //@}
  //! Print service methods; defines behavior of ostream << operator.
  virtual void Print(ostream& os) const;
 protected:

  void CopyMat(double * A, int LDA, int NumRows, int NumCols, double * B, int LDB);
  void DeleteArrays(void);

  int M_;
  int N_;
  int LDA_;
  bool A_Copied_;
  double * A_;


};

#endif /* _PETRA_BLAS_DGE_MATRIX_H_ */
