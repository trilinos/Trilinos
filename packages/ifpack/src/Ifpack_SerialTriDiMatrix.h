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

#ifndef IFPACK_SERIALTRIDIMATRIX_H
#define IFPACK_SERIALTRIDIMATRIX_H

#if defined(Ifpack_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Ifpack package is deprecated"
#endif
#endif

#include "Epetra_ConfigDefs.h"
#include "Epetra_Object.h"
#include "Epetra_CompObject.h"
#include "Epetra_BLAS.h"

class Epetra_VbrMatrix;

//! Ifpack_SerialTriDiMatrix: A class for constructing and using real double precision general TriDi matrices.

/*! The Ifpack_SerialTriDiMatrix class enables the construction and use of real-valued, general,
    double-precision TriDi matrices.  It is built on the BLAS, and derives from the Epetra_BLAS.

The Ifpack_SerialTriDiMatrix class is intended to provide very basic support for TriDiagonal matrices.


<b>Constructing Ifpack_SerialTriDiMatrix Objects</b>

There are four Ifpack_SerialTriDiMatrix constructors.  The first constructs a zero-sized object which should be made
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

<b>Extracting Data from Ifpack_SerialTriDiMatrix Objects</b>

Once a Ifpack_SerialTriDiMatrix is constructed, it is possible to view the data via access functions.

\warning Use of these access functions cam be \e extremely dangerous from a data hiding perspective.


<b>Vector and Utility Functions</b>

Once a Ifpack_SerialTriDiMatrix is constructed, several mathematical functions can be applied to
the object.  Specifically:
<ul>
  <li> Multiplication.
  <li> Norms.
</ul>

<b>Counting floating point operations </b>
The Ifpack_SerialTriDiMatrix class has Epetra_CompObject as a base class.  Thus, floating point operations
are counted and accumulated in the Epetra_Flop object (if any) that was set using the SetFlopCounter()
method in the Epetra_CompObject base class.

*/


//=========================================================================
class Ifpack_SerialTriDiMatrix : public Epetra_CompObject, public Epetra_Object, public Epetra_BLAS {

  public:

    //! @name Constructor/Destructor Methods
  //@{
  //! Default constructor; defines a zero size object.
  /*!
    Ifpack_SerialTriDiMatrix objects defined by the default constructor should be sized with the
    Shape() or Reshape functions.
    Values should be defined by using the () operators.
   */
  Ifpack_SerialTriDiMatrix(bool set_object_label=true);

  //! Shaped constructor; defines a variable-sized object
  /*!
    \param In
           NumRowCol - Number of rows and columns in object.

    Ifpack_SerialTriDiMatrix objects defined by the shaped constructor are already shaped to the
                dimensions given as a parameters. All values are initialized to 0. Calling this constructor
                is equivalent to using the default constructor, and then calling the Shape function on it.
    Values should be defined by using the [] or () operators.
   */
  Ifpack_SerialTriDiMatrix(int NumRowCol, bool set_object_label=true);

  //! Set object values from two-dimensional array.
  /*!
    \param In
           Epetra_DataAccess - Enumerated type set to Copy or View.
    \param In
           A - Pointer to an array of double precision numbers.
           The
    \param In
           NumRows - Number of rows and columns in object.

           See Detailed Description section for further discussion.
  */
  Ifpack_SerialTriDiMatrix(Epetra_DataAccess CV, double* A_in, int NumRowCol,
                           bool set_object_label=true);

  //! Ifpack_SerialTriDiMatrix copy constructor.

  Ifpack_SerialTriDiMatrix(const Ifpack_SerialTriDiMatrix& Source);

  //! Ifpack_SerialTriDiMatrix destructor.
  virtual ~Ifpack_SerialTriDiMatrix ();
  //@}

  //! @name Shaping/sizing Methods
  //@{
  //! Set dimensions of a Ifpack_SerialTriDiMatrix object; init values to zero.
  /*!
    \param In
           NumRowCol - Number of rows and columns in object.

           Allows user to define the dimensions of a Ifpack_SerialTriDiMatrix at any point. This function can
           be called at any point after construction.  Any values that were previously in this object are
           destroyed and the resized matrix starts off with all zero values.

    \return Integer error code, set to 0 if successful.
  */
  int Shape(int NumRowCol);

  int Reshape(int, int);

  //! @name Mathematical methods
  //@{

  //! Matrix-Matrix multiplication, \e this = ScalarThis*\e this + ScalarAB*A*B.
  /*! This function performs a variety of matrix-matrix multiply operations.

  \param In
         TransA - Operate with the transpose of A if = 'T', else no transpose if = 'N'.
  \param In
         TransB - Operate with the transpose of B if = 'T', else no transpose if = 'N'.

  \param In
         ScalarAB - Scalar to multiply with A*B.
  \param In
         A - TriDi Matrix.
  \param In
         B - TriDi Matrix.
  \param In
         ScalarThis - Scalar to multiply with \e this.

    \return Integer error code, set to 0 if successful.

  */
  int Multiply(char TransA, char TransB, double ScalarAB,
               const Ifpack_SerialTriDiMatrix& A,
               const Ifpack_SerialTriDiMatrix& B,
               double ScalarThis);

  //! Matrix-Vector multiplication, y = A*x, where 'this' == A.
  /* This method is intended to imitate the semantics of the matrix-vector
    multiplication provided by Epetra's sparse matrices. The 'vector' arguments
    are actually matrices; this method will return an error if the
    dimensions of 'x' are not compatible. 'y' will be reshaped if necessary.
  */
  /* int Multiply(bool transA, */
  /*              const Ifpack_SerialTriDiMatrix& x, */
  /*              Ifpack_SerialTriDiMatrix& y); */

  //! Matrix-Matrix multiplication with a symmetric matrix A.
  /*! If SideA = 'L', compute \e this = ScalarThis*\e this + ScalarAB*A*B.
      If SideA = 'R', compute \e this = ScalarThis*\e this + ScalarAB*B*A.

This function performs a variety of matrix-matrix multiply operations.

  \param In
         SideA - Specifies order of A relative to B.

  \param In
         ScalarAB - Scalar to multiply with A*B.
  \param In
         A - Symmetric TriDi Matrix, either upper or lower triangle will be used depending on
         value of A.Upper().
  \param In
         B - TriDi Matrix.
  \param In
         ScalarThis - Scalar to multiply with \e this.

    \return Integer error code, set to 0 if successful.



  \param ScalarA (In) Scalar to multiply with A.

   \return Integer error code, set to 0 if successful.

  */
  int Scale(double ScalarA);

  //! Computes the 1-Norm of the \e this matrix.
  /*!
    \return Integer error code, set to 0 if successful.
  */
  virtual double NormOne() const;

  //! Computes the Infinity-Norm of the \e this matrix.
  virtual double NormInf() const;

  //@}

  //! @name Data Accessor methods
  //@{

  //! Value copy from one matrix to another.
  /*!
    The operator= allows one to copy the values from one existing SerialTriDiMatrix to another, as
    long as there is enough room in the target to hold the source.

    \return Values of the left hand side matrix are modified by the values of the right hand side matrix.
  */
    Ifpack_SerialTriDiMatrix & operator = (const Ifpack_SerialTriDiMatrix& Source);

    //! Comparison operator.
    /*! operator== compares two Ifpack_SerialTriDiMatrix objects, returns false if sizes are different,
      or if any coefficients differ by an amount greater than Epetra_MinDouble.
    */
    bool operator==(const Ifpack_SerialTriDiMatrix& rhs) const;

    //! Inequality operator
    /*! operator!= simply returns the negation of operator==.
     */
    bool operator!=(const Ifpack_SerialTriDiMatrix& rhs) const
    { return !(*this == rhs); }

  //! Add one matrix to another.
  /*!
    The operator+= allows one to add the values from one existin SerialTriDiMatrix to another, as
    long as there is enough room in the target to hold the source.

    \return Values of the left hand side matrix are modified by the addition
    of the values of the right hand side matrix.
  */
    Ifpack_SerialTriDiMatrix & operator += (const Ifpack_SerialTriDiMatrix& Source);

  //! Element access function.
  /*!
    The parentheses operator returns the element in the ith row and jth column if A(i,j) is
    specified,

    \return Element from the specified row and column.

                \warning No bounds checking is done unless Epetra is compiled with HAVE_EPETRA_ARRAY_BOUNDS_CHECK.
  */
            double& operator () (int RowIndex, int ColIndex);

  //! Element access function.
  /*!
    The parentheses operator returns the element in the ith row and jth column if A(i,j) is
    specified, the expression A[j][i] (note that i and j are reversed) will return the same element.
    Thus, A(i,j) = A[j][i] for all valid i and j.

    \return Element from the specified row and column.

                \warning No bounds checking is done unless Epetra is compiled with HAVE_EPETRA_ARRAY_BOUNDS_CHECK.
  */
    const double& operator () (int RowIndex, int ColIndex) const;

  //! Column access function.
  /*!
    The parentheses operator returns the element in the ith row and jth column if A(i,j) is
    specified, the expression A[j][i] (note that i and j are reversed) will return the same element.
    Thus, A(i,j) = A[j][i] for all valid i and j.

    \return Pointer to address of specified column.

    \warning No bounds checking can be done for the index i in the expression A[j][i].
                \warning No bounds checking is done unless Epetra is compiled with HAVE_EPETRA_ARRAY_BOUNDS_CHECK.
  */
    //    double* operator [] (int ColIndex);

  //! Column access function.
  /*!
    The parentheses operator returns the element in the ith row and jth column if A(i,j) is
    specified, the expression A[j][i] (note that i and j are reversed) will return the same element.
    Thus, A(i,j) = A[j][i] for all valid i and j.

    \return Pointer to address of specified column.

    \warning No bounds checking can be done for the index i in the expression A[j][i].
                \warning No bounds checking is done unless Epetra is compiled with HAVE_EPETRA_ARRAY_BOUNDS_CHECK.
  */
    //    const double* operator [] (int ColIndex) const;

  //! Set matrix values to random numbers.
  /*!
                SerialTriDiMatrix uses the random number generator provided by Epetra_Util.
                The matrix values will be set to random values on the interval (-1.0, 1.0).

                \return Integer error code, set to 0 if successful.
  */
  int Random();

  //! Returns column dimension of system.
  int N() const {return(N_);};

  int LDA() const {return(LDA_);};

  //! Returns pointer to the \e this matrix.
  double* A() const {return(A_);};

  //! Returns pointer to the \e this matrix.
  //  double* A() {return(A_);};

  double* DL() { return DL_;};
  double* DL() const { return DL_;};
  double* D() { return D_;};
  double* D() const { return D_;};
  double* DU() { return DU_;};
  double* DU() const { return DU_;};
  double* DU2() { return DU2_;};
  double* DU2() const { return DU2_;};

  //! Returns the data access mode of the \e this matrix.
  Epetra_DataAccess CV() const {return(CV_);};
  //@}

  //! @name I/O methods
  //@{
  //! Print service methods; defines behavior of ostream << operator.
  virtual void Print(std::ostream& os) const;
  //@}

  //! @name Deprecated methods (will be removed in later versions of this class)
  //@{

  //! Computes the 1-Norm of the \e this matrix (identical to NormOne() method).
  /*!
    \return Integer error code, set to 0 if successful.
  */
  virtual double OneNorm() const {return(NormOne());};

  //! Computes the Infinity-Norm of the \e this matrix (identical to NormInf() method).
  virtual double InfNorm() const {return(NormInf());};
  //@}

  //! @name Additional methods to support Ifpack_SerialTriDiOperator interface
  //@{

    //! If set true, transpose of this operator will be applied.
    /*! This flag allows the transpose of the given operator to be used implicitly.  Setting this flag
        affects only the Apply() and ApplyInverse() methods.  If the implementation of this interface
        does not support transpose use, this method should return a value of -1.

    \param In
           UseTranspose -If true, multiply by the transpose of operator, otherwise just use operator.

    \return Integer error code, set to 0 if successful.  Set to -1 if this implementation does not support transpose.
  */
    virtual int SetUseTranspose(bool UseTranspose_in) { UseTranspose_ = UseTranspose_in; return (0); }

    //! Returns the result of a Ifpack_SerialTriDiOperator inverse applied to an Ifpack_SerialTriDiMatrix X in Y.
    /*!
    \param In
           X - A Ifpack_SerialTriDiMatrix to solve for.
    \param Out
           Y -A Ifpack_SerialTriDiMatrix containing result.

    \return Integer error code, set to 0 if successful.

  */
    virtual int ApplyInverse(const Ifpack_SerialTriDiMatrix & X, Ifpack_SerialTriDiMatrix & Y)
    {
      (void)X;//prevents unused variable compiler warning
      (void)Y;
      return (-1);
    }

    //! Returns a character string describing the operator
    virtual const char * Label() const { return Epetra_Object::Label(); }

    //! Returns the current UseTranspose setting.
    virtual bool UseTranspose() const { return UseTranspose_; }

    //! Returns true if the \e this object can provide an approximate Inf-norm, false otherwise.
    virtual bool HasNormInf() const { return true; }

    //! Returns the column dimension of operator
    virtual int RowColDim() const { return N(); }
  //@}

 protected:

  void CopyMat(const double* Source, int NumRowCol,
               double* Target, int NRC2, bool add=false);
  void CleanupData();

  int N_;
  int LDA_;
  bool A_Copied_;
  Epetra_DataAccess CV_;

  //For performance reasons, it's better if Epetra_VbrMatrix can access the
  //A_ members of this class directly without going through an
  //accessor method. Rather than making them public members, we'll make
  //Epetra_VbrMatrix a friend class.

  friend class Epetra_VbrMatrix;

  double* A_;
  double* DL_;
  double* D_;
  double* DU_;
  double* DU2_;

  bool UseTranspose_;
};

// inlined definitions of op() and op[]
//=========================================================================
inline double& Ifpack_SerialTriDiMatrix::operator () (int RowIndex, int ColIndex) {

 int diff = ColIndex - RowIndex;

#ifdef HAVE_EPETRA_ARRAY_BOUNDS_CHECK
 if (ColIndex >= N_ || ColIndex < 0)
                throw ReportError("Column index = " +toString(ColIndex) +
                                  " Out of Range 0 - " + toString(N_-1),-2);
 if (RowIndex >= N_ || RowIndex < 0)
                throw ReportError("Row index = " +toString(RowIndex) +
                                  " Out of Range 0 - " + toString(N_-1),-2);

 if ( diff > 1 || diff < -1 )
   throw ReportError("Row index = " +toString(RowIndex) + " differs from Col_Index " + toString(ColIndex) +
                     " Out of Range -1 to 1",-2);
#endif

 switch (diff) {
 case -1:
   // DL
   return DL_[ColIndex];
   // break; // unreachable
 case 0:
   return D_[ColIndex];
   // break; // unreachable
 case 1:
   return DU_[RowIndex];
   // break; // unreachable
 default:
   throw ReportError("Row index = " +toString(RowIndex) + " differs from Col_Index " + toString(ColIndex) +" Out of Range -1 to 1",1);
   // return D_[0]; // unreachable
 }
 //throw ReportError("Row index = " +toString(RowIndex) + " differs from Col_Index " + toString(ColIndex) + " Out of Range -1 to 1",1); // unreachable
 // return D_[0]; // unreachable
}
//=========================================================================
inline const double& Ifpack_SerialTriDiMatrix::operator () (int RowIndex, int ColIndex) const {
 int diff = ColIndex - RowIndex;

#ifdef HAVE_EPETRA_ARRAY_BOUNDS_CHECK
 if (ColIndex >= N_ || ColIndex < 0)
                throw ReportError("Column index = " +toString(ColIndex) +
                                  " Out of Range 0 - " + toString(N_-1),-2);
 if (RowIndex >= N_ || RowIndex < 0)
                throw ReportError("Row index = " +toString(RowIndex) +
                                  " Out of Range 0 - " + toString(N_-1),-2);
 if ( diff > 1 || diff < -1 )
   throw ReportError("Row index = " +toString(RowIndex) + " differs from Col_Index " + toString(ColIndex) + " Out of Range -1 to 1",-2);
#endif
 switch (diff) {
 case -1:
   // DL
   return DL_[ColIndex];
   // break; // unreachable
 case 0:
   return D_[ColIndex];
   // break; // unreachable
 case 1:
   return DU_[RowIndex];
   // break; // unreachable
 default:
   throw ReportError("Row index = " +toString(RowIndex) + " differs from Col_Index " + toString(ColIndex) + " Out of Range -1 to 1",-2);
   // return D_[0]; // unreachable
 }
 // return D_[0]; // unreachable
}


#endif /* EPETRA_SERIALTRIDIMATRIX_H */
