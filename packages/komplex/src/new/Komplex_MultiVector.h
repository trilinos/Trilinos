//@HEADER
/*
************************************************************************

              Komplex: Complex Linear Solver Package 
                Copyright (2002) Sandia Corporation

Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
license for use of this work by or on behalf of the U.S. Government.

This library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as
published by the Free Software Foundation; either version 2.1 of the
License, or (at your option) any later version.
 
This library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.
 
You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
USA
Questions? Contact Michael A. Heroux (maherou@sandia.gov) 

************************************************************************
*/
//@HEADER

//! Komplex_MultiVector: A class for constructing and using dense complex multi-vectors, vectors 
    and matrices, stored in equivalent real formulation, in parallel.

/*! The Komplex_MultiVector class enables the construction and use of complex-valued, 
  double-precision dense vectors, multi-vectors,
  and matrices in a distributed memory environment.  The dimensions and distribution of the dense
  multi-vectors is determined in part by a Epetra_Comm object, a Epetra_Map (or Epetra_LocalMap
  or Epetra_BlockMap) and the number of vectors passed to the constructors described below.

  There are several concepts that are important for understanding the Komplex_MultiVector class:

  <ul>
  <li>  Multi-vectors, Vectors and Matrices.  
  <ul>
  <li> Vector - A list of real-valued, double-precision numbers.  Also a multi-vector with one vector.
  <li> Multi-Vector - A collection of one or more vectors, all having the same length and distribution.
  <li> (Dense) Matrix - A special form of multi-vector such that stride in memory between any 
  two consecutive vectors in the multi-vector is the same for all vectors.  This is identical
  to a two-dimensional array in Fortran and plays an important part in high performance
  computations.
  </ul>
  <li> Distributed Global vs. Replicated Local.
  <ul>
  <li> Distributed Global Multi-vectors - In most instances, a multi-vector will be partitioned
  across multiple memory images associated with multiple processors.  In this case, there is 
  a unique copy of each element and elements are spread across all processors specified by 
  the Epetra_Comm communicator.
  <li> Replicated Local Multi-vectors - Some algorithms use multi-vectors that are too small to
  be distributed across all processors, for instance, the Hessenberg matrix in a GMRES
  computation.  In other cases, such as with block iterative methods,  block dot product 
  functions produce small dense matrices that are required by all processors.  Replicated 
  local multi-vectors handle these types of situation.
  </ul>
  <li> Multi-vector Functions vs. Dense Matrix Functions.
  <ul>
  <li> Multi-vector functions - These functions operate simultaneously but independently
  on each vector in the multi-vector and produce individual results for each vector.
  <li> Dense matrix functions - These functions operate on the multi-vector as a matrix, 
  providing access to selected dense BLAS and LAPACK operations.
  </ul>
  </ul>

  <b>Constructing Komplex_MultiVectors</b>

  Except for the basic constructor, general constructor, and copy constructor, Komplex_MultiVector constructors
  have two data access modes:
  <ol>
  <li> Copy mode - Allocates memory and makes a copy of the user-provided data. In this case, the
  user data is not needed after construction.
  <li> View mode - Creates a "view" of the user data. In this case, the
  user data is required to remain intact for the life of the multi-vector.
  </ol>

  \warning View mode is \e extremely dangerous from a data hiding perspective.
  Therefore, we strongly encourage users to develop code using Copy mode first and 
  only use the View mode in a secondary optimization phase.

  All Komplex_MultiVector constructors require a map argument that describes the layout of elements
  on the parallel machine.  Specifically, 
  \c map is a Epetra_Map, Epetra_LocalMap or Epetra_BlockMap object describing the desired
  memory layout for the multi-vector.

  There are five different Komplex_MultiVector constructors: 
  <ul>
  <li> Default - All values are zero.
  <li> General - Create a new Komplex_MultiVector from two Epetra_MultiVectors.
  <li> Copy - Copy an existing Komplex_MultiVector.
  <li> Copy or make view of a list of vectors from another Komplex_MultiVector object.
  <li> Copy or make view of a range of vectors from another Komplex_MultiVector object.
  </ul>

  <b>Vector, Matrix and Utility Functions</b>

  Once a Komplex_MultiVector is constructed, a variety of mathematical functions can be applied to
  the individual vectors.  Specifically:
  <ul>
  <li> Vector Updates.
  <li> \e p Norms.
  </ul>

  In addition, a matrix-matrix multiply function supports a variety of operations on any viable
  combination of global distributed and local replicated multi-vectors using calls to DGEMM, a
  high performance kernel for matrix operations.  In the near future we will add support for calls
  to other selected BLAS and LAPACK functions.

  \warning A Epetra_Map, Epetra_LocalMap or Epetra_BlockMap object is required for all 
  Komplex_MultiVector constructors.
*/

//==========================================================================
class Komplex_MultiVector: public Epetra_MultiVector {

 public:

  //@{ \name Constructors/destructors.
  //! Basic Komplex_MultiVector constuctor.
  /*! Creates a Komplex_MultiVector object and, by default, fills with zero values.  

  \param In 
         Map - A Epetra_LocalMap, Epetra_Map or Epetra_BlockMap.
  \warning Note that, because Epetra_LocalMap
  derives from Epetra_Map and Epetra_Map derives from Epetra_BlockMap, this constructor works
  for all three types of Epetra map classes.
  \param In 
         NumVectors - Number of vectors in multi-vector.
  \param In
         zeroOut - If <tt>true</tt> then the allocated memory will be zeroed
                   out initially.  If <tt>false</tt> then this memory will not
                   be touched which can be significantly faster.
  \return Pointer to a Komplex_MultiVector.
  */
  Komplex_MultiVector(const Epetra_BlockMap & Map, int NumVectors, bool zeroOut = true);

 //! General Komplex_MultiVector constructor.
 /*!
	\param In
		 Map - A Epetra_LocalMap, Epetra_Map, or Epetra_BlockMap
	\param In
		 Br - A Epetra_MultiVector containing the real parts of the complex multi-vector
	\param In
		 Bi - A Epetra_MultiVector containing the imaginary parts of the complex multi-vector
  */
  Komplex_MultiVector(const Epetra_BlockMap & Map, const Epetra_MultiVector & Br, 
			    const Epetra_MultiVector & Bi);

  //! Komplex_MultiVector copy constructor.
  Komplex_MultiVector(const Komplex_MultiVector& Source);
  
  //! Set multi-vector values from list of vectors in an existing Komplex_MultiVector.
  /*!
    \param In 
	     Epetra_DataAccess - Enumerated type set to Copy or View.
    \param In
    	     Source - An existing fully constructed Komplex_MultiVector.
    \param In
	     Indices - Integer list of the vectors to copy.  
    \param In 
	     NumVectors - Number of vectors in multi-vector.
    \return Integer error code, set to 0 if successful.
    See Detailed Description section for further discussion.
  */
  Komplex_MultiVector(Epetra_DataAccess CV, const Komplex_MultiVector & Source, 
			    int * Indices, int NumVectors);

  //! Set multi-vector values from range of vectors in an existing Komplex_MultiVector.
  /*!
    \param In 
	     Epetra_DataAccess - Enumerated type set to Copy or View.
    \param In
	     Source - An existing fully constructed Komplex_MultiVector.
    \param In
	     StartIndex - First of the vectors to copy.  
    \param In 
	     NumVectors - Number of vectors in multi-vector.
    \return Integer error code, set to 0 if successful.
    See Detailed Description section for further discussion.
  */
  Komplex_MultiVector(Epetra_DataAccess CV, const Komplex_MultiVector& Source, 
			    int StartIndex, int NumVectors);
  
  //! Komplex_MultiVector destructor.  
  virtual ~Komplex_MultiVector();
  //@}

  //@{ \name Post-construction modification routines.

  //! Replace current value at the specified (GlobalRow, VectorIndex) location with ScalarValue.
  /*!
    Replaces the  existing value for a single entry in the multivector.  The
    specified global row must correspond to a GID owned by the map of the multivector on the
    calling processor.  In other words, this method does not perform cross-processor communication.

    If the map associated with this multivector is an Epetra_BlockMap, only the first point entry associated
    with the global row will be modified.  To modify a different point entry, use the other version of
    this method.

    \param In
           GlobalRow - Row of Multivector to modify in global index space.
    \param In
	     VectorIndex - Vector within MultiVector to modify.
    \param In
	     ScalarValue - Value to replace the existing value.
    \return Integer error code, set to 0 if successful, set to 1 if GlobalRow not associated with calling processor
    set to -1 if VectorIndex >= NumVectors().
  */
  int ReplaceGlobalValue(int GlobalRow, int VectorIndex, double ScalarValue);

  //! Replace current value at the specified (GlobalBlockRow, BlockRowOffset, VectorIndex) location with ScalarValue.
  /*!
    Replaces the existing value for a single entry in the multivector.  The
    specified global block row and block row offset 
    must correspond to a GID owned by the map of the multivector on the
    calling processor.  In other words, this method does not perform cross-processor communication.

    \param In
	     GlobalBlockRow - BlockRow of Multivector to modify in global index space.
    \param In
	     BlockRowOffset - Offset into BlockRow of Multivector to modify in global index space.
    \param In
	     VectorIndex - Vector within MultiVector to modify.
    \param In
	     ScalarValue - Value to replace the existing value.
    \return Integer error code, set to 0 if successful, set to 1 if GlobalRow not associated with calling processor
    set to -1 if VectorIndex >= NumVectors(), set to -2 if BlockRowOffset is out-of-range.
  */
  int ReplaceGlobalValue(int GlobalBlockRow, int BlockRowOffset, int VectorIndex, double ScalarValue);

  //! Adds ScalarValue to existing value at the specified (GlobalRow, VectorIndex) location.
  /*!
    Sums the given value into the existing value for a single entry in the multivector.  The
    specified global row must correspond to a GID owned by the map of the multivector on the
    calling processor.  In other words, this method does not perform cross-processor communication.

    If the map associated with this multivector is an Epetra_BlockMap, only the first point entry associated
    with the global row will be modified.  To modify a different point entry, use the other version of
    this method.

    \param In
	     GlobalRow - Row of Multivector to modify in global index space.
    \param In
	     VectorIndex - Vector within MultiVector to modify.
    \param In
	     ScalarValue - Value to add to existing value.
    \return Integer error code, set to 0 if successful, set to 1 if GlobalRow not associated with calling processor
    set to -1 if VectorIndex >= NumVectors().
  */
  int SumIntoGlobalValue(int GlobalRow, int VectorIndex, double ScalarValue);

  //! Adds ScalarValue to existing value at the specified (GlobalBlockRow, BlockRowOffset, VectorIndex) location.
  /*!
    Sums the given value into the existing value for a single entry in the multivector.  The
    specified global block row and block row offset 
    must correspond to a GID owned by the map of the multivector on the
    calling processor.  In other words, this method does not perform cross-processor communication.

    \param In
	     GlobalBlockRow - BlockRow of Multivector to modify in global index space.
    \param In
	     BlockRowOffset - Offset into BlockRow of Multivector to modify in global index space.
    \param In
	     VectorIndex - Vector within MultiVector to modify.
    \param In
	     ScalarValue - Value to add to existing value.
    \return Integer error code, set to 0 if successful, set to 1 if GlobalRow not associated with calling processor
    set to -1 if VectorIndex >= NumVectors(), set to -2 if BlockRowOffset is out-of-range.
  */
  int SumIntoGlobalValue(int GlobalBlockRow, int BlockRowOffset, int VectorIndex, double ScalarValue);

  //! Replace current value  at the specified (MyRow, VectorIndex) location with ScalarValue.
  /*!
    Replaces the existing value for a single entry in the multivector.  The
    specified local row must correspond to a GID owned by the map of the multivector on the
    calling processor.  In other words, this method does not perform cross-processor communication.

    This method is intended for use with vectors based on an Epetra_Map.  If used 
    on a vector based on a non-trivial Epetra_BlockMap, this will update only block 
    row 0, i.e. 

    Komplex_MultiVector::ReplaceMyValue  (  MyRow,  VectorIndex,  ScalarValue )  is 
    equivalent to:  
    Komplex_MultiVector::ReplaceMyValue  (  0, MyRow,  VectorIndex,  ScalarValue )


    \param In
 	     MyRow - Row of Multivector to modify in local index space.
    \param In
	     VectorIndex - Vector within MultiVector to modify.
    \param In
	     ScalarValue - Value to replace the existing value.
    \return Integer error code, set to 0 if successful, set to 1 if MyRow not associated with calling processor
    set to -1 if VectorIndex >= NumVectors().
  */
  int ReplaceMyValue(int MyRow, int VectorIndex, double ScalarValue);

  //! Replace current value at the specified (MyBlockRow, BlockRowOffset, VectorIndex) location with ScalarValue.
  /*!
    Replaces the existing value for a single entry in the multivector.  The
    specified local block row and block row offset 
    must correspond to a GID owned by the map of the multivector on the
    calling processor.  In other words, this method does not perform cross-processor communication.

    \param In
	     MyBlockRow - BlockRow of Multivector to modify in local index space.
    \param In
	     BlockRowOffset - Offset into BlockRow of Multivector to modify in local index space.
    \param In
	     VectorIndex - Vector within MultiVector to modify.
    \param In
	     ScalarValue - Value to replace the existing value.
    \return Integer error code, set to 0 if successful, set to 1 if MyRow not associated with calling processor
    set to -1 if VectorIndex >= NumVectors(), set to -2 if BlockRowOffset is out-of-range.
  */
  int ReplaceMyValue(int MyBlockRow, int BlockRowOffset, int VectorIndex, double ScalarValue);

  //! Adds ScalarValue to existing value at the specified (MyRow, VectorIndex) location.
  /*!
    Sums the given value into the existing value for a single entry in the multivector.  The
    specified local row must correspond to a GID owned by the map of the multivector on the
    calling processor.  In other words, this method does not perform cross-processor communication.

    If the map associated with this multivector is an Epetra_BlockMap, only the first point entry associated
    with the local row will be modified.  To modify a different point entry, use the other version of
    this method.

    \param In
	     MyRow - Row of Multivector to modify in local index space.
    \param In
	     VectorIndex - Vector within MultiVector to modify.
    \param In
	     ScalarValue - Value to add to existing value.
    \return Integer error code, set to 0 if successful, set to 1 if MyRow not associated with calling processor
    set to -1 if VectorIndex >= NumVectors().
  */
  int SumIntoMyValue(int MyRow, int VectorIndex, double ScalarValue);

  //! Adds ScalarValue to existing value at the specified (MyBlockRow, BlockRowOffset, VectorIndex) location.
  /*!
    Sums the given value into the existing value for a single entry in the multivector.  The
    specified local block row and block row offset 
    must correspond to a GID owned by the map of the multivector on the
    calling processor.  In other words, this method does not perform cross-processor communication.

    \param In
	     MyBlockRow - BlockRow of Multivector to modify in local index space.
    \param In
	     BlockRowOffset - Offset into BlockRow of Multivector to modify in local index space.
    \param In
	     VectorIndex - Vector within MultiVector to modify.
    \param In
	     ScalarValue - Value to add to existing value.
    \return Integer error code, set to 0 if successful, set to 1 if MyRow not associated with calling processor
    set to -1 if VectorIndex >= NumVectors(), set to -2 if BlockRowOffset is out-of-range.
  */
  int SumIntoMyValue(int MyBlockRow, int BlockRowOffset, int VectorIndex, double ScalarValue);
  //@}

  //@{ \name Mathematical methods.

  //! Scale the current values of a multi-vector, \e this = ScalarValue*\e this.
  /*!
    \param In
	     ScalarValue - Scale value.
    \param Out
	     \e This - Multi-vector with scaled values.
    \return Integer error code, set to 0 if successful.
  */
  int Scale(double ScalarValue);

  //! Replace multi-vector values with scaled values of A, \e this = ScalarA*A.
  /*!
    \param In
	     ScalarA - Scale value.
    \param In
	     A - Multi-vector to copy.
    \param Out
	     \e This - Multi-vector with values overwritten by scaled values of A.

    \return Integer error code, set to 0 if successful.
  */
  int Scale(double ScalarA, const Komplex_MultiVector & A);

  //! Compute 1-norm of each vector in multi-vector.
  /*!
    \param Out
	     Result - Result[i] contains 1-norm of ith vector.
    \return Integer error code, set to 0 if successful.
  */
  int Norm1(double * Result) const;

  //! Compute 2-norm of each vector in multi-vector.
  /*!
    \param Out
	     Result - Result[i] contains 2-norm of ith vector.
    \return Integer error code, set to 0 if successful.
  */
  int Norm2(double * Result) const;

  //! Compute Inf-norm of each vector in multi-vector.
  /*!
    \param Out
	     Result - Result[i] contains Inf-norm of ith vector.
    \return Integer error code, set to 0 if successful.
  */
  int NormInf(double * Result) const;

  //! Matrix-Matrix multiplication, \e this = ScalarThis*\e this + ScalarAB*A*B.
  /*! This function performs a variety of matrix-matrix multiply operations, interpreting
    the Komplex_MultiVectors (\e this-aka C , A and B) as 2D matrices.  Variations are due to
    the fact that A, B and C can be local replicated or global distributed
    Komplex_MultiVectors and that we may or may not operate with the transpose of 
    A and B.  Possible cases are:
    \verbatim

    Total of 32 case (2^5).
    Num
    OPERATIONS                        case  Notes
    1) C(local) = A^X(local) * B^X(local)  4   (X=Transpose or Not, No comm needed) 
    2) C(local) = A^T(distr) * B  (distr)  1   (2D dot product, replicate C)
    3) C(distr) = A  (distr) * B^X(local)  2   (2D vector update, no comm needed)

    Note that the following operations are not meaningful for 
    1D distributions:

    1) C(local) = A^T(distr) * B^T(distr)  1
    2) C(local) = A  (distr) * B^X(distr)  2
    3) C(distr) = A^X(local) * B^X(local)  4
    4) C(distr) = A^X(local) * B^X(distr)  4
    5) C(distr) = A^T(distr) * B^X(local)  2
    6) C(local) = A^X(distr) * B^X(local)  4
    7) C(distr) = A^X(distr) * B^X(local)  4
    8) C(local) = A^X(local) * B^X(distr)  4

    \endverbatim

    \param In
	     TransA - Operate with the transpose of A if = 'T', else no transpose if = 'N'.
    \param In
	     TransB - Operate with the transpose of B if = 'T', else no transpose if = 'N'.
    \param In
	     ScalarAB - Scalar to multiply with A*B.
    \param In
	     A - Multi-vector.
    \param In
	     B - Multi-vector.
    \param In
	     ScalarThis - Scalar to multiply with \e this.
    \return Integer error code, set to 0 if successful.

    \warning {Each multi-vector A, B and \e this is checked if it has constant stride using the
    ConstantStride() query function.  If it does not have constant stride, a temporary
    copy is made and used for the computation.  This activity is transparent to the user,
    except that there is memory and computation overhead.  All temporary space is deleted
    prior to exit.}	 
  */
  int Multiply(char TransA, char TransB, double ScalarAB, 
	       const Komplex_MultiVector & A, const Komplex_MultiVector & B,
	       double ScalarThis);
  
  //! Multiply a Komplex_MultiVector with another, element-by-element.
  /*! This function supports diagonal matrix multiply.  A is usually a single vector
    while B and \e this may have one or more columns.  Note that B and \e this must
    have the same shape.  A can be one vector or have the same shape as B.  The actual
    computation is \e this = ScalarThis * \e this + ScalarAB * B @ A where @ denotes element-wise
    multiplication.
  */
  int Multiply(double ScalarAB, const Komplex_MultiVector & A, const Komplex_MultiVector & B,
	       double ScalarThis);
  //@}

  //@{ \name Overloaded operators

  //! = Operator.
  /*!
    \param In
	     A - Komplex_MultiVector to copy.
    \return Komplex_MultiVector.
  */
  Komplex_MultiVector & operator = (const Komplex_MultiVector & Source);

  //! Vector access function.
  /*!
    \return Pointer to the array of doubles containing the local values of the ith vector in the multi-vector.
  */
  double*& operator [] (int i);

  //! Vector access function.
  /*!
    \return Pointer to the array of doubles containing the local values of the ith vector in the multi-vector.
  */
  double * const & operator [] (int i) const;

  //! Vector access function.
  /*!
    \return A Komplex_Vector pointer to the ith vector in the multi-vector.
  */
  Komplex_Vector * & operator () (int i);

  //! Vector access function.
  /*!
    \return A Komplex_Vector pointer to the ith vector in the multi-vector.
  */
  const Komplex_Vector * & operator () (int i) const;
  //@}

  //@{ \name Attribute access functions
  
  //! Returns the number of vectors in the multi-vector.
  int NumVectors() const;

  //! Returns the local vector length on the calling processor of vectors in the multi-vector.
  int MyLength() const;

  //! Returns the global vector length of vectors in the multi-vector.
  int GlobalLength() const;

  //! Returns the stride between  vectors in the multi-vector (only meaningful if ConstantStride() is true).
  int Stride() const;
  
  //! Returns true if this multi-vector has constant stride between vectors.
  bool ConstantStride() const;
  //@}

  //! Replace map, only if new map has same point-structure as current map.
  /*  Return 0 if map is replaced, -1 if not.
  */
  int ReplaceMap(const Epetra_BlockMap & map);

  //@{ \name I/O methods

  //! Print method
  void Print(ostream & os) const;
  //@}

 protected:

 private:
 bool ConstantStride_;
 int Stride_;
 int GlobalLength_;
 int MyLength_;
 int NumVectors_;
 int IndexBase_;
};