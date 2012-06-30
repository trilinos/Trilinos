//@HEADER
// ***********************************************************************
// 
//                Komplex: Complex Linear Solver Package
//                 Copyright (2002) Sandia Corporation
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
//@HEADER

#ifndef KOMPLEX_MULTIVECTOR_H
#define KOMPLEX_MULTIVECTOR_H

class Epetra_BlockMap;
class Epetra_Vector;

#include "Komplex_KForms.hpp"
#include "Epetra_MultiVector.h"
#include "Epetra_DataAccess.h"
#include "Komplex_Ordering.hpp"
#include "Komplex_KForms.hpp"


//! Komplex_MultiVector: A class for constructing and using equivalent real formulations of dense complex multivectors, vectors and matrices in parallel.

/*! The Komplex_MultiVector class enables the construction and use of equivalent real formulations of complex-valued,
  double-precision dense vectors, multivectors,
  and matrices in a distributed memory environment.  The dimensions and distribution of the dense
  multivectors is determined by the Komplex_MultiVector object(s) as described below.
  
  There are several concepts that important for understanding the Komplex_MultiVector class:
  
  <ul>
  <li>  Vectors, Multivectors and Matrices.  
  <ul>
  <li> Vector - A list of real-valued, double-precision numbers.  Also a multivector with one vector.
  <li> Multivector - A collection of one or more vectors, all having the same length and distribution.
  <li> (Dense) Matrix - A special form of multivector such that stride in memory between any 
  two consecutive vectors in the multivector is the same for all vectors.  This is identical
  to a two-dimensional array in Fortran and plays an important part in high performance
  computations.
  </ul>
  <li> Distributed Global vs. Replicated Local.
  <ul>
  <li> Distributed Global Multivectors - In most instances, a multivector will be partitioned
  across multiple memory images associated with multiple processors.  In this case, there is 
  a unique copy of each element, and elements are spread across all processors specified by 
  the Epetra_Comm communicator.
  <li> Replicated Local Multivectors - Some algorithms use multivectors that are too small to
  be distributed across all processors; for instance, the Hessenberg matrix in a GMRES
  computation.  In other cases, such as with block iterative methods,  block dot product 
  functions produce small
  dense matrices that are required by all processors.  Replicated local multivectors handle
  these types of situation.
  </ul>
  <li> Multivector Functions vs. Dense Matrix Functions.
  <ul>
  <li> Multivector Functions - These functions operate simultaneously but independently
  on each vector in the multivector and produce individual results for each vector.
  <li> Dense matrix Functions - These functions operate on the multivector as a matrix, 
  providing access to selected dense BLAS and LAPACK operations.
  </ul>
  </ul>
  
  <b>Constructing Komplex_MultiVectors</b>
  
  Except for the basic constructor and copy constructor, Komplex_MultiVector constructors
  have two data access modes:
  <ol>
  <li> Copy mode - Allocates memory and makes a copy of the user-provided data. In this case, the
  user data is not needed after construction.
  <li> View mode - Creates a "view" of the user data. In this case, the
  user data is required to remain intact for the life of the multivector.
  </ol>
  
  \warning View mode is \e extremely dangerous from a data hiding perspective.
  Therefore, we strongly encourage users to develop code using Copy mode first and 
  only use the View mode in a secondary optimization phase.
  
  There are ten different Komplex_MultiVector constructors:
  <ul>
  <li> Basic - All values are zero.
  <li> Copy - Copy an existing multivector.
  <li> Copy from or make view of two-dimensional Fortran style array, with the real and imaginary parts interleaved.
  <li> Copy from or make view of two two-dimensional Fortran style arrays, representing the real and imaginary parts.
  <li> Copy from or make view of an array of pointers, with the real and imaginary parts interleaved.
  <li> Copy from or make view of two arrays of pointers, representing the real and imaginary parts.
  <li> Copy or make view of a list of vectors from another Komplex_MultiVector object.
  <li> Copy or make view of a range of vectors from another Komplex_MultiVector object.
  <li> Copy of make view of a Epetra_MultiVector object, with the real and imaginary parts interleaved.
  <li> Copy or make view of two Epetra_MultiVector objects, representing the real and imaginary parts.
  </ul>
  
  <b>Extracting Data from Komplex_MultiVectors</b>
  
  Once a Komplex_MultiVector is constructed, it is possible to view it as an Epetra_MultiVector.
  
  <b>Vector, Matrix and Utility Functions</b>
  
  Once a Komplex_MultiVector is constructed, a variety of mathematical functions can be applied to
  the individual vectors.  Specifically:
  <ul>
  <li> Dot Products.
  <li> Vector Updates.
  <li> \e p Norms, treating the vector either as real-valued or complex-valued.
  <li> Weighted Norms, treating the vector either as real-valued or complex-valued.
  <li> Minimum, Maximum and Average Values.
  </ul>
  
  In addition, a matrix-matrix multiply function supports a variety of operations on any viable
  combination of global distributed and local replicated multivectors using calls to DGEMM, a
  high performance kernel for matrix operations.  In the near future we will add support for calls
  to other selected BLAS and LAPACK functions.
  
  <b> Counting Floating Point Operations </b>
  
  Each Komplex_MultiVector object keeps track of the number
  of \e serial floating point operations performed using the specified object as the \e this argument
  to the function.  The Flops() function returns this number as a double precision number.  Using this 
  information, in conjunction with the Epetra_Time class, one can get accurate parallel performance
  numbers.  The ResetFlops() function resets the floating point counter.
  
*/

//==========================================================================
class Komplex_MultiVector {
  
public:
  
  //@{ \name Constructors/destructors.
  //! Basic Komplex_MultiVector constuctor.
  /*! Creates a Komplex_MultiVector object and, by default, fills with zero values.     
    \param Map (In) A Epetra_LocalMap, Epetra_Map or Epetra_BlockMap.
    \param NumVectors (In) Number of vectors in multivector.
    \param RHS (In) If true, then the multivector will be treated as a right-hand side multivector; if false, as a left-hand side multivector.
    \param zeroOut (In) If <tt>true</tt> then the allocated memory will be zeroed
    out initially.  If <tt>false</tt> then this memory will not
    be touched, which can be significantly faster.
    \param KForm (In) The Komplex_KForms to use for this multivector; by default, it is set to K1.

    \return Pointer to a Komplex_MultiVector.

    \warning Note that, because Epetra_LocalMap
    derives from Epetra_Map and Epetra_Map derives from Epetra_BlockMap, this constructor works
    for all three types of Epetra map classes. 
  */
  Komplex_MultiVector(const Epetra_BlockMap& Map, int NumVectors, bool RHS, bool zeroOut = true, Komplex_KForms KForm = K1);
  
  //! Komplex_MultiVector copy constructor.  
  Komplex_MultiVector(const Komplex_MultiVector& Source);
  
  //! Set multivector values from two-dimensional array.
  /*!
    \param Epetra_DataAccess (In) Enumerated type set to Copy or View.
    \param Map (In) A Epetra_LocalMap, Epetra_Map or Epetra_BlockMap.
    \param A (In) Pointer to an array of double precision numbers.  The first vector starts at A.
    The second vector starts at A+MyLDA, the third at A+2*MyLDA, and so on.
    \param MyLDA (In) The "Leading Dimension" or stride between vectors in memory.

    \warning This value refers to the stride on the calling processor.  Thus it is a
    local quantity, not a global quantity.
    
    \param NumVectors (In) Number of vectors in multivector.
    \param RHS (In) If true, then the multivector will be treated as a right-hand side multivector; if false, as a left-hand side multivector.
    \param KForm (In) The Komplex_KForms to use for this multivector; by default, it is set to K1.
    
    \return Integer error code, set to 0 if successful.
    
    See Detailed Description section for further discussion.
  */
  Komplex_MultiVector(Epetra_DataAccess CV, const Epetra_BlockMap& Map, 
		          double* A, int MyLDA, int NumVectors, bool RHS, Komplex_KForms KForm = K1);
  
  //! Set multivector values from two two-dimensional arrays.
  /*!
    \param Epetra_DataAccess (In) Enumerated type set to Copy or View.
    \param Map (In) A Epetra_LocalMap, Epetra_Map or Epetra_BlockMap.
    \param Real (In) Pointer to an array of double precision numbers representing the real parts.
    The first vector starts at Real.  The second vector starts at Real+MyLDA, the third at Real+2*MyLDA, and so on.
    \param Imag (In) Pointer to an array of double precision numbers representing the imaginary parts.
    The first vector starts at Imag.  The second vector starts at Imag+MyLDA, the third at Imag+2*MyLDA, and so on.
    \param MyLDA (In) The "Leading Dimension" or stride between vectors in memory.

    \warning This value refers to the stride on the calling processor.  Thus it is a 
    local quantity, not a global quantity.

    \param NumVectors (In) Number of vectors in multivector.
    \param RHS (In) If true, then the multivector will be treated as a right-hand side multivector; if false, as a left-hand side multivector.
    \param KForm (In) The Komplex_KForms to use for this multivector; by default, it is set to K1.
    
    \return Integer error code, set to 0 if successful.
    
    See Detailed Description section for further discussion.  (YES OR NO????????)
  */
  Komplex_MultiVector(Epetra_DataAccess CV, const Epetra_BlockMap& Map, double* Real, 
			    double* Imag, int MyLDA, int NumVectors, bool RHS, Komplex_KForms KForm = K1);
   
  //! Set multivector values from array of pointers.
  /*!
    \param Epetra_DataAccess (In) Enumerated type set to Copy or View.
    \param Map (In) A Epetra_LocalMap, Epetra_Map or Epetra_BlockMap.
    \param ArrayOfPointers (In) An array of pointers such that ArrayOfPointers[i] points to the memory
    location containing the ith vector to be copied.
    \param NumVectors (In) Number of vectors in multivector.
    \param RHS (In) If true, then the multivector will be treated as a right-hand side multivector; if false, as a left-hand side multivector.
    \param KForm (In) The Komplex_KForms to use for this multivector; by default, it is set to K1.
    
    \return Integer error code, set to 0 if successful.
    
    See Detailed Description section for further discussion.
  */
  Komplex_MultiVector(Epetra_DataAccess CV, const Epetra_BlockMap& Map, double** ArrayOfPointers, 
			    int NumVectors, bool RHS, Komplex_KForms KForm = K1);
  
  //! Set multivector values from two arrays of pointers, representing the real and imaginary parts.
  /*!
    \param Epetra_DataAccess (In) Enumerated type set to Copy or View.
    \param Map (In) A Epetra_LocalMap, Epetra_Map or Epetra_BlockMap.
    \param AOPReal (In) An array of pointers such that AOPReal[i] points to the memory
    location containing the ith real vector to be copied.
    \param AOPImag (In) An array of pointers such that AOPImag[i] points to the memory
    location containing the ith imaginary vector to be copied.
    \param NumVectors (In) Number of vectors in multivector.
    \param RHS (In) If true, then the multivector will be treated as a right-hand side multivector; if false, as a left-hand side multivector.
    \param KForm (In) The Komplex_KForms to use for this multivector; by default, it is set to K1.
    
    \return Integer error code, set to 0 if successful.
    
    See Detailed Description section for further discussion.
  */
  Komplex_MultiVector(Epetra_DataAccess CV, const Epetra_BlockMap& Map, double** AOPReal, 
			    double** AOPImag, int NumVectors, bool RHS, Komplex_KForms KForm = K1);
  
  //! Set multivector values from list of vectors in an existing Epetra_MultiVector.
  /*!
    \param Epetra_DataAccess (In) Enumerated type set to Copy or View.
    \param Source (In) An existing fully constructed Epetra_MultiVector.
    \param Indices (In) Integer list of the vectors to copy.  
    \param NumVectors (In) Number of vectors in multivector.
    \param RHS (In) If true, then the multivector will be treated as a right-hand side multivector; if false, as a left-hand side multivector.
    \param KForm (In) The Komplex_KForms to use for this multivector; by default, it is set to K1.
    
    \return Integer error code, set to 0 if successful.
    
    See Detailed Description section for further discussion.
  */
  Komplex_MultiVector(Epetra_DataAccess CV, const Epetra_MultiVector& Source, 
		          int* Indices, int NumVectors, bool RHS, Komplex_KForms KForm = K1);
  
  //! Set multivector values from range of vectors in an existing Epetra_MultiVector.
  /*!
    \param Epetra_DataAccess (In) Enumerated type set to Copy or View.
    \param Source (In) An existing fully constructed Epetra_MultiVector.
    \param StartIndex (In) First of the vectors to copy.  
    \param NumVectors (In) Number of vectors in multivector.
    \param RHS (In) If true, then the multivector will be treated as a right-hand side multivector; if false, as a left-hand side multivector.
    \param KForm (In) The Komplex_KForms to use for this multivector; by default, it is set to K1.
    
    \return Integer error code, set to 0 if successful.
    
    See Detailed Description section for further discussion.
  */
  Komplex_MultiVector(Epetra_DataAccess CV, const Epetra_MultiVector& Source, int StartIndex, 
		          int NumVectors, bool RHS, Komplex_KForms KForm = K1);
  
  //! Set multivector values from an existing Epetra_MultiVector, with the real and imaginary parts interleaved.
  /*!
    \param Epetra_DataAccess (In) Enumerated type set to Copy of View.
    \param Source (In) An existing fully constructed Epetra_MultiVector.
    \param RHS (In) If true, then the multivector will be treated as a right-hand side multivector; if false, as a left-hand side multivector.
    \param KForm (In) The Komplex_KForms to use for this multivector; by default, it is set to K1.
    
    \return Integer error code, set to 0 if successful.
    
    See Detailed Description secton for further discussion.
  */
  Komplex_MultiVector(Epetra_DataAccess CV, const Epetra_MultiVector& Source, bool RHS, Komplex_KForms KForm = K1);
  
  //! Set multivector values from two Epetra_MultiVectors, one representing the real and the other the imaginary values.
  /*!
    \param Epetra_DataAccess (In) Enumerated type set to Copy or View.
    \param Real (In) An existing fully constructed Epetra_MultiVector representing the real values.
    \param Imag (In) An existing fully constructed Epetra_MultiVector representing the imaginary values.
    \param RHS (In) If true, then the multivector will be treated as a right-hand side multivector; if false, as a left-hand side multivector.
    \param KForm (In) The Komplex_KForms to use for this multivector; by default, it is set to K1.
    
    \return Integer error code, set to 0 if successful.
    
    See Detailed Description section for further discussion.

    \warning Real and Imag must have the same Epetra_Map
  */
  Komplex_MultiVector(Epetra_DataAccess CV, const Epetra_MultiVector& Real, const Epetra_MultiVector& Imag,
                      bool RHS, Komplex_KForms KForm = K1);
  
  //! Komplex_MultiVector destructor.  
  virtual ~Komplex_MultiVector();
  //@}
  
  //@{ \name Post-construction modification routines.
  
  //! Replace current value  at the specified (GlobalRow, VectorIndex) location with ScalarValue.
  /*!
    Replaces the  existing value for a single entry in the multivector.  The
    specified global row must correspond to a GID owned by the map of the multivector on the
    calling processor.  In other words, this method does not perform cross-processor communication.
    
    If the map associated with this multivector is an Epetra_BlockMap, only the first point entry associated
    with the global row will be modified.  To modify a different point entry, use the other version of
    this method
    
    \param GlobalRow (In) Row of multivector to modify in global index space.
    \param VectorIndex (In) Vector within multivector to modify.
    \param ScalarValue (In) Value to add to existing value.
    
    \return Integer error code, set to 0 if successful, set to 1 if GlobalRow is not associated with calling processor,
    set to -1 if VectorIndex >= NumVectors().
  */
  int ReplaceGlobalValue(int GlobalRow, int VectorIndex, double ScalarValue);
    
  //! Add ScalarValue to existing value at the specified (GlobalRow, VectorIndex) location.
  /*!
    Sums the given value into the existing value for a single entry in the multivector.  The
    specified global row must correspond to a GID owned by the map of the multivector on the
    calling processor.  In other words, this method does not perform cross-processor communication.
    
    If the map associated with this multivector is an Epetra_BlockMap, only the first point entry associated
    with the global row will be modified.  To modify a different point entry, use the other version of
    this method
    
    \param GlobalRow (In) Row of multivector to modify in global index space.
    \param VectorIndex (In) Vector within multivector to modify.
    \param ScalarValue (In) Value to add to existing value.
    
    \return Integer error code, set to 0 if successful, set to 1 if GlobalRow is not associated with calling processor,
    set to -1 if VectorIndex >= NumVectors().
  */
  int SumIntoGlobalValue(int GlobalRow, int VectorIndex, double ScalarValue);
    
  //! Replace current value  at the specified (MyRow, VectorIndex) location with ScalarValue.
  /*!
    Replaces the existing value for a single entry in the multivector.  The
    specified local row must correspond to a GID owned by the map of the multivector on the
    calling processor.  In other words, this method does not perform cross-processor communication.
    
    This method is intended for use with vectors based on an Epetra_Map.  If used 
    on a vector based on a non-trivial Epetra_BlockMap, this will update only block 
    row 0, i.e. 
    
    Komplex_MultiVector::ReplaceMyValue  (MyRow, VectorIndex, ScalarValue)  is 
    equivalent to:  
    Komplex_MultiVector::ReplaceMyValue  (0, MyRow, VectorIndex, ScalarValue)
    
    
    \param MyRow (In) Row of multivector to modify in local index space.
    \param VectorIndex (In) Vector within multivector to modify.
    \param ScalarValue (In) Value to add to existing value.
    
    \return Integer error code, set to 0 if successful, set to 1 if MyRow is not associated with calling processor,
    set to -1 if VectorIndex >= NumVectors().
  */
  int ReplaceMyValue(int MyRow, int VectorIndex, double ScalarValue);
  
  //! Add ScalarValue to existing value at the specified (MyRow, VectorIndex) location.
  /*!
    Sums the given value into the existing value for a single entry in the multivector.  The
    specified local row must correspond to a GID owned by the map of the multivector on the
    calling processor.  In other words, this method does not perform cross-processor communication.
    
    If the map associated with this multivector is an Epetra_BlockMap, only the first point entry associated
    with the local row will be modified.  To modify a different point entry, use the other version of
    this method
    
    \param MyRow (In) Row of multivector to modify in local index space.
    \param VectorIndex (In) Vector within multivector to modify.
    \param ScalarValue (In) Value to add to existing value.
    
    \return Integer error code, set to 0 if successful, set to 1 if MyRow is not associated with calling processor,
    set to -1 if VectorIndex >= NumVectors().
  */
  int SumIntoMyValue(int MyRow, int VectorIndex, double ScalarValue);
    
  //! Initialize all values in a multivector with constant value.
  /*!
    \param ScalarConstant (In) Value to use.
    
    \return Integer error code, set to 0 if successful.
  */
  int PutScalar (double ScalarConstant);
  
  //! Set multivector values to random numbers.
  /*! This uses the random number generator provided by Epetra_Util. 
    The multivector values will be set to random values on the interval (-1.0, 1.0).
    
    \return Integer error code, set to 0 if successful.   
  */
  int Random();

  //! Creates a map one-half or twice the size of the existing map, allowing for return of the real parts, the imaginary parts, or an interleaved multivector when given the opposite in the constructor.
  void CreateOtherMap();
  
  //@}
  
  //@{ \name Mathematical methods.
  
  //! Computes dot product of each corresponding pair of vectors.
  /*!
    \param A (In) Multivector to be used with the \e this multivector.
    \param Result (Out) Result[i] will contain the ith dot product result.
    
    \return Integer error code, set to 0 if successful.
  */
  int Dot(const Komplex_MultiVector& A, double* Result) const;
  
  //! Puts element-wise absolute values of input multivector in target.
  /*!
    \param A (In) Input multivector.

    \post On return, \e this will contain the absolute values of the entries of A.
    
    \return Integer error code, set to 0 if successful.
    
    Note:  It is possible to use the same argument for A and \e this.
  */
  int Abs(const Komplex_MultiVector& A);
  
  //! Puts element-wise reciprocal values of input multivector in target.
  /*!
    \param A (In) Input multivector.

    \post On return, \e this will contain the element-wise reciprocal values of the entries of A.
    
    \return Integer error code, set to 0 if successful.  Returns 2 if some entry
    is too small, but not zero.  Returns 1 if some entry is zero.
    
    Note:  It is possible to use the same argument for A and \e this.  Also, 
    if a given value of A is smaller than Epetra_DoubleMin (defined in Epetra_Epetra.h),
    but nonzero, then the return code is 2.  If an entry is zero, the return code
    is 1.  However, in all cases the reciprocal value is still used, even
    if a NaN is the result.
  */
  int Reciprocal(const Komplex_MultiVector& A);
  
  //! Scale the current values of a multivector, \e this = ScalarValue*\e this.
  /*!
    \param ScalarValue (In) Scale value.

    \post On return, \e this will contain the scaled values.
    
    \return Integer error code, set to 0 if successful.
  */
  int Scale(double ScalarValue);
  
  //! Replace multivector values with scaled values of A, \e this = ScalarA*A.
  /*!
    \param ScalarA (In) Scale value.
    \param A (In) Multivector to copy.

    \post On return, \e this will contain the scaled values of A.
    
    \return Integer error code, set to 0 if successful.
  */
  int Scale(double ScalarA, const Komplex_MultiVector& A);
  
  //! Update multivector values with scaled values of A, \e this = ScalarThis*\e this + ScalarA*A.
  /*!
    \param ScalarA (In) Scale value for A.
    \param A (In) Multivector to add.
    \param ScalarThis (In) Scale value for \e this.

    \post On return, \e this will contain the updated values.
    
    \return Integer error code, set to 0 if successful.
  */
  int Update(double ScalarA, const Komplex_MultiVector& A, double ScalarThis);
  
  //! Update multivector with scaled values of A and B, \e this = ScalarThis*\e this + ScalarA*A + ScalarB*B.
  /*!
    \param ScalarA (In) Scale value for A.
    \param A (In) Multivector to add.
    \param ScalarB (In) Scale value for B.
    \param B (In) Multivector to add.
    \param ScalarThis (In) Scale value for \e this.

    \post On return, \e this will contain the updated values.
    
    \return Integer error code, set to 0 if successful.
  */
  int Update(double ScalarA, const Komplex_MultiVector& A, 
	     double ScalarB, const Komplex_MultiVector& B, double ScalarThis);
  
  //! Compute the 1-norm of each vector in multivector.
  /*!
    \param Result (Out) Result[i] contains the 1-norm of the ith vector.
    
    \return Integer error code, set to 0 if successful.
  */
  int Norm1(double* Result) const;
  
  //! Compute the 1-norm of each vector, regarded as a complex vector, in multivector.
  /*!
    \param Result (Out) Result[i] contains the 1-norm of the ith vector.
    
    \return Integer error code, set to 0 if successful.
  */
  int ComplexNorm1(double* Result) const;
  
  //! Compute the 2-norm of each vector in multivector.
  /*!
    \param Result (Out) Result[i] contains the 2-norm of the ith vector.
    
    \return Integer error code, set to 0 if successful.
  */
  int Norm2(double* Result) const;
  
  //! Compute the 2-norm of each vector, regarded as a complex vector, in multivector.
  /*!
    \param Result (Out) Result[i] contains the 2-norm of the ith vector.
    
    \return Integer error code, set to 0 if successful.
  */
  int ComplexNorm2(double* Result) const;
  
  //! Compute the Inf-norm of each vector in multivector.
  /*!
    \param Result (Out) Result[i] contains the Inf-norm of the ith vector.
    
    \return Integer error code, set to 0 if successful.
  */
  int NormInf(double* Result) const;
  
  //! Compute the Inf-norm of each vector, regarded as a comnplex vector, in multivector.
  /*!
    \param Result (Out) Result[i] contains the Inf-norm of the ith vector.
    
    \return Integer error code, set to 0 if successful.
  */
  int ComplexNormInf(double* Result) const;
  
  //! Compute the Weighted 2-norm (RMS Norm) of each vector in multivector.
  /*!
    \param Weights (In) Multivector of weights.  If Weights contains a single vector,
    that vector will be used as the weights for all vectors of \e this.  Otherwise,
    Weights should have the same number of vectors as \e this.
    \param Result (Out) Result[i] contains the weighted 2-norm of ith vector.  Specifically
    if we denote the ith vector in the multivector by \f$x\f$, and the ith weight
    vector by \f$w\f$ and let j represent the jth entry of each vector, on return
    Result[i] will contain the following result:
    \f[\sqrt{(1/n)\sum_{j=1}^n(x_j/w_j)^2}\f],
    where \f$n\f$ is the global length of the vectors.
    
    \return Integer error code, set to 0 if successful.
  */
  int NormWeighted(const Epetra_MultiVector& Weights, double* Result) const;
    
  //! Compute minimum value of each vector in multivector.
  /*!
    \param Result (Out) Result[i] contains the minimum value of the ith vector.
    
    \return Integer error code, set to 0 if successful.
  */
  int MinValue(double* Result) const;
  
  //! Compute maximum value of each vector in multivector.
  /*!
    \param Result (Out) Result[i] contains the maximum value of the ith vector.
    
    \return Integer error code, set to 0 if successful.
  */
  int MaxValue(double* Result) const;
  
  //! Compute mean (average) value of each vector in multivector.
  /*!
    \param Result (Out) Result[i] contains the mean value of the ith vector.
    
    \return Integer error code, set to 0 if successful.
  */
  int MeanValue(double* Result) const;
  
  //@}
  
  //@{ \name Random number utilities
  
  //! Set seed for Random function.
  /*!
    \param Seed (In) Should be an integer on the interval (0, 2^31-1).
    
    \return Integer error code, set to 0 if successful.
  */
  int SetSeed(unsigned int Seed);
  
  //! Get seed from Random function.
  /*!
    \return Current random number seed.
  */
  unsigned int Seed() const;
  
  //@}
  
  //@{ \name Overloaded operators
  
  //! = Operator.
  /*!
    \param A (In) Komplex_MultiVector to copy.
    
    \return Komplex_MultiVector.
  */
  Komplex_MultiVector & operator = (const Komplex_MultiVector& Source);
  
  // Local element access functions
  
  //! Vector access function.
  /*!
    \return Pointer to the array of doubles containing the local values of the ith vector in the multivector.
  */
  double*& operator [] (int i);
  //! Vector access function.
  /*!
    \return Pointer to the array of doubles containing the local values of the ith vector in the multivector.
  */
  double * const & operator [] (int i) const;
  
  //! Vector access function.
  /*!
    \return A Komplex_Vector pointer to the ith vector in the multivector.
  */
  //######### Komplex_Vector * & operator () (int i);
  //! Vector access function.
  /*!
    \return A Komplex_Vector pointer to the ith vector in the multivector.
  */
  //######### const Komplex_Vector * & operator () (int i) const;
  
  //! Conversion to Epetra_MultiVector.
  /*!
    \return A Epetra_MultiVector pointer to the multivector.
  */
  Epetra_MultiVector * EpetraMultiVector() const;
  
  //! Conversion of real parts to Epetra_MultiVector.
  /*!
    \return A Epetra_MultiVector with the real parts.

  */
  Epetra_MultiVector * RealMultiVector() const;

  //! Conversion of imaginary parts to Epetra_MultiVector.
  /*!
    \return A Epetra_MultiVector with the imaginary parts.
  */
  Epetra_MultiVector * ImagMultiVector() const;

  //! Single vector conversion to Epetra_Vector.
  /*!
    \return A Epetra_Vector pointer to the ith vector in the multivector.
  */
  Epetra_Vector * EpetraVector(int index) const;

  //! Single vector conversion to Epetra_Vector, including only the real values.
  /*!
    \return A Epetra_Vector containing the real values of the ith vector.

  */
  Epetra_Vector * RealVector(int index) const;

  //! Single vector conversion to Epetra_Vector, including only the imaginary values.
  /*! 
    \return A Epetra_Vector containing the imaginary values of the ith vector.
  */
  Epetra_Vector * ImagVector(int index) const;

  //! Vector access function.
  /*!
    \return Pointer to the array of doubles containing the real parts of the local values of the ith vector in the multivector.
  */
  double * & RealValues(int i) const;

  //! Vector access function.
  /*!
    \return Pointer to the array of doubles containing the imaginary parts of the local values of the ith vector in the multivector.
  */
  double * & ImagValues(int i) const;

  //@}
  
  //@{ \name Attribute access functions
  
  //! Returns the number of vectors in the multivector.
  int NumVectors() const;
  
  //! Returns the local vector length on the calling processor of vectors in the multivector.
  int MyLength() const;
  
  //! Returns the global vector length of vectors in the multivector.
  int GlobalLength() const;
  
  //! Returns the current K form.
  Komplex_KForms KForm() const;
  
  //! Returns true if this is a right-hand side multivector, false otherwise.
  bool RHS() const;

  //! Switches the current K form.
  /*!
    \return Integer error code, set to 0 if successful.
  */
  int SwitchKForm(Komplex_KForms NewKForm);
  //@}
  
  //! Replace map, only if new map has same point-structure as current map.
  /*!  
    \return 0 if map is replaced, -1 if not.
  */
  int ReplaceMap(const Epetra_BlockMap& map);
  
  //@{ \name I/O methods
  
  //! Print method
  virtual void Print(ostream& os) const;
  //@}

  void PrintReal(ostream& os);
  void PrintImag(ostream& os);
  
  
protected:
  //if IsOneMatrix_ is true, then this is called for creating the "half" map
  void CreateHalfMap(); 
  //if IsOneMatrix_ is false, then this is called for creating the "double" map
  void CreateDoubleMap(); 
  
private:
  
  Epetra_MultiVector* Real_;
  Epetra_MultiVector* Imag_;
  Komplex_Ordering* Ordering_;
  bool IsOneObject_; 
  bool RHS_;
  Epetra_BlockMap* OtherMap_;
  mutable double* TempDouble_;
  
};

#endif /* KOMPLEX_MULTIVECTOR_H */
