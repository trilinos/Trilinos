
/* Copyright (2001) Sandia Corportation. Under the terms of Contract 
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this 
 * work by or on behalf of the U.S. Government.  Export of this program
 * may require a license from the United States Government. */


/* NOTICE:  The United States Government is granted for itself and others
 * acting on its behalf a paid-up, nonexclusive, irrevocable worldwide
 * license in ths data to reproduce, prepare derivative works, and
 * perform publicly and display publicly.  Beginning five (5) years from
 * July 25, 2001, the United States Government is granted for itself and
 * others acting on its behalf a paid-up, nonexclusive, irrevocable
 * worldwide license in this data to reproduce, prepare derivative works,
 * distribute copies to the public, perform publicly and display
 * publicly, and to permit others to do so.
 * 
 * NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED STATES DEPARTMENT
 * OF ENERGY, NOR SANDIA CORPORATION, NOR ANY OF THEIR EMPLOYEES, MAKES
 * ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR
 * RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY
 * INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS
 * THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS. */

#ifndef _EPETRA_MULTIVECTOR_H_
#define _EPETRA_MULTIVECTOR_H_

class Epetra_Comm;
class Epetra_BlockMap;
class Epetra_Map;
class Epetra_Import;
class Epetra_Export;
class Epetra_Distributor;
#include "Epetra_DistObject.h"
#include "Epetra_CompObject.h"
#include "Epetra_BLAS.h"

//! Epetra_MultiVector: A class for constructing and using dense multi-vectors, vectors and matrices in parallel.

/*! The Epetra_MultiVector class enables the construction and use of real-valued, 
    double-precision dense vectors, multi-vectors,
    and matrices in a distributed memory environment.  The dimensions and distribution of the dense
    multi-vectors is determined in part by a Epetra_Comm object, a Epetra_Map (or Epetra_LocalMap
    or Epetra_BlockMap) and the number of vectors passed to the constructors described below.

    There are several concepts that important for understanding the Epetra_MultiVector class:

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
       be distributed across all processors, the Hessenberg matrix in a GMRES
       computation.  In other cases, such as with block iterative methods,  block dot product 
       functions produce small
       dense matrices that are required by all processors.  Replicated local multi-vectors handle
       these types of situation.
</ul>
<li> Multi-vector Functions vs. Dense Matrix Functions.
<ul>
  <li> Multi-vector functions - These functions operate simultaneously but independently
       on each vector in the multi-vector and produce individual results for each vector.
  <li> Dense matrix functions - These functions operate on the multi-vector as a matrix, 
       providing access to selected dense BLAS and LAPACK operations.
</ul>
</ul>

<b>Constructing Epetra_MultiVectors</b>

Except for the basic constructor and copy constructor, Epetra_MultiVector constructors
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

All Epetra_MultiVector constructors require a map argument that describes the layout of elements
on the parallel machine.  Specifically, 
\c map is a Epetra_Map, Epetra_LocalMap or Epetra_BlockMap object describing the desired
memory layout for the multi-vector.

There are six different Epetra_MultiVector constructors:
<ul>
  <li> Basic - All values are zero.
  <li> Copy - Copy an existing multi-vector.
  <li> Copy from or make view of two-dimensional Fortran style array.
  <li> Copy from or make view of an array of pointers.
  <li> Copy or make view of a list of vectors from another Epetra_MultiVector object.
  <li> Copy or make view of a range of vectors from another Epetra_MultiVector object.
</ul>

<b>Extracting Data from Epetra_MultiVectors</b>

Once a Epetra_MultiVector is constructed, it is possible to extract a copy of the values or create
a view of them.

\warning ExtractView functions are \e extremely dangerous from a data hiding perspective.
For both ExtractView fuctions, there is a corresponding ExtractCopy function.  We
strongly encourage users to develop code using ExtractCopy functions first and 
only use the ExtractView functions in a secondary optimization phase.

There are four Extract functions:
<ul>
  <li> ExtractCopy - Copy values into a user-provided two-dimensional array.
  <li> ExtractCopy - Copy values into a user-provided array of pointers.
  <li> ExtractView - Set user-provided two-dimensional array parameters 
       to point to Epetra_MultiVector data.
  <li> ExtractView - Set user-provided array of pointer parameters 
       to point to Epetra_MultiVector data.
</ul>

<b>Vector, Matrix and Utility Functions</b>

Once a Epetra_MultiVector is constructed, a variety of mathematical functions can be applied to
the individual vectors.  Specifically:
<ul>
  <li> Dot Products.
  <li> Vector Updates.
  <li> \e p Norms.
  <li> Weighted Norms.
  <li> Minimum, Maximum and Average Values.
</ul>

In addition, a matrix-matrix multiply function supports a variety of operations on any viable
combination of global distributed and local replicated multi-vectors using calls to DGEMM, a
high performance kernel for matrix operations.  In the near future we will add support for calls
to other selected BLAS and LAPACK functions.

<b> Counting Floating Point Operations </b>

Each Epetra_MultiVector object keep track of the number
of \e serial floating point operations performed using the specified object as the \e this argument
to the function.  The Flops() function returns this number as a double precision number.  Using this 
information, in conjunction with the Epetra_Time class, one can get accurate parallel performance
numbers.  The ResetFlops() function resets the floating point counter.

\warning A Epetra_Map, Epetra_LocalMap or Epetra_BlockMap object is required for all 
  Epetra_MultiVector constructors.

*/

//==========================================================================
class Epetra_MultiVector: public Epetra_DistObject, public Epetra_CompObject, public Epetra_BLAS {

  public:

  //! Basic Epetra_MultiVector constuctor.
  /*! Creates a Epetra_MultiVector object and fills with zero values.  

    \param In 
           Map - A Epetra_LocalMap, Epetra_Map or Epetra_BlockMap.

	   \warning Note that, because Epetra_LocalMap
	   derives from Epetra_Map and Epetra_Map derives from Epetra_BlockMap, this constructor works
	   for all three types of Epetra map classes.
    \param In 
           NumVectors - Number of vectors in multi-vector.

    \return Pointer to a Epetra_MultiVector.

  */
  Epetra_MultiVector(const Epetra_BlockMap& Map, int NumVectors);

  //! Epetra_MultiVector copy constructor.
  
  Epetra_MultiVector(const Epetra_MultiVector& Source);
  
  //! Set multi-vector values from two-dimensional array.
  /*!
    \param In 
           Epetra_DataAccess - Enumerated type set to Copy or View.
    \param In 
           Map - A Epetra_LocalMap, Epetra_Map or Epetra_BlockMap.
    \param In
           A - Pointer to an array of double precision numbers.  The first vector starts at A.
	   The second vector starts at A+MyLDA, the third at A+2*MyLDA, and so on.
    \param In
           MyLDA - The "Leading Dimension", or stride between vectors in memory.
	   \warning This value refers to the stride on the calling processor.  Thus it is a
	   local quantity, not a global quantity.
    \param In 
           NumVectors - Number of vectors in multi-vector.

    \return Integer error code, set to 0 if successful.

	   See Detailed Description section for further discussion.
  */
  Epetra_MultiVector(Epetra_DataAccess CV, const Epetra_BlockMap& Map, 
			double *A, int MyLDA, int NumVectors);

  //! Set multi-vector values from array of pointers.
  /*!
    \param In 
           Epetra_DataAccess - Enumerated type set to Copy or View.
    \param In 
           Map - A Epetra_LocalMap, Epetra_Map or Epetra_BlockMap.
    \param In
           ArrayOfPointers - An array of pointers such that ArrayOfPointers[i] points to the memory
	   location containing ith vector to be copied.
    \param In 
           NumVectors - Number of vectors in multi-vector.

    \return Integer error code, set to 0 if successful.

	   See Detailed Description section for further discussion.
  */
  Epetra_MultiVector(Epetra_DataAccess CV, const Epetra_BlockMap& Map, 
			 double **ArrayOfPointers, int NumVectors);

  //! Set multi-vector values from list of vectors in an existing Epetra_MultiVector.
  /*!
    \param In 
           Epetra_DataAccess - Enumerated type set to Copy or View.
    \param In
           Source - An existing fully constructed Epetra_MultiVector.
    \param In
           Indices - Integer list of the vectors to copy.  
    \param In 
           NumVectors - Number of vectors in multi-vector.

    \return Integer error code, set to 0 if successful.

	   See Detailed Description section for further discussion.
  */
  Epetra_MultiVector(Epetra_DataAccess CV,  
			const Epetra_MultiVector& Source, int *Indices, int NumVectors);

  //! Set multi-vector values from range of vectors in an existing Epetra_MultiVector.
  /*!
    \param In 
           Epetra_DataAccess - Enumerated type set to Copy or View.
    \param In
           Source - An existing fully constructed Epetra_MultiVector.
    \param In
           StartIndex - First of the vectors to copy.  
    \param In 
           NumVectors - Number of vectors in multi-vector.

    \return Integer error code, set to 0 if successful.

	   See Detailed Description section for further discussion.
  */
  Epetra_MultiVector(Epetra_DataAccess CV, 
			const Epetra_MultiVector& Source, int StartIndex, 
			int NumVectors);
  
  //! Epetra_MultiVector destructor.  
  virtual ~Epetra_MultiVector();
  
  //! Set multi-vector values to random numbers.
  /*!
    \return Integer error code, set to 0 if successful.

  */
   int Random();

  //! Put multi-vector values into user-provided two-dimensional array.
  /*!
    \param Out
           A - Pointer to memory space that will contain the multi-vector values.  
	   The first vector will be copied to the memory pointed to by A.
	   The second vector starts at A+MyLDA, the third at A+2*MyLDA, and so on.
    \param In
           MyLDA - The "Leading Dimension", or stride between vectors in memory.
	   \warning This value refers to the stride on the calling processor.  Thus it is a
	   local quantity, not a global quantity.

    \return Integer error code, set to 0 if successful.

	   See Detailed Description section for further discussion.
  */
  int ExtractCopy(double *A, int MyLDA) const;

  //! Put multi-vector values into user-provided array of pointers.
  /*!
    \param Out
           ArrayOfPointers - An array of pointers to memory space that will contain the 
	   multi-vector values, such that ArrayOfPointers[i] points to the memory
	   location where the ith vector to be copied.

    \return Integer error code, set to 0 if successful.

	   See Detailed Description section for further discussion.
  */
  int ExtractCopy(double **ArrayOfPointers) const;

  // ExtractView functions

  
  //! Set user-provided addresses of A and MyLDA.
  /*!
    \param Out
           A - Address of a pointer to that will be set to point to the values of the multi-vector.  
	   The first vector will be at the memory pointed to by A.
	   The second vector starts at A+MyLDA, the third at A+2*MyLDA, and so on.
    \param In
           MyLDA - Address of the "Leading Dimension", or stride between vectors in memory.
	   \warning This value refers to the stride on the calling processor.  Thus it is a
	   local quantity, not a global quantity.

    \return Integer error code, set to 0 if successful.

	   See Detailed Description section for further discussion.
  */
  int ExtractView(double **A, int *MyLDA) const;

  //! Set user-provided addresses of ArrayOfPointers.
  /*!
    \param Out
           ArrayOfPointers - Address of array of pointers to memory space that will set to the
	   multi-vector array of pointers, such that ArrayOfPointers[i] points to the memory
	   location where the ith vector is located.

    \return Integer error code, set to 0 if successful.

	   See Detailed Description section for further discussion.
  */
  int ExtractView(double ***ArrayOfPointers) const;



  //! Initialize all values in a multi-vector with constant value.
  /*!
    \param In
           Scalar - Value to use.

    \return Integer error code, set to 0 if successful.
  */
  int PutScalar (double Scalar);
  // Mathematical functions.

  //! Computes dot product of each corresponding pair of vectors.
  /*!
    \param In
           A - Multi-vector to be used with the "\e this" multivector.
    \param Out
           Result - Result[i] will contain the ith dot product result.

    \return Integer error code, set to 0 if successful.
  */
  int Dot(const Epetra_MultiVector& A, double *Result) const;

  //! Puts element-wise absolute values of input Multi-vector in target.
  /*!
    \param In
           A - Input Multi-vector.
    \param Out
           \e this will contain the absolute values of the entries of A.

    \return Integer error code, set to 0 if successful.
    
    Note:  It is possible to use the same argument for A and \e this.
  */
  int Abs(const Epetra_MultiVector& A);

  //! Puts element-wise reciprocal values of input Multi-vector in target.
  /*!
    \param In
           A - Input Multi-vector.
    \param Out
           \e this will contain the element-wise reciprocal values of the entries of A.

    \return Integer error code, set to 0 if successful.  Returns 2 if some entry
            is too small, but not zero.  Returns 1 if some entry is zero.
    
    Note:  It is possible to use the same argument for A and \e this.  Also, 
    if a given value of A is smaller than Epetra_DoubleMin (defined in Epetra_Epetra.h),
    but nonzero, then the return code is 2.  If an entry is zero, the return code
    is 1.  However, in all cases the reciprocal value is still used, even
    if a NaN is the result.
  */
  int Reciprocal(const Epetra_MultiVector& A);

  //! Scale the current values of a multi-vector, \e this = Scalar*\e this.
  /*!
    \param In
           Scalar - Scale value.
    \param Out
           \e This - Multi-vector with scaled values.

    \return Integer error code, set to 0 if successful.
  */
  int Scale(double Scalar);

  //! Replace multi-vector values with scaled values of A, \e this = Scalar*A.
  /*!
    \param In
           ScalarA - Scale value.
    \param In
           A - Multi-vector to copy.
    \param Out
           \e This - Multi-vector with values overwritten by scaled values of A.

    \return Integer error code, set to 0 if successful.
  */
  int Scale(double ScalarA, const Epetra_MultiVector& A);

  //! Update multi-vector values with scaled values of A, \e this = Scalar*\e this + ScalarA*A.
  /*!
    \param In
           ScalarA - Scale value for A.
    \param In
           A - Multi-vector to add.
    \param In
           Scalar - Scale value for \e this.
    \param Out
           \e This - Multi-vector with updatede values.

    \return Integer error code, set to 0 if successful.
  */
  int Update(double ScalarA, const Epetra_MultiVector& A, double Scalar);

  //! Update multi-vector with scaled values of A and B, \e this = Scalar*\e this + ScalarA*A + ScalarB*B.
  /*!
    \param In
           ScalarA - Scale value for A.
    \param In
           A - Multi-vector to add.
    \param In
           ScalarB - Scale value for B.
    \param In
           B - Multi-vector to add.
    \param In
           Scalar - Scale value for \e this.
    \param Out
           \e This - Multi-vector with updatede values.

    \return Integer error code, set to 0 if successful.
  */
  int Update(double ScalarA, const Epetra_MultiVector& A, 
		     double ScalarB, const Epetra_MultiVector& B, double Scalar);

  //! Compute 1-norm of each vector in multi-vector.
  /*!
    \param Out
           Result - Result[i] contains 1-norm of ith vector.

    \return Integer error code, set to 0 if successful.
  */
  int Norm1   (double * Result) const;

  //! Compute 2-norm of each vector in multi-vector.
  /*!
    \param Out
           Result - Result[i] contains 2-norm of ith vector.

    \return Integer error code, set to 0 if successful.
  */
  int Norm2   (double * Result) const;

  //! Compute Inf-norm of each vector in multi-vector.
  /*!
    \param Out
           Result - Result[i] contains Inf-norm of ith vector.

    \return Integer error code, set to 0 if successful.
  */
  int NormInf (double * Result) const;

  //! Compute Weighted 2-norm (RMS Norm) of each vector in multi-vector.
  /*!
    \param In
           Weights - Multi-vector of weights.  If Weights contains a single vector,
           that vector will be used as the weights for all vectors of \e this.  Otherwise,
           Weights should have the same number of vectors as \e this.
    \param Out
           Result - Result[i] contains the weighted 2-norm of ith vector.  Specifically
           if we denote the ith vector in the multivector by \f$x\f$, and the ith weight
           vector by \f$w\f$ and let j represent the jth entry of each vector, on return
           Result[i] will contain the following result:
           \f[\sqrt{(1/n)\sum_{j=1}^n(w_jx_j)^2}\f],
           where \f$n\f$ is the global length of the vectors.

    \return Integer error code, set to 0 if successful.
  */
  int NormWeighted   (const Epetra_MultiVector& Weights, double * Result) const;

  //! Compute minimum value of each vector in multi-vector.
  /*!
    \param Out
           Result - Result[i] contains minimum value of ith vector.

    \return Integer error code, set to 0 if successful.
  */
  int MinValue  (double * Result) const;

  //! Compute maximum value of each vector in multi-vector.
  /*!
    \param Out
           Result - Result[i] contains maximum value of ith vector.

    \return Integer error code, set to 0 if successful.
  */
  int MaxValue  (double * Result) const;

  //! Compute mean (average) value of each vector in multi-vector.
  /*!
    \param Out
           Result - Result[i] contains mean value of ith vector.

    \return Integer error code, set to 0 if successful.
  */
  int MeanValue (double * Result) const;
  
  
  //! Matrix-Matrix multiplication, \e this = Scalar*\e this + ScalarAB*A*B.
  /*! This function performs a variety of matrix-matrix multiply operations, interpreting
      the Epetra_MultiVectors (\e this-aka C , A and B) as 2D matrices.  Variations are due to
      the fact that A, B and C can be local replicated or global distributed
      Epetra_MultiVectors and that we may or may not operate with the transpose of 
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
         Scalar - Scalar to multiply with \e this.

    \return Integer error code, set to 0 if successful.

\warning {Each multi-vector A, B and \e this is checked if it has constant stride using the
         ConstantStride() query function.  If it does not have constant stride, a temporary
	 copy is made and used for the computation.  This activity is transparent to the user,
	 except that there is memory and computation overhead.  All temporary space is deleted
	 prior to exit.}
	 
  */
  int Multiply(char TransA, char TransB, double ScalarAB, 
		       const Epetra_MultiVector& A, const Epetra_MultiVector& B,
		       double Scalar );
  


  //! Multiply a Epetra_MultiVector with another, element-by-element.
  /*! This function supports diagonal matrix multiply.  A is usually a single vector
      while B and \e this may have one or more columns.  Note that B and \e this must
      have the same shape.  A can be one vector or have the same shape as B.  The actual
      computation is \e this = Scalar * \e this + ScalarAB * B @ A where @ denotes element-wise
      multiplication.
  */
  int Multiply(double ScalarAB, const Epetra_MultiVector& A, const Epetra_MultiVector& B,
		       double Scalar );


  //! Multiply a Epetra_MultiVector by the reciprocal of another, element-by-element.
  /*! This function supports diagonal matrix scaling.  A is usually a single vector
      while B and \e this may have one or more columns.  Note that B and \e this must
      have the same shape.  A can be one vector or have the same shape as B. The actual
      computation is \e this = Scalar * \e this + ScalarAB * B @ A where @ denotes element-wise
      division.
  */
  int ReciprocalMultiply(double ScalarAB, const Epetra_MultiVector& A, const Epetra_MultiVector& B,
		       double Scalar );


  // Random number utilities


  //! Set seed for Random function.
  /*!
    \param In
           Seed - Should be an odd positive integer value (stored as double).

    \return Integer error code, set to 0 if successful.
  */
  int SetSeed(double Seed){Seed_=Seed; return(0);};

  //! Get seed from Random function.
  /*!
    \return Current random number seed.
  */
  double Seed(){return(Seed_);};

  //! = Operator.
  /*!
    \param In
           A - Epetra_MultiVector to copy.

    \return Epetra_MultiVector.
  */
  Epetra_MultiVector& operator = (const Epetra_MultiVector& Source);
  
  // Local element access functions

  // 

  //! Vector access function.
  /*!
    \return Pointer to vector in multi-vector specified by Index.
  */
  double*& operator [] (int Index);
  //! Vector access function.
  /*!
    \return Pointer to vector in multi-vector specified by Index.
  */
  const double*& operator [] (int Index) const;

  // Attribute access functions
  
  //! Returns the number of vectors in the multi-vector.
  int NumVectors() const {return(NumVectors_);};

  //! Returns the local vector length on the calling processor of vectors in the multi-vector.
  int MyLength() const {return(MyLength_);};

  //! Returns the global vector length of vectors in the multi-vector.
  int GlobalLength() const {return(GlobalLength_);};

  //! Returns the stride between  vectors in the multi-vector (only meaningful if ConstantStride() is true).
  int Stride() const {return(Stride_);};
  
  //! Returns true if this multi-vector has constant stride between vectors.
  bool ConstantStride() const {return(ConstantStride_);};

  //! Print method
  virtual void Print(ostream & os) const;
  
 protected:

  //! Get pointer to MultiVector values
  double* Values() const {return Values_;};

  //! Get pointer to individual vector pointers
  double** Pointers() const {return Pointers_;};

  // Internal utilities
  void Assign(const Epetra_MultiVector& rhs);
  int CheckInput();

  int IndexBase_;
  double *Values_;    // local MultiVector coefficients
  double **Pointers_;        // Pointers to each vector;
  
  friend class Epetra_CrsMatrix;
  friend class Epetra_VbrMatrix;
 private:


  // Internal utilities

  int Reduce(); 
  int AllocateForCopy(void);
  int DoCopy(void);
  int AllocateForView(void);
  int DoView(void);
  
  int CheckSizes(const Epetra_DistObject& A);
  int CopyAndPermute(const Epetra_DistObject & Source, int NumSameIDs, 
			 int NumPermuteIDs, int * PermuteToLIDs, int * PermuteFromLIDs);

  int PackAndPrepare(const Epetra_DistObject & Source, int NumExportIDs, int * ExportLIDs,
				      int Nsend, int Nrecv,
				      int & LenExports, char * & Exports, int & LenImports, 
				      char * & Imports, 
				      int & SizeOfPacket, Epetra_Distributor & Distor);
  
  int UnpackAndCombine(const Epetra_DistObject & Source,
		       int NumImportIDs, int * ImportLIDs, 
		       char * Imports, int & SizeOfPacket, 
		       Epetra_Distributor & Distor, Epetra_CombineMode CombineMode );


  int MyLength_;
  int GlobalLength_;
  int NumVectors_;
  bool UserAllocated_;
  bool ConstantStride_;
  int Stride_;
  bool Allocated_;
  double Seed_;
  double * DoubleTemp_;
  int * IntTemp_;

};

#endif /* _EPETRA_MULTIVECTOR_H_ */
