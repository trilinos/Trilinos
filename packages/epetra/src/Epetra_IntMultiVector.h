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

#ifndef EPETRA_INTMULTIVECTOR_H
#define EPETRA_INTMULTIVECTOR_H

#if defined(Epetra_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Epetra package is deprecated"
#endif
#endif



class Epetra_Comm;
class Epetra_BlockMap;
class Epetra_Map;
class Epetra_Import;
class Epetra_Export;
class Epetra_Distributor;
class Epetra_IntVector;

#include "Epetra_ConfigDefs.h"
#include "Epetra_DistObject.h"
#include "Epetra_CompObject.h"
#include "Epetra_BLAS.h"
#include "Epetra_Util.h"

//! Epetra_IntMultiVector: A class for constructing and using dense multi-vectors, vectors and matrices in parallel.

/*! The Epetra_IntMultiVector class enables the construction and use of ordinal-valued,
  dense vectors, multi-vectors,
  and matrices in a distributed memory environment.  The dimensions and distribution of the dense
  multi-vectors is determined in part by a Epetra_Comm object, a Epetra_Map (or Epetra_LocalMap
  or Epetra_BlockMap) and the number of vectors passed to the constructors described below.

  There are several concepts that important for understanding the Epetra_IntMultiVector class:

  <ul>
  <li>  IntMulti-vectors, IntVectors and IntMatrices.
  <ul>
  <li> IntVector - A list of ordinal-valued numbers.  Also an IntMultiVector with one vector.
  <li> IntMulti-Vector - A collection of one or more vectors, all having the same length and distribution.
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

  <b>Constructing Epetra_IntMultiVectors</b>

  Except for the basic constructor and copy constructor, Epetra_IntMultiVector constructors
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

  All Epetra_IntMultiVector constructors require a map argument that describes the layout of elements
  on the parallel machine.  Specifically,
  \c map is a Epetra_Map, Epetra_LocalMap or Epetra_BlockMap object describing the desired
  memory layout for the multi-vector.

  There are six different Epetra_IntMultiVector constructors:
  <ul>
  <li> Basic - All values are zero.
  <li> Copy - Copy an existing multi-vector.
  <li> Copy from or make view of two-dimensional Fortran style array.
  <li> Copy from or make view of an array of pointers.
  <li> Copy or make view of a list of vectors from another Epetra_IntMultiVector object.
  <li> Copy or make view of a range of vectors from another Epetra_IntMultiVector object.
  </ul>

  <b>Extracting Data from Epetra_IntMultiVectors</b>

  Once a Epetra_IntMultiVector is constructed, it is possible to extract a copy of the values
  or create a view of them.

  \warning ExtractView functions are \e extremely dangerous from a data hiding perspective.
  For both ExtractView fuctions, there is a corresponding ExtractCopy function.  We
  strongly encourage users to develop code using ExtractCopy functions first and
  only use the ExtractView functions in a secondary optimization phase.

  There are four Extract functions:
  <ul>
  <li> ExtractCopy - Copy values into a user-provided two-dimensional array.
  <li> ExtractCopy - Copy values into a user-provided array of pointers.
  <li> ExtractView - Set user-provided two-dimensional array parameters
  to point to Epetra_IntMultiVector data.
  <li> ExtractView - Set user-provided array of pointer parameters
  to point to Epetra_IntMultiVector data.
  </ul>

  <b>Vector, Matrix and Utility Functions</b>

  Once a Epetra_IntMultiVector is constructed, a variety of mathematical functions can be applied to
  the individual vectors.  Specifically:
  <ul>
  <li> Dot Products.
  <li> Vector Updates.
  <li> \e p Norms.
  <li> Weighted Norms.
  <li> Minimum, Maximum and Average Values.
  </ul>

  <b> Counting Floating Point Operations </b>

  Each Epetra_IntMultiVector object keep track of the number
  of \e serial floating point operations performed using the specified object as the \e this argument
  to the function.  The Flops() function returns this number as a double precision number.  Using this
  information, in conjunction with the Epetra_Time class, one can get accurate parallel performance
  numbers.  The ResetFlops() function resets the floating point counter.

  \warning A Epetra_Map, Epetra_LocalMap or Epetra_BlockMap object is required for all
  Epetra_IntMultiVector constructors.

*/

//==========================================================================
class EPETRA_LIB_DLL_EXPORT Epetra_IntMultiVector: public Epetra_DistObject, public Epetra_CompObject, public Epetra_BLAS {

 public:

   //! @name Constructors/destructors
  //@{
  //! Basic Epetra_IntMultiVector constuctor.
  /*! Creates a Epetra_IntMultiVector object and, by default, fills with zero values.

  \param In
  Map - A Epetra_LocalMap, Epetra_Map or Epetra_BlockMap.

  \warning Note that, because Epetra_LocalMap
  derives from Epetra_Map and Epetra_Map derives from Epetra_BlockMap, this constructor works
  for all three types of Epetra map classes.
  \param In
  NumVectors - Number of vectors in multi-vector.
  \param In
  zeroOut - If <tt>true</tt> then the allocated memory will be zeroed
            out initialy.  If <tt>false</tt> then this memory will not
            be touched which can be significantly faster.
  \return Pointer to a Epetra_IntMultiVector.

  */
  Epetra_IntMultiVector(const Epetra_BlockMap& Map, int NumVectors, bool zeroOut = true);

  //! Epetra_MultiVector copy constructor.

  Epetra_IntMultiVector(const Epetra_IntMultiVector& Source);

  //! Set multi-vector values from two-dimensional array.
  /*!
    \param In
    Epetra_DataAccess - Enumerated type set to Copy or View.
    \param In
    Map - A Epetra_LocalMap, Epetra_Map or Epetra_BlockMap.
    \param In
    A - Pointer to an array of ordinal numbers.  The first vector starts at A.
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
  Epetra_IntMultiVector(Epetra_DataAccess CV, const Epetra_BlockMap& Map,
         int *A, int MyLDA, int NumVectors);

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
  Epetra_IntMultiVector(Epetra_DataAccess CV, const Epetra_BlockMap& Map,
         int **ArrayOfPointers, int NumVectors);

  //! Set multi-vector values from list of vectors in an existing Epetra_IntMultiVector.
  /*!
    \param In
    Epetra_DataAccess - Enumerated type set to Copy or View.
    \param In
    Source - An existing fully constructed Epetra_IntMultiVector.
    \param In
    Indices - Integer list of the vectors to copy.
    \param In
    NumVectors - Number of vectors in multi-vector.

    \return Integer error code, set to 0 if successful.

    See Detailed Description section for further discussion.
  */
  Epetra_IntMultiVector(Epetra_DataAccess CV,
         const Epetra_IntMultiVector& Source, int *Indices, int NumVectors);

  //! Set multi-vector values from range of vectors in an existing Epetra_IntMultiVector.
  /*!
    \param In
    Epetra_DataAccess - Enumerated type set to Copy or View.
    \param In
    Source - An existing fully constructed Epetra_IntMultiVector.
    \param In
    StartIndex - First of the vectors to copy.
    \param In
    NumVectors - Number of vectors in multi-vector.

    \return Integer error code, set to 0 if successful.

    See Detailed Description section for further discussion.
  */
  Epetra_IntMultiVector(Epetra_DataAccess CV,
         const Epetra_IntMultiVector& Source, int StartIndex,
         int NumVectors);

  //! Epetra_MultiVector destructor.
  virtual ~Epetra_IntMultiVector();
  //@}

  //! @name Post-construction modification routines
  //@{

  //! Replace current value  at the specified (GlobalRow, VectorIndex) location with OrdinalValue.
  /*!
    Replaces the  existing value for a single entry in the multivector.  The
    specified global row must correspond to a GID owned by the map of the multivector on the
    calling processor.  In other words, this method does not perform cross-processor communication.

    If the map associated with this multivector is an Epetra_BlockMap, only the first point entry associated
    with the global row will be modified.  To modify a different point entry, use the other version of
    this method

    \param In
    GlobalRow - Row of Multivector to modify in global index space.
    \param In
    VectorIndex - Vector within MultiVector that should to modify.
    \param In
    OrdinalValue - Value to add to existing value.

    \return Integer error code, set to 0 if successful, set to 1 if GlobalRow not associated with calling processor
    set to -1 if VectorIndex >= NumVectors().
  */
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
  int ReplaceGlobalValue(int GlobalRow, int VectorIndex, int OrdinalValue);
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
  int ReplaceGlobalValue(long long GlobalRow, int VectorIndex, int OrdinalValue);
#endif


  //! Replace current value at the specified (GlobalBlockRow, BlockRowOffset, VectorIndex) location with OrdinalValue.
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
    VectorIndex - Vector within MultiVector that should to modify.
    \param In
    OrdinalValue - Value to add to existing value.

    \return Integer error code, set to 0 if successful, set to 1 if GlobalRow not associated with calling processor
    set to -1 if VectorIndex >= NumVectors(), set to -2 if BlockRowOffset is out-of-range.
  */
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
  int ReplaceGlobalValue(int GlobalBlockRow, int BlockRowOffset, int VectorIndex, int OrdinalValue);
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
  int ReplaceGlobalValue(long long GlobalBlockRow, int BlockRowOffset, int VectorIndex, int OrdinalValue);
#endif


  //! Adds OrdinalValue to existing value at the specified (GlobalRow, VectorIndex) location.
  /*!
    Sums the given value into the existing value for a single entry in the multivector.  The
    specified global row must correspond to a GID owned by the map of the multivector on the
    calling processor.  In other words, this method does not perform cross-processor communication.

    If the map associated with this multivector is an Epetra_BlockMap, only the first point entry associated
    with the global row will be modified.  To modify a different point entry, use the other version of
    this method

    \param In
    GlobalRow - Row of Multivector to modify in global index space.
    \param In
    VectorIndex - Vector within MultiVector that should to modify.
    \param In
    OrdinalValue - Value to add to existing value.

    \return Integer error code, set to 0 if successful, set to 1 if GlobalRow not associated with calling processor
    set to -1 if VectorIndex >= NumVectors().
  */
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
  int SumIntoGlobalValue(int GlobalRow, int VectorIndex, int OrdinalValue);
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
  int SumIntoGlobalValue(long long GlobalRow, int VectorIndex, int OrdinalValue);
#endif


  //! Adds OrdinalValue to existing value at the specified (GlobalBlockRow, BlockRowOffset, VectorIndex) location.
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
    VectorIndex - Vector within MultiVector that should to modify.
    \param In
    ScalarValue - Value to add to existing value.

    \return Integer error code, set to 0 if successful, set to 1 if GlobalRow not associated with calling processor
    set to -1 if VectorIndex >= NumVectors(), set to -2 if BlockRowOffset is out-of-range.
  */
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
  int SumIntoGlobalValue(int GlobalBlockRow, int BlockRowOffset, int VectorIndex, int OrdinalValue);
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
  int SumIntoGlobalValue(long long GlobalBlockRow, int BlockRowOffset, int VectorIndex, int OrdinalValue);
#endif

  //! Replace current value  at the specified (MyRow, VectorIndex) location with OrdinalValue.
  /*!
    Replaces the existing value for a single entry in the multivector.  The
    specified local row must correspond to a GID owned by the map of the multivector on the
    calling processor.  In other words, this method does not perform cross-processor communication.

    This method is intended for use with vectors based on an Epetra_Map.  If used
    on a vector based on a non-trivial Epetra_BlockMap, this will update only block
    row 0, i.e.

    Epetra_IntMultiVector::ReplaceMyValue  (  MyRow,  VectorIndex,  OrdinalValue )  is
    equivalent to:
    Epetra_IntMultiVector::ReplaceMyValue  (  0, MyRow,  VectorIndex,  OrdinalValue )


    \param In
    MyRow - Row of IntMultivector to modify in local index space.
    \param In
    VectorIndex - Vector within IntMultiVector that should to modify.
    \param In
    OrdinalValue - Value to add to existing value.

    \return Integer error code, set to 0 if successful, set to 1 if MyRow not associated with calling processor
    set to -1 if VectorIndex >= NumVectors().
  */
  int ReplaceMyValue(int MyRow, int VectorIndex, int OrdinalValue);


  //! Replace current value at the specified (MyBlockRow, BlockRowOffset, VectorIndex) location with OrdinalValue.
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
    VectorIndex - Vector within MultiVector that should to modify.
    \param In
    OrdinalValue - Value to add to existing value.

    \return Integer error code, set to 0 if successful, set to 1 if MyRow not associated with calling processor
    set to -1 if VectorIndex >= NumVectors(), set to -2 if BlockRowOffset is out-of-range.
  */
  int ReplaceMyValue(int MyBlockRow, int BlockRowOffset, int VectorIndex, int OrdinalValue);


  //! Adds OrdinalValue to existing value at the specified (MyRow, VectorIndex) location.
  /*!
    Sums the given value into the existing value for a single entry in the multivector.  The
    specified local row must correspond to a GID owned by the map of the multivector on the
    calling processor.  In other words, this method does not perform cross-processor communication.

    If the map associated with this multivector is an Epetra_BlockMap, only the first point entry associated
    with the local row will be modified.  To modify a different point entry, use the other version of
    this method

    \param In
    MyRow - Row of Multivector to modify in local index space.
    \param In
    VectorIndex - Vector within MultiVector that should to modify.
    \param In
    OrdinalValue - Value to add to existing value.

    \return Integer error code, set to 0 if successful, set to 1 if MyRow not associated with calling processor
    set to -1 if VectorIndex >= NumVectors().
  */
  int SumIntoMyValue(int MyRow, int VectorIndex, int OrdinalValue);


  //! Adds OrdinalValue to existing value at the specified (MyBlockRow, BlockRowOffset, VectorIndex) location.
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
    VectorIndex - Vector within MultiVector that should to modify.
    \param In
    OrdinalValue - Value to add to existing value.

    \return Integer error code, set to 0 if successful, set to 1 if MyRow not associated with calling processor
    set to -1 if VectorIndex >= NumVectors(), set to -2 if BlockRowOffset is out-of-range.
  */
  int SumIntoMyValue(int MyBlockRow, int BlockRowOffset, int VectorIndex, int OrdinalValue);

  //! Initialize all values in a multi-vector with constant value.
  /*!
    \param In
    OrdinalConstant - Value to use.

    \return Integer error code, set to 0 if successful.
  */
  int PutScalar (int OrdinalConstant);

  //@}

  //! @name Extraction methods
  //@{

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
  int ExtractCopy(int *A, int MyLDA) const;

  //! Put multi-vector values into user-provided array of pointers.
  /*!
    \param Out
    ArrayOfPointers - An array of pointers to memory space that will contain the
    multi-vector values, such that ArrayOfPointers[i] points to the memory
    location where the ith vector to be copied.

    \return Integer error code, set to 0 if successful.

    See Detailed Description section for further discussion.
  */
  int ExtractCopy(int **ArrayOfPointers) const;

  // ExtractView functions


  //! Set user-provided addresses of A and MyLDA.
  /*!
    \param
    A (Out) - Address of a pointer to that will be set to point to the values of the multi-vector.
    The first vector will be at the memory pointed to by A.
    The second vector starts at A+MyLDA, the third at A+2*MyLDA, and so on.
    \param
    MyLDA (Out) - Address of the "Leading Dimension", or stride between vectors in memory.
    \warning This value refers to the stride on the calling processor.  Thus it is a
    local quantity, not a global quantity.

    \return Integer error code, set to 0 if successful.

    See Detailed Description section for further discussion.
  */
  int ExtractView(int **A, int *MyLDA) const;

  //! Set user-provided addresses of ArrayOfPointers.
  /*!
    \param
    ArrayOfPointers (Out) - Address of array of pointers to memory space that will set to the
    multi-vector array of pointers, such that ArrayOfPointers[i] points to the memory
    location where the ith vector is located.

    \return Integer error code, set to 0 if successful.

    See Detailed Description section for further discussion.
  */
  int ExtractView(int ***ArrayOfPointers) const;

  //@}

  //! @name Mathematical methods
  //@{

  //! Compute minimum value of each vector in multi-vector.
  /*! Note that the vector contents must be already initialized for this
      function to compute a well-defined result. The length of the
      vector need not be greater than zero on all processors. If length is
      greater than zero on any processor then a valid result will be computed.
    \param Out
    Result - Result[i] contains minimum value of ith vector.

    \return Integer error code, set to 0 if successful.
  */
  int MinValue  (int * Result) const;

  //! Compute maximum value of each vector in multi-vector.
  /*! Note that the vector contents must be already initialized for this
      function to compute a well-defined result. The length of the
      vector need not be greater than zero on all processors. If length is
      greater than zero on any processor then a valid result will be computed.
    \param Out
    Result - Result[i] contains maximum value of ith vector.

    \return Integer error code, set to 0 if successful.
  */
  int MaxValue  (int * Result) const;

  //@}

  //! @name Overloaded operators
  //@{

  //! = Operator.
  /*!
    \param In
    A - Epetra_IntMultiVector to copy.

    \return Epetra_IntMultiVector.
  */
  Epetra_IntMultiVector& operator = (const Epetra_IntMultiVector& Source);

  // Local element access functions

  //

  //! Vector access function.
  /*!
    \return Pointer to the array of ordinals containing the local values of the ith vector in the multi-vector.
  */
  int*& operator [] (int i) { return Pointers_[i]; }
  //! Vector access function.
  /*!
    \return Pointer to the array of ordinal containing the local values of the ith vector in the multi-vector.
  */
  //  const int*& operator [] (int i) const;
  int * const & operator [] (int i) const { return Pointers_[i]; }

  //! Vector access function.
  /*!
    \return An Epetra_IntVector pointer to the ith vector in the multi-vector.
  */
  Epetra_IntVector * & operator () (int i);
  //! Vector access function.
  /*!
    \return An Epetra_IntVector pointer to the ith vector in the multi-vector.
  */
  const Epetra_IntVector * & operator () (int i) const;

  //@}

  //! @name Attribute access functions
  //@{

  //! Returns the number of vectors in the multi-vector.
  int NumVectors() const {return(NumVectors_);}

  //! Returns the local vector length on the calling processor of vectors in the multi-vector.
  int MyLength() const {return(MyLength_);}

  //! Returns the global vector length of vectors in the multi-vector.
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
  int GlobalLength() const {
    if(Map().GlobalIndicesInt())
      return (int) GlobalLength_;
    throw "Epetra_MultiVector::GlobalLength: GlobalIndices not int.";
  }
#endif
  long long GlobalLength64() const {return(GlobalLength_);}

  //! Returns the stride between  vectors in the multi-vector (only meaningful if ConstantStride() is true).
  int Stride() const {return(Stride_);}

  //! Returns true if this multi-vector has constant stride between vectors.
  bool ConstantStride() const {return(ConstantStride_);}
  //@}

  /** Replace map, only if new map has same point-structure as current map.
      return 0 if map is replaced, -1 if not.
   */
  int ReplaceMap(const Epetra_BlockMap& map);

  //! @name I/O methods
  //@{

  //! Print method
  virtual void Print(std::ostream & os) const;
  //@}

  //! @name Expert-only unsupported methods
  //@{

  //! Reset the view of an existing multivector to point to new user data.
  /*! Allows the (very) light-weight replacement of multivector values for an
    existing multivector that was constructed using an Epetra_DataAccess mode of View.
    No checking is performed to see if the array of values passed in contains valid
    data.  It is assumed that the user has verified the integrity of data before calling
    this method. This method is useful for situations where a multivector is needed
    for use with an Epetra operator or matrix and the user is not passing in a multivector,
    or the multivector is being passed in with another map that is not exactly compatible
    with the operator, but has the correct number of entries.

    This method is used by AztecOO and Ifpack in the matvec, and solve methods to improve
    performance and reduce repeated calls to constructors and destructors.

    @param ArrayOfPointers Contains the array of pointers containing the multivector data.

    \return Integer error code, set to 0 if successful, -1 if the multivector was not created as a View.

    \warning This method is extremely dangerous and should only be used by experts.
  */

  int ResetView(int ** ArrayOfPointers);

  //! Get pointer to MultiVector values
  int* Values() const {return Values_;}

  //! Get pointer to individual vector pointers
  int** Pointers() const {return Pointers_;}
  //@}

  // Expert-only function
  //  SuperLU defines Reduce to be a macro in util.h
#ifdef Reduce
#undef Reduce
#endif
  int Reduce();

 protected:

  // Internal utilities
  void Assign(const Epetra_IntMultiVector& rhs);
  int CheckInput();

  int *Values_;    // local MultiVector coefficients

 private:


  // Internal utilities

  int AllocateForCopy(void);
  int DoCopy(void);

  inline void UpdateOrdinalTemp() const
  {if (OrdinalTemp_==0) OrdinalTemp_=new int[NumVectors_+1]; return;}

  inline void UpdateIntVectors()  const {if (IntVectors_==0) { IntVectors_ = new Epetra_IntVector *[NumVectors_];
    for (int i=0; i<NumVectors_; i++) IntVectors_[i] = 0;}
    return;
  }

  int AllocateForView(void);
  int DoView(void);
  template<typename int_type>
  int ChangeGlobalValue(int_type GlobalBlockRow,
                        int BlockRowOffset,
                        int VectorIndex,
                        int OrdinalValue,
                        bool SumInto);
  int ChangeMyValue(int MyBlockRow,
                    int BlockRowOffset,
                    int VectorIndex,
                    int OrdinalValue,
                    bool SumInto);

  int CheckSizes(const Epetra_SrcDistObject& A);

  int CopyAndPermute(const Epetra_SrcDistObject & Source,
                     int NumSameIDs,
                     int NumPermuteIDs,
                     int * PermuteToLIDs,
                     int * PermuteFromLIDs,
                     const Epetra_OffsetIndex * Indexor,
                     Epetra_CombineMode CombineMode = Zero);

  int PackAndPrepare(const Epetra_SrcDistObject & Source,
                     int NumExportIDs,
                     int * ExportLIDs,
                     int & LenExports,
                     char * & Exports,
                     int & SizeOfPacket,
                     int * Sizes,
                     bool & VarSizes,
                     Epetra_Distributor & Distor);

  int UnpackAndCombine(const Epetra_SrcDistObject & Source,
                       int NumImportIDs,
                       int * ImportLIDs,
                       int LenImports,
                       char * Imports,
                       int & SizeOfPacket,
                       Epetra_Distributor & Distor,
                       Epetra_CombineMode CombineMode,
                       const Epetra_OffsetIndex * Indexor );

  int **Pointers_;        // Pointers to each vector;

  int MyLength_;
  long long GlobalLength_;
  int NumVectors_;
  bool UserAllocated_;
  bool ConstantStride_;
  int Stride_;
  bool Allocated_;
  mutable int * OrdinalTemp_;
  mutable Epetra_IntVector ** IntVectors_;
  Epetra_Util Util_;

};

#endif /* EPETRA_INTMULTIVECTOR_H */
