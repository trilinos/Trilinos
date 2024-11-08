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

#ifndef EPETRA_LONGLONGVECTOR_H
#define EPETRA_LONGLONGVECTOR_H

#if defined(Epetra_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Epetra package is deprecated"
#endif
#endif



#include "Epetra_ConfigDefs.h"
#include "Epetra_DistObject.h"
#include "Epetra_BlockMap.h"
#include "Epetra_Distributor.h"
class Epetra_Map;

#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES

//! Epetra_LongLongVector: A class for constructing and using dense integer vectors on a parallel computer.

/*! The Epetra_LongLongVector class enables the construction and use of integer
     dense vectors in a distributed memory environment.  The distribution of the dense
    vector is determined in part by a Epetra_Comm object and a Epetra_Map (or Epetra_LocalMap
    or Epetra_BlockMap).


<b> Distributed Global vs. Replicated Local</b>
<ul>
  <li> Distributed Global Vectors - In most instances, a multi-vector will be partitioned
       across multiple memory images associated with multiple processors.  In this case, there is
       a unique copy of each element and elements are spread across all processors specified by
       the Epetra_Comm communicator.
  <li> Replicated Local Vectors - Some algorithms use vectors that are too small to
       be distributed across all processors.  Replicated local vectors handle
       these types of situation.
</ul>

<b>Constructing Epetra_LongLongVectors</b>

There are four Epetra_LongLongVector constructors.  The first is a basic constructor that allocates
space and sets all values to zero, the second is a
copy constructor. The third and fourth constructors work with user data.  These constructors have
two data access modes:
<ol>
  <li> Copy mode - Allocates memory and makes a copy of the user-provided data. In this case, the
       user data is not needed after construction.
  <li> View mode - Creates a "view" of the user data. In this case, the
       user data is required to remain intact for the life of the vector.
</ol>

\warning View mode is \e extremely dangerous from a data hiding perspective.
Therefore, we strongly encourage users to develop code using Copy mode first and
only use the View mode in a secondary optimization phase.

All Epetra_LongLongVector constructors require a map argument that describes the layout of elements
on the parallel machine.  Specifically,
\c map is a Epetra_Map, Epetra_LocalMap or Epetra_BlockMap object describing the desired
memory layout for the vector.

There are four different Epetra_LongLongVector constructors:
<ul>
  <li> Basic - All values are zero.
  <li> Copy - Copy an existing vector.
  <li> Copy from or make view of user int array.
</ul>

<b>Extracting Data from Epetra_LongLongVectors</b>

Once a Epetra_LongLongVector is constructed, it is possible to extract a copy of the values or create
a view of them.

\warning ExtractView functions are \e extremely dangerous from a data hiding perspective.
For both ExtractView fuctions, there is a corresponding ExtractCopy function.  We
strongly encourage users to develop code using ExtractCopy functions first and
only use the ExtractView functions in a secondary optimization phase.

There are two Extract functions:
<ul>
  <li> ExtractCopy - Copy values into a user-provided array.
  <li> ExtractView - Set user-provided array to point to Epetra_LongLongVector data.
</ul>


\warning A Epetra_Map, Epetra_LocalMap or Epetra_BlockMap object is required for all
  Epetra_LongLongVector constructors.

*/

//=========================================================================
class EPETRA_LIB_DLL_EXPORT Epetra_LongLongVector : public Epetra_DistObject {

  public:

    //! @name Constructors/destructors
  //@{
  //! Basic Epetra_LongLongVector constuctor.
  /*! Creates a Epetra_LongLongVector object and, by default, fills with zero values.

    \param In
           Map - A Epetra_LocalMap, Epetra_Map or Epetra_BlockMap.

     \warning Note that, because Epetra_LocalMap
     derives from Epetra_Map and Epetra_Map derives from Epetra_BlockMap, this constructor works
     for all three types of Epetra map classes.
  \param In
  zeroOut - If <tt>true</tt> then the allocated memory will be zeroed
            out initialy.  If <tt>false</tt> then this memory will not
            be touched which can be significantly faster.

    \return Pointer to a Epetra_LongLongVector.

  */
  Epetra_LongLongVector(const Epetra_BlockMap& Map, bool zeroOut = true);

  //! Epetra_LongLongVector copy constructor.

  Epetra_LongLongVector(const Epetra_LongLongVector& Source);

  //! Set vector values from user array.
  /*!
    \param In
           Epetra_DataAccess - Enumerated type set to Copy or View.
    \param In
           Map - A Epetra_LocalMap, Epetra_Map or Epetra_BlockMap.
    \param In
           V - Pointer to an array of long long numbers..

    \return Integer error code, set to 0 if successful.

     See Detailed Description section for further discussion.
  */
  Epetra_LongLongVector(Epetra_DataAccess CV, const Epetra_BlockMap& Map, long long *V);

  //! Epetra_LongLongVector destructor.
  virtual ~Epetra_LongLongVector ();
  //@}


  //! @name Post-construction modification methods
  //@{
  //! Set all elements of the vector to Value
  int PutValue(long long Value);
  //@}


  //! @name Extraction methods
  //@{


  //! Put vector values into user-provided array.
  /*!
    \param Out
           V - Pointer to memory space that will contain the vector values.

    \return Integer error code, set to 0 if successful.
  */
  int ExtractCopy(long long *V) const;

  //! Set user-provided address of V.
  /*!
    \param Out
           V - Address of a pointer to that will be set to point to the values of the vector.

    \return Integer error code, set to 0 if successful.
  */
  int ExtractView(long long **V) const;
  //@}

  //! @name Mathematical methods
  //@{
  //! Find maximum value
  /*!
    \return Maximum value across all processors.
  */
  long long MaxValue();

  //! Find minimum value
  /*!
    \return Minimum value across all processors.
  */
  long long MinValue();

  //@}

  //! @name Overloaded operators
  //@{

  //! = Operator.
  /*!
    \param In
           A - Epetra_LongLongVector to copy.

    \return Epetra_LongLongVector.
  */
  Epetra_LongLongVector& operator = (const Epetra_LongLongVector& Source);

  //! Element access function.
  /*!
    \return V[Index].
  */
    long long& operator [] (int index) { return Values_[index]; }
  //! Element access function.
  /*!
    \return V[Index].
  */
    const long long& operator [] (int index) const { return Values_[index]; }
    //@}

    //! @name Attribute access functions
  //@{

  //! Returns a pointer to an array containing the values of this vector.
  long long * Values() const {return(Values_);};

  //! Returns the local vector length on the calling processor of vectors in the multi-vector.
  int MyLength() const {return(Map().NumMyPoints());};

  //! Returns the global vector length of vectors in the multi-vector.
  long long GlobalLength64() const {return(Map().NumGlobalPoints64());};
  //@}

  //! @name I/O methods
  //@{

  //! Print method
  virtual void Print(std::ostream & os) const;
    //@}
 private:

  int AllocateForCopy();
  int DoCopy(long long * V);
  int AllocateForView();
  int DoView(long long * V);

   // Routines to implement Epetra_DistObject virtual methods
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
                     bool& VarSizes,
                     Epetra_Distributor & Distor);

  int UnpackAndCombine(const Epetra_SrcDistObject & Source,
                       int NumImportIDs,
                       int * ImportLIDs,
                       int LenImports,
                       char * Imports,
                       int & SizeOfPacket,
                       Epetra_Distributor & Distor,
                       Epetra_CombineMode CombineMode,
                       const Epetra_OffsetIndex * Indexor);

  long long * Values_;
  bool UserAllocated_;
  bool Allocated_;
};

#endif // EPETRA_NO_64BIT_GLOBAL_INDICES

#endif /* EPETRA_LONGLONGVECTOR_H */
