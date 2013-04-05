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

#ifndef EPETRA_MAP_H
#define EPETRA_MAP_H

//! Epetra_Map: A class for partitioning vectors and matrices.

/*! It is often the case that multiple matrix and vector objects have an identical distribution 
  of elements on a parallel machine. The Epetra_Map class keep information that describes 
  this distribution for matrices and vectors.  

  Epetra_Map allows the storage and retrieval of the following information.  Depending on the
  constructor that is used, some of the information is defined by the user and some is 
  determined by the constructor.  Once a Epetra_Map is constructed any of the following attributes can 
  be obtained
  by calling a query function that has the name as the attribute, e.g. to get the
  value of NumGlobalElements, you can call a function NumGlobalElements().  For attributes that
  are lists, the query functions return the list values in a user allocated array.

  <ul>
  <li> NumGlobalElements - The total number of elements across all processors. If this parameter and
       NumMyElements are both passed into the constructor, one of the three cases will apply: 
       <ol> 
       <li> If NumGlobalElements = NumMyElements (and not equal to zero)
            the map is defined to be a local replicated
            map.  In this case, objects constructed using this map will be identically replicated across
      all processors in the communicator.
       <li> If NumGlobalElements = -1 and NumMyElements is passed in then NumGlobalElements will
            be computed as the sum of NumMyElements across all processors.
       <li> If neither of the above is true, NumGlobalElements will be checked against the sum of 
            NumMyElements across all processors.  An error is issued if the comparison is not equal.
       </ol>
  <li> NumMyElements - The number of elements owned by the calling processor.
  <li> MyGlobalElements - A list of length NumMyElements that contains the global element IDs
       of the elements owned by the calling processor.
  <li> IndexBase - The base integer value for indexed array references.  Typically this is 0
       for C/C++ and 1 for Fortran, but it can be set to any integer value.
  <li> Comm - The Epetra_Comm communicator.  This communicator can in turn be queried for
       processor rank and size information.
  </ul>


  In addition to the information above that is passed in to or created by the Epetra_Map constructor,
  the following attributes are computed and available via query to the user using the same scheme
  as above, e.g., use NumGlobalPoints() to get the value of NumGlobalPoints.

  <ul>
  <li> NumGlobalPoints - The total number of points across all processors.
  <li> NumMyPoints - The number of points on the calling processor.
  <li> MinAllGID - The minimum global index value across all processors.
  <li> MaxAllGID - The maximum global index value across all processors.
  <li> MinMyGID - The minimum global index value on the calling processor.
  <li> MaxMyGID - The maximum global index value on the calling processor.
  <li> MinLID - The minimum local index value on the calling processor.
  <li> MaxLID - The maximum local index value on the calling processor.
  </ul>

  The following functions allow boolean tests for certain properties.    

  <ul>
  <li> LinearMap() - Returns true if the elements are distributed linear across processors, i.e.,
       processor 0 gets the first n/p elements, processor 1 gets the next n/p elements, etc. where
       n is the number of elements and p is the number of processors.
  <li> DistributedGlobal() - Returns true if the element space of the map spans more than one processor.
       This will be true in most cases, but will be false in serial cases and for objects
       that are created via the derived Epetra_LocalMap class.
  </ul>

  \warning An Epetra_Comm object is required for all Epetra_Map constructors.

  \note In the current implementation, Epetra_BlockMap is the base class for Epetra_Map.

*/

#include "Epetra_ConfigDefs.h"
#include "Epetra_BlockMap.h"

class EPETRA_LIB_DLL_EXPORT Epetra_Map : public Epetra_BlockMap {
    
  public:

  //! Epetra_Map constructor for a Epetra-defined uniform linear distribution of elements.
  /*! Creates a map that distributes NumGlobalElements elements evenly across all processors in the
      Epetra_Comm communicator. If NumGlobalElements does not divide exactly into the number of processors,
      the first processors in the communicator get one extra element until the remainder is gone.

    \param In
            NumGlobalElements - Number of elements to distribute.
    
    \param In
            IndexBase - Minimum index value used for arrays that use this map.  Typically 0 for
      C/C++ and 1 for Fortran.
      
    \param In
            Comm - Epetra_Comm communicator containing information on the number of
      processors.

    \return Pointer to a Epetra_Map object.

  */ 
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
  Epetra_Map(int NumGlobalElements, int IndexBase, const Epetra_Comm& Comm);
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
  Epetra_Map(long long NumGlobalElements, int IndexBase, const Epetra_Comm& Comm);
  Epetra_Map(long long NumGlobalElements, long long IndexBase, const Epetra_Comm& Comm);
#endif





  //! Epetra_Map constructor for a user-defined linear distribution of elements.
  /*! Creates a map that puts NumMyElements on the calling processor. If
        NumGlobalElements=-1, the number of global elements will be
        the computed sum of NumMyElements across all processors in the
        Epetra_Comm communicator.

    \param In
            NumGlobalElements - Number of elements to distribute.  Must be
     either -1 or equal to the computed sum of NumMyElements across all
     processors in the Epetra_Comm communicator.

    \param In
            NumMyElements - Number of elements owned by the calling processor.
    
    \param In
            IndexBase - Minimum index value used for arrays that use this map.  Typically 0 for
      C/C++ and 1 for Fortran.
      
    \param In
            Comm - Epetra_Comm communicator containing information on the number of
      processors.

    \return Pointer to a Epetra_Map object.

  */ 
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
  Epetra_Map(int NumGlobalElements, int NumMyElements, int IndexBase, const Epetra_Comm& Comm);
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
  Epetra_Map(long long NumGlobalElements, int NumMyElements, int IndexBase, const Epetra_Comm& Comm);
  Epetra_Map(long long NumGlobalElements, int NumMyElements, long long IndexBase, const Epetra_Comm& Comm);
#endif





  //! Epetra_Map constructor for a user-defined arbitrary distribution of elements.
  /*! Creates a map that puts NumMyElements on the calling processor. The indices of the elements
      are determined from the list MyGlobalElements.  If 
       NumGlobalElements=-1, the number of global elements will be 
       the computed sum of NumMyElements across all processors in the
       Epetra_Comm communicator.

    \param In
            NumGlobalElements - Number of elements to distribute.  Must be
     either -1 or equal to the computed sum of NumMyElements across all
     processors in the Epetra_Comm communicator.

    \param In
            NumMyElements - Number of elements owned by the calling processor.
    
    \param In
            MyGlobalElements - Integer array of length NumMyElements.  The ith entry contains the
      global index value of the ith element on this processor.  Index values are not required to
      be contiguous on a processor, or to be within the range of 0 to NumGlobalElements.  As
      long as the index values are consistently defined and used, any set of NumGlobalElements
      distinct integer values is acceptable.

    \param In
            IndexBase - Minimum index value used for arrays that use this map.  Typically 0 for
      C/C++ and 1 for Fortran.
      
    \param In
            Comm - Epetra_Comm communicator containing information on the number of
      processors.

    \return Pointer to a Epetra_Map object.

  */ 
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
  Epetra_Map(int NumGlobalElements, int NumMyElements,
             const int *MyGlobalElements,
             int IndexBase, const Epetra_Comm& Comm);
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
  Epetra_Map(long long NumGlobalElements, int NumMyElements,
             const long long *MyGlobalElements,
             int IndexBase, const Epetra_Comm& Comm);
  Epetra_Map(long long NumGlobalElements, int NumMyElements,
             const long long *MyGlobalElements,
             long long IndexBase, const Epetra_Comm& Comm);
#endif
  
#if defined(EPETRA_NO_32BIT_GLOBAL_INDICES) && defined(EPETRA_NO_64BIT_GLOBAL_INDICES)
  // default implementation so that no compiler/linker error in case neither 32 nor 64
  // bit indices present.
  Epetra_Map() {}
#endif

 //! Epetra_Map constructor for a user-defined arbitrary distribution of elements, where the user provides all the globals.
  /*! \warning This method is intended for expert developer use only, and should never be called by user code.
   */
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
  Epetra_Map(long long NumGlobal_Elements, int NumMy_Elements,
	     const long long * myGlobalElements, 
	     int indexBase,
	     const Epetra_Comm& comm,
	     bool UserIsDistributedGlobal,
	     long long UserMinAllGID, long long UserMaxAllGID);
#endif
  
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
  Epetra_Map(int NumGlobal_Elements, int NumMy_Elements,
	     const int * myGlobalElements, 
	     int indexBase,
	     const Epetra_Comm& comm,
	     bool UserIsDistributedGlobal,
	     int UserMinAllGID, int UserMaxAllGID);
#endif

 
  //! Epetra_Map copy constructor.
  Epetra_Map(const Epetra_Map& map);
  
  //! Epetra_Map destructor.
  virtual ~Epetra_Map(void);

  //! Assignment Operator
  Epetra_Map & operator=(const Epetra_Map & map);
  
  /// \brief Return a new BlockMap with processes with zero elements removed.
  ///
  /// \warning This method is only for expert users.  Understanding
  ///   how to use this method correctly requires some familiarity
  ///   with semantics of MPI communicators.
  ///
  /// \warning We make no promises of backwards compatibility for
  ///   this method.  It may go away or change at any time.
  ///
  /// This method first computes a new communicator, which contains
  /// only those processes in this Map's communicator (the "original
  /// communicator") that have a nonzero number of elements in this
  /// BlockMap (the "original BlockMap").  It then returns a new BlockMap
  /// distributed over the new communicator.  The new BlockMap represents
  /// the same distribution as the original BlockMap, except that
  /// processes containing zero elements are not included in the new
  /// BlockMap or its communicator.  On processes not included in the new
  /// BlockMap or communicator, this method returns NULL.
  ///
  /// The returned BlockMap always has a distinct communicator from this
  /// BlockMap's original communicator.  The new communicator contains a
  /// subset of processes from the original communicator.  Even if
  /// the number of processes in the new communicator equals the
  /// number of processes in the original communicator, the new
  /// communicator is distinct.  (In an MPI implementation, the new
  /// communicator is created using MPI_Comm_split.)
  ///
  /// This method must be called collectively on the original
  /// communicator.  It leaves the original Map and communicator
  /// unchanged.
  ///
  /// This method was intended for applications such as algebraic
  /// multigrid or other multilevel preconditioners.  Construction
  /// of each level of the multilevel preconditioner typically
  /// requires constructing sparse matrices, which in turn requires
  /// all-reduces over all participating processes at that level.
  /// Matrix sizes at successively coarser levels shrink
  /// geometrically.  At the coarsest levels, some processes might
  /// be left with zero rows of the matrix, or the multigrid
  /// implementation might "rebalance" (redistribute the matrix) and
  /// intentionally leave some processes with zero rows.  Removing
  /// processes with zero rows makes the all-reduces and other
  /// communication operations cheaper.
  Epetra_Map * RemoveEmptyProcesses() const;

  /// \brief Replace this Map's communicator with a subset communicator.
  ///
  /// \warning This method is only for expert users.  Understanding
  ///   how to use this method correctly requires some familiarity
  ///   with semantics of MPI communicators.
  ///
  /// \warning We make no promises of backwards compatibility for
  ///   this method.  It may go away or change at any time.
  ///
  /// \pre The input communicator's processes are a subset of this
  ///   Map's current communicator's processes.
  /// \pre On processes which are not included in the input
  ///   communicator, the input communicator is null.
  ///
  /// This method must be called collectively on the original
  /// communicator.  It leaves the original BlockMap and communicator
  /// unchanged.
  ///
  /// \note This method differs from removeEmptyProcesses(), in that
  ///   it does not assume that excluded processes have zero
  ///   entries.  For example, one might wish to remove empty
  ///   processes from the row Map of a CrsGraph using
  ///   removeEmptyProcesses(), and then apply the resulting subset
  ///   communicator to the column, domain, and range Maps of the
  ///   same graph.  For the latter three Maps, one would in general
  ///   use this method instead of removeEmptyProcesses(), giving
  ///   the new row Map's communicator to this method.
  Epetra_Map* ReplaceCommWithSubset(const Epetra_Comm * Comm) const;


};

#endif /* EPETRA_MAP_H */
