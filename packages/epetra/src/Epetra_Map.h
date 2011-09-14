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
  Epetra_Map(int NumGlobalElements, int IndexBase, const Epetra_Comm& Comm);





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
  Epetra_Map(int NumGlobalElements, int NumMyElements, int IndexBase, const Epetra_Comm& Comm);





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
  Epetra_Map(int NumGlobalElements, int NumMyElements,
             const int *MyGlobalElements,
             int IndexBase, const Epetra_Comm& Comm);
  
	//! Epetra_Map copy constructor.
  Epetra_Map(const Epetra_Map& map);
  
  //! Epetra_Map destructor.
  virtual ~Epetra_Map(void);

	//! Assignment Operator
	Epetra_Map & operator=(const Epetra_Map & map);
  
};

#endif /* EPETRA_MAP_H */
