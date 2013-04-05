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

#ifndef EPETRA_LOCALMAP_H
#define EPETRA_LOCALMAP_H

//! Epetra_LocalMap: A class for replicating vectors and matrices across multiple processors.

/*! Small matrix and vector objects are often replicated on distributed memory
  parallel machines. The Epetra_LocalMap class allows construction of these replicated
  local objects and keeps information that describes 
  this distribution.  

  Epetra_LocalMap allows the storage and retrieval of the following information.  
  Once a Epetra_Map is constructed any of the following attributes can 
  be obtained
  by calling a query function that has the name as the attribute, e.g. to get the
  value of NumGlobalPoints, you can call a function NumGlobalElements().
  For attributes that
  are lists, the query functions return the list values in a user allocated array.  


  <ul>
  <li> NumMyElements - The number of elements owned by the calling processor.
  <li> IndexBase - The base integer value for indexed array references.  Typically this is 0
       for C/C++ and 1 for Fortran, but it can be set to any integer value.
  <li> Comm - The Epetra_Comm communicator.  This communicator can in turn be queried for
       processor rank and size information.
  </ul>

  The Epetra_LocalMap class is actually a derived class of Epetra_Map.  Epetra_Map is in turn derived
  from Epetra_BlockMap.  As such,  Epetra_LocalMap has full access to all the functions in these other
  map classes.

  In particular, the following function allows a boolean test:    

  <ul>
  <li> DistributedGlobal() - Returns false for a Epetra_LocalMap object.
  </ul>

  \warning A Epetra_Comm object is required for all Epetra_LocalMap constructors.

  \internal In the current implementation, Epetra_Map is the base class for Epetra_LocalMap.

*/
#include "Epetra_ConfigDefs.h"
#include "Epetra_Map.h"

class EPETRA_LIB_DLL_EXPORT Epetra_LocalMap : public Epetra_Map {
    
  public:
  //! Epetra_LocalMap constructor for a user-defined replicate distribution of elements.
  /*! Creates a map that puts NumMyElements on the calling processor. Each processor should
      pass in the same value for NumMyElements.

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
	Epetra_LocalMap(int NumMyElements, int IndexBase, const Epetra_Comm& Comm);
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
	Epetra_LocalMap(long long NumMyElements, int IndexBase, const Epetra_Comm& Comm);
	Epetra_LocalMap(long long NumMyElements, long long IndexBase, const Epetra_Comm& Comm);
#endif
	
  //! Epetra_LocalMap copy constructor.
  
	Epetra_LocalMap(const Epetra_LocalMap& map);
  
  //! Epetra_LocalMap destructor.
	
	virtual ~Epetra_LocalMap();
	
	//! Assignment Operator
	Epetra_LocalMap & operator=(const Epetra_LocalMap & map);
	
 private:
	
	int CheckInput();
	
};
#endif /* EPETRA_LOCALMAP_H */
