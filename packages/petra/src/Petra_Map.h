#ifndef _PETRA_MAP_H_
#define _PETRA_MAP_H_

//! Petra_Map: A class for partitioning vectors and matrices.

/*! It is often the case that multiple matrix and vector objects have an identical distribution 
  of elements on a parallel machine. The Petra_Map class keep information that describes 
  this distribution for matrices and vectors.  

  Petra_Map allows the storage and retrieval of the following information.  Depending on the
  constructor that is used, some of the information is defined by the user and some is 
  determined by the constructor.  Once a Petra_Map is constructed any of the following attributes can 
  be obtained
  by calling a query function that has the name as the attribute, e.g. to get the
  value of NumGlobalEquations, you can call a function NumGlobalElements().  For attributes that
  are lists, the query functions return the list values in a user allocated array.

  <ul>
  <li> NumGlobalElements - The total number of elements across all processors. If this parameter and
       NumMyElements are both passed in to the constructor, one of the three cases will apply: 
       <ol> 
       <li> If NumGlobalElements = NumMyElements (and not equal to zero)
            the map is defined to be a local replicated
            map.  In this case, objects constructed using this map will be identically replicated across
	    all processors in the communicator.
       <li> If NumGlobalElements = 0 and NumMyElements is passed in then NumGlobalElements will
            be computed as the sum of NumMyElements across all processors.
       <li> If neither of the above is true, NumGlobalElements will be checked against the sum of 
            NumMyElements across all processors.  An error is issued if the comparison is not equal.
       </ol>
  <li> NumMyElements - The number of elements owned by the calling processor.
  <li> MyGlobalElements - A list of length NumMyElements that contains the global element IDs
       of the elements owned by the calling processor.
  <li> IndexBase - The base integer value for indexed array references.  Typically this is 0
       for C/C++ and 1 for Fortran, but it can be set to any integer value.
  <li> Comm - The Petra_Comm communicator.  This communicator can in turn be queried for
       processor rank and size information.
  </ul>


  In addition to the information above that is passed in to or created by the Petra_Map constructor,
  the following attributes are computed and available via query to the user using the same scheme
  as above, e.g., use NumGlobalEquations() to get the value of NumGlobalEquations.

  <ul>
  <li> NumGlobalEquations - The total number of equations across all processors.
  <li> NumMyEquations - The number of equations on the calling processor.
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
       This will be true in most cases, but will be false on in serial and for objects
       that are created via the derived Petra_LocalMap class.
  </ul>

  \warning A Petra_Comm object is required for all Petra_Map constructors.

  \internal In the current implementation, Petra_BlockMap is the base class for Petra_Map.

*/

#include "Petra_Petra.h"
#include "Petra_Comm.h"
#include "Petra_BlockMap.h"

class Petra_Map : public Petra_BlockMap {
    
  public:

  //! Petra_Map constructor for a Petra-defined uniform linear distribution of elements.
  /*! Creates a map that distributes NumGlobalElements elements evenly across all processors in the
      Petra_Comm communicator. If NumGlobalElements does not divide exactly into the number of processors,
      the first processors in the communicator get one extra element until the remainder is gone.

    \param In
            NumGlobalElements - Number of elements to distribute.
    
    \param In
            IndexBase - Minimum index value used for arrays that use this map.  Typically 0 for
	    C/C++ and 1 for Fortran.
	    
    \param In
            Comm - Petra_Comm communicator containing information on the number of
	    processors.

    \return Pointer to a Petra_Map object.

  */ 
  Petra_Map(int NumGlobalElements, int IndexBase, const Petra_Comm& Comm);





  //! Petra_Map constructor for a user-defined linear distribution of elements.
  /*! Creates a map that puts NumMyElements on the calling processor. NumGlobalElements will be the
      computed as the sum of NumMyElements across all processors in the Petra_Comm communicator.

    \param In
            NumGlobalElements - Number of elements to distribute.
    \param In
            NumMyElements - Number of elements owned by the calling processor.
    
    \param In
            IndexBase - Minimum index value used for arrays that use this map.  Typically 0 for
	    C/C++ and 1 for Fortran.
	    
    \param In
            Comm - Petra_Comm communicator containing information on the number of
	    processors.

    \return Pointer to a Petra_Map object.

  */ 
  Petra_Map(int NumGlobalElements, int NumMyElements, int IndexBase, const Petra_Comm& Comm);





  //! Petra_Map constructor for a user-defined arbitrary distribution of elements.
  /*! Creates a map that puts NumMyElements on the calling processor. The indices of the elements
      are determined from the list MyGlobalElements.  NumGlobalElements will be the
      computed as the sum of NumMyElements across all processors in the Petra_Comm communicator.

    \param In
            NumGlobalElements - Number of elements to distribute.
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
            Comm - Petra_Comm communicator containing information on the number of
	    processors.

    \return Pointer to a Petra_Map object.

  */ 
  Petra_Map(int NumGlobalElements, int NumMyElements, int *MyGlobalElements,  
	    int IndexBase, const Petra_Comm& Comm);
  //! Petra_Map copy constructor.
  
  Petra_Map(const Petra_Map& map);
  
  //! Petra_Map destructor.
  
  virtual ~Petra_Map(void);
  
};

#endif /* _PETRA_MAP_H_ */
