#ifndef _PETRA_LOCALMAP_H_
#define _PETRA_LOCALMAP_H_

//! Petra_LocalMap: A class for replicating vectors and matrices across multiple processors.

/*! Small matrix and vector objects are often replicated on distributed memory
  parallel machines. The Petra_LocalMap class allows construction of these replicated
  local objects and keeps information that describes 
  this distribution.  

  Petra_LocalMap allows the storage and retrieval of the following information.  
  Once a Petra_Map is constructed any of the following attributes can 
  be obtained
  by calling a query function that has the name as the attribute, e.g. to get the
  value of NumGlobalEquations, you can call a function NumGlobalElements().
  For attributes that
  are lists, the query functions return the list values in a user allocated array.  


  <ul>
  <li> NumMyElements - The number of elements owned by the calling processor.
  <li> IndexBase - The base integer value for indexed array references.  Typically this is 0
       for C/C++ and 1 for Fortran, but it can be set to any integer value.
  <li> Comm - The Petra_Comm communicator.  This communicator can in turn be queried for
       processor rank and size information.
  </ul>

  The Petra_LocalMap class is actually a derived class of Petra_Map.  Petra_Map is in turn derived
  from Petra_BlockMap.  As such,  Petra_LocalMap has full access to all the functions in these other
  map classes.

  In particular, the following function allows a boolean test:    

  <ul>
  <li> DistributedGlobal() - Returns false for a Petra_LocalMap object.
  </ul>

  \warning A Petra_Comm object is required for all Petra_LocalMap constructors.

  \internal In the current implementation, Petra_Map is the base class for Petra_LocalMap.

*/
#include "Petra_Map.h"

class Petra_LocalMap : public Petra_Map {
    
  public:
  //! Petra_LocalMap constructor for a user-defined replicate distribution of elements.
  /*! Creates a map that puts NumMyElements on the calling processor. Each processor should
      pass in the same value for NumMyElements.

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
    Petra_LocalMap(int NumMyElements, int IndexBase, const Petra_Comm& Comm);

  //! Petra_LocalMap copy constructor.
  
    Petra_LocalMap(const Petra_LocalMap& map);
  
  //! Petra_LocalMap destructor.

    virtual ~Petra_LocalMap();

  private:

    int CheckInput();

};
#endif /* _PETRA_LOCALMAP_H_ */
