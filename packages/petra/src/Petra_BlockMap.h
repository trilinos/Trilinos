#ifndef _PETRA_BLOCKMAP_H_
#define _PETRA_BLOCKMAP_H_

//! Petra_BlockMap: A class for partitioning block element vectors and matrices.

/*! It is often the case that multiple matrix and vector objects have an identical distribution 
  of elements on a parallel machine. The Petra_BlockMap class keep information that describes 
  this distribution for matrices and vectors that have block elements.  The definition of an 
  element can vary depending on the situation.  For
  vectors (and multi-vectors), an element is a span of one or more contiguous entries.
  For matrices, it is a span of one or more matrix rows.

  This class has a variety of constructors that can be separated into two categories:
  <ol>
  <li> Fixed block size constructors:
       All map elements have an identical block size.
       This corresponds to a block partitioning of matrices and vectors where the block
       size is the same for all blocks. A common example is multiple degrees of freedom
       per mesh node in finite element computations where the number of degrees of
       freedom is the same for all nodes.
  <li> Variable block size constructor:
       Map element sizes may vary and are individually defined via a list of element sizes.
       This is the most general case and corresponds to a variable block partitioning of the
       matrices and vectors. A common example is 
       multiple degrees of freedom per mesh node in finite element computations where the
       number of degrees of freedom varies.  This happens, for example, if regions have differing
       material types or there are chemical reactions in the simulation.
  </ol>

  Petra_BlockMap allows the storage and retrieval of the following information.  Depending on the
  constructor that is used, some of the information is defined by the user and some is 
  determined by the constructor.  Once a Petra_BlockMap is constructed any of the following can 
  be obtained
  by calling a query function that has the same name as the attribute, e.g. to get the
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
  <li> ElementSize - The size of elements if the size of all elements is the same.
       This will be the case if the query function ConstantElementSize() returns true.
       Otherwise this value will be set to zero.
  <li> ElementSizeList - A list of the element sizes for elements owned by the calling
       processor.  This list is always accessible, even if the element sizes are all one
       or of constant value.  However, in these cases, the ElementSizeList will not be 
       generated unless a query for the list is called.
  <li> IndexBase - The base integer value for indexed array references.  Typically this is 0
       for C/C++ and 1 for Fortran, but it can be set to any integer value.
  <li> Comm - The Petra_Comm communicator.  This communicator can in turn be queried for
       processor rank and size information.
  </ul>


  In addition to the information above that is passed in to or created by the Petra_BlockMap constructor,
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
  <li> MinElementSize - The minimum element size across all processors.
  <li> MaxElementSize - The maximum element size across all processors.
  </ul>

  The following functions allow boolean tests for certain properties.    

  <ul>
  <li> ConstantElementSize() - Returns true if the element size for this map is the same 
       for all elements.
  <li> LinearMap() - Returns true if the elements are distributed linear across processors, i.e.,
       processor 0 gets the first n/p elements, processor 1 gets the next n/p elements, etc. where
       n is the number of elements and p is the number of processors.
  <li> DistributedGlobal() - Returns true if the element space of the map spans more than one processor.
       This will be true in most cases, but will be false on in serial and for objects
       that are created via the derived Petra_LocalMap class.
  </ul>

  \warning A Petra_Comm object is required for all Petra_BlockMap constructors.

  \internal 
    {
    In the current implementation, Petra_BlockMap is the base class for:
    <ul>
    <li> Petra_Map.
    <li> Petra_LocalBlockMap.
    </ul>
    }

*/

#include "Petra_Petra.h"
#include "Petra_Comm.h"
#include "Petra_Directory.h"

class Petra_BlockMap {
    
  // Give ostream << function some access to private and protected data/functions.

  friend ostream& operator << (ostream& os, const Petra_BlockMap & Map);
  public:

  //! Petra_BlockMap constructor for a Petra-defined uniform linear distribution of constant block size elements.
  /*! Creates a map that distributes NumGlobalElements elements evenly across all processors in the
      Petra_Comm communicator. If NumGlobalElements does not divide exactly into the number of processors,
      the first processors in the communicator get one extra element until the remainder is gone.

      The elements are defined to have a constant fixed block size specified by ElementSize.

    \param In
            NumGlobalElements - Number of elements to distribute.
    
    \param In
            ElementSize - Number of equations or vector entries per element.

    \param In
            IndexBase - Minimum index value used for arrays that use this map.  Typically 0 for
	    C/C++ and 1 for Fortran.
	    
    \param In
            Comm - Petra_Comm communicator containing information on the number of
	    processors.

    \return Pointer to a Petra_BlockMap object.

  */ 
  Petra_BlockMap(int NumGlobalElements, int ElementSize, int IndexBase, const Petra_Comm& Comm);





  //! Petra_BlockMap constructor for a user-defined linear distribution of constant block size elements.
  /*! Creates a map that puts NumMyElements on the calling processor. NumGlobalElements will be the
      computed as the sum of NumMyElements across all processors in the Petra_Comm communicator.

      The elements are defined to have a constant fixed block size specified by ElementSize.

    \param In
            NumGlobalElements - Number of elements to distribute.
    \param In
            NumMyElements - Number of elements owned by the calling processor.
    
    \param In
            ElementSize - Number of equations or vector entries per element.

    \param In
            IndexBase - Minimum index value used for arrays that use this map.  Typically 0 for
	    C/C++ and 1 for Fortran.
	    
    \param In
            Comm - Petra_Comm communicator containing information on the number of
	    processors.

    \return Pointer to a Petra_BlockMap object.

  */ 
  Petra_BlockMap(int NumGlobalElements, int NumMyElements, 
	    int ElementSize, int IndexBase, const Petra_Comm& Comm);





  //! Petra_BlockMap constructor for a user-defined arbitrary distribution of constant block size elements.
  /*! Creates a map that puts NumMyElements on the calling processor. The indices of the elements
      are determined from the list MyGlobalElements.  NumGlobalElements will be the
      computed as the sum of NumMyElements across all processors in the Petra_Comm communicator.

      The elements are defined to have a constant fixed block size specified by ElementSize.

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
            ElementSize - Number of equations or vector entries per element.

    \param In
            IndexBase - Minimum index value used for arrays that use this map.  Typically 0 for
	    C/C++ and 1 for Fortran.
	    
    \param In
            Comm - Petra_Comm communicator containing information on the number of
	    processors.

    \return Pointer to a Petra_BlockMap object.

  */ 
  Petra_BlockMap(int NumGlobalElements, int NumMyElements, int *MyGlobalElements,  
	    int ElementSize, int IndexBase, const Petra_Comm& Comm);



  //! Petra_BlockMap constructor for a user-defined arbitrary distribution of variable block size elements.
  /*! Creates a map that puts NumMyElements on the calling processor. NumGlobalElements will be the
      computed as the sum of NumMyElements across all processors in the Petra_Comm communicator.

      The elements are defined to have a variable block size defined by ElementSizeList.

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
            ElementSizeList - A list of the element sizes for elements owned by the calling
	    processor. The ith entry contains the element size of the ith element on this processor.

    \param In
            IndexBase - Minimum index value used for arrays that use this map.  Typically 0 for
	    C/C++ and 1 for Fortran.
	    
    \param In
            Comm - Petra_Comm communicator containing information on the number of
	    processors.

    \return Pointer to a Petra_BlockMap object.

  */ 
  Petra_BlockMap(int NumGlobalElements, int NumMyElements, int *MyGlobalElements, 
	    int *ElementSizeList, int IndexBase, const Petra_Comm& Comm);
  
  //! Petra_BlockMap copy constructor.
  
  Petra_BlockMap(const Petra_BlockMap& map);
  
  //! Petra_BlockMap destructor.
  
  virtual ~Petra_BlockMap(void);

  
  //! Returns the processor IDs and corresponding local index value for a given list of global indices
  /*! For each element (GID) of a given a list of global element numbers (stored in GIDList) of length NumIDs,
      this function returns (in PIDList) the with processor that owns the GID for this map and returns the
      local index (in LIDList) of the GID on that processor.
  */
  int RemoteIDList(int NumIDs, const int * GIDList, int * PIDList, int * LIDList) const {
    return(RemoteIDList(NumIDs, GIDList, PIDList, LIDList, 0));};

  //! Returns the processor IDs, corresponding local index value, and element size for a given list of global indices
  /*! For each element (GID) of a given a list of global element numbers (stored in GIDList) of length NumIDs,
      this function returns (in PIDList) the with processor that owns the GID for this map and returns the
      local index (in LIDList) of the GID on that processor.  Finally it returns the element sizes in
      SizeList.
  */
  int RemoteIDList(int NumIDs, const int * GIDList, int * PIDList, int * LIDList, int * SizeList) const;

  //! Returns local ID of global ID, return -1 if not found on this processor.
  int  LID(int GID) const;
  
  //! Returns global ID of local ID, return IndexBase-1 if not found on this processor.
  int  GID(int LID) const; 
  
  //! Returns the LID of the block that contains the given EquationID, and the Offset of the equation in that block.
  int FindLocalBlockID(int EquationID, int & BlockID, int & BlockOffset)  const;

  //! Returns true if the GID passed in belongs to the calling processor in this map, otherwise returns false.
  bool  MyGID(int GID) const {return(LID(GID)!=-1);};
   
  //! Returns true if the LID passed in belongs to the calling processor in this map, otherwise returns false.
  bool  MyLID(int LID) const {return(GID(LID)!=IndexBase_-1);};

  //! Access function for NumGlobalElements.
  int  NumGlobalElements() const {return(NumGlobalElements_);};
  
  //! Access function for NumMyElements.
  int  NumMyElements() const {return(NumMyElements_);};
  
  //! Access function for MyGlobalElements.
  int MyGlobalElements(int * MyGlobalElements) const;
  
  //! Access function for ElementSize.
  int  ElementSize() const {return(ElementSize_);};
    
  //! Size of element for specified LID.
  int  ElementSize(int LID) const;

  //! Access function for ElementSizeList.
  int ElementSizeList(int * ElementSizeList)const;
  
  //! Access function for FirstElementEntryList.
  int FirstElementEntryList(int * FirstElementEntryList)const;
  
  //! Access function for IndexBase.
  int  IndexBase() const {return(IndexBase_);};
  
  //! Access function for NumGlobalEquations.
  int  NumGlobalEquations() const {return(NumGlobalEquations_);};
  
  //! Access function for NumMyEquations.
  int  NumMyEquations() const {return(NumMyEquations_);};
  
  //! Access function for MinAllGID.
  int  MinAllGID() const {return(MinAllGID_);};
  
  //! Access function for MaxAllGID.
  int  MaxAllGID() const {return(MaxAllGID_);};
  
  //! Access function for MinMyGID.
  int  MinMyGID() const {return(MinMyGID_);};
  
  //! Access function for MaxMyGID.
  int  MaxMyGID() const {return(MaxMyGID_);};
  
  //! Access function for MinLID.
  int  MinLID() const {return(MinLID_);};
  
  //! Access function for MaxLID.
  int  MaxLID() const {return(MaxLID_);};
  
  //! Access function for MinMyElementSize.
  int  MinMyElementSize() const {return(MinMyElementSize_);};
  
  //! Access function for MaxMyElementSize.
  int  MaxMyElementSize() const {return(MaxMyElementSize_);};
  
  //! Access function for MinElementSize.
  int  MinElementSize() const {return(MinElementSize_);};
  
  //! Access function for MaxElementSize.
  int  MaxElementSize() const {return(MaxElementSize_);};
  
  //! Access function for ConstantElementSize.
  bool  ConstantElementSize() const {return(ConstantElementSize_);};

  //! Returns true if \e this and Map are identical maps
  bool SameAs(const Petra_BlockMap & Map) const;
  
  //! Access function for LinearMap.
  bool  LinearMap() const {return(LinearMap_);};

  //! Access function for DistributedGlobal.
  bool  DistributedGlobal() const {return(DistributedGlobal_);};

  //! Access function for Petra_Comm communicator.
  const Petra_Comm& Comm() const {return(*Comm_);};
  
  // Access Function for MyGlobalElements
  int * MyGlobalElements() const;

  // Access Function for FirstElementEntryList
  int * FirstElementEntryList() const;

  // Access Function for ElementSizeList
  int * ElementSizeList() const;

  // Copy EquationToBlockList into given array
  int EquationToBlockList(int * EquationToBlockList) const;

  // Access Function for EquationToBlockList
  int * EquationToBlockList() const;
  

  friend class Petra_Directory;
  friend class Petra_LocalMap;
  
 private: // These need to be accessible to derived map classes.
  
  void GlobalToLocalSetup();
  int *LID_;
  int NumGlobalElements_;
  int NumMyElements_;
  int* MyGlobalElements_;
  int* FirstElementEntryList_;
  int ElementSize_;
  int* ElementSizeList_;
  int* EquationToBlockList_;
  int IndexBase_;
  const Petra_Comm * Comm_;
  Petra_Directory * Directory_;
  
  int NumGlobalEquations_;
  int NumMyEquations_;
  int MinAllGID_;
  int MaxAllGID_;
  int MinMyGID_;
  int MaxMyGID_;
  int MinLID_;
  int MaxLID_;
  int MinMyElementSize_;
  int MaxMyElementSize_;
  int MinElementSize_;
  int MaxElementSize_;
  
  bool ConstantElementSize_;
  bool LinearMap_;
  bool DistributedGlobal_;
};

//! << operator will work for Petra_RDP_MultiVectors.
ostream& operator << (ostream& os, const Petra_BlockMap & Map);
#endif /* _PETRA_BLOCKMAP_H_ */
