
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

#ifndef _EPETRA_BLOCKMAP_H_
#define _EPETRA_BLOCKMAP_H_

#include "Epetra_Object.h"

class Epetra_Comm;
class Epetra_Directory;

class Epetra_HashTable;


//! Epetra_BlockMap: A class for partitioning block element vectors and matrices.

/*! It is often the case that multiple matrix and vector objects have an identical distribution 
  of elements on a parallel machine. The Epetra_BlockMap class keeps information that describes 
  this distribution for matrices and vectors that have block elements.  The definition of an 
  element can vary depending on the situation.  For
  vectors (and multi-vectors), an element is a span of one or more contiguous entries.
  For matrices, it is a span of one or more matrix rows.  More generally, an element in the BlockMap class is an ordered list of points.  (NOTE: Points do not have global ID's.)  Two additional definitions useful in understanding the BlockMap class follow:
  <ul>
  <li> BlockMap - A distributed ordered list of elements.
  <li> First Point - First ordered point in an element
  </ul>

  This class has a variety of constructors that can be separated into two categories:
  <ol>
  <li> Fixed element size constructors:
       All map elements have an identical size.
       This corresponds to a block partitioning of matrices and vectors where the element
       size is the same for all elements. A common example is multiple degrees of freedom
       per mesh node in finite element computations where the number of degrees of
       freedom is the same for all nodes.
  <li> Variable element size constructor:
       Map element sizes may vary and are individually defined via a list of element sizes.
       This is the most general case and corresponds to a variable block partitioning of the
       matrices and vectors. A common example is 
       multiple degrees of freedom per mesh node in finite element computations where the
       number of degrees of freedom varies.  This happens, for example, if regions have differing
       material types or there are chemical reactions in the simulation.
  </ol>

  Epetra_BlockMap allows the storage and retrieval of the following information.  Depending on the
  constructor that is used, some of the information is defined by the user and some is 
  determined by the constructor.  Once an Epetra_BlockMap is constructed any of the following can 
  be obtained
  by calling a query function that has the same name as the attribute, e.g. to get the
  value of NumGlobalElements, you can call a function NumGlobalElements().  For attributes that
  are lists, the query functions return the list values in a user allocated array.

  <ul>
  <li> NumGlobalElements - The total number of elements across all processors. If this parameter and
       NumMyElements are both passed in to the constructor, one of the three cases will apply: 
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
  <li> ElementSize - The size of elements if the size of all elements is the same.
       This will be the case if the query function ConstantElementSize() returns true.
       Otherwise this value will be set to zero.
  <li> ElementSizeList - A list of the element sizes for elements owned by the calling
       processor.  This list is always accessible, even if the element sizes are all one
       or of constant value.  However, in these cases, the ElementSizeList will not be 
       generated unless a query for the list is called.
  <li> IndexBase - The base integer value for indexed array references.  Typically this is 0
       for C/C++ and 1 for Fortran, but it can be set to any integer value.
  <li> Comm - The Epetra_Comm communicator.  This communicator can in turn be queried for
       processor rank and size information.
  </ul>


  In addition to the information above that is passed in to or created by the Epetra_BlockMap constructor,
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
       that are created via the derived Epetra_LocalMap class.
  </ul>

  \warning A Epetra_Comm object is required for all Epetra_BlockMap constructors.

  \bf {error handling}

  Most methods in Epetra_BlockMap return an integer error code.  If the error code is 0, then no error occurred.
  If > 0 then a warning error occurred.  If < 0 then a fatal error occurred.

  Epetra_BlockMap constructors will throw an exception of an error occurrs.  These exceptions will alway be negative integer values
  as follows:
  <ol>
  <li> -1  NumGlobalElements < -1.  Should be >= -1 (Should be >= 0 for first BlockMap constructor).
  <li> -2  NumMyElements < 0.  Should be >= 0.
  <li> -3  ElementSize <= 0. Should be > 0.
  <li> -4  Invalid NumGlobalElements.  Should equal sum of MyGlobalElements, or set to -1 to compute automatically.
  <li> -5  Minimum global element index is less than index base.
  <li> -99 Internal Epetra_BlockMap error.  Contact developer.
  </ol>

  For robust code, Epetra_BlockMap constructor calls should be caught using the try {...} catch {...} mechanism.  For example:

\verbatim
  try {

    Epetra_BlockMap * map = new Epetra_BlockMap(NumGlobalElements, ElementSize, IndexBase, Comm);
  }
  catch (int Error) {
    if (Error==-1) { // handle error }
    if (Error==-2) ...
\endverbatim
 

  \note
    {
    In the current implementation, Epetra_BlockMap is the base class for:
    <ul>
    <li> Epetra_Map.
    <li> Epetra_LocalBlockMap.
    </ul>
    }

*/

class Epetra_BlockMap: public Epetra_Object {
    
  public:
  //@{ \name Constructors/destructors.
  //! Epetra_BlockMap constructor for a Epetra-defined uniform linear distribution of constant size elements.
  /*! Creates a map that distributes NumGlobalElements elements evenly across all processors in the
      Epetra_Comm communicator. If NumGlobalElements does not divide exactly into the number of processors,
      the first processors in the communicator get one extra element until the remainder is gone.

      The elements are defined to have a constant fixed size specified by ElementSize.

    \param In
            NumGlobalElements - Number of elements to distribute.
    
    \param In
            ElementSize - Number of points or vector entries per element.

    \param In
            IndexBase - Minimum index value used for arrays that use this map.  Typically 0 for
	    C/C++ and 1 for Fortran.
	    
    \param In
            Comm - Epetra_Comm communicator containing information on the number of
	    processors.

    \return Pointer to a Epetra_BlockMap object.

  */ 
  Epetra_BlockMap(int NumGlobalElements, int ElementSize, int IndexBase, const Epetra_Comm& Comm);

  //! Epetra_BlockMap constructor for a user-defined linear distribution of constant size elements.
  /*! Creates a map that puts NumMyElements on the calling processor. NumGlobalElements will be the
      computed sum of NumMyElements across all processors in the Epetra_Comm communicator.

      The elements are defined to have a constant fixed size specified by ElementSize.

    \param In
            NumGlobalElements - Number of elements to distribute.
    \param In
            NumMyElements - Number of elements owned by the calling processor.
    
    \param In
            ElementSize - Number of points or vector entries per element.

    \param In
            IndexBase - Minimum index value used for arrays that use this map.  Typically 0 for
	    C/C++ and 1 for Fortran.
	    
    \param In
            Comm - Epetra_Comm communicator containing information on the number of
	    processors.

    \return Pointer to a Epetra_BlockMap object.

  */ 
  Epetra_BlockMap(int NumGlobalElements, int NumMyElements, 
	    int ElementSize, int IndexBase, const Epetra_Comm& Comm);





  //! Epetra_BlockMap constructor for a user-defined arbitrary distribution of constant size elements.
  /*! Creates a map that puts NumMyElements on the calling processor. The indices of the elements
      are determined from the list MyGlobalElements.  NumGlobalElements will be the
      computed sum of NumMyElements across all processors in the Epetra_Comm communicator.

      The elements are defined to have a constant fixed size specified by ElementSize.

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
            ElementSize - Number of points or vector entries per element.

    \param In
            IndexBase - Minimum index value used for arrays that use this map.  Typically 0 for
	    C/C++ and 1 for Fortran.
	    
    \param In
            Comm - Epetra_Comm communicator containing information on the number of
	    processors.

    \return Pointer to a Epetra_BlockMap object.

  */ 
  Epetra_BlockMap(int NumGlobalElements, int NumMyElements, int *MyGlobalElements,  
	    int ElementSize, int IndexBase, const Epetra_Comm& Comm);



  //! Epetra_BlockMap constructor for a user-defined arbitrary distribution of variable size elements.
  /*! Creates a map that puts NumMyElements on the calling processor. NumGlobalElements will be the
      computed sum of NumMyElements across all processors in the Epetra_Comm communicator.

      The elements are defined to have a variable size defined by ElementSizeList.

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
            Comm - Epetra_Comm communicator containing information on the number of
	    processors.

    \return Pointer to a Epetra_BlockMap object.

  */ 
  Epetra_BlockMap(int NumGlobalElements, int NumMyElements, int *MyGlobalElements, 
	    int *ElementSizeList, int IndexBase, const Epetra_Comm& Comm);
  
  //! Epetra_BlockMap copy constructor.
  
  Epetra_BlockMap(const Epetra_BlockMap& map);
  
  //! Epetra_BlockMap destructor.
  
  virtual ~Epetra_BlockMap(void);
  //@}
  
  //@{ \name Local/Global ID accessor methods.
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
  
  //! Returns the LID of the element that contains the given local PointID, and the Offset of the point in that element.
  int FindLocalElementID(int PointID, int & ElementID, int & ElementOffset)  const;

  //! Returns true if the GID passed in belongs to the calling processor in this map, otherwise returns false.
  bool  MyGID(int GID) const {return(LID(GID)!=-1);};
   
  //! Returns true if the LID passed in belongs to the calling processor in this map, otherwise returns false.
  bool  MyLID(int LID) const {return(GID(LID)!=IndexBase_-1);};
  
  //!Returns the minimum global ID across the entire map.
  int  MinAllGID() const {return(MinAllGID_);};
  
  //! Returns the maximum global ID across the entire map.
  int  MaxAllGID() const {return(MaxAllGID_);};
  
  //! Returns the maximum global ID owned by this processor.
  int  MinMyGID() const {return(MinMyGID_);};
  
  //! Returns the maximum global ID owned by this processor.
  int  MaxMyGID() const {return(MaxMyGID_);};
  
  //!  The minimum local index value on the calling processor.
  int  MinLID() const {return(MinLID_);};
  
  //! The maximum local index value on the calling processor.
  int  MaxLID() const {return(MaxLID_);};
  //@}

  //@{ \name Size and dimension accessor functions
  //! Number of elements across all processors.
  int  NumGlobalElements() const {return(NumGlobalElements_);};
  
  //! Number of elements on the calling processor.
  int  NumMyElements() const {return(NumMyElements_);};
  
  //! Puts list of global elements on this processor into the user-provided array.
  int MyGlobalElements(int * MyGlobalElementList) const;
  
  //! Returns the size of elements in the map; only valid if map has constant element size.
  int  ElementSize() const {return(ElementSize_);};
    
  //! Size of element for specified LID.
  int  ElementSize(int LID) const;
    
  //! Returns the requested entry in the FirstPointInElementList; see FirstPointInElementList() for details.
  /*! This function provides similar functionality to FirstPointInElementList(), but for simple maps may avoid
      the explicit construction of the FirstPointInElementList array.  Returns -1 if LID is out-of-range.
  */
  int  FirstPointInElement(int LID) const;
  
  //! Index base for this map.
  int  IndexBase() const {return(IndexBase_);};
  
  //! Number of global points for this map; equals the sum of all element sizes across all processors.
  int  NumGlobalPoints() const {return(NumGlobalPoints_);};
  
  //! Number of global points for this map; equals the sum of all element sizes on the calling processor.
  int  NumMyPoints() const {return(NumMyPoints_);};
  
  //! Minimum element size on the calling processor.
  int  MinMyElementSize() const {return(MinMyElementSize_);};
  
  //! Maximum element size on the calling processor.
  int  MaxMyElementSize() const {return(MaxMyElementSize_);};
  
  //! Minimum element size across all processors.
  int  MinElementSize() const {return(MinElementSize_);};
  
  //! Maximum element size across all processors.
  int  MaxElementSize() const {return(MaxElementSize_);};
  //@}

  //@{ \name Miscellaneous boolean tests
  //! Returns true if map has constant element size.
  bool  ConstantElementSize() const {return(ConstantElementSize_);};

  //! Returns true if \e this and Map are identical maps
  bool SameAs(const Epetra_BlockMap & Map) const;
  
  //! Returns true if the global ID space is contiguously divided (but not necessarily uniformly) across all processors.
  bool  LinearMap() const {return(LinearMap_);};

  //! Returns true if map is defined across more than one processor.
  bool  DistributedGlobal() const {return(DistributedGlobal_);};
  //@}

  //@{ \name Array accessor functions

  //! Pointer to internal array containing list of global IDs assigned to the calling processor.
  int * MyGlobalElements() const;

  //! Pointer to internal array containing a mapping between the local elements and the first local point number in each element.
  /*! This array is a scan sum of the ElementSizeList such that the ith entry in FirstPointInElementList is the sum of the first
      i-1 entries of ElementSizeList().
  */
  int * FirstPointInElementList() const;

  //! List of the element sizes corresponding to the array MyGlobalElements().
  int * ElementSizeList() const;

  //! For each local point, indicates the local element ID that the point belongs to.
  int * PointToElementList() const;

  //! Same as ElementSizeList() except it fills the user array that is passed in.
  int ElementSizeList(int * ElementSizeList)const;
  
  //! Same as FirstPointInElementList() except it fills the user array that is passed in.
  int FirstPointInElementList(int * FirstPointInElementList)const;

  //! Same as PointToElementList() except it fills the user array that is passed in.
  int PointToElementList(int * PointToElementList) const;

  //@}

  //@{ \name Miscellaneous

  //! Print object to an output stream
  virtual void Print(ostream & os) const;

  //! Access function for Epetra_Comm communicator.
  const Epetra_Comm& Comm() const {return(*Comm_);};
  //@}

  //@{ \name Expert Users and Developers Only
  Epetra_BlockMap * NewReference() const { NumReferences_++; return((Epetra_BlockMap *) this);};
  //@}

  friend class Epetra_Directory;
  friend class Epetra_LocalMap;
  
 private: // These need to be accessible to derived map classes.
  
  void GlobalToLocalSetup();
  bool IsDistributedGlobal(int NumGlobalElements, int NumMyElements) const;
  int *LID_;
  int NumGlobalElements_;
  int NumMyElements_;
  int* MyGlobalElements_;
  int* FirstPointInElementList_;
  int ElementSize_;
  int* ElementSizeList_;
  int* PointToElementList_;
  int IndexBase_;
  const Epetra_Comm * Comm_;
  Epetra_Directory * Directory_;
  
  int NumGlobalPoints_;
  int NumMyPoints_;
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

  int LastContiguousGIDLoc_;
  Epetra_HashTable * LIDHash_;
  mutable int NumReferences_;

};

#endif /* _EPETRA_BLOCKMAP_H_ */
