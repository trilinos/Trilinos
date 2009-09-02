/*
  An abstract interface for a list of elements. 
 Hasn't been debated on
  Jaideep Ray, SNL, Livermore, 08/20/02
*/

#ifndef EntityListHSeen
#define  EntityListHSeen

namespace LoadPartitionerSpace
{
  class EntityList
  {
    public :

    EntityList() {}

    virtual ~EntityList() {}

    /// Returns how long the list is (N)
    virtual unsigned int GetListLength() = 0 ;
    
    /// Methods to access mesh entities
    /** Mesh entities (points/elements etc) are each identified by a 
	Q-element long 1D array of integers. Thus to access the ith
	element in an array a returned by methods in this interface
	one considers all  from a[ i * Q ] to a[ (i+1)*Q -1], both
	inclusive, 0 <= i < N. Elements have a Global ID, unique across all processors
	and a Local ID, unique to the processor where the entity currently
	is. Each element's Local ID can be a 1D array of unsigned longs, R-elements
	in size, Q != R. 
    */
    //@{
    /**
       Returns how long the array containing GIDs is. This number is 
       equal to length of the list (as returned by GetListLength() above and
       the number of unsigned ints needed to represent an GID. 
    */
    virtual int GetGIDListLength() = 0 ;

    /// Returns a pointer to memory containing a list of GIDs
    virtual unsigned int *GetAllGIDs() = 0 ;

    /**  
	 Returns how long the array containing LIDs is. This number is 
	 equal to length of the list (as returned by GetListLength() above and
	 the number of unsigned ints needed to represent an LID. 
    */
    virtual int GetLIDListLength() = 0 ;
    
    /// Returns a pointer to memory containing a list of LIDs
    virtual unsigned int *GetAllLIDs() = 0 ;
    
    /// Returns which processors the N elements live on - integer array N-elements long
    virtual int *GetResidentProcsList() = 0 ;

    /// Return the length of this array
    virtual int GetExtraDataSize() = 0 ;

    /// Return an array of unsigned ints  containing any other data that might be 
    /// needed
    virtual unsigned int *GetExtraData() = 0 ;

    //@}

    /// Reference counting 
    /** Do a addRef() if you are passing a pointer to this object outside your control
	Do a deleteRef() if you are done with this and don't need it any more.
    */
    virtual void addRef() = 0 ;
    virtual void deleteRef() = 0 ;

  } ; // End of EntityList
} ; // End of namespace
#endif
