#ifndef __ESI_VectorReplaceAccess_h
#define __ESI_VectorReplaceAccess_h

namespace esi {

/** The esi::VectorReplaceAccess class extends the base 'real' 
    vector class to provide direct pointer replacement access to 
    the data.  Contiguous memory layout is required.

    The intent is to allow a F77/C client to create a vector 
    object and then to loop over a {setArrayPointer;doSomething(v);}
    This has all sorts of inherent safety problems, but those 
    high performance programmers never make mistakes.

    \verbatim
    Change log:

      10/29/2001 RLC Added comments about the behavior of zero-length 
                 arrays (pointers) which are returned as output args, 
                 per Ben Allan's suggestion. Cleaned up the 
                 documentation for consistency.
    \endverbatim
*/
template<class Scalar, class Ordinal>
class VectorReplaceAccess : public virtual Vector<Scalar, Ordinal>
{
 public:
  
  /** Default destructor.  */
  virtual ~VectorReplaceAccess( void ) {};
  
  /** Set pointer to the vector's data and indicate the array length.
      This function will return an error if it is called on a Vector
      that is locked by any of the getCoef*Lock functions.
      @param array The data to replace the current data in the vector.
      @param length The length of the given array. length must match
      the partition local size obtainable from the IndexSpace
      associated with this Vector.
  
      Even if the supplied data has 0 length, data pointer arguments to 
      this function shall not be 0/NULL and shall be properly aligned 
      for the data type.
      
      \verbatim
      Implementation Notes:
      - Here no memory copy required. The length is supplied strictly 
        as a consistency check; use of this operation must not cause 
        the vector to change size!
      - The length is determined at the time the Vector and IndexSpace
        are created.
      - This CAN BE called twice on the same vector.
      - Arrays given via this function will not be deleted when the 
        vector is destroyed.
      - The array given must not be deleted elsewhere before the vector
        is destroyed or another array is substituted by a subsequent 
        call to this function.
      - Start looking for misuse of this function if you leak memory.
      \endverbatim
  */
  virtual ErrorCode setArrayPointer( Scalar * array, 
                                     Ordinal length ) = 0;

}; // esi::VectorReplaceAccess class
}; // esi namespace
#endif // __ESI_VectorReplaceAccess_h
