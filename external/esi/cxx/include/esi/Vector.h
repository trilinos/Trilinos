#ifndef __ESI_Vector_h
#define __ESI_Vector_h

namespace esi {

/** The esi::Vector class is the base class for 'real' ESI vector objects.

    Unless otherwise noted, the notation here is such that you 
    are 'in' the object defined by <em>this</em>, as defined in the 
    C++ standard (i.e., in C++, <em>this</em> = 'pointer to calling object').

    \verbatim
    Developer notes: 
    This class will probably be morphed/extended into one that 
    supports 'complex' types, as some of our core apps need that
    capability.

    \verbatim
    Change log:

      10/29/2001 RLC Added comments about the behavior of zero-length arrays
                 (pointers) which are returned as output args, per Ben Allan's
                 suggestion. Cleaned up the documentation for consistency.

      10/25/2001 RLC Changed the documentation to use syntax suggested
                 by Barry Smith (email 10/23/2001) whereby we ditched 
                 the 'this' = &y notation, in favor of simply using 
                 'this' whereever appropriate.  The new way seems more 
                 straight-forward and explicit.  I also edited and
                 extended the content of the documentation, trying to 
                 get things more explicit and consistent.  Removed some
                 buggy doxygen annotations which were corrupting the 
                 doxygen parsing and hence output.  Added  'return 0 
                 if successful' comments for all relevant methods.

      8/31/2001  RLC Added getLocalSize method. Tuned up the 
                 documentation toward more consistent overall format.
    \endverbatim
*/
template<class Scalar, class Ordinal>
class Vector : public virtual Object
{
 public:
  
  /** Default destructor. */
  virtual ~Vector( void ) {};
  
  /** Magnitude type definition */
  typedef TYPENAME scalarTraits<Scalar>::magnitude_type magnitude_type;
  
  //@{ \name Easy functions
  
  /** Clone <em>this</em> vector to create new vector <em>x</em>.
      
      This method is expected to duplicate the current esi::Vector 
      object by creating a new esi::Vector object which is returned 
      in the pointer <em>x</em>. Note that this function does 
      <em<not</em> copy the vector contents into the new
      vector. It simply clones the size/structure/layout, creating a 
      new vector that is compatible with <em>this</em> vector.
      Deletion of the memory for this new object is left to the 'owner' 
      of <em>x</em>.

      @return ErrorCode = 0 if successful.
  */
  virtual ErrorCode clone( Vector * & x ) = 0;

  /** Get the Hilbert-space global size.

      This is the mathematical dimension of the vector, regardless of 
      data-processor distribution.

      @return ErrorCode = 0 if successful.
  */
  virtual ErrorCode getGlobalSize( Ordinal & dim ) = 0;

  /** Query the <em>local</em> vector-element array size, where 
      by 'local' we mean processor locality in an MPP 
      (distributed-memory) data distribution environment.

      When there are no overlaps, the sum of the <em>localSize</em> 
      values on all processors should equal the <em>globalSize</em> 
      value.  Unless otherwise stated, the default is to assume no 
      overlaps.

      @param localSize   OUTPUT: The number of local vector elements.
      @return ErrorCode = 0 if successful.
  */
  virtual ErrorCode getLocalSize( Ordinal & localSize ) = 0;

  /** Get (and cache) a pointer to an esi::IndexSpace object. 

      @return ErrorCode = 0 if successful.
  */
  virtual ErrorCode getIndexSpace( IndexSpace<Ordinal> * & indexSpace ) = 0;

  /** <em>this</em>[i] = <em>scalar</em> for all i. 

      @return ErrorCode = 0 if successful.
  */
  virtual ErrorCode put( Scalar scalar ) = 0;

  /** scaling <em>this</em> <- <em>this</em> * <em>scalar</em>. 

      @return ErrorCode = 0 if successful.
  */
  virtual ErrorCode scale( Scalar scalar ) = 0;

  /** sum(abs(<em>this</em>[i])). 

      @return ErrorCode = 0 if successful.
  */
  virtual ErrorCode norm1( magnitude_type & norm ) = 0;

  /** sqrt(sum(sqr(<em>this</em>[i]))). 

      @return ErrorCode = 0 if successful.
  */
  virtual ErrorCode norm2( magnitude_type & norm ) = 0;

  /** sum(sqr(<em>this</em>[i])). 

      @return ErrorCode = 0 if successful.
  */
  virtual ErrorCode norm2squared( magnitude_type & norm ) = 0;

  /** max(fabs(<em>this</em>[i])). 

      @return ErrorCode = 0 if successful.
  */
  virtual ErrorCode normInfinity( magnitude_type & norm ) = 0;

  /** scalar <- min(abs(<em>this</em>[i])).
      
      Non-norm but useful in some scaling algorithms.

      @return ErrorCode = 0 if successful.
  */
  virtual ErrorCode minAbsCoef( magnitude_type & scalar ) = 0;
  //@} // end easy
    
  //@{ \name HARD

  /** scale diagonal: <em>this</em>[i] <- <em>this</em>[i] * <em>x</em>[i]. 

      \verbatim
      Developer notes:

      When the matrix/operator/preconditioner/solver folks get to 
      it, this function may move to an esi::OperatorDiagonal 
      class.
      \endverbatim 

      @return ErrorCode = 0 if successful.
  */
  virtual ErrorCode scaleDiagonal( Vector & x ) = 0;
  
  /** copy: <em>this</em> <-- <em>x</em>. 

      @return ErrorCode = 0 if successful.
  */
  virtual ErrorCode copy( Vector & x ) = 0;

  /** dot product: <em>product</em> <-- <em>this</em>.<em>x</em>. 

      @return ErrorCode = 0 if successful.
  */
  virtual ErrorCode dot( Vector & x, 
                         Scalar & product ) = 0;
  
  /** axpy: <em>this</em> <-- <em>this</em> + <em>scalar</em> * <em>x</em>. 

      @return ErrorCode = 0 if successful.
  */
  virtual ErrorCode axpy( Vector & x, 
                          Scalar scalar ) = 0;
  
  /** aypx: <em>this</em> <-- <em>this</em> * <em>scalar</em> +<em>x</em>. 

      @return ErrorCode = 0 if successful.
  */
  virtual ErrorCode aypx( Scalar scalar, 
                          Vector & x ) = 0;

  /** This method computes <em>this</em> = <em>a*x</em> + <em>b*y</em>, where 
      esi::Vector objects <em>x</em> and <em>y</em> are different and conforming
      to <em>this</em> vector.  
        
      In a quality implementation several special cases are handled separately:
      Either of (a,b) in (-1,0,1), <em>x=y, x=this, y=this, x=y=this</em>.
      <em>x, y</em> = NULL (NULL objects) are not expected special cases.
      The other functions in this class that handle 2 vectors overwrite
      results on one of the numerically contributing arguments.
      It is desired to have a function that does not overwrite
      the inputs. <em>this</em> = <em>a*x</em> + <em>b*y</em> is such a function.

      @return ErrorCode = 0 if successful.
  */
  virtual ErrorCode axpby( Scalar a, 
                           Vector & x, 
                           Scalar b, 
                           Vector & y ) = 0 ;
  //@} // end hard

  //@{ \name Common direct memory access.
  
  /**  Direct pointer access to data, whether contiguous storage 
       is 'native' form or not.  The pointer returned is 
       laid out according the result of getIndexSpace().
  */

  /** Get read-only access to an array pointer.

      Depending on the implementation this may be the native data.  
      Multiple simultaneous reads are allowed.
      No one else can use the getCoefPtrReadWriteLock function of 
      the Vector until releaseCoefPtrLock is called for each Read 
      outstanding.
      Modification of the data array returned from Read
      will have unpredictable results and should not be done.
      The address coefPtr (not just *coefPtr) may be stored by the
      implementation and must be the same at the time
      of the corresponding releaseCoefPtrLock call.

      As with ANSI malloc, each data pointer retrieved with this function
      shall not be 0/NULL and shall be properly aligned for the data type
      even if the data has 0 length.

      @return ErrorCode = 0 if successful.
  */
  virtual ErrorCode getCoefPtrReadLock( Scalar * & coefPtr ) = 0;
  
  /** Get read/write access to an array pointer.
      
      Depending on the implementation this may be the native data.  
      Only one ReadWrite is allowed at a time, and getting it will
      fail if there is a Read outstanding.
      No one else can use getCoefPtrReadLock function of the Vector
      until releaseCoefPtrLock is called. 
      The Vector values are updated from the array during 
      releaseCoefPtrLock or sooner, depending on the implementation.
      The address coefPtr (not just *coefPtr) may be stored by the
      implementation and must be the same at the time of the 
      corresponding releaseCoefPtrLock call.

      As with ANSI malloc, each data pointer retrieved with this function
      shall not be 0/NULL and shall be properly aligned for the data type
      even if the data has 0 length.

      @return ErrorCode = 0 if successful.
  */
  virtual ErrorCode getCoefPtrReadWriteLock( Scalar * & coefPtr ) = 0;
  
  /** Release the Vector from a previous getCoefPtr*Lock call,
      ensuring that write-back to the vector has occured if
      the accessArray call was Write or ReadWrite. The pointer 
      given must be the one given to the previous getCoefPtrLock 
      call and *coefPtr will be reset to NULL/0 during the 
      releaseCoefPtrLock call.

      @return ErrorCode = 0 if successful.
  */
  virtual ErrorCode releaseCoefPtrLock( Scalar * & coefPtr ) = 0;

  //@} // end common direct memory access

}; // esi::Vector class
}; // esi namespace
#endif //__ESI_Vector_h
