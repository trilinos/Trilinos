#ifndef __ESI_PreconditionerTranspose_h
#define __ESI_PreconditionerTranspose_h

namespace esi {

/** The ESI PreconditionerTranspose class.

    \verbatim
    The esi::PreconditionerTranspose class is the basis for deriving ESI
    transpose preconditioner implementations.

    The functionality provided is: solveLeftT, solveRightT, solveT, and
    applyBT -- i.e., the tranposes of the operations in the 
    esi::Preconditioner class.
    
    Notes:
    1. It inherits from esi::OperatorTranspose for interoperability.
       For this class, the 'applyTranpose' method is mapped to 'solveT'.
    \endverbatim
*/
template<class Scalar, class Ordinal>
class PreconditionerTranspose : public virtual Preconditioner<Scalar, Ordinal>,
                                public virtual OperatorTranspose<Scalar, Ordinal>
{
 public:

  /** Default destructor. */
  virtual ~PreconditionerTranspose( void ) {};

  //@{ \name Transpose Preconditioning Methods
  /** z = M1^(-T) y */
  virtual ErrorCode solveLeftT( Vector<Scalar, Ordinal> & y, 
                                Vector<Scalar, Ordinal> & z ) = 0;

  /** z = M2^(-T) y */
  virtual ErrorCode solveRightT( Vector<Scalar, Ordinal> & y, 
                                 Vector<Scalar, Ordinal> & z ) = 0;

  /** z = M^(-T) y */
  virtual ErrorCode solveT( Vector<Scalar, Ordinal> & y, 
                            Vector<Scalar, Ordinal> & z ) = 0;

  /** z = B^(T) y */
  virtual ErrorCode applyBT( Vector<Scalar, Ordinal> & y, 
                             Vector<Scalar, Ordinal> & z ) = 0;
  //@}

};     // esi::PreconditionTranspose class
};     // esi namespace
#endif // __ESI_PreconditionerTranspose_h
