#ifndef __ESI_OperatorTranspose_h
#define __ESI_OperatorTranspose_h

namespace esi {

/** The ESI OperatorTranspose class.

    The esi::OperatorTranspose class is designed to support interchangeability 
    of matrices, preconditioners, and solvers when viewed conceptually as just 
    linear operators.
*/
template<class Scalar, class Ordinal>
class OperatorTranspose : public virtual Operator<Scalar, Ordinal>
{
 public:

  /** Default destructor. */  
  virtual ~OperatorTranspose( void ) {};

  /** Function for applying the transpose of this esi::Operator to an
      esi::Vector (<em>x</em>, and producing the result in another 
      esi::Vector (<em>y</em>.  If the operator is <em>this</em>, then the
      apply operation is y = <em>this^T</em> * x (e.g., y = A^Tx 
      matrix-transpose/vector multiply).

      \param x   INPUT: esi::Vector.
      \param y   OUTPUT: esi::Vector.
  */
  virtual ErrorCode applyTranspose( Vector<Scalar, Ordinal> & x, 
                                    Vector<Scalar, Ordinal> & y ) = 0;

};     // esi::OperatorTranspose class
};     // esi namespace
#endif //__ESI_OperatorTranspose_h
