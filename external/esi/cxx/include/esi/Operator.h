#ifndef __ESI_Operator_h
#define __ESI_Operator_h

namespace esi {

/** The ESI Operator class.

    The esi::Operator class is designed to support interchangeability of 
    matrices, preconditioners, and solvers when viewed conceptually as just 
    linear operators.
*/
template<class Scalar, class Ordinal>
class Operator : public virtual Object 
{
 public:

  /** Default destructor. */  
  virtual ~Operator( void ) {};

  /** Function for performing initial calculations (e.g., factorization). */
  virtual ErrorCode setup( void ) = 0;

  /** Function for applying this operator to a Vector, and producing the 
      result in another Vector.  If the operator is <em>this</em>, then the
      apply operation is y = <em>this</em> * x (e.g., y = Ax matrix/vector
      multiply).

      \param x   INPUT: esi::Vector.
      \param y   OUTPUT: esi::Vector.
  */
  virtual ErrorCode apply( Vector<Scalar, Ordinal> & x,
                           Vector<Scalar, Ordinal> & y ) = 0;

};     // esi::Operator class
};     // esi namespace
#endif //__ESI_Operator_h
