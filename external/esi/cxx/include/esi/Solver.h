#ifndef __ESI_Solver_h
#define __ESI_Solver_h

namespace esi {

/** The ESI Solver class.

    \verbatim
    The esi::Solver class is the basis for deriving other ESI solver 
    classes.
    
    It provides a solve(x,b) and a standard way of passing string 
    parameters to a solver via parameters().
    
    Notes:
    1. It inherits from the esi::Operator base class so that it can be 
       used as a member of that class for interoperability.  For this 
       class, the 'apply' method is mapped to 'solve'.
   \endverbatim
*/
template<class Scalar, class Ordinal>
class Solver : public virtual Operator<Scalar, Ordinal>
{
 public:

  virtual ~Solver(){}

  /** Solve function: x = A^(-1) b. 
      \param b The right-hand-side ESI vector.
      \param x The solution ESI vector.
  */
  virtual ErrorCode solve( Vector<Scalar, Ordinal> & b, 
                           Vector<Scalar, Ordinal> & x ) = 0;

  /** Set the input control parameters. */
  virtual ErrorCode parameters( int numParams,
                                char** paramStrings ) = 0;

};     // esi::Solver class
};     // esi namespace
#endif //__ESI_Solver_h
