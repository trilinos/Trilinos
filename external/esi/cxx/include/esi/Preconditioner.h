#ifndef __ESI_Preconditioner_h
#define __ESI_Preconditioner_h

namespace esi {

/** ESI preconditioner side terms */
typedef enum { PRECONDITIONER_LEFT,
               PRECONDITIONER_RIGHT,
               PRECONDITIONER_TWO_SIDED } PreconditionerSide;

/** The ESI Preconditioner class.

    \verbatim
    The esi::Preconditioner class is the basis for deriving ESI preconditioner 
    implementations.

    The preconditioner functionality provided is: solveLeft, solveRight,
    solve, and applyB. The meaning of the solve functions is as follows.
    Use M to denote the preconditioner, and M = M1 * M2 is a splitting
    (where applicable), and B = M1^(-1) A M2^(-1). Then solveLeft() solves the
    system M1*z = y for z, solveRight() solves M2*z = y for z, and solve()
    solves M*z = y for z. Note that for many preconditioners, either M1 or
    M2 is the identity matrix.

    Check esi/cxx/include/esi/basicTypes.h for definitions of preconditioner 
    'sides'.  The following enum is defined there:

    typedef enum { PRECONDITIONER_LEFT, 
                   PRECONDITIONER_RIGHT,
                   PRECONDITIONER_TWO_SIDED } PreconditionerSide;

    Developer Notes:

    1. It inherits from esi::Operator so that it can be used as a member 
       of that class for interoperability.  For this class, the 'apply' 
       method is mapped to 'solve'.

    \endverbatim
*/
template<class Scalar, class Ordinal>
class Preconditioner : public virtual Operator<Scalar, Ordinal>
{
 public:

  /** Default destructor. */
  virtual ~Preconditioner( void ) {};

  /** Input control parameters. */
  virtual ErrorCode parameters( int numParams,
                                char** paramStrings ) = 0;

  //@{ \name Preconditioning Methods
  /** z = M1^(-1) y */
  virtual ErrorCode solveLeft( Vector<Scalar, Ordinal> & y, 
                               Vector<Scalar, Ordinal> & z ) = 0;
  /** z = M2^(-1) y */
  virtual ErrorCode solveRight( Vector<Scalar, Ordinal> & y, 
                                Vector<Scalar, Ordinal> & z ) = 0;

  /** z = M^(-1) y */
  virtual ErrorCode solve( Vector<Scalar, Ordinal> & y, 
                           Vector<Scalar, Ordinal> & z ) = 0;
  
  /** z = B y */
  virtual ErrorCode applyB( Vector<Scalar, Ordinal> & y, 
                            Vector<Scalar, Ordinal> & z ) = 0;
  //@}

  /** Get the preconditioning side. */
  virtual ErrorCode getPreconditionerSide( PreconditionerSide & side ) = 0;

  /** Set the preconditioning side. */
  virtual ErrorCode setPreconditionerSide( PreconditionerSide side ) = 0;

  /** Set the operator A (the one being solved for, i.e., the one that
      appears in the "B" operator described at the top of this file). 
  */
  virtual ErrorCode setOperator(Operator<Scalar, Ordinal> & A) = 0;

};     // esi::Preconditioner class
};     // esi namespace
#endif // __ESI_Preconditioner_h
