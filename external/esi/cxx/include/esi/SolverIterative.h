#ifndef __ESI_SolverIterative_h
#define __ESI_SolverIterative_h

namespace esi {

/** The ESI SolverIterative class.

    The esi::SolverIterative class is the basis for deriving ESI
    iterative solver implementations.
*/
template<class Scalar, class Ordinal>
class SolverIterative : public virtual Solver<Scalar, Ordinal>
{
 public:

  /** Default destructor. */
  virtual ~SolverIterative( void ) {};

  typedef TYPENAME scalarTraits<Scalar>::magnitude_type magnitude_type;

  /** Get the operator. */
  virtual ErrorCode getOperator( Operator<Scalar, Ordinal> * & A ) = 0;

  /** Set the operator. */
  virtual ErrorCode setOperator( Operator<Scalar, Ordinal> & A ) = 0;

  /** Get the preconditioner. */
  virtual ErrorCode getPreconditioner( Preconditioner<Scalar, Ordinal> * & pc ) = 0;

  /** Set the preconditioner. */
  virtual ErrorCode setPreconditioner( Preconditioner<Scalar, Ordinal> & pc ) = 0;

  /** Get the convergence tolerance. */
  virtual ErrorCode getTolerance( magnitude_type & tol ) = 0;

  /** Set the convergence tolerance. */
  virtual ErrorCode setTolerance( magnitude_type tol ) = 0;

  /** Get the maximum number of iterations. */
  virtual ErrorCode getMaxIterations( Ordinal & maxIterations ) = 0;

  /** Set the maximum number of iterations. */
  virtual ErrorCode setMaxIterations( Ordinal maxIterations ) = 0;

  /** Query the number of iterations that were taken during the previous solve.
  */
  virtual ErrorCode getNumIterationsTaken(Ordinal& itersTaken) = 0;

};     // esi::SolverIterative class
};     // esi namespace
#endif //__ESI_SolverIterative_h
