#ifndef IFPACK_JACOBI_H
#define IFPACK_JACOBI_H

#include "Ifpack_ConfigDefs.h"
#include "Ifpack_Preconditioner.h"
#include "Ifpack_PointPreconditioner.h"
class Epetra_MultiVector;
class Epetra_Vector;
class Epetra_Map;
class Epetra_Comm;
class Epetra_Vector;
class Epetra_RowMatrix;

//! Ifpack_Jacobi: a class to define point-Jacobi preconditioners of Epetra_RowMatrix's.

/*! The Ifpack_Jacobi class enables the construction of point-Jacobi
 * preconditioners of an Epetra_RowMatrix. Ifpack_Jacobi is derived from 
 * the Ifpack_Preconditioner class, which is derived from Epetra_Operator.
 * Therefore this object can be used as preconditioner everywhere an
 * ApplyInverse() method is required in the preconditioning step.
 * 
 * Consider a linear system of type
 * \f[
 * A x = b,
 * \f]
 * where \f$A\f$ is a square, real matrix, and \f$x, b\f$ are two real
 * vectors. We begin with the decomposition
 * \f[
 * A = D - E - F
 * \f]
 * where \f$D\f$ is the diagonal of A, \f$-E\f$ is the strict lower part, and
 * \f$-F\f$ is the strict upper part. It is assumed that the diagonal entries
 * of \f$A\f$ are different from zero.
 * 
 * Given an starting solution \f$x_0\f$, an iteration of the Jacobi 
 * method can be written in matrix form as follows:
 * \f[
 * x_{k+1} = \omega D^{-1}(E + F) x_k + D{-1}b,
 * \f]
 * for \f$k < k_{max}\f$, and \f$\omega \f$ a damping parameter. 
 * 
 * Using Ifpack_Jacobi, the user can apply the specified number of sweeps
 * (\f$k_{max}\f$), and the damping parameter. If only one sweep is used, then
 * the class simply applies the inverse of the diagonal of A to the input
 * vector.
 *
 * An example of use of Ifpack_Jacobi is the following.
 \code
 #include "Ifpack_Jacobi.h"

 ...
 Epetra_RowMatrix* A;
 Epetra_MultiVector* X, *Y;
 ...                                  // fill the elements of A, X, Y

 Teuchos::ParameterList List;
 List.set("point: sweeps", 1555);            // maximum number of Jacobi sweeps
 List.set("point: omega", 0.67);             // damping parameter
 List.set("point: print frequency", 5);      // prints out convergence info
                                             // every 5 sweeps

 Ifpack_Jacobi APrec(A);              // create an object

 APrec.Compute(List);                 // sets up the preconditioner

 assert (APrec.IsComputed() == true); // verify that it was created

 APrec.ApplyInverse(X,Y);             // solve the linear system with A
                                      // and X as rhs, put the result in Y			  
 \endcode

 To use Jacobi as an AztecOO preconditioner is as simple done as in the 
 following example:
 \code
 AztecOO AztecOOSolver;               // creates the AztecOO solver

 Ifpack_Jacobi APrec(A);
 ...                                  // sets up the preconditioner

 AztecOOSolver.SetPrec(APrec);        // Now Jacobi is the preconditioner
 \endcode

 * \note This class requires the storage of the diagonal of the input matrix.
 * If more than one sweep of Jacobi is used, two additional vectors are
 * allocated.
 *
 * \note It is supposed that RowMatrixRowMap() == OperatorDomainMap().
 * This has to be fixed in the future.
 *
 * \author Marzio Sala, SNL 9214
 *
 * \date Sep-04
 * 
 */
class Ifpack_Jacobi : public Ifpack_PointPreconditioner {

public:

  //@{ \name Constructors/Destructors
  //! Ifpack_Jacobi constructor with given Epetra_RowMatrix.
  /*! Creates an Ifpack_Jacobi preconditioner. 
   *
   * \param In
   * Matrix - Pointer to matrix. It is supposed that the values on the
   *          diagonal of \c Matrix will not change during the 
   *          application of the preconditioner.
   */
  Ifpack_Jacobi(const Epetra_RowMatrix* Matrix) :
    Ifpack_PointPreconditioner(Matrix),
    FirstTime_(true)
  {
    SetLabel();
  }

  //! Destructor.
  ~Ifpack_Jacobi()
  {}

  //! Applies the Jacobi preconditioner to X, returns the result in Y.
  virtual int ApplyInverse(const Epetra_MultiVector& X, 
			   Epetra_MultiVector& Y) const;

private:

  virtual void SetLabel()
  {
    sprintf(Label_, "Ifpack_Jacobi, sweeps = %d, damping = %e",
	    NumSweeps(), DampingFactor());
  }


  mutable bool FirstTime_;
};

#endif // IFPACK_JACOBI_H
