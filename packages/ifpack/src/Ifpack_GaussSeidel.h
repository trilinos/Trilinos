#ifndef IFPACK_GAUSSSEIDEL_H
#define IFPACK_GAUSSSEIDEL_H

#include "Ifpack_ConfigDefs.h"
#include "Ifpack_Preconditioner.h"
#include "Ifpack_PointPreconditioner.h"
class Epetra_MultiVector;
class Epetra_Vector;
class Epetra_Map;
class Epetra_Comm;
class Epetra_Vector;
class Epetra_RowMatrix;

//! Ifpack_GaussSeidel: a class to define point Gauss-Seidel preconditioners of Epetra_RowMatrix's.

/*! The Ifpack_GaussSeidel class enables the construction of point Gauss-Seidel
 * preconditioners of an Epetra_RowMatrix. Ifpack_GaussSeidel is derived from 
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
 * Given an starting solution \f$x_0\f$, an iteration of the GaussSeidel 
 * method can be written in matrix form as follows:
 * \f[
 * (D - E) x_{k+1} = \omega F x_k + b,
 * \f]
 * for \f$k < k_{max}\f$, and \f$\omega \f$ a damping parameter. Equivalently,
 * the Gauss-Seidel preconditioner can be defined as
 * \f[
 * P_{GS}^{-1} = (D - E)^{-1}.
 * \f]
 * 
 * Clearly, the role of E and F can be interchanged. However,
 * Ifpack_GaussSeidel does not consider backward Gauss-Seidel methods.
 * 
 * Using Ifpack_GaussSeidel, the user can apply the specified number of sweeps
 * (\f$k_{max}\f$), and the damping parameter. If only one sweep is used, then
 * the class simply applies the inverse of the diagonal of A to the input
 * vector.
 *
 * An example of use of Ifpack_GaussSeidel is the following.
 \code
 #include "Ifpack_GaussSeidel.h"

 ...
 Epetra_RowMatrix* A;
 Epetra_MultiVector* X, *Y;
 ...                                  // fill the elements of A, X, Y

 Teuchos::ParameterList List;
 List.set("point: sweeps", 1555);     // maximum number of GaussSeidel sweeps
 List.set("point: omega", 0.67);      // damping parameter
 List.set("point: print frequency", true); // prints out convergence info
                                           // every 5 sweeps

 Ifpack_GaussSeidel APrec(A);         // create an object

 APrec.Compute(List);                 // set up the preconditioner

 assert (APrec.IsComputed() == true); // verify that it was created

 APrec.ApplyInverse(X,Y);             // solve the linear system with A
                                      // and X as rhs, put the result in Y			  
 \endcode

 * \note 
 * If more than one sweep of GaussSeidel is used, two additional vectors are
 * allocated.
 *
 * \note It is supposed that RowMatrixRowMap() == OperatorDomainMap().
 * This has to be fixed in the future.
 *
 * \date Sep-04
 * 
 */
class Ifpack_GaussSeidel : public Ifpack_PointPreconditioner {

public:

  //@{ \name Constructor/Destructor
  //! Constructs a Gauss-Seidel preconditioner object for the input Epetra_RowMatrix.
  Ifpack_GaussSeidel(const Epetra_RowMatrix* Matrix) :
    Ifpack_PointPreconditioner(Matrix),
    FirstTime_(true)
  {}

  //! Destructor.  
  ~Ifpack_GaussSeidel()
  {}

  virtual int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

  virtual void SetLabel();

private:

  mutable bool FirstTime_;
};

#endif // IFPACK_GAUSSSEIDEL_H
