#ifndef IFPACK_H
#define IFPACK_H

#include "Ifpack_ConfigDefs.h"
#include "Ifpack_Preconditioner.h"

//! Ifpack: a function class to define Ifpack preconditioners.
/*!
Class Ifpack is a function class, that contains just one method:
Create(). Using Create(), users can easily define a variety of 
IFPACK preconditioners. 

Create requires 3 arguments:
- a string, indicating the preconditioner to be built;
- a pointer to an Epetra_RowMatrix, representing the matrix
  to be used to define the preconditioner;
- an interger (defaulted to 0), that specifies the amoung of
  overlap among the processes.

The first argument can assume the following values:
- \c "Jacobi"
- \c "Gauss-Seidel"
- \c "symmetric Gauss-Seidel"
- \c "SOR"
- \c "SSOR"
- \c "block Jacobi": block Jacobi preconditioner,
  using LAPACK to apply the inverse of each  block.
- \c "block Gauss-Seidel": block Gauss-Seidel preconditioner,
  using LAPACK to apply the inverse of each  block.
- \c "block symmetric Gauss-Seidel": block 
  symmetric Gauss-Seidel preconditioner,
  using LAPACK to apply the inverse of each  block.
- \c "block Jacobi (Amesos)": block Jacobi, using
  Amesos to apply the inverse of each block
  (require \c --enable-amesos)
- \c "block Gauss-Seidel (Amesos)": block Gauss-Seidel, using
  Amesos to apply the inverse of each block
  (require \c --enable-amesos)
- \c "block symmetric Gauss-Seidel (Amesos)": 
  block symmetric Gauss-Seidel using
  Amesos to apply the inverse of each block
  (require \c --enable-amesos)
- \c "Amesos"
- \c "ICT": incomplete Cholesky factorization
- \c "RILUK": RILU(K) factorization
- otherwise, Create() returns 0.

<P> The following fragment of code shows the
basic usage of this class.
\code
#include "Ifpack.h"

...

Ifpack Factory;

Epetra_RowMatrix* A; // A is FillComplete()'d.
string PrecType = "ICT"; // use incomplete Cholesky on each process
int OverlapLevel = 1; // one row of overlap among the processes
Ifpack_Preconditioner* Prec = Factory.Create(PrecType, A, OverlapLevel);
assert (Prec != 0);

Teuchos::ParameterList List;
List.set("level-of-fill", 5); // use ICT(5)
List.set("schwarz: use RCM reordering", false); // no reordering

Prec->SetParameters(List);
Prec->Initialize();
Prec->Compute();

// now Prec can be used as AztecOO preconditioner
// like for instance
AztecOO AztecOOSolver(*Problem);

// specify solver
AztecOOSolver.SetAztecOption(AZ_solver,AZ_gmres);
AztecOOSolver.SetAztecOption(AZ_output,32);

// Set Prec as preconditioning operator
AztecOOSolver.SetPrecOperator(Prec);

// Call the solver
AztecOOSolver.Iterate(1550,1e-8);

// delete the preconditioner
delete Prec;
\endcode

\author Marzio Sala, SNL 9214

\date Started Oct-04, last update Oct-04
*/

class Ifpack {

public:
  //! Creates an instance of Ifpack_Preconditioner.
  /*! Creates an Ifpack_Preconditioner.
   * \param PrecType (In) - type of preconditioner to be created. 
   *
   * \param Matrix (In) - Matrix used to define the preconditioner
   *
   * \param overlap (In) - specified overlap, defaulted to 0.
   */
  Ifpack_Preconditioner* Create(const string PrecType,
				Epetra_RowMatrix* Matrix,
				const int overlap = 0);


};

#endif
