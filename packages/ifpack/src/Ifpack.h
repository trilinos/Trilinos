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
- an interger (defaulted to 0), that specifies the amount of
  overlap among the processes.

The first argument can assume the following values:
- \c "point relaxation" : returns an instance of Ifpack_AdditiveSchwarz<Ifpack_PointRelaxation>
- \c "point relaxation stand-alone" : returns an instance of Ifpack_PointRelaxation (value of overlap is ignored).
- \c "block relaxation" : returns an instance of Ifpack_AdditiveSchwarz<Ifpack_BlockRelaxation>
- \c "block relaxation stand-alone)" : returns an instance of Ifpack_BlockRelaxation.
- \c "Amesos" : returns an instance of Ifpack_AdditiveSchwarz<Ifpack_Amesos>.
- \c "Amesos" : returns an instance of Ifpack_Amesos.
- \c "IC" : returns an instance of Ifpack_AdditiveSchwarz<Ifpack_IC>.
- \c "IC stand-alone" : returns an instance of Ifpack_AdditiveSchwarz<Ifpack_IC>.
- \c "ICT" : returns an instance of Ifpack_AdditiveSchwarz<Ifpack_ICT>.
- \c "ICT stand-alone" : returns an instance of Ifpack_ICT.
- \c "ILU" : returns an instance of Ifpack_AdditiveSchwarz<Ifpack_ILU>.
- \c "ILU stand-alone" : returns an instance of Ifpack_ILU.
- \c "ILUT" : returns an instance of Ifpack_AdditiveSchwarz<Ifpack_ILUT>.
- \c "ILUT stand-alone" : returns an instance of Ifpack_ILUT.
- otherwise, Create() returns 0.

\note Objects in stand-alone mode cannot use reordering, variable overlap, and singleton filters.
However, their construction can be slightly faster than the non stand-alone counterpart. 

<P> The following fragment of code shows the
basic usage of this class.
\code
#include "Ifpack.h"

...

Ifpack Factory;

Epetra_RowMatrix* A; // A is FillComplete()'d.
string PrecType = "ILU"; // use incomplete LU on each process
int OverlapLevel = 1; // one row of overlap among the processes
Ifpack_Preconditioner* Prec = Factory.Create(PrecType, A, OverlapLevel);
assert (Prec != 0);

Teuchos::ParameterList List;
List.set("fact: level-of-fill", 5); // use ILU(5)

IFPACK_CHK_ERR(Prec->SetParameters(List));
IFPACK_CHK_ERR(Prec->Initialize());
IFPACK_CHK_ERR(Prec->Compute());

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

// print information on stdout
cout << *Prec;

// delete the preconditioner
delete Prec;
\endcode

\author Marzio Sala, SNL 9214

\date Last updated on 25-Jan-05.
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

  //! Sets the options in List from the command line
  int SetParameters(int argc, char* argv[],
                    Teuchos::ParameterList& List, string& PrecType,
                    int& Overlap);

};

#endif
