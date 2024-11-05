/*@HEADER
// ***********************************************************************
//
//       Ifpack: Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2002) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER
*/

#ifndef IFPACK_H
#define IFPACK_H

#if defined(Ifpack_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Ifpack package is deprecated"
#endif
#endif


#include "Ifpack_ConfigDefs.h"
#include "Ifpack_Preconditioner.h"
#include "Teuchos_iostream_helpers.hpp"


#ifdef HAVE_HYPRE
#include "Ifpack_Hypre.h"
#endif


//! Ifpack: a function class to define Ifpack preconditioners.
/*!
Class Ifpack is a function class, that contains just one method:
Create(). Using Create(), users can easily define a variety of
IFPACK preconditioners.

Create requires 3 arguments:
- a std::string, indicating the preconditioner to be built;
- a pointer to an Epetra_RowMatrix, representing the matrix
  to be used to define the preconditioner;
- an interger (defaulted to 0), that specifies the amount of
  overlap among the processes.

The first argument can assume the following values:
- \c "point relaxation" : returns an instance of Ifpack_AdditiveSchwarz<Ifpack_PointRelaxation> (no Additive Schwarz in serial)
- \c "point relaxation stand-alone" : returns an instance of Ifpack_PointRelaxation (value of overlap is ignored).
- \c "block relaxation" : returns an instance of Ifpack_AdditiveSchwarz<Ifpack_BlockRelaxation> (no Additive Schwarz in serial)
- \c "block relaxation stand-alone)" : returns an instance of Ifpack_BlockRelaxation.
- \c "Amesos" : returns an instance of Ifpack_AdditiveSchwarz<Ifpack_Amesos> (no Additive Schwarz in serial)
- \c "Amesos stand-alone" : returns an instance of Ifpack_Amesos.
- \c "IC" : returns an instance of Ifpack_AdditiveSchwarz<Ifpack_IC> (no Additive Schwarz in serial)
- \c "IC stand-alone" : returns an instance of Ifpack_IC.
- \c "ICT" : returns an instance of Ifpack_AdditiveSchwarz<Ifpack_ICT> (no Additive Schwarz in serial)
- \c "ICT stand-alone" : returns an instance of Ifpack_ICT.
- \c "ILU" : returns an instance of Ifpack_AdditiveSchwarz<Ifpack_ILU> (no Additive Schwarz in serial)
- \c "ILU stand-alone" : returns an instance of Ifpack_ILU.
- \c "ILUT" : returns an instance of Ifpack_AdditiveSchwarz<Ifpack_ILUT> (no Additive Schwarz in serial)
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
std::string PrecType = "ILU"; // use incomplete LU on each process
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

\author Marzio Sala, (formally) SNL org. 1414

\date Last updated on 25-Jan-05.
*/

class Ifpack {
public:

  /** \brief Enum for the type of preconditioner. */
  enum EPrecType {
    POINT_RELAXATION
    ,POINT_RELAXATION_STAND_ALONE
    ,BLOCK_RELAXATION
    ,BLOCK_RELAXATION_STAND_ALONE
    ,BLOCK_RELAXATION_STAND_ALONE_ILU
    ,BLOCK_RELAXATION_STAND_ALONE_ILUT
    ,BLOCK_RELAXATION_STAND_ALONE_IC
#ifdef HAVE_IFPACK_SUPERLU
    ,BLOCK_RELAXATION_STAND_ALONE_SILU
#endif
#ifdef HAVE_IFPACK_AMESOS
    ,BLOCK_RELAXATION_STAND_ALONE_AMESOS
    ,BLOCK_RELAXATION_AMESOS
    ,AMESOS
    ,AMESOS_STAND_ALONE
#endif // HAVE_IFPACK_AMESOS
    ,IC
    ,IC_STAND_ALONE
    ,ICT
    ,ICT_STAND_ALONE
    ,ILU
    ,ILU_STAND_ALONE
    ,ILUT
    ,ILUT_STAND_ALONE
#ifdef HAVE_IFPACK_SPARSKIT
    ,SPARSKIT
#endif // HAVE_IFPACK_SPARSKIT
#ifdef HAVE_IFPACK_HIPS
  ,HIPS
#endif
#ifdef HAVE_HYPRE
  ,HYPRE
#endif
#ifdef HAVE_IFPACK_SUPERLU
  ,SILU
#endif
#if defined (HAVE_IFPACK_SUPPORTGRAPH) && defined (HAVE_IFPACK_AMESOS)
    ,MSF_AMESOS
#endif
#ifdef HAVE_IFPACK_SUPPORTGRAPH
    ,MSF_IC
#endif
    ,CHEBYSHEV
    ,POLYNOMIAL
    ,KRYLOV
    ,IHSS
    ,SORA
    ,TRIDI_RELAXATION
    ,TRIDI_RELAXATION_STAND_ALONE
  };

  /** \brief . */
  static const int numPrecTypes =
    +7
#ifdef HAVE_IFPACK_AMESOS
    +4
#endif
    +8
#ifdef HAVE_IFPACK_SPARSKIT
    +1
#endif
#ifdef HAVE_IFPACK_HIPS
    +1
#endif
#ifdef HAVE_HYPRE
    +1
#endif
#ifdef HAVE_IFPACK_SUPERLU
    +2
#endif
#if defined (HAVE_IFPACK_SUPPORTGRAPH) && defined (HAVE_IFPACK_AMESOS)
    +1
#endif
#ifdef HAVE_IFPACK_SUPPORTGRAPH
    +1
#endif
    +7
    ;

  /** \brief List of the preconditioner types as enum values . */
  static const EPrecType precTypeValues[numPrecTypes];

  /** \brief List of preconditioner types as std::string values. */
  static const char* precTypeNames[numPrecTypes];

  /** \brief List of bools that determines if the preconditioner type supports
   * unsymmetric matrices. */
  static const bool supportsUnsymmetric[numPrecTypes];

  /** \brief Function that gives the std::string name for preconditioner given its
   * enumerication value. */
  static const char* toString(const EPrecType precType)
      { return precTypeNames[precType]; }

  /** \brief Creates an instance of Ifpack_Preconditioner given the enum value
   * of the preconditioner type (can not fail, no bad input possible).
   *
   * \param PrecType (In) - Enum value of preconditioner type to be created.
   *
   * \param Matrix (In) - Matrix used to define the preconditioner
   *
   * \param overlap (In) - specified overlap, defaulted to 0.
   */
  static Ifpack_Preconditioner* Create(
    EPrecType PrecType, Epetra_RowMatrix* Matrix, const int overlap = 0, bool overrideSerialDefault = false
    );

  /** \brief Creates an instance of Ifpack_Preconditioner given the std::string
   * name of the preconditioner type (can fail with bad input).
   *
   * \param PrecType (In) - String name of preconditioner type to be created.
   *
   * \param Matrix (In) - Matrix used to define the preconditioner
   *
   * \param overlap (In) - specified overlap, defaulted to 0.
   *
   * Returns <tt>0</tt> if the preconditioner with that input name does not
   * exist.  Otherwise, return a newly created preconditioner object.  Note
   * that the client is responsible for calling <tt>delete</tt> on the
   * returned object once it is finished using it!
   */
  Ifpack_Preconditioner* Create(const std::string PrecType,
                                Epetra_RowMatrix* Matrix,
                                const int overlap = 0,
                                bool overrideSerialDefault = false);

  /** \brief Sets the options in List from the command line.
   *
   * Note: If you want full support for all parameters, consider reading in a
   * parameter list from an XML file as supported by the Teuchos helper
   * function <tt>Teuchos::updateParametersFromXmlFile()</tt> or
   * <tt>Teuchos::updateParametersFromXmlStream()</tt>.
   */
  int SetParameters(int argc, char* argv[],
                    Teuchos::ParameterList& List, std::string& PrecType,
                    int& Overlap);

};


TEUCHOS_ENUM_INPUT_STREAM_OPERATOR(Ifpack::EPrecType)


#endif
