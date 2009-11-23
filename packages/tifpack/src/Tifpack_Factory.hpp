/*@HEADER
// ***********************************************************************
//
//       Tifpack: Tempated Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER
*/

#ifndef TIFPACK_FACTORY_HPP
#define TIFPACK_FACTORY_HPP

#include "Tifpack_ConfigDefs.hpp"
#include "Tifpack_Preconditioner.hpp"
#include "Tifpack.hpp"
#include "Tifpack_Preconditioner.hpp"
#include "Tifpack_PointRelaxation.hpp"
#include "Tifpack_BlockRelaxation.hpp"
#include "Tifpack_IC.hpp"
#include "Tifpack_ICT.hpp"
#include "Tifpack_ILU.hpp"
#include "Tifpack_ILUT.hpp"
#include "Tifpack_SPARSKIT.hpp"
#include "Tifpack_AdditiveSchwarz.hpp"
#include "Tifpack_DenseContainer.hpp"
#include "Tifpack_SparseContainer.hpp"
#ifdef HAVE_TIFPACK_AMESOS
#include "Tifpack_Amesos.hpp"
#endif
#include "Tifpack_Chebyshev.hpp"

#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_StringToIntMap.hpp"

namespace {

const Teuchos::StringToIntMap
precTypeNameToIntMap(
  "parameter \"Prec Type\"", Tifpack::numPrecTypes, Tifpack::precTypeNames
  );

} // namespace <anonymous>

namespace Tifpack {

//! Tifpack::Factory a factory class to create Tifpack preconditioners.
/*!
Class Tifpack::Factory is a function class, that contains just one method:
Create(). Using Create(), users can easily define a variety of 
Tifpack preconditioners. 

Create requires 3 arguments:
- a string, indicating the preconditioner to be built;
- a pointer to a Tpetra::RowMatrix, representing the matrix
  to be used to define the preconditioner;
- an integer (defaulted to 0), that specifies the amount of
  overlap among the processes.

The first argument can assume the following values:
- \c "point relaxation" : returns an instance of Tifpack_AdditiveSchwarz<Tifpack_PointRelaxation>
- \c "point relaxation stand-alone" : returns an instance of Tifpack_PointRelaxation (value of overlap is ignored).
- \c "block relaxation" : returns an instance of Tifpack_AdditiveSchwarz<Tifpack_BlockRelaxation>
- \c "block relaxation stand-alone)" : returns an instance of Tifpack_BlockRelaxation.
- \c "Amesos" : returns an instance of Tifpack_AdditiveSchwarz<Tifpack_Amesos>.
- \c "Amesos" : returns an instance of Tifpack_Amesos.
- \c "IC" : returns an instance of Tifpack_AdditiveSchwarz<Tifpack_IC>.
- \c "IC stand-alone" : returns an instance of Tifpack_AdditiveSchwarz<Tifpack_IC>.
- \c "ICT" : returns an instance of Tifpack_AdditiveSchwarz<Tifpack_ICT>.
- \c "ICT stand-alone" : returns an instance of Tifpack_ICT.
- \c "ILU" : returns an instance of Tifpack_AdditiveSchwarz<Tifpack_ILU>.
- \c "ILU stand-alone" : returns an instance of Tifpack_ILU.
- \c "ILUT" : returns an instance of Tifpack_AdditiveSchwarz<Tifpack_ILUT>.
- \c "ILUT stand-alone" : returns an instance of Tifpack_ILUT.
- otherwise, Create() returns 0.

\note Objects in stand-alone mode cannot use reordering, variable overlap, and singleton filters.
However, their construction can be slightly faster than the non stand-alone counterpart. 

<P> The following fragment of code shows the
basic usage of this class.
\code
#include "Tifpack.hpp"

...

Tifpack::Factory Factory;

Tpetra_RowMatrix* A; // A is FillComplete()'d.
string PrecType = "ILU"; // use incomplete LU on each process
int OverlapLevel = 1; // one row of overlap among the processes
Tifpack_Preconditioner* Prec = Factory.Create(PrecType, A, OverlapLevel);
assert (Prec != 0);

Teuchos::ParameterList List;
List.set("fact: level-of-fill", 5); // use ILU(5)

TIFPACK_CHK_ERR(Prec->SetParameters(List));
TIFPACK_CHK_ERR(Prec->Initialize());
TIFPACK_CHK_ERR(Prec->Compute());

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

\author Michael Heroux, (formally) SNL org. 1416

\date Last updated on 25-Jan-05.
*/

class Factory {
public:

  /** \brief Enum for the type of preconditioner. */
  enum EPrecType {
    POINT_RELAXATION
    ,POINT_RELAXATION_STAND_ALONE
    ,BLOCK_RELAXATION
    ,BLOCK_RELAXATION_STAND_ALONE
    ,BLOCK_RELAXATION_STAND_ALONE_ILU
#ifdef HAVE_TIFPACK_AMESOS
    ,BLOCK_RELAXATION_STAND_ALONE_AMESOS
    ,BLOCK_RELAXATION_AMESOS
    ,AMESOS
    ,AMESOS_STAND_ALONE
#endif // HAVE_TIFPACK_AMESOS
    ,IC
    ,IC_STAND_ALONE
    ,ICT
    ,ICT_STAND_ALONE
    ,ILU
    ,ILU_STAND_ALONE
    ,ILUT
    ,ILUT_STAND_ALONE
#ifdef HAVE_TIFPACK_SPARSKIT
    ,SPARSKIT
#endif // HAVE_TIFPACK_SPARSKIT
    ,CHEBYSHEV
  };

  /** \brief . */
  static const int numPrecTypes =
    +5
#ifdef HAVE_TIFPACK_AMESOS
    +4
#endif
    +8
#ifdef HAVE_TIFPACK_SPARSKIT
    +1
#endif
    +1
    ;

  /** \brief List of the preconditioner types as enum values . */
  static const EPrecType precTypeValues[numPrecTypes];

  /** \brief List of preconditioner types as string values. */
  static const char* precTypeNames[numPrecTypes];

  /** \brief List of bools that determines if the precondtioner type supports
   * unsymmetric matrices. */
  static const bool supportsUnsymmetric[numPrecTypes];

  /** \brief Function that gives the string name for preconditioner given its
   * enumerication value. */
  static const char* toString(const EPrecType precType)
      { return precTypeNames[precType]; }

  /** \brief Creates an instance of Tifpack_Preconditioner given the enum value
   * of the preconditioner type (can not fail, no bad input possible).
   *
   * \param PrecType (In) - Enum value of preconditioner type to be created. 
   *
   * \param Matrix (In) - Matrix used to define the preconditioner
   *
   * \param overlap (In) - specified overlap, defaulted to 0.
   */
  static Tifpack_Preconditioner* Create(
    EPrecType PrecType, Tpetra_RowMatrix* Matrix, const int overlap = 0
    );

  /** \brief Creates an instance of Tifpack_Preconditioner given the string
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
  Tifpack_Preconditioner* Create(const string PrecType,
				Tpetra_RowMatrix* Matrix,
				const int overlap = 0);

  /** \brief Sets the options in List from the command line.
   *
   * Note: If you want full support for all parameters, consider reading in a
   * parameter list from an XML file as supported by the Teuchos helper
   * function <tt>Teuchos::updateParametersFromXmlFile()</tt> or
   * <tt>Teuchos::updateParametersFromXmlStream()</tt>.
   */
  int SetParameters(int argc, char* argv[],
                    Teuchos::ParameterList& List, string& PrecType,
                    int& Overlap);

};

/////////////////////////////////////////////
/////////////////////////////////////////////

//==============================================================================
const Tifpack::EPrecType Tifpack::precTypeValues[Tifpack::numPrecTypes] =
{
  POINT_RELAXATION
  ,POINT_RELAXATION_STAND_ALONE
  ,BLOCK_RELAXATION
  ,BLOCK_RELAXATION_STAND_ALONE
  ,BLOCK_RELAXATION_STAND_ALONE_ILU
#ifdef HAVE_TIFPACK_AMESOS
  ,BLOCK_RELAXATION_STAND_ALONE_AMESOS
  ,BLOCK_RELAXATION_AMESOS
  ,AMESOS
  ,AMESOS_STAND_ALONE
#endif // HAVE_TIFPACK_AMESOS
  ,IC
  ,IC_STAND_ALONE
  ,ICT
  ,ICT_STAND_ALONE
  ,ILU
  ,ILU_STAND_ALONE
  ,ILUT
  ,ILUT_STAND_ALONE
#ifdef HAVE_TIFPACK_SPARSKIT
  ,SPARSKIT
#endif // HAVE_TIFPACK_SPARSKIT
  ,CHEBYSHEV
};

//==============================================================================
const char* Tifpack::precTypeNames[Tifpack::numPrecTypes] =
{
  "point relaxation"
  ,"point relaxation stand-alone"
  ,"block relaxation"
  ,"block relaxation stand-alone"
  ,"block relaxation stand-alone (ILU)"
#ifdef HAVE_TIFPACK_AMESOS
  ,"block relaxation stand-alone (Amesos)"
  ,"block relaxation (Amesos)"
  ,"Amesos"
  ,"Amesos stand-alone"
#endif
  ,"IC"
  ,"IC stand-alone"
  ,"ICT"
  ,"ICT stand-alone"
  ,"ILU"
  ,"ILU stand-alone"
  ,"ILUT"
  ,"ILUT stand-alone"
#ifdef HAVE_TIFPACK_SPARSKIT
  ,"SPARSKIT"
#endif
  ,"Chebyshev"
};

//==============================================================================
const bool Tifpack::supportsUnsymmetric[Tifpack::numPrecTypes] =
{
  true // point relaxation
  ,true // point relaxation stand-alone
  ,true // block relaxation
  ,true // block relaxation stand-alone
  ,true // block relaxation stand-alone (ILU)
#ifdef HAVE_TIFPACK_AMESOS
  ,true // block relaxation stand-alone (Amesos)
  ,true // block relaxation (Amesos)
  ,true // Amesos
  ,true // Amesos stand-alone 
#endif
  ,false // IC
  ,false // IC stand-alone
  ,false // ICT
  ,false // ICT stand-alone
  ,true // ILU
  ,true // ILU stand-alone
  ,true // ILUT
  ,true // ILUT stand-alone
#ifdef HAVE_TIFPACK_SPARSKIT
  ,true // SPARSKIT
#endif
  ,false // CHEBYSHEV
};

//==============================================================================
Tifpack_Preconditioner* Tifpack::Create(EPrecType PrecType,
                                      Tpetra_RowMatrix* Matrix,
                                      const int Overlap)
{
  switch(PrecType) {
    case POINT_RELAXATION:
      return(new Tifpack_AdditiveSchwarz<Tifpack_PointRelaxation>(Matrix, Overlap));
    case POINT_RELAXATION_STAND_ALONE:
      return(new Tifpack_PointRelaxation(Matrix));
    case BLOCK_RELAXATION:
      return(new Tifpack_AdditiveSchwarz<
             Tifpack_BlockRelaxation<Tifpack_DenseContainer> >(Matrix,Overlap));
    case BLOCK_RELAXATION_STAND_ALONE:
      return(new Tifpack_BlockRelaxation<Tifpack_DenseContainer>(Matrix));
    case BLOCK_RELAXATION_STAND_ALONE_ILU:
      return(new Tifpack_BlockRelaxation<Tifpack_SparseContainer<Tifpack_ILU> >(Matrix));
#ifdef HAVE_TIFPACK_AMESOS
    case BLOCK_RELAXATION_STAND_ALONE_AMESOS:
      return(new Tifpack_BlockRelaxation<Tifpack_SparseContainer<Tifpack_Amesos> >(Matrix));
    case BLOCK_RELAXATION_AMESOS:
      return(new Tifpack_AdditiveSchwarz<
             Tifpack_BlockRelaxation<Tifpack_SparseContainer<Tifpack_Amesos> > >(Matrix,Overlap));
    case AMESOS:
      return(new Tifpack_AdditiveSchwarz<Tifpack_Amesos>(Matrix,Overlap));
    case AMESOS_STAND_ALONE:
      return(new Tifpack_Amesos(Matrix));
#endif
    case IC:
      return(new Tifpack_AdditiveSchwarz<Tifpack_IC>(Matrix,Overlap));
    case IC_STAND_ALONE:
      return(new Tifpack_IC(Matrix));
    case ICT:
      return(new Tifpack_AdditiveSchwarz<Tifpack_ICT>(Matrix,Overlap));
    case ICT_STAND_ALONE:
      return(new Tifpack_ICT(Matrix));
    case ILU:
      return(new Tifpack_AdditiveSchwarz<Tifpack_ILU>(Matrix,Overlap));
    case ILU_STAND_ALONE:
      return(new Tifpack_ILU(Matrix));
    case ILUT:
      return(new Tifpack_AdditiveSchwarz<Tifpack_ILUT>(Matrix,Overlap));
    case ILUT_STAND_ALONE:
      return(new Tifpack_ILUT(Matrix));
#ifdef HAVE_TIFPACK_SPARSKIT
    case SPARSKIT:
      return(new Tifpack_SPARSKIT(Matrix));
#endif
    case CHEBYSHEV:
      return(new Tifpack_Chebyshev(Matrix));
    default:
      TEST_FOR_EXCEPT(true);
      // The only way to get here is if some code developer does a cast like
      // (EPrecType)(anyNumber).  You can never get here by passing in a
      // string value for the preconditioner!
  } // end switch
  return 0; // This will never ever be executed!
}

//==============================================================================
Tifpack_Preconditioner* Tifpack::Create(const string PrecType,
                                      Tpetra_RowMatrix* Matrix,
                                      const int Overlap)
{
  try {
    return Tifpack::Create(Teuchos::get<EPrecType>(::precTypeNameToIntMap,PrecType),Matrix,Overlap);
  }
  catch( const Teuchos::StringToIntMap::DoesNotExist &excpt ) {
    // The old implementation of this function just silently returned a NULL
    // when a preconditiner type name was not recognized.  If you like this
    // behavior then you should use this function.  If you do not like this
    // behavior, then consider using the Tifpack/Thyra adapter
    // Thyra::TifpackPreconditionerFactory or better yet the Stratimikos
    // wrapper class Stratimikos::DefaultLinearSolverBuilder.
  }
  return 0;
}

// ======================================================================
int Tifpack::SetParameters(int argc, char* argv[],
                          Teuchos::ParameterList& List, string& PrecType,
                          int& Overlap)
{
  // THIS FUNCTION IS VERY INCOMPLETE...

  Teuchos::CommandLineProcessor CLP;

  // prec type
  string ifp_prec_type = "ILU";
  CLP.setOption("ifp-prec-type",&ifp_prec_type,"Preconditioner type");
  // overlap among the processors
  int ifp_overlap = 0;
  CLP.setOption("ifp-overlap",&ifp_overlap,"Overlap among processors");
  // relaxation type
  string ifp_relax_type = "Jacobi";
  CLP.setOption("ifp-relax-type",&ifp_relax_type,"Relaxation type");
  // sweeps (for relax only)
  int ifp_relax_sweeps = 1;
  CLP.setOption("ifp-relax-sweeps",
                &ifp_relax_sweeps,"Number of sweeps for relaxation");
  // damping (for relax only)
  double ifp_relax_damping = 1.0;
  CLP.setOption("ifp-relax-damping",
                &ifp_relax_damping,"Damping for relaxation");
  // partitioner type (for block relaxation only)
  string ifp_part_type = "greedy";
  CLP.setOption("ifp-part-type",&ifp_part_type,"Partitioner type");
  // number of local parts (for block relaxation only)
  int ifp_part_local = 1;
  CLP.setOption("ifp-part-local",&ifp_part_local,"number of local partitions");

  // allow users to specify other options for other packages
  CLP.recogniseAllOptions(false);
  CLP.throwExceptions(false);
  CLP.parse(argc,argv);

  // I cannot really set those in the List, I pass them back to the user
  PrecType = ifp_prec_type;
  Overlap = ifp_overlap;

  // set the list here
  List.set("relaxation: type", ifp_relax_type);
  List.set("relaxation: sweeps", ifp_relax_sweeps);
  List.set("relaxation: damping factor", ifp_relax_damping);
  List.set("partitioner: type", ifp_part_type);
  List.set("partitioner: local parts", ifp_part_local);

  return(0);
}
#endif
