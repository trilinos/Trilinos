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
#include "Tifpack_PointRelaxation.hpp"
#include "Tifpack_ILUT.hpp"

namespace Tifpack {

/** \brief Return true if the specified precondtioner type supports
 * unsymmetric matrices. */
bool supportsUnsymmetric(const std::string& prec_type);

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
#include "Tifpack_Factory.hpp"

...

Tifpack::Factory Factory;

Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>* A; // A is FillComplete()'d.
string PrecType = "ILU"; // use incomplete LU on each process
int OverlapLevel = 1; // one row of overlap among the processes
Teuchos::RCP<Tifpack::Preconditioner<Scalar,LocalOrdinal,GlobalOrdinal,Node>*> Prec =
     Factory.Create(PrecType, A, OverlapLevel);
assert (Prec != Teuchos::null);

Teuchos::ParameterList List;
List.set("fact: level-of-fill", 5); // use ILU(5)

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

// print information on stdout
std::cout << *Prec;
\endcode

\author Michael Heroux, (formally) SNL org. 1416

\date Last updated on Dec-01-2009.
*/
template<class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node>
class Factory {
public:

  /** \brief Creates an instance of Tifpack_Preconditioner given the string
   * name of the preconditioner type (throws exception if given unrecognized name).
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
  static
  Teuchos::RCP<Tifpack::Preconditioner<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
  Create(const std::string& prec_type,
         const Teuchos::RCP<const Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >& Matrix,
         const int overlap = 0);

  /** \brief Sets the options in List from the command line.
   *
   * Note: If you want full support for all parameters, consider reading in a
   * parameter list from an XML file as supported by the Teuchos helper
   * function <tt>Teuchos::updateParametersFromXmlFile()</tt> or
   * <tt>Teuchos::updateParametersFromXmlStream()</tt>.
   */
//  int SetParameters(int argc, char* argv[],
//                    Teuchos::ParameterList& List, string& PrecType,
//                    int& Overlap);

};

/////////////////////////////////////////////
/////////////////////////////////////////////

template<class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node>
Teuchos::RCP<Tifpack::Preconditioner<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
Factory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Create(const std::string& prec_type,
         const Teuchos::RCP<const Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >& Matrix,
         const int /* overlap */)
{
  Teuchos::RCP<Tifpack::Preconditioner<Scalar,LocalOrdinal,GlobalOrdinal,Node> > prec;

  if (prec_type == "POINT_RELAXATION") {
    prec = Teuchos::rcp(new Tifpack::PointRelaxation<Scalar,LocalOrdinal,GlobalOrdinal,Node>(Matrix));
  }
  else if (prec_type == "ILUT") {
    prec = Teuchos::rcp(new Tifpack::ILUT<Scalar,LocalOrdinal,GlobalOrdinal,Node>(Matrix));
  }
  else {
    std::ostringstream os;
    os << "Tifpack::Factory::Create ERROR, invalid preconditioner type ("
       << prec_type;
    std::string str = os.str();
    throw std::runtime_error(str);
  }
  return prec;
}

} //namespace Tifpack

#endif
