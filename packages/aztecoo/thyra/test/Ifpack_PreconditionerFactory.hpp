/*@HEADER
// ***********************************************************************
// 
//       Ifpack: Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2002) Sandia Corporation
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

#ifndef IFPACK_PRECONDITIONER_FACTORY_HPP
#define IFPACK_PRECONDITIONER_FACTORY_HPP

#include "Teuchos_StandardMemberCompositionMacros.hpp"
#include "Teuchos_RefCountPtr.hpp"

class Epetra_Operator;

namespace Ifpack {

/** \brief Utility class that automates the creation of preconditions from
 * Epetra_Operator objects..
 *
 * ToDo: Finish documentation!
 *
 * Note, the default copy constructor and assignment operators functions are
 * allowed since they have the correct semantics.  Copies a
 * <tt>PreconditionerFactory</tt> effectively just copy the options (see the
 * constructor <tt>PreconditionerFactory()</tt>).
 *
 * Note, this class maintains no specific state data except for the options
 * (see the constructor <tt>PreconditionerFactory()</tt>).  Therefore, the
 * same <tt>PreconditionerFactory</tt> object can be used to generate
 * preconditioner objects for as many different types and kinds of operators
 * as the client would like.
 *
 * Also, once a preconditioner is created, it is a self contained entity that
 * exists indepenent of the operator it was created form or from
 * <tt>this</tt>.  This is the best that you could ever hope.
 */
class PreconditionerFactory {
public:

  /** \brief Sets all of the adjustable options to default values.
   *
   * Note, these where the default values used in NOX as of 2004/1/16
   *
   * ToDo: Finish documentation!
   */
  PreconditionerFactory(
    const int        levelFill     = 1
    ,const int       levelOverlap  = 0
    ,const double    absThreshold  = 0.0
    ,const double    relThreshold  = 1.0
		,const bool      calcCondEst   = false
    );

  /** \brief . */
	STANDARD_MEMBER_COMPOSITION_MEMBERS( int, levelFill )
  /** \brief . */
	STANDARD_MEMBER_COMPOSITION_MEMBERS( int, levelOverlap )
  /** \brief . */
	STANDARD_MEMBER_COMPOSITION_MEMBERS( double, absThreshold )
  /** \brief . */
	STANDARD_MEMBER_COMPOSITION_MEMBERS( double, relThreshold )
  /** \brief . */
	STANDARD_MEMBER_COMPOSITION_MEMBERS( bool, calcCondEst )

  /** \brief Setup (and create if not done so yet) a preconditioner.
   *
   * ToDo: Finish documentation!
   */
  void setupPrec(
    const Epetra_Operator                         &Op
    ,Teuchos::RefCountPtr<Epetra_Operator>        *Prec
		,std::ostream                                 *out            = NULL
    ,const std::string                            &leadingIndent  = ""
    ,const std::string                            &indentSpacer   = "  "
    ) const;

}; // class PreconditionerFactory

} // namespace Ifpack

#endif // IFPACK_PRECONDITIONER_FACTORY_HPP
