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

#ifndef THYRA_IFPACK_PRECONDITIONER_FACTORY_DECL_HPP
#define THYRA_IFPACK_PRECONDITIONER_FACTORY_DECL_HPP

#include "Thyra_Ifpack_Types.hpp"
#include "Thyra_PreconditionerFactoryBase.hpp"
#include "Thyra_EpetraOperatorViewExtractorBase.hpp"
#include "Teuchos_StandardCompositionMacros.hpp"

namespace Thyra {

/** \brief Concrete preconditioner factory subclass based on Ifpack.
 *
 * ToDo: Finish documentation!
 */
class IfpackPreconditionerFactory : public PreconditionerFactoryBase<double> {
public:

  /** @name Constructors/initializers/accessors */
  //@{

  /** \brief . */
  IfpackPreconditionerFactory();
    
  /** \brief Set the strategy object used to extract an
   * <tt>Epetra_Operator</tt> view of an input forward operator.
   *
   * This view will then be dynamically casted to <tt>Epetra_RowMatrix</tt>
   * before it is used.
   *
   * The default implementation used is <tt>EpetraOperatorViewExtractorBase</tt>.
   */
  STANDARD_COMPOSITION_MEMBERS( EpetraOperatorViewExtractorBase, epetraFwdOpViewExtractor )

  //@}

  /** @name Overridden from PreconditionerFactoryBase */
  //@{

  /** \brief . */
  bool isCompatible( const LinearOpSourceBase<double> &fwdOpSrc ) const;
  /** \brief . */
  bool applySupportsConj(EConj conj) const;
  /** \brief . */
  bool applyTransposeSupportsConj(EConj conj) const;
  /** \brief . */
  Teuchos::RefCountPtr<PreconditionerBase<double> > createPrec() const;
  /** \brief . */
  void initializePrec(
    const Teuchos::RefCountPtr<const LinearOpSourceBase<double> >    &fwdOpSrc
    ,PreconditionerBase<double>                                      *prec
    ,const ESupportSolveUse                                           supportSolveUse
    ) const;
  /** \brief . */
  void uninitializePrec(
    PreconditionerBase<double>                                *prec
    ,Teuchos::RefCountPtr<const LinearOpSourceBase<double> >  *fwdOpSrc
    ,ESupportSolveUse                                         *supportSolveUse
    ) const;

  //@}

  /** @name Overridden from ParameterListAcceptor */
  //@{

  /** \brief . */
  void setParameterList(Teuchos::RefCountPtr<Teuchos::ParameterList> const& paramList);
  /** \brief . */
  Teuchos::RefCountPtr<Teuchos::ParameterList> getParameterList();
  /** \brief . */
  Teuchos::RefCountPtr<Teuchos::ParameterList> unsetParameterList();
  /** \brief . */
  Teuchos::RefCountPtr<const Teuchos::ParameterList> getParameterList() const;
  /** \brief . */
  Teuchos::RefCountPtr<const Teuchos::ParameterList> getValidParameters() const;

  //@}

  /** \name Public functions overridden from Teuchos::Describable. */
  //@{

  /** \brief . */
  std::string description() const;

  //@}

private:

  // ////////////////////////////////
  // Private data members

  Ifpack::EPrecType                                  precType_;
  int                                                overlap_;
  Teuchos::RefCountPtr<Teuchos::ParameterList>       paramList_;

  // ////////////////////////////////
  // Private member functions

  static void initializeTimers();

  static Teuchos::RefCountPtr<const Teuchos::ParameterList> generateAndGetValidParameters();

};

} // namespace Thyra

#endif // THYRA_IFPACK_PRECONDITIONER_FACTORY_DECL_HPP
