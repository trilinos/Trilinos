/*
// @HEADER
// ***********************************************************************
// 
//                Amesos: Direct Sparse Solver Package
//                 Copyright (2004) Sandia Corporation
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
// @HEADER
*/

#ifndef THYRA_AMESOS_LINEAR_OP_WITH_SOLVE_FACTORY_HPP
#define THYRA_AMESOS_LINEAR_OP_WITH_SOLVE_FACTORY_HPP

#include "Thyra_LinearOpWithSolveFactoryBase.hpp"
#include "Thyra_AmesosTypes.hpp"
#include "Amesos_BaseSolver.h"
#include "Thyra_EpetraOperatorViewExtractorBase.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"
#include "Teuchos_StandardCompositionMacros.hpp"

namespace Thyra {

/** \brief Concrete <tt>LinearOpWithSolveFactoryBase</tt> adapter subclass that uses
 * Amesos direct solvers.
 *
 * This class creates objects of type <tt>AmesosLinearOpWithSolve</tt>
 * (through the <tt>LinearOpWithSolveBase</tt> interface) which can then be
 * used to solve for linear systems.  The <tt>%AmesosLinearOpWithSolve</tt>
 * objects created an initialized by this object are completely indpendent
 * from <tt>*this</tt>.  This allows for multiple
 * <tt>%AmesosLinearOpWithSolve</tt> objects to be created and maintained
 * simultaneously and for <tt>*this</tt> factory object to be destroyed
 * without affecting the created <tt>%AmesosLinearOpWithSolve</tt> objects.
 *
 * ToDo: Mention parameter list usage.
 *
 * <b>Development notes:</b> This class has been designed to allow for "smart"
 * <tt>EpetraLinearOpBase</tt> subclasses that can create an
 * <tt>Epetra_Operator</tt> view on command where the underlying storage may
 * not be an <tt>Epetra</tt> object.  However, the current implementation of
 * at least some of the <tt>Amesos_BaseSolver</tt> subclasses do not allow the
 * <tt>%Epetra_Operator</tt> object to change after construction.  Therefore,
 * this current flawed implementation requires that every call to the
 * <tt>EpetraLinearOpBase::getEpetraOpView()</tt> function return the same
 * <tt>%Epetra_Operator</tt> object.
 *
 * \ingroup Amesos_Thyra_adapters_grp
 */
class AmesosLinearOpWithSolveFactory : public LinearOpWithSolveFactoryBase<double> {
public:

  /** \name Parameter names for Paramter List */
  //@{

  /** \brief . */
  static const std::string SolverType_name;
  /** \brief . */
  static const std::string RefactorizationPolicy_name;
  /** \brief . */
  static const std::string ThrowOnPreconditionerInput_name;
  /** \brief . */
  static const std::string Amesos_Settings_name;

  //@}

  /** \name Constructors/initializers/accessors */
  //@{

  /** \brief . */
  ~AmesosLinearOpWithSolveFactory();

  /** \brief Constructor which sets the defaults.
   */
  AmesosLinearOpWithSolveFactory(
    const Amesos::ESolverType                 solverType
#ifdef HAVE_AMESOS_KLU                        
                                                                     = Amesos::KLU
#else                                         
                                                                     = Amesos::LAPACK
#endif                                        
    ,const Amesos::ERefactorizationPolicy     refactorizationPolicy  = Amesos::REPIVOT_ON_REFACTORIZATION
    ,const bool                               throwOnPrecInput       = true
    );
    
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

  /** @name Overridden public functions from LinearOpWithSolveFactoryBase */
  //@{

  /** \brief Returns true if <tt>dynamic_cast<const
   * EpetraLinearOpBase*>(fwdOpSrc)!=NULL</tt> .
   */
  bool isCompatible( const LinearOpSourceBase<double> &fwdOpSrc ) const;

  /** \brief . */
  Teuchos::RefCountPtr<LinearOpWithSolveBase<double> > createOp() const;

  /** \brief . */
  void initializeOp(
    const Teuchos::RefCountPtr<const LinearOpSourceBase<double> >    &fwdOpSrc
    ,LinearOpWithSolveBase<double>                                   *Op
    ,const ESupportSolveUse                                          supportSolveUse
    ) const;

  /** \brief Returns <tt>false</tt> . */
  bool supportsPreconditionerInputType(const EPreconditionerInputType precOpType) const;

  /** \brief Throws exception if <tt>this->throwOnPrecInput()==true</tt> and
   * calls <tt>this->initializeOp(fwdOpSrc,Op)</tt> otherwise
   */
  void initializePreconditionedOp(
    const Teuchos::RefCountPtr<const LinearOpSourceBase<double> >       &fwdOpSrc
    ,const Teuchos::RefCountPtr<const PreconditionerBase<double> >      &prec
    ,LinearOpWithSolveBase<double>                                      *Op
    ,const ESupportSolveUse                                             supportSolveUse
    ) const;

  /** \brief Throws exception if <tt>this->throwOnPrecInput()==true</tt> and
   * calls <tt>this->initializeOp(fwdOpSrc,Op)</tt> otherwise
   */
  void initializePreconditionedOp(
    const Teuchos::RefCountPtr<const LinearOpSourceBase<double> >       &fwdOpSrc
    ,const Teuchos::RefCountPtr<const LinearOpSourceBase<double> >      &approxFwdOpSrc
    ,LinearOpWithSolveBase<double>                                      *Op
    ,const ESupportSolveUse                                             supportSolveUse
    ) const;

  /** \brief . */
  void uninitializeOp(
    LinearOpWithSolveBase<double>                               *Op
    ,Teuchos::RefCountPtr<const LinearOpSourceBase<double> >    *fwdOpSrc
    ,Teuchos::RefCountPtr<const PreconditionerBase<double> >    *prec
    ,Teuchos::RefCountPtr<const LinearOpSourceBase<double> >    *approxFwdOpSrc
    ,ESupportSolveUse                                           *supportSolveUse
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

  // /////////////////////////
  // Private data members

  Amesos::ESolverType                             solverType_;
  Amesos::ERefactorizationPolicy                  refactorizationPolicy_;
  bool                                            throwOnPrecInput_;
  Teuchos::RefCountPtr<Teuchos::ParameterList>    paramList_;

  // /////////////////////////
  // Private member functions

  static void initializeTimers();

  static Teuchos::RefCountPtr<const Teuchos::ParameterList>
  generateAndGetValidParameters();

};

} // namespace Thyra

#endif // THYRA_AMESOS_LINEAR_OP_WITH_SOLVE_FACTORY_HPP
