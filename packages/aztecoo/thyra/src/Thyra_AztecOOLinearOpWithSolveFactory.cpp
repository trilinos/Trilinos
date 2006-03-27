/*@Header
// ***********************************************************************
// 
//        AztecOO: An Object-Oriented Aztec Linear Solver Package 
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

#ifndef __sun

#include "Thyra_AztecOOLinearOpWithSolveFactory.hpp"
#include "Thyra_AztecOOLinearOpWithSolve.hpp"
#include "Thyra_EpetraOperatorViewExtractorStd.hpp"
#include "Thyra_ScaledAdjointLinearOpBase.hpp"
#include "Thyra_EpetraLinearOpBase.hpp"
#include "EpetraExt_ProductOperator.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_dyn_cast.hpp"

namespace {
const std::string AOOLOWSF_epetraPrecOp_str = "AOOLOWSF::epetraPrecOp";
const std::string AOOLOWSF_aztec_epetra_epetraFwdOp_str = "AOOLOWSF::aztec_epetra_epetraFwdOp";
const std::string AOOLOWSF_aztec_epetra_epetraAdjOp_str = "AOOLOWSF::aztec_epetra_epetraAdjOp";
const std::string AOOLOWSF_rowmatrix_epetraFwdOp_str = "AOOLOWSF::rowmatrix_epetraFwdOp";
const std::string AOOLOWSF_rowmatrix_epetraPrecOp_str = "AOOLOWSF::rowmatrix_epetraPrecOp";
const std::string AOOLOWSF_aztec_fwd_epetra_epetraPrecOp_str = "AOOLOWSF::aztec_fwd_epetra_epetraPrecOp";
const std::string AOOLOWSF_aztec_adj_epetra_epetraPrecOp_str = "AOOLOWSF::aztec_adj_epetra_epetraPrecOp";
const std::string AOOLOWSF_setPrecondtionerOperator_str = "AOOLOWSF::setPrecondtionerOperator";
const std::string AOOLOWSF_constructedAztecPreconditoner_str = "AOOLOWSF::constructedAztecPreconditoner";
} // namespace

namespace Thyra {

// Constructors/initializers/accessors

AztecOOLinearOpWithSolveFactory::AztecOOLinearOpWithSolveFactory(
  const int       fwdDefaultMaxIterations
  ,const double   fwdDefaultTol
  ,const int      adjDefaultMaxIterations
  ,const double   adjDefaultTol
  )
  :fwdDefaultMaxIterations_(fwdDefaultMaxIterations)
  ,fwdDefaultTol_(fwdDefaultTol)
  ,adjDefaultMaxIterations_(adjDefaultMaxIterations)
  ,adjDefaultTol_(adjDefaultTol)
  ,epetraFwdOpViewExtractor_(Teuchos::rcp(new EpetraOperatorViewExtractorStd()))
{}

// Overridden from LinearOpWithSolveFactoryBase

bool AztecOOLinearOpWithSolveFactory::isCompatible(
  const LinearOpBase<double> &fwdOp
  ) const
{
  return epetraFwdOpViewExtractor_->isCompatible(fwdOp);
}

Teuchos::RefCountPtr<LinearOpWithSolveBase<double> >
AztecOOLinearOpWithSolveFactory::createOp() const
{
  return Teuchos::rcp(new AztecOOLinearOpWithSolve());
}

void AztecOOLinearOpWithSolveFactory::initializeOp(
  const Teuchos::RefCountPtr<const LinearOpBase<double> >    &fwdOp
  ,LinearOpWithSolveBase<double>                             *Op
  ,const ESupportSolveUse                                     supportSolveUse
  ) const
{
  this->initializeOp_impl(fwdOp,Teuchos::null,Teuchos::null,false,Op);
}

void AztecOOLinearOpWithSolveFactory::initializeAndReuseOp(
  const Teuchos::RefCountPtr<const LinearOpBase<double> >    &fwdOp
  ,LinearOpWithSolveBase<double>                             *Op
  ) const
{
  this->initializeOp_impl(fwdOp,Teuchos::null,Teuchos::null,true,Op);
}

bool AztecOOLinearOpWithSolveFactory::supportsPreconditionerInputType(const EPreconditionerInputType precOpType) const
{
  const_cast<bool&>(useAztecPrec_) = (
    paramList_.get() && paramList_->sublist("Forward Solve").get("AZ_precond","none")!="none"
    );
  switch(precOpType) {
    case PRECONDITIONER_INPUT_TYPE_AS_OPERATOR:
      return true;
      break;
    case PRECONDITIONER_INPUT_TYPE_AS_MATRIX:
      return useAztecPrec_;
      break;
    default:
      TEST_FOR_EXCEPT(true);
  }
  return PRECONDITIONER_INPUT_TYPE_AS_OPERATOR; // Should never be called!
}

void AztecOOLinearOpWithSolveFactory::initializePreconditionedOp(
  const Teuchos::RefCountPtr<const LinearOpBase<double> >             &fwdOp
  ,const Teuchos::RefCountPtr<const PreconditionerBase<double> >      &prec
  ,LinearOpWithSolveBase<double>                                      *Op
  ,const ESupportSolveUse                                             supportSolveUse
  ) const
{
  TEST_FOR_EXCEPT(prec.get()==NULL);
  this->initializeOp_impl(fwdOp,prec,Teuchos::null,false,Op);
}

void AztecOOLinearOpWithSolveFactory::initializeApproxPreconditionedOp(
  const Teuchos::RefCountPtr<const LinearOpBase<double> >             &fwdOp
  ,const Teuchos::RefCountPtr<const LinearOpBase<double> >            &approxFwdOp
  ,LinearOpWithSolveBase<double>                                      *Op
  ,const ESupportSolveUse                                             supportSolveUse
  ) const
{
  TEST_FOR_EXCEPT(approxFwdOp.get()==NULL);
  this->initializeOp_impl(fwdOp,Teuchos::null,approxFwdOp,false,Op);
}

void AztecOOLinearOpWithSolveFactory::uninitializeOp(
  LinearOpWithSolveBase<double>                               *Op
  ,Teuchos::RefCountPtr<const LinearOpBase<double> >          *fwdOp
  ,Teuchos::RefCountPtr<const PreconditionerBase<double> >    *prec
  ,Teuchos::RefCountPtr<const LinearOpBase<double> >          *approxFwdOp
  ,ESupportSolveUse                                           *supportSolveUse
  ) const
{
#ifdef _DEBUG
  TEST_FOR_EXCEPT(Op==NULL);
#endif
  AztecOOLinearOpWithSolve
    *aztecOp = &Teuchos::dyn_cast<AztecOOLinearOpWithSolve>(*Op);
  // Extract and unset the fwdOP and approxFwdOp objects
  Teuchos::RefCountPtr<const LinearOpBase<double> >
    _fwdOp = aztecOp->extract_fwdOp(),             // Will be null if not initialized!
    _approxFwdOp = aztecOp->extract_approxFwdOp(); // Will be null if no approxFwdOp set
  if(fwdOp) *fwdOp = _fwdOp;
  if(approxFwdOp) *approxFwdOp = _approxFwdOp;
  // Only extract and uset the prec object if it is external.  If it is
  // internal, then we need to hold on to this so that we can reinitialize it
  // later.
  if(aztecOp->isExternalPrec()) {
    Teuchos::RefCountPtr<const PreconditionerBase<double> >
      _prec = aztecOp->extract_prec(); // Will be null if not external preconditioner was set
    if(prec) *prec = _prec;
  }
  // ToDo: Extract the Epetra_Operator views what where used to initialize the
  // forward and adjoint solvers!  This is needed to make this totally
  // stateless.
}

// Overridden from ParameterListAcceptor

void AztecOOLinearOpWithSolveFactory::setParameterList(Teuchos::RefCountPtr<Teuchos::ParameterList> const& paramList)
{
  TEST_FOR_EXCEPT(paramList.get()==NULL);
  paramList->validateParameters(*this->getValidParameters(),0); // Only validate very top level for now!
  paramList_ = paramList;
}

Teuchos::RefCountPtr<Teuchos::ParameterList>
AztecOOLinearOpWithSolveFactory::getParameterList()
{
  return paramList_;
}

Teuchos::RefCountPtr<Teuchos::ParameterList>
AztecOOLinearOpWithSolveFactory::unsetParameterList()
{
  Teuchos::RefCountPtr<Teuchos::ParameterList> _paramList = paramList_;
  paramList_ = Teuchos::null;
  return _paramList;
}

Teuchos::RefCountPtr<const Teuchos::ParameterList>
AztecOOLinearOpWithSolveFactory::getParameterList() const
{
  return paramList_;
}

Teuchos::RefCountPtr<const Teuchos::ParameterList>
AztecOOLinearOpWithSolveFactory::getValidParameters() const
{
  return generateAndGetValidParameters();
}

// Public functions overridden from Teuchos::Describable

std::string AztecOOLinearOpWithSolveFactory::description() const
{
  std::ostringstream oss;
  oss << "Thyra::AztecOOLinearOpWithSolveFactory";
  return oss.str();
}

// private

Teuchos::RefCountPtr<const Teuchos::ParameterList>
AztecOOLinearOpWithSolveFactory::generateAndGetValidParameters()
{
  static Teuchos::RefCountPtr<Teuchos::ParameterList> validParamList;
  if(validParamList.get()==NULL) {
    validParamList = Teuchos::rcp(new Teuchos::ParameterList("Thyra::AztecOOLinearOpWithSolveFactory"));
    Teuchos::ParameterList
      fwdSolvePL = validParamList->sublist("Forward Solve");
    // ToDo: Figure out how to add parameters for this!
    Teuchos::ParameterList
      adjSolvePL = validParamList->sublist("Adjoint Solve");
    // ToDo: Figure out how to add parameters for this!
  }
  return validParamList;
}

void AztecOOLinearOpWithSolveFactory::initializeOp_impl(
  const Teuchos::RefCountPtr<const LinearOpBase<double> >             &fwdOp
  ,const Teuchos::RefCountPtr<const PreconditionerBase<double> >      &prec
  ,const Teuchos::RefCountPtr<const LinearOpBase<double> >            &approxFwdOp
  ,const bool                                                         reusePrec
  ,LinearOpWithSolveBase<double>                                      *Op
  ) const
{
  using Teuchos::RefCountPtr;
  using Teuchos::null;
  using Teuchos::rcp;
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::rcp_const_cast;
  using Teuchos::set_extra_data;
  using Teuchos::get_optional_extra_data;
  typedef EpetraExt::ProductOperator PO;
#ifdef _DEBUG
  TEST_FOR_EXCEPT(Op==NULL);
#endif
  //
  // Unwrap and get the forward operator or matrix
  //
  Teuchos::RefCountPtr<const Epetra_Operator> epetra_epetraFwdOp;
  ETransp                                     epetra_epetraFwdOpTransp;
  EApplyEpetraOpAs                            epetra_epetraFwdOpApplyAs;
  EAdjointEpetraOp                            epetra_epetraFwdOpAdjointSupport;
  double                                      epetra_epetraFwdOpScalar;
  epetraFwdOpViewExtractor_->getEpetraOpView(
    fwdOp
    ,&epetra_epetraFwdOp,&epetra_epetraFwdOpTransp,&epetra_epetraFwdOpApplyAs
    ,&epetra_epetraFwdOpAdjointSupport,&epetra_epetraFwdOpScalar
    );
  TEST_FOR_EXCEPTION(
    epetra_epetraFwdOp.get()==NULL, std::logic_error
    ,"Error, The input fwdOp object must be fully initialized before calling this function!"
    );
  //
  // Unwrap and get the external preconditioner operator
  //
  RefCountPtr<const LinearOpBase<double> >     rightPrecOp;
  if(prec.get()) {
    RefCountPtr<const LinearOpBase<double> > unspecified = prec->getUnspecifiedPrecOp();
    RefCountPtr<const LinearOpBase<double> > left        = prec->getLeftPrecOp();
    RefCountPtr<const LinearOpBase<double> > right       = prec->getRightPrecOp();
    TEST_FOR_EXCEPTION(
      !( left.get() || right.get() || unspecified.get() ), std::logic_error
      ,"Error, at least one preconditoner linear operator objects must be set!"
      );
    if(unspecified.get()) {
      rightPrecOp = unspecified;
    }
    else {
      // Set a left, right or split preconditioner
      TEST_FOR_EXCEPTION(
        left.get(),std::logic_error
        ,"Error, we can not currently handle a left preconditioner with the Thyra/AztecOO adapters!"
        );
      rightPrecOp = right;
    }
  }
  double                                       wrappedPrecOpScalar = 0.0;
  ETransp                                      wrappedPrecOpTransp = NOTRANS;
  RefCountPtr<const LinearOpBase<double> >     wrappedPrecOp = null;
  RefCountPtr<const EpetraLinearOpBase>        epetraPrecOp;
  Teuchos::RefCountPtr<const Epetra_Operator>  epetra_epetraPrecOp;
  ETransp                                      epetra_epetraPrecOpTransp;
  EApplyEpetraOpAs                             epetra_epetraPrecOpApplyAs;
  EAdjointEpetraOp                             epetra_epetraPrecOpAdjointSupport;
  ETransp                                      overall_epetra_epetraPrecOpTransp;
  if(rightPrecOp.get()) {
    unwrap(rightPrecOp,&wrappedPrecOpScalar,&wrappedPrecOpTransp,&wrappedPrecOp);
    epetraPrecOp = rcp_dynamic_cast<const EpetraLinearOpBase>(wrappedPrecOp,true);
    epetraPrecOp->getEpetraOpView(&epetra_epetraPrecOp,&epetra_epetraPrecOpTransp,&epetra_epetraPrecOpApplyAs,&epetra_epetraPrecOpAdjointSupport);
    TEST_FOR_EXCEPTION(
      epetra_epetraPrecOp.get()==NULL,std::logic_error
      ,"Error, The input prec object and its embedded precondtioner operator must be fully initialized before calling this function!"
      );
    set_extra_data(epetraPrecOp,AOOLOWSF_epetraPrecOp_str,&epetra_epetraPrecOp,Teuchos::POST_DESTROY,false);
    overall_epetra_epetraPrecOpTransp
      = trans_trans(real_trans(wrappedPrecOpTransp),real_trans(epetra_epetraPrecOpTransp));
  }
  //
  // Unwrap and get the approximate forward operator to be used to generate a preconditioner
  //
  if(approxFwdOp.get()) {
    // Note, here we just use the same members data that would be set for an
    // extenral preconditioner operator since it is not getting used.
    unwrap(approxFwdOp,&wrappedPrecOpScalar,&wrappedPrecOpTransp,&wrappedPrecOp);
    epetraPrecOp = rcp_dynamic_cast<const EpetraLinearOpBase>(wrappedPrecOp,true);
    epetraPrecOp->getEpetraOpView(&epetra_epetraPrecOp,&epetra_epetraPrecOpTransp,&epetra_epetraPrecOpApplyAs,&epetra_epetraPrecOpAdjointSupport);
    TEST_FOR_EXCEPTION(
      epetra_epetraPrecOp.get()==NULL,std::logic_error
      ,"Error, The input approxFwdOp object must be fully initialized before calling this function!"
      );
    set_extra_data(epetraPrecOp,AOOLOWSF_epetraPrecOp_str,&epetra_epetraPrecOp,Teuchos::POST_DESTROY,false);
    overall_epetra_epetraPrecOpTransp
      = trans_trans(real_trans(wrappedPrecOpTransp),real_trans(epetra_epetraPrecOpTransp));
  }
  //
  // Get the AztecOOLinearOpWithSolve object
  //
  AztecOOLinearOpWithSolve
    *aztecOp = &Teuchos::dyn_cast<AztecOOLinearOpWithSolve>(*Op);
  //
  // Determine if the forward and preconditioner operators are a row matrices or not
  //
  RefCountPtr<const Epetra_RowMatrix>
    rowmatrix_epetraFwdOp  = rcp_dynamic_cast<const Epetra_RowMatrix>(epetra_epetraFwdOp),
    rowmatrix_epetraPrecOp = rcp_dynamic_cast<const Epetra_RowMatrix>(epetra_epetraPrecOp);
  //
  // Determine the type of preconditoner
  //
  this->supportsPreconditionerInputType(PRECONDITIONER_INPUT_TYPE_AS_MATRIX); // Updates useAztecPrec_, input value does not matter
  enum ELocalPrecType { PT_NONE, PT_AZTEC_FROM_OP, PT_AZTEC_FROM_APPROX_FWD_MATRIX, PT_FROM_PREC_OP };
  ELocalPrecType localPrecType;
  if( prec.get()==NULL && approxFwdOp.get()==NULL && !useAztecPrec_ ) {
    // No preconditioning at all!
    localPrecType = PT_NONE;
  }
  else if( prec.get()==NULL && approxFwdOp.get()==NULL && useAztecPrec_ ) {
    // We are using the forward matrix for the preconditioner using aztec preconditioners
    localPrecType = PT_AZTEC_FROM_OP;
  }
  else if( approxFwdOp.get() && useAztecPrec_ ) {
    // The preconditioner comes from the input as a matrix and we are using aztec preconditioners
    localPrecType = PT_AZTEC_FROM_APPROX_FWD_MATRIX;
  }
  else if( prec.get() ) {
    // The preconditioner comes as an external operator so let's use it as such
    localPrecType = PT_FROM_PREC_OP;
  }
  //
  // Determine if aztecOp already contains solvers and if we need to reinitialize or not
  //
  RefCountPtr<AztecOO> aztecFwdSolver, aztecAdjSolver;
  bool startingOver;
  if(1){
    // Let's assume that fwdOp, prec and/or approxFwdOp are compatible with
    // the already created AztecOO objects.  If they are not, then the client
    // should have created a new LOWSB object from scratch!
    Teuchos::RefCountPtr<const LinearOpBase<double> >        old_fwdOp;
    Teuchos::RefCountPtr<const PreconditionerBase<double> >  old_prec;
    bool                                                     old_isExternalPrec;
    Teuchos::RefCountPtr<const LinearOpBase<double> >        old_approxFwdOp;
    Teuchos::RefCountPtr<AztecOO>                            old_aztecFwdSolver;
    Teuchos::RefCountPtr<AztecOO>                            old_aztecAdjSolver;
    double                                                   old_aztecSolverScalar;
    aztecOp->uninitialize(
      &old_fwdOp
      ,&old_prec
      ,&old_isExternalPrec
      ,&old_approxFwdOp
      ,&old_aztecFwdSolver
      ,NULL
      ,&old_aztecAdjSolver
      ,NULL
      ,&old_aztecSolverScalar
      );
    if( old_aztecFwdSolver.get()==NULL ) {
      // This has never been initialized before
      startingOver = true;
    }
    else {
      // Let's assume that fwdOp, prec and/or approxFwdOp are compatible with
      // the already created AztecOO objects.  If they are not, then the
      // client should have created a new LOWSB object from scratch!
      aztecFwdSolver = old_aztecFwdSolver;
      aztecAdjSolver = old_aztecAdjSolver;
      startingOver = false;
      // We must wipe out the old preconditoner if we are not reusing the preconditioner
      bool *constructedAztecPreconditioner = NULL;
      if(
        !reusePrec
        && ( constructedAztecPreconditioner = get_optional_extra_data<bool>(aztecFwdSolver,"AOOLOWSF::constructedAztecPreconditoner") )
        && *constructedAztecPreconditioner
        )
      {
        aztecFwdSolver->DestroyPreconditioner();
        *constructedAztecPreconditioner = false;
      }
      // We must see if we set an external preconditioner but will not do so
      // again in which case we must blow away AztecOO and start over again!
      bool *setPreconditionerOperator = NULL;
      if(
        localPrecType != PT_FROM_PREC_OP
        && ( setPreconditionerOperator = get_optional_extra_data<bool>(aztecFwdSolver,"AOOLOWSF::setPreconditonerOperator") )
        && *setPreconditionerOperator
        )
      {
        // We must start over again since there is no way to unset an external preconditioner!
        startingOver = true;
      }
    }
  }
  //
  // Create the AztecOO solvers if we are starting over
  //
  startingOver = true; // ToDo: Remove this and figure out why this is not working!
  if(startingOver) {
    // Forward solver
    aztecFwdSolver = rcp(new AztecOO());
    aztecFwdSolver->SetAztecOption(AZ_output,AZ_none); // Don't mess up output
    aztecFwdSolver->SetAztecOption(AZ_conv,AZ_rhs);    // Specified by this interface (may change)
    aztecFwdSolver->SetAztecOption(AZ_diagnostics,AZ_none); // This was turned off in NOX?
    aztecFwdSolver->SetAztecOption(AZ_keep_info, 1);
    // Adjoint solver (if supported)
    if(
      epetra_epetraFwdOpAdjointSupport==EPETRA_OP_ADJOINT_SUPPORTED
      && localPrecType!=PT_AZTEC_FROM_OP && localPrecType!=PT_AZTEC_FROM_APPROX_FWD_MATRIX
      //&& (epetra_epetraPrecOp.get()==NULL ||epetra_epetraPrecOpAdjointSupport==EPETRA_OP_ADJOINT_SUPPORTED)
      )
    {
      aztecAdjSolver = rcp(new AztecOO());
      aztecAdjSolver->SetAztecOption(AZ_output,AZ_none);
      aztecAdjSolver->SetAztecOption(AZ_conv,AZ_rhs);
      aztecAdjSolver->SetAztecOption(AZ_diagnostics,AZ_none);
      //aztecAdjSolver->SetAztecOption(AZ_keep_info, 1);
    }
  }
  //
  // Set the options on the AztecOO solvers
  //
  if( startingOver ) {
    if(paramList_.get())
      aztecFwdSolver->SetParameters(paramList_->sublist("Forward Solve"),true);
    if(aztecAdjSolver.get()) {
      if(paramList_.get())
        aztecAdjSolver->SetParameters(paramList_->sublist("Adjoint Solve"),true);
    }
    // ToDo: Once you can get a full list of valid parameters then change the
    // above "cerr_warn_if_unused" from true to false since we will compare
    // against the valid parameter list ahead of time.
  }
  //
  // Process the forward operator
  //
  RefCountPtr<const Epetra_Operator>
    aztec_epetra_epetraFwdOp,
    aztec_epetra_epetraAdjOp;
  // Forward solve
  RefCountPtr<const Epetra_Operator>
    epetraOps[]
    = { epetra_epetraFwdOp };
  Teuchos::ETransp
    epetraOpsTransp[]
    = { epetra_epetraFwdOpTransp==NOTRANS ? Teuchos::NO_TRANS : Teuchos::TRANS };
  PO::EApplyMode
    epetraOpsApplyMode[]
    = { epetra_epetraFwdOpApplyAs==EPETRA_OP_APPLY_APPLY ? PO::APPLY_MODE_APPLY : PO::APPLY_MODE_APPLY_INVERSE };
  if( epetraOpsTransp[0] == Teuchos::NO_TRANS && epetraOpsApplyMode[0] == PO::APPLY_MODE_APPLY )
    aztec_epetra_epetraFwdOp = epetra_epetraFwdOp;
  else
    aztec_epetra_epetraFwdOp = rcp(new PO(1,epetraOps,epetraOpsTransp,epetraOpsApplyMode));
  if( startingOver || aztec_epetra_epetraFwdOp.get() != aztecFwdSolver->GetUserOperator() ) {
    // Here we will be careful not to reset the forward operator in fears that
    // it will blow out the internally created stuff.
    aztecFwdSolver->SetUserOperator(const_cast<Epetra_Operator*>(&*aztec_epetra_epetraFwdOp));
    set_extra_data(aztec_epetra_epetraFwdOp,AOOLOWSF_aztec_epetra_epetraFwdOp_str,&aztecFwdSolver,Teuchos::POST_DESTROY,false);
  }
  // Adjoint solve
  if( aztecAdjSolver.get() ) {
    epetraOpsTransp[0] = ( epetra_epetraFwdOpTransp==NOTRANS ? Teuchos::TRANS : Teuchos::NO_TRANS );
    if( epetraOpsTransp[0] == Teuchos::NO_TRANS && epetraOpsApplyMode[0] == PO::APPLY_MODE_APPLY )
      aztec_epetra_epetraAdjOp = epetra_epetraFwdOp;
    else
      aztec_epetra_epetraAdjOp = rcp(new PO(1,epetraOps,epetraOpsTransp,epetraOpsApplyMode));
    aztecAdjSolver->SetUserOperator(const_cast<Epetra_Operator*>(&*aztec_epetra_epetraAdjOp));
    set_extra_data(aztec_epetra_epetraAdjOp,AOOLOWSF_aztec_epetra_epetraAdjOp_str,&aztecAdjSolver,Teuchos::POST_DESTROY,false);
  }
  //
  // Process the preconditioner
  //
  RefCountPtr<const Epetra_Operator>
    aztec_fwd_epetra_epetraPrecOp,
    aztec_adj_epetra_epetraPrecOp;
  bool setAztecPreconditioner = false;
  switch(localPrecType) {
    case PT_NONE: {
      //
      // No preconditioning at all!
      //
      break;
    }
    case PT_AZTEC_FROM_OP: {
      //
      // We are using the forward matrix for the preconditioner using aztec preconditioners
      //
      if( startingOver || !reusePrec ) {
        TEST_FOR_EXCEPTION(
          rowmatrix_epetraFwdOp.get()==NULL, std::logic_error
          ,"AztecOOLinearOpWithSolveFactor::initializeOp_impl(...): Error, There is no preconditioner given by client, but the client "
          "passed in an Epetra_Operator for the forward operator of type \'" <<typeid(*epetra_epetraFwdOp).name()<<"\' that does not "
          "support the Epetra_RowMatrix interface!"
          );
        TEST_FOR_EXCEPTION(
          epetra_epetraFwdOpTransp!=NOTRANS, std::logic_error
          ,"AztecOOLinearOpWithSolveFactor::initializeOp_impl(...): Error, There is no preconditioner given by client and the client "
          "passed in an Epetra_RowMatrix for the forward operator but the overall transpose is not NOTRANS and therefore we can can just "
          "hand this over to aztec without making a copy which is not supported here!"
          );
        aztecFwdSolver->SetPrecMatrix(const_cast<Epetra_RowMatrix*>(&*rowmatrix_epetraFwdOp));
        set_extra_data(rowmatrix_epetraFwdOp,AOOLOWSF_rowmatrix_epetraFwdOp_str,&aztecFwdSolver,Teuchos::POST_DESTROY,false);
      }
      setAztecPreconditioner = true;
      break;
    }
    case PT_AZTEC_FROM_APPROX_FWD_MATRIX: {
      //
      // The preconditioner comes from the input as a matrix and we are using aztec preconditioners
      //
      if( startingOver || !reusePrec ) {
        TEST_FOR_EXCEPTION(
          rowmatrix_epetraPrecOp.get()==NULL, std::logic_error
          ,"AztecOOLinearOpWithSolveFactor::initializeOp_impl(...): The client "
          "passed in an Epetra_Operator for the preconditioner matrix of type \'" <<typeid(*epetra_epetraPrecOp).name()<<"\' that does not "
          "support the Epetra_RowMatrix interface!"
          );
        TEST_FOR_EXCEPTION(
          overall_epetra_epetraPrecOpTransp!=NOTRANS, std::logic_error
          ,"AztecOOLinearOpWithSolveFactor::initializeOp_impl(...): Error, The client "
          "passed in an Epetra_RowMatrix for the preconditoner matrix but the overall transpose is not NOTRANS and therefore we can can just "
          "hand this over to aztec without making a copy which is not supported here!"
          );
        aztecFwdSolver->SetPrecMatrix(const_cast<Epetra_RowMatrix*>(&*rowmatrix_epetraPrecOp));
        set_extra_data(rowmatrix_epetraPrecOp,AOOLOWSF_rowmatrix_epetraPrecOp_str,&aztecFwdSolver,Teuchos::POST_DESTROY,false);
      }
      setAztecPreconditioner = true;
      break;
    }
    case PT_FROM_PREC_OP: {
      //
      // The preconditioner comes as an operator so let's use it as such
      //
      // Forawrd solve
      RefCountPtr<const Epetra_Operator>
        epetraOps[]
        = { epetra_epetraPrecOp };
      Teuchos::ETransp
        epetraOpsTransp[]
        = { overall_epetra_epetraPrecOpTransp==NOTRANS ? Teuchos::NO_TRANS : Teuchos::TRANS };
      PO::EApplyMode
        epetraOpsApplyMode[] // Here we must toggle the apply mode since aztecoo applies the preconditioner using ApplyInverse(...)
        = { epetra_epetraPrecOpApplyAs==EPETRA_OP_APPLY_APPLY ? PO::APPLY_MODE_APPLY_INVERSE : PO::APPLY_MODE_APPLY };
      if( epetraOpsTransp[0] == Teuchos::NO_TRANS && epetra_epetraPrecOpApplyAs==EPETRA_OP_APPLY_APPLY_INVERSE )
        aztec_fwd_epetra_epetraPrecOp = epetra_epetraPrecOp;
      else
        aztec_fwd_epetra_epetraPrecOp = rcp(new PO(1,epetraOps,epetraOpsTransp,epetraOpsApplyMode));
      aztecFwdSolver->SetPrecOperator(const_cast<Epetra_Operator*>(&*aztec_fwd_epetra_epetraPrecOp));
      set_extra_data(aztec_fwd_epetra_epetraPrecOp,AOOLOWSF_aztec_fwd_epetra_epetraPrecOp_str,&aztecFwdSolver,Teuchos::POST_DESTROY,false);
      // Adjoint solve
      if( aztecAdjSolver.get() && epetra_epetraPrecOpAdjointSupport == EPETRA_OP_ADJOINT_SUPPORTED ) {
        epetraOpsTransp[0] = ( overall_epetra_epetraPrecOpTransp==NOTRANS ? Teuchos::TRANS : Teuchos::NO_TRANS );
        if( epetraOpsTransp[0] == Teuchos::NO_TRANS && epetra_epetraPrecOpApplyAs==EPETRA_OP_APPLY_APPLY_INVERSE )
          aztec_adj_epetra_epetraPrecOp = epetra_epetraPrecOp;
        else
          aztec_adj_epetra_epetraPrecOp = rcp(new PO(1,epetraOps,epetraOpsTransp,epetraOpsApplyMode));
        aztecAdjSolver->SetPrecOperator(const_cast<Epetra_Operator*>(&*aztec_adj_epetra_epetraPrecOp));
        set_extra_data(aztec_adj_epetra_epetraPrecOp,AOOLOWSF_aztec_adj_epetra_epetraPrecOp_str,&aztecAdjSolver,Teuchos::POST_DESTROY,false);
        set_extra_data<bool>(true,AOOLOWSF_setPrecondtionerOperator_str,&aztecFwdSolver,Teuchos::POST_DESTROY,false);
      }
      break;
    }
    default:
      TEST_FOR_EXCEPT(true);
  }
  //
  // Initialize the aztec preconditioner
  //
  if(setAztecPreconditioner) {
    if( startingOver || !reusePrec ) {
      double condNumEst = -1.0;
      TEST_FOR_EXCEPT(0!=aztecFwdSolver->ConstructPreconditioner(condNumEst));
      //aztecFwdSolver->SetAztecOption(AZ_pre_calc, AZ_calc);
      set_extra_data<bool>(true,AOOLOWSF_constructedAztecPreconditoner_str,&aztecFwdSolver,Teuchos::POST_DESTROY,false);
    }
    else {
      //aztecFwdSolver->SetAztecOption(AZ_pre_calc, AZ_reuse);
    }
  }
  //
  // Initialize the AztecOOLinearOpWithSolve object and set its options
  //
  if(aztecAdjSolver.get()) {
    aztecOp->initialize(
      fwdOp,prec,true,approxFwdOp
      ,aztecFwdSolver,true,aztecAdjSolver,true,epetra_epetraFwdOpScalar
      );
  }
  else {
    aztecOp->initialize(
      fwdOp,prec,true,approxFwdOp
      ,aztecFwdSolver,true,null,false,epetra_epetraFwdOpScalar
      );
  }
  aztecOp->fwdDefaultMaxIterations(fwdDefaultMaxIterations());
  aztecOp->fwdDefaultTol(fwdDefaultTol());
  aztecOp->adjDefaultMaxIterations(adjDefaultMaxIterations());
  aztecOp->adjDefaultTol(adjDefaultTol());
#ifdef _DEBUG
  paramList_->validateParameters(*this->getValidParameters(),0); // Only validate very top level for now!
#endif
}

} // namespace Thyra

#endif // __sun
