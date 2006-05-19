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

#include "Thyra_IfpackPreconditionerFactory.hpp"
#include "Thyra_EpetraOperatorViewExtractorStd.hpp"
#include "Thyra_EpetraLinearOp.hpp"
#include "Thyra_DefaultPreconditioner.hpp"
#include "Ifpack_ValidParameters.h"
#include "Ifpack_Preconditioner.h"
#include "Ifpack.h"
#include "Epetra_RowMatrix.h"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_dyn_cast.hpp"

namespace {

Teuchos::RefCountPtr<Teuchos::Time> overallTimer, creationTimer, factorizationTimer;

} // namespace

namespace Thyra {

// Constructors/initializers/accessors

IfpackPreconditionerFactory::IfpackPreconditionerFactory()
  :epetraFwdOpViewExtractor_(Teuchos::rcp(new EpetraOperatorViewExtractorStd()))
   ,precType_(Ifpack::ILU)
   ,overlap_(0)
{
  initializeTimers();
}

// Overridden from PreconditionerFactoryBase

bool IfpackPreconditionerFactory::isCompatible( const LinearOpBase<double> &fwdOp ) const
{
  Teuchos::RefCountPtr<const Epetra_Operator> epetraFwdOp;
  ETransp                                     epetraFwdOpTransp;
  EApplyEpetraOpAs                            epetraFwdOpApplyAs;
  EAdjointEpetraOp                            epetraFwdOpAdjointSupport;
  double                                      epetraFwdOpScalar;
  epetraFwdOpViewExtractor_->getEpetraOpView(
    Teuchos::rcp(&fwdOp,false)
    ,&epetraFwdOp,&epetraFwdOpTransp,&epetraFwdOpApplyAs,&epetraFwdOpAdjointSupport,&epetraFwdOpScalar
    );
  if( !dynamic_cast<const Epetra_RowMatrix*>(&*epetraFwdOp) )
    return false;
  return true;
}

bool IfpackPreconditionerFactory::applySupportsConj(EConj conj) const
{
  return true;
}

bool IfpackPreconditionerFactory::applyTransposeSupportsConj(EConj conj) const
{
  return false; // See comment below
}

Teuchos::RefCountPtr<PreconditionerBase<double> >
IfpackPreconditionerFactory::createPrec() const
{
  return Teuchos::rcp(new DefaultPreconditioner<double>());
}

void IfpackPreconditionerFactory::initializePrec(
  const Teuchos::RefCountPtr<const LinearOpBase<double> >    &fwdOp
  ,PreconditionerBase<double>                                *prec
  ,const ESupportSolveUse                                    supportSolveUse
  ) const
{
  using Teuchos::OSTab;
  using Teuchos::dyn_cast;
  using Teuchos::RefCountPtr;
  using Teuchos::null;
  using Teuchos::rcp;
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::rcp_const_cast;
  using Teuchos::set_extra_data;
  using Teuchos::get_optional_extra_data;
  Teuchos::Time totalTimer(""), timer("");
  totalTimer.start(true);
  Teuchos::TimeMonitor overallTimeMonitor(*overallTimer);
  const Teuchos::RefCountPtr<Teuchos::FancyOStream> out       = this->getOStream();
  const Teuchos::EVerbosityLevel                    verbLevel = this->getVerbLevel();
  Teuchos::OSTab tab(out);
  if(out.get() && static_cast<int>(verbLevel) >= static_cast<int>(Teuchos::VERB_LOW))
    *out << "\nEntering Thyra::IfpackPreconditionerFactory::initializePrec(...) ...\n";
#ifdef _DEBUG
  TEST_FOR_EXCEPT(fwdOp.get()==NULL);
  TEST_FOR_EXCEPT(prec==NULL);
#endif
  //
  // Unwrap and get the forward Epetra_Operator object
  //
  Teuchos::RefCountPtr<const Epetra_Operator> epetraFwdOp;
  ETransp                                     epetraFwdOpTransp;
  EApplyEpetraOpAs                            epetraFwdOpApplyAs;
  EAdjointEpetraOp                            epetraFwdOpAdjointSupport;
  double                                      epetraFwdOpScalar;
  epetraFwdOpViewExtractor_->getEpetraOpView(
    fwdOp,&epetraFwdOp,&epetraFwdOpTransp,&epetraFwdOpApplyAs,&epetraFwdOpAdjointSupport,&epetraFwdOpScalar
    );
  // Validate what we get is what we need
  RefCountPtr<const Epetra_RowMatrix>
    epetraFwdRowMat = rcp_dynamic_cast<const Epetra_RowMatrix>(epetraFwdOp,true);
  TEST_FOR_EXCEPTION(
    epetraFwdOpApplyAs != EPETRA_OP_APPLY_APPLY, std::logic_error
    ,"Error, incorrect apply mode for an Epetra_RowMatrix"
    );
  //
  // Get the concrete precondtioner object
  //
  DefaultPreconditioner<double>
    *defaultPrec = &Teuchos::dyn_cast<DefaultPreconditioner<double> >(*prec);
  //
  // Get the EpetraLinearOp object that is used to implement the preconditoner linear op
  //
  RefCountPtr<EpetraLinearOp>
    epetra_precOp = rcp_dynamic_cast<EpetraLinearOp>(defaultPrec->getNonconstUnspecifiedPrecOp(),true);
  //
  // Get the embedded Ifpack_Preconditioner object if it exists
  //
  Teuchos::RefCountPtr<Ifpack_Preconditioner>
    ifpack_precOp;
  if(epetra_precOp.get())
    ifpack_precOp = rcp_dynamic_cast<Ifpack_Preconditioner>(epetra_precOp->epetra_op(),true);
  //
  // Get the attached forward operator if it exists and make sure that it matches
  //
  if(ifpack_precOp.get()) {
    // ToDo: Get the forward operator and make sure that it matches what is
    // already being used!
  }
  //
  // Permform initialization if needed
  //
  const bool startingOver = (ifpack_precOp.get() == NULL);
  if(startingOver) {
    if(out.get() && static_cast<int>(verbLevel) >= static_cast<int>(Teuchos::VERB_LOW))
      *out << "\nCreating the initial Ifpack_Preconditioner object of type \'"<<toString(precType_)<<"\' ...\n";
    timer.start(true);
    Teuchos::TimeMonitor creationTimeMonitor(*creationTimer);
    // Create the initial preconditioner
    ::Ifpack ifpackFcty; // Should be a static function!
    const std::string &precTypeName = toString(precType_);
    ifpack_precOp = rcp(
      ifpackFcty.Create(
        precTypeName
        ,const_cast<Epetra_RowMatrix*>(&*epetraFwdRowMat)
        ,overlap_
        )
      );
    timer.stop();
    if(out.get() && static_cast<int>(verbLevel) >= static_cast<int>(Teuchos::VERB_LOW))
      *OSTab(out).getOStream() <<"\n=> Creation time = "<<timer.totalElapsedTime()<<" sec\n";
    // RAB: Above, I am just passing a string to Ifpack::Create(...) in order
    // get this code written.  However, in the future, it would be good to
    // copy the contents of what is in Ifpack::Create(...) into a local
    // function and then use switch(...) to create the initial
    // Ifpack_Preconditioner object.  This would result in better validation
    // and faster code.
    TEST_FOR_EXCEPTION(
      ifpack_precOp.get()==NULL, std::logic_error
      ,"Error, Ifpack::Create(precTypeName,...) returned NULL for precType = \""<<precTypeName<<"\"!"
      );
    // Set parameters if the list exists
    if(paramList_.get())
      TEST_FOR_EXCEPT(0!=ifpack_precOp->SetParameters(paramList_->sublist("Ifpack"))); // This will create new sublist if it does not exist!
    // Initailize the structure for the preconditioner
    TEST_FOR_EXCEPT(0!=ifpack_precOp->Initialize());
  }
  //
  // Attach the epetraFwdOp to the ifpack_precOp to guarantee that it will not go away
  //
  set_extra_data(epetraFwdOp,"IFPF::epetraFwdOp",&ifpack_precOp,Teuchos::POST_DESTROY,false);
  //
  // Update the factorization
  //
  if(1){
    if(out.get() && static_cast<int>(verbLevel) >= static_cast<int>(Teuchos::VERB_LOW))
      *out << "\nComputing the factorization of the preconditioner ...\n";
    Teuchos::TimeMonitor factorizationTimeMonitor(*factorizationTimer);
    timer.start(true);
    TEST_FOR_EXCEPT(0!=ifpack_precOp->Compute());
    timer.stop();
    if(out.get() && static_cast<int>(verbLevel) >= static_cast<int>(Teuchos::VERB_LOW))
      *OSTab(out).getOStream() <<"\n=> Factorization time = "<<timer.totalElapsedTime()<<" sec\n";
  }
  //
  // Compute the conditioner number estimate if asked
  //

  // ToDo: Implement

  //
  // Attach fwdOp to the ifpack_precOp
  //
  set_extra_data(fwdOp,"IFPF::fwdOp",&ifpack_precOp,Teuchos::POST_DESTROY,false);
  //
  // Initialize the output EpetraLinearOp
  //
  if(startingOver) {
    epetra_precOp = rcp(new EpetraLinearOp);
  }
  epetra_precOp->initialize(
    ifpack_precOp
    ,epetraFwdOpTransp
    ,EPETRA_OP_APPLY_APPLY_INVERSE
    ,EPETRA_OP_ADJOINT_UNSUPPORTED  // ToDo: Look into adjoints again.
    );
  //
  // Initialize the preconditioner
  //
  defaultPrec->initializeUnspecified(
    Teuchos::rcp_implicit_cast<LinearOpBase<double> >(epetra_precOp)
    );
  totalTimer.stop();
  if(out.get() && static_cast<int>(verbLevel) >= static_cast<int>(Teuchos::VERB_LOW))
    *out
      << "\nTotal time = "<<totalTimer.totalElapsedTime()<<" sec\n"
      << "\nLeaving Thyra::IfpackPreconditionerFactory::initializePrec(...) ...\n";
}

void IfpackPreconditionerFactory::uninitializePrec(
  PreconditionerBase<double>                          *prec
  ,Teuchos::RefCountPtr<const LinearOpBase<double> >  *fwdOp
  ,ESupportSolveUse                                   *supportSolveUse
  ) const
{
  TEST_FOR_EXCEPT(true);
}

// Overridden from ParameterListAcceptor

void IfpackPreconditionerFactory::setParameterList(Teuchos::RefCountPtr<Teuchos::ParameterList> const& paramList)
{
  TEST_FOR_EXCEPT(paramList.get()==NULL);
  paramList->validateParameters(*this->getValidParameters(),1);
  paramList_ = paramList;
  overlap_ = paramList_->get("Overlap",overlap_);
  std::ostringstream oss;
  oss << "(sub)list \""<<paramList->name()<<"\"parameter \"Prec Type\"";
  precType_ = Ifpack::precTypeNameToEnum(
    paramList_->get("Prec Type",toString(precType_))
    ,oss.str()
    );
}

Teuchos::RefCountPtr<Teuchos::ParameterList>
IfpackPreconditionerFactory::getParameterList()
{
  return paramList_;
}

Teuchos::RefCountPtr<Teuchos::ParameterList>
IfpackPreconditionerFactory::unsetParameterList()
{
  Teuchos::RefCountPtr<Teuchos::ParameterList> _paramList = paramList_;
  paramList_ = Teuchos::null;
  return _paramList;
}

Teuchos::RefCountPtr<const Teuchos::ParameterList>
IfpackPreconditionerFactory::getParameterList() const
{
  return paramList_;
}

Teuchos::RefCountPtr<const Teuchos::ParameterList>
IfpackPreconditionerFactory::getValidParameters() const
{
  return generateAndGetValidParameters();
}

// Public functions overridden from Teuchos::Describable

std::string IfpackPreconditionerFactory::description() const
{
  std::ostringstream oss;
  oss << "Thyra::IfpackPreconditionerFactory{";
  oss << "precType=\"" << toString(precType_) << "\"";
  oss << ",overlap=" << overlap_;
  oss << "}";
  return oss.str();
}

// private

void IfpackPreconditionerFactory::initializeTimers()
{
  if(!overallTimer.get()) {
    overallTimer       = Teuchos::TimeMonitor::getNewTimer("IfpackPF");
    creationTimer      = Teuchos::TimeMonitor::getNewTimer("IfpackPF:Creation");
    factorizationTimer = Teuchos::TimeMonitor::getNewTimer("IfpackPF:Factorization");
  }
}

Teuchos::RefCountPtr<const Teuchos::ParameterList>
IfpackPreconditionerFactory::generateAndGetValidParameters()
{
  static Teuchos::RefCountPtr<Teuchos::ParameterList> validParamList;
  if(validParamList.get()==NULL) {
    validParamList = Teuchos::rcp(new Teuchos::ParameterList("IfpackPreconditionerFactory"));
    validParamList->set("Prec Type","ILU");
    validParamList->set("Overlap",0);
    validParamList->sublist("Ifpack").setParameters(Ifpack_GetValidParameters());
  }
  return validParamList;
}

} // namespace Thyra
