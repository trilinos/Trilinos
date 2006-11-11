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

#include "Thyra_MLPreconditionerFactory.hpp"

#include "Thyra_EpetraOperatorViewExtractorStd.hpp"
#include "Thyra_EpetraLinearOp.hpp"
#include "Thyra_DefaultPreconditioner.hpp"
#include "ml_MultiLevelPreconditioner.h"
#include "ml_MultiLevelOperator.h"
#include "Epetra_RowMatrix.h"
#include "Teuchos_dyn_cast.hpp"
#include "Teuchos_TimeMonitor.hpp"

using namespace Thyra ;


// Constructors/initializers/accessors
  
MLPreconditionerFactory
::MLPreconditionerFactory()
  :epetraFwdOpViewExtractor_(Teuchos::rcp(new EpetraOperatorViewExtractorStd())),
   paramList_(defaultParameters(ML_DomainDecomposition))
{}
  
MLPreconditionerFactory
::MLPreconditionerFactory(const RefCountPtr<ParameterList>& params)
  :epetraFwdOpViewExtractor_(Teuchos::rcp(new EpetraOperatorViewExtractorStd())),
   paramList_(params)
{}
  
MLPreconditionerFactory
::MLPreconditionerFactory(const EMLProblemType& probType,
                          const ParameterList& revisions)
  :epetraFwdOpViewExtractor_(Teuchos::rcp(new EpetraOperatorViewExtractorStd())),
   paramList_(reviseDefaultList(*defaultParameters(probType), revisions))
{}
  
MLPreconditionerFactory
::MLPreconditionerFactory(const std::string& probType,
                          const ParameterList& revisions)
  : epetraFwdOpViewExtractor_(Teuchos::rcp(new EpetraOperatorViewExtractorStd())),
    paramList_(reviseDefaultList(*defaultParameters(probType), revisions))
{}
  
std::string MLPreconditionerFactory
::probToString(const EMLProblemType& probType) const
{
  switch(probType)
    {
    case ML_Maxwell:
      return "maxwell";
    case ML_DomainDecomposition:
      return "DD";
    case ML_DomainDecompositionML:
      return "DD-ML";
    default:
      return "SA";
    }
}
  
RefCountPtr<ParameterList> MLPreconditionerFactory
::defaultParameters(const EMLProblemType& probType) const 
{
  return defaultParameters(probToString(probType));
}

  
RefCountPtr<ParameterList> MLPreconditionerFactory
::defaultParameters(const string& probType) const 
{
  RefCountPtr<ParameterList> rtn = rcp(new ParameterList());
  int err = ML_Epetra::SetDefaults(probType, *rtn);
  TEST_FOR_EXCEPTION(err != 0, runtime_error,
                     "unable to find default parameters for problem type "
                     << probType);
  return rtn;
}


// Overridden from PreconditionerFactoryBase

bool MLPreconditionerFactory::isCompatible( const LinearOpSourceBase<double> &fwdOpSrc ) const
{
  Teuchos::RefCountPtr<const Epetra_Operator> epetraFwdOp;
  ETransp                                     epetraFwdOpTransp;
  EApplyEpetraOpAs                            epetraFwdOpApplyAs;
  EAdjointEpetraOp                            epetraFwdOpAdjointSupport;
  double                                      epetraFwdOpScalar;
  Teuchos::RefCountPtr<const LinearOpBase<double> > fwdOp = fwdOpSrc.getOp();
  epetraFwdOpViewExtractor_->getEpetraOpView(
                                             fwdOp,
                                             &epetraFwdOp,&epetraFwdOpTransp,
                                             &epetraFwdOpApplyAs,
                                             &epetraFwdOpAdjointSupport,
                                             &epetraFwdOpScalar
                                             );
  if( !dynamic_cast<const Epetra_RowMatrix*>(&*epetraFwdOp) )
    return false;
  return true;
}

bool MLPreconditionerFactory::applySupportsConj(EConj conj) const
{
  return true;
}

bool MLPreconditionerFactory::applyTransposeSupportsConj(EConj conj) const
{
  return false; // See comment below
}

Teuchos::RefCountPtr<PreconditionerBase<double> >
MLPreconditionerFactory::createPrec() const
{
  return Teuchos::rcp(new DefaultPreconditioner<double>());
}

void MLPreconditionerFactory::initializePrec(
                                             const Teuchos::RefCountPtr<const LinearOpSourceBase<double> >    &fwdOpSrc
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
  const Teuchos::RefCountPtr<Teuchos::FancyOStream> out       = this->getOStream();
  const Teuchos::EVerbosityLevel                    verbLevel = this->getVerbLevel();
  Teuchos::OSTab tab(out);
  if(out.get() && static_cast<int>(verbLevel) >= static_cast<int>(Teuchos::VERB_LOW))
    *out << "\nEntering Thyra::MLPreconditionerFactory::initializePrec(...) ...\n";

  Teuchos::RefCountPtr<const LinearOpBase<double> > fwdOp = fwdOpSrc->getOp();
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
  // Get the embedded ML_Epetra::MultiLevelPreconditioner object if it exists
  //
  Teuchos::RefCountPtr<ML_Epetra::MultiLevelPreconditioner> ml_precOp;
  if(epetra_precOp.get())
    ml_precOp = rcp_dynamic_cast<ML_Epetra::MultiLevelPreconditioner>(epetra_precOp->epetra_op(),true);
  //
  // Get the attached forward operator if it exists and make sure that it matches
  //
  if(ml_precOp.get()) {
    // ToDo: Get the forward operator and make sure that it matches what is
    // already being used!
  }
  //
  // Permform initialization if needed
  //
  const bool startingOver = (ml_precOp.get() == NULL);
  if(startingOver) 
    {
      if(out.get() && static_cast<int>(verbLevel) >= static_cast<int>(Teuchos::VERB_LOW))
        *out << "\nCreating the initial ML_Epetra::MultiLevelPreconditioner object...\n";
      timer.start(true);
      // Create the initial preconditioner
      ml_precOp = rcp(new ML_Epetra::MultiLevelPreconditioner(*epetraFwdRowMat,
                                                              *paramList_));
      
      timer.stop();
      if(out.get() && static_cast<int>(verbLevel) >= static_cast<int>(Teuchos::VERB_LOW))
        OSTab(out).o() <<"\n=> Creation time = "<<timer.totalElapsedTime()<<" sec\n";
      // RAB: Above, I am just passing a string to ML::Create(...) in order
      // get this code written.  However, in the future, it would be good to
      // copy the contents of what is in ML::Create(...) into a local
      // function and then use switch(...) to create the initial
      // ML_Epetra::MultiLevelPreconditioner object.  This would result in better validation
      // and faster code.
      // Set parameters if the list exists
      if(paramList_.get())
        TEST_FOR_EXCEPT(0!=ml_precOp->SetParameterList(*paramList_)); // This will create new sublist if it does not exist!
      // Initailize the structure for the preconditioner
      //      TEST_FOR_EXCEPT(0!=ml_precOp->Initialize());
    }
  //
  // Attach the epetraFwdOp to the ml_precOp to guarantee that it will not go away
  //
  set_extra_data(epetraFwdOp,"IFPF::epetraFwdOp",&ml_precOp,Teuchos::POST_DESTROY,false);
  //
  // Update the factorization
  //
  if(out.get() && static_cast<int>(verbLevel) >= static_cast<int>(Teuchos::VERB_LOW))
    *out << "\nComputing the factorization of the preconditioner ...\n";
  timer.start(true);
  TEST_FOR_EXCEPT(0!=ml_precOp->ComputePreconditioner());
  timer.stop();
  if(out.get() && static_cast<int>(verbLevel) >= static_cast<int>(Teuchos::VERB_LOW))
    OSTab(out).o() <<"\n=> Factorization time = "<<timer.totalElapsedTime()<<" sec\n";
  //
  // Compute the conditioner number estimate if asked
  //

  // ToDo: Implement

  //
  // Attach fwdOp to the ml_precOp
  //
  set_extra_data(fwdOp,"IFPF::fwdOp",&ml_precOp,Teuchos::POST_DESTROY,false);
  //
  // Initialize the output EpetraLinearOp
  //
  if(startingOver) {
    epetra_precOp = rcp(new EpetraLinearOp);
  }
  epetra_precOp->initialize(
                            ml_precOp
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
      << "\nLeaving Thyra::MLPreconditionerFactory::initializePrec(...) ...\n";
}

void MLPreconditionerFactory::uninitializePrec(
                                               PreconditionerBase<double>                          *prec
                                               ,Teuchos::RefCountPtr<const LinearOpSourceBase<double> >  *fwdOp
                                               ,ESupportSolveUse                                   *supportSolveUse
                                               ) const
{
  TEST_FOR_EXCEPT(true);
}

// Overridden from ParameterListAcceptor

void MLPreconditionerFactory::setParameterList(Teuchos::RefCountPtr<Teuchos::ParameterList> const& paramList)
{
  TEST_FOR_EXCEPT(paramList.get()==NULL);
  // Don't know how to validate an ML list
  //  paramList->validateParameters(*this->getValidParameters(),1);
  paramList_ = paramList;
}

Teuchos::RefCountPtr<Teuchos::ParameterList>
MLPreconditionerFactory::getParameterList()
{
  return paramList_;
}

Teuchos::RefCountPtr<Teuchos::ParameterList>
MLPreconditionerFactory::unsetParameterList()
{
  Teuchos::RefCountPtr<Teuchos::ParameterList> _paramList = paramList_;
  paramList_ = Teuchos::null;
  return _paramList;
}

Teuchos::RefCountPtr<const Teuchos::ParameterList>
MLPreconditionerFactory::getParameterList() const
{
  return paramList_;
}

Teuchos::RefCountPtr<const Teuchos::ParameterList>
MLPreconditionerFactory::getValidParameters() const
{
  if(!validPL_.get()) {
    validPL_ = defaultParameters(ML_DomainDecomposition);
    // Todo: the above is not really all of the valid parameters.  We need to
    // get ML to generate the entire list!
  }
  return validPL_;
}

// Public functions overridden from Teuchos::Describable

std::string MLPreconditionerFactory::description() const
{
  std::ostringstream oss;
  oss << "Thyra::MLPreconditionerFactory";
  return oss.str();
}

RefCountPtr<ParameterList> MLPreconditionerFactory
::reviseDefaultList(const ParameterList& defaults, const ParameterList& revisions) const
{
  RefCountPtr<ParameterList> rtn = rcp(new ParameterList(defaults));

  for (ParameterList::ConstIterator i=revisions.begin(); i != revisions.end(); i++)
    {
      rtn->setEntry(revisions.name(i), revisions.entry(i));
    }
  return rtn;
}

